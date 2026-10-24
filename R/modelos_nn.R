# Redes distribucionais em torch. Todas produzem a FUNÇÃO QUANTIL Q(tau | x)
# e são treinadas pela mesma perda: CRPS escrito como integral do escore
# quantílico, com peso opcional v(tau) nos quantis altos:
#
#   qwCRPS(F, y) = 2 * integral QS_tau(Q(tau), y) v(tau) dtau,
#   v(tau) = 1 + gama * 1{tau > tau_peso}
#
# qwCRPS é um escore próprio (Gneiting & Ranjan 2011) para qualquer v >= 0,
# ao contrário das perdas de "regressão desbalanceada" (EVL, DenseLoss, ...).
#
# Cabeças (forma de Q):
#   "drn"  Q = mu + sigma * Phi^-1(tau)                       (Rasp & Lerch 2018)
#   "bqn"  Q = polinômio de Bernstein crescente em tau        (Bremnes 2020)
#   "bgp"  corpo Bernstein em [0, tau0] + cauda GPD condicional acima de
#          u(x) = Q(tau0 | x):
#          Q(tau) = u + s(x) * (exp(xi(x) * L) - 1) / xi(x),  L = log((1-tau0)/(1-tau))
#          => P(Y > u(x) | x) = 1 - tau0 por construção, e acima de u(x) os
#             excessos são GPD(s(x), xi(x)). xi é global ou depende de x com
#             encolhimento lambda * (xi(x) - média)^2 (análogo à penalização do ERF).
#
# Várias estações num modelo só: n_estacoes > 0 acrescenta um embedding de
# estação (dimensão dim_emb) concatenado às entradas (Rasp & Lerch 2018); o
# índice da estação (1..n_estacoes) entra em forward(x, tau, id).

library(torch)

.bernstein_base <- function(tau, grau) {
  sapply(0:grau, function(k) stats::dbinom(k, grau, tau))   # length(tau) x (grau+1)
}

.coef_crescentes <- function(raw) {
  # raw: n x (grau+1) -> coeficientes não decrescentes (quantis não cruzam)
  inc <- nnf_softplus(raw[, 2:raw$size(2)])
  torch_cat(list(raw[, 1:1], raw[, 1:1] + torch_cumsum(inc, dim = 2)), dim = 2)
}

rede_distribucional <- nn_module(
  "rede_distribucional",
  initialize = function(n_in, cabeca, ocultas = c(64, 32), grau = 12, tau0 = 0.9,
                        xi_modo = "global", xi_lim = c(-0.6, 0.4), dropout = 0.1,
                        n_estacoes = 0, dim_emb = 4) {
    self$cabeca <- cabeca; self$grau <- grau; self$tau0 <- tau0
    self$xi_modo <- xi_modo; self$xi_lim <- xi_lim
    self$usa_emb <- n_estacoes > 0
    # n_estacoes + 1 linhas: seguro tanto para índices 1..n quanto 0..n-1
    if (self$usa_emb) self$emb <- nn_embedding(n_estacoes + 1, dim_emb)
    camadas <- list(); d <- n_in + if (self$usa_emb) dim_emb else 0
    for (h in ocultas) {
      camadas <- c(camadas, nn_linear(d, h), nn_elu(), nn_dropout(dropout)); d <- h
    }
    self$corpo <- do.call(nn_sequential, camadas)
    n_out <- switch(cabeca, drn = 2, bqn = grau + 1,
                    bgp = grau + 1 + 1 + (xi_modo == "x"))
    self$saida <- nn_linear(d, n_out)
    if (cabeca == "bgp" && xi_modo == "global") self$xi_raw <- nn_parameter(torch_zeros(1))
  },

  xi = function(z) {
    raw <- if (self$xi_modo == "global") self$xi_raw$expand(c(z$size(1), 1)) else z[, -1, drop = FALSE]
    self$xi_lim[1] + (self$xi_lim[2] - self$xi_lim[1]) * torch_sigmoid(raw)
  },

  # devolve os quantis n x length(tau) e (para bgp) xi n x 1
  forward = function(x, tau, id = NULL) {
    if (self$usa_emb) x <- torch_cat(list(x, self$emb(id)), dim = 2)
    z <- self$saida(self$corpo(x))
    tau_t <- torch_tensor(tau, dtype = torch_float())
    if (self$cabeca == "drn") {
      mu <- z[, 1:1]; sig <- nnf_softplus(z[, 2:2]) + 1e-3
      zt <- torch_tensor(stats::qnorm(tau), dtype = torch_float())$unsqueeze(1)
      return(list(q = mu + sig * zt, xi = NULL))
    }
    if (self$cabeca == "bqn") {
      a <- .coef_crescentes(z)
      B <- torch_tensor(.bernstein_base(tau, self$grau), dtype = torch_float())
      return(list(q = torch_matmul(a, B$t()), xi = NULL))
    }
    # bgp
    g <- self$grau
    a <- .coef_crescentes(z[, 1:(g + 1)])
    s <- nnf_softplus(z[, (g + 2):(g + 2)]) + 1e-3
    xi <- self$xi(z)
    corpo <- tau <= self$tau0
    q_corpo <- torch_matmul(a, torch_tensor(.bernstein_base(pmin(tau / self$tau0, 1), g),
                                            dtype = torch_float())$t())
    u <- a[, (g + 1):(g + 1)]                                # Q(tau0) = último coeficiente
    L <- torch_tensor(log((1 - self$tau0) / (1 - pmax(tau, self$tau0))),
                      dtype = torch_float())$unsqueeze(1)
    xiL <- xi * L
    # (exp(xi L) - 1)/xi estável perto de xi = 0
    fator <- torch_where(torch_abs(xi) < 1e-4, L * (1 + xiL / 2), torch_expm1(xiL) / xi)
    q_cauda <- u + s * fator
    m <- torch_tensor(corpo, dtype = torch_bool())$unsqueeze(1)
    list(q = torch_where(m, q_corpo, q_cauda), xi = xi, u = u, s = s)
  }
)

# Log-verossimilhança negativa GPD dos excessos acima de u(x), somada nos casos
# com y > u e dividida pelo tamanho do lote (contribuição por observação).
# Fora do suporte (xi < 0 e excesso além do ponto final) aplica-se uma
# penalidade linear, que mantém gradiente para trazer o ponto final para fora.
.perda_gpd_nll <- function(u, s, xi, y) {
  # u entra sem gradiente: o limiar é definido pela perda quantílica (é o
  # quantil tau0); senão a rede poderia subir u só para ter menos excessos
  z <- y$unsqueeze(2) - u$detach()
  exc <- (z > 0)$to(dtype = torch_float())
  z <- torch_clamp(z, min = 0)
  t <- 1 + xi * z / s
  t_seg <- torch_clamp(t, min = 1e-4)
  nll_xi <- torch_log(s) + (1 + 1 / xi) * torch_log(t_seg)
  nll_exp <- torch_log(s) + z / s
  nll <- torch_where(torch_abs(xi) < 1e-4, nll_exp, nll_xi) + 100 * nnf_relu(1e-4 - t)
  torch_sum(nll * exc) / y$size(1)
}

# w: pesos por amostra (opcional), normalizados no lote
.perda_qw <- function(q, y, tau, v, w = NULL) {
  tau_t <- torch_tensor(tau, dtype = torch_float())$unsqueeze(1)
  v_t <- torch_tensor(v, dtype = torch_float())$unsqueeze(1)
  d <- q - y$unsqueeze(2)
  qs <- ((d > 0)$to(dtype = torch_float()) - tau_t) * d
  por_caso <- torch_mean(qs * v_t, dim = 2)
  if (is.null(w)) return(2 * torch_mean(por_caso))
  2 * torch_sum(por_caso * w) / torch_sum(w)
}

# Adam / AdamW em C++ (optim_ignite_*): mesmo algoritmo, ~2,5x mais rápido que a versão
# em R, que faz um laço em R sobre os tensores de parâmetros a cada passo
otimizador <- function(rede, lr, wd) {
  if (wd > 0) optim_ignite_adamw(rede$parameters, lr = lr, weight_decay = wd)
  else optim_ignite_adam(rede$parameters, lr = lr)
}

# Grade de treino: passo 0.01 no corpo, 0.001 acima de 0.9 (cauda)
TAU_TREINO <- c((1:90 - 0.5) / 100, 0.9 + (1:100 - 0.5) / 1000)

# Largura de integração de cada nível (regra do ponto médio em grade não
# uniforme sobre [0, fim])
.larguras <- function(tau, fim = 1) {
  b <- c(0, (tau[-1] + tau[-length(tau)]) / 2, fim)
  diff(b)
}

# Treina UMA rede com early stopping na validação (mesma perda do treino).
#
# perda = "crps":    qwCRPS em toda a grade (corpo e cauda), como descrito acima.
# perda = "hibrida": só para "bgp". Escore quantílico integrado no CORPO
#   (tau <= tau0, o que inclui u(x) = Q(tau0 | x)) + lambda_nll x
#   verossimilhança GPD dos excessos acima de u(x). Limiar e cauda são
#   estimados juntos, num único passo (o EQRN faz em dois passos).
#
# pesos: pesos por amostra de treino (ex.: dar mais peso aos anos recentes).
# Xva = NULL: sem validação; treina exatamente `epocas` épocas.
# wd > 0: AdamW com esse weight decay (Loshchilov & Hutter 2019) em vez de Adam.
treinar_rede <- function(Xtr, ytr, Xva, yva, cabeca, gama = 0, tau_peso = 0.9,
                         lambda_xi = 0, seed = 1, epocas = 300, paciencia = 25,
                         lote = 128, lr = 2e-3, tau_treino = TAU_TREINO,
                         perda = c("crps", "hibrida"), lambda_nll = 1, pesos = NULL,
                         id_tr = NULL, id_va = NULL, wd = 0, ...) {
  perda <- match.arg(perda)
  torch_manual_seed(seed); set.seed(seed)
  fim <- 1
  if (perda == "hibrida") {
    stopifnot(cabeca == "bgp")
    tau0 <- if (is.null(list(...)$tau0)) 0.9 else list(...)$tau0
    tau_treino <- c(tau_treino[tau_treino < tau0], tau0)
    fim <- tau0
  }
  # pesos = largura de integração x v(tau); normalizados para que, com gama = 0,
  # mean(QS * v) aproxime (1/fim) integral QS dtau
  v <- .larguras(tau_treino, fim) * (1 + gama * (tau_treino > tau_peso))
  v <- v / sum(v) * length(tau_treino)
  funcao_perda <- function(out, y, w = NULL) {
    p <- .perda_qw(out$q, y, tau_treino, v, w)
    if (perda == "hibrida") p <- p + lambda_nll * .perda_gpd_nll(out$u, out$s, out$xi, y)
    p
  }
  rede <- rede_distribucional(ncol(Xtr), cabeca, ...)
  opt <- otimizador(rede, lr, wd)
  # alvo padronizado com média/desvio do TREINO (xi é invariante a isso)
  ym <- mean(ytr); ys <- stats::sd(ytr)
  ytr <- (ytr - ym) / ys
  xtr <- torch_tensor(Xtr, dtype = torch_float()); ytr_t <- torch_tensor(ytr, dtype = torch_float())
  w_t <- if (is.null(pesos)) NULL else torch_tensor(pesos / mean(pesos), dtype = torch_float())
  idt <- function(id) if (is.null(id)) NULL else torch_tensor(as.integer(id), dtype = torch_long())
  idtr_t <- idt(id_tr)
  com_validacao <- !is.null(Xva)
  if (com_validacao) {
    xva <- torch_tensor(Xva, dtype = torch_float())
    yva_t <- torch_tensor((yva - ym) / ys, dtype = torch_float())
    idva_t <- idt(id_va)
  }
  melhor <- Inf; melhor_estado <- NULL; melhor_epoca <- NA_integer_
  sem_melhora <- 0; hist <- numeric()
  n <- nrow(Xtr)
  for (ep in seq_len(epocas)) {
    rede$train()
    for (idx in split(sample.int(n), ceiling(seq_len(n) / lote))) {
      opt$zero_grad()
      out <- rede(xtr[idx, , drop = FALSE], tau_treino,
                  if (is.null(idtr_t)) NULL else idtr_t[idx, drop = FALSE])
      # drop = FALSE: um lote de 1 caso não pode virar escalar
      valor <- funcao_perda(out, ytr_t[idx, drop = FALSE], if (is.null(w_t)) NULL else w_t[idx, drop = FALSE])
      if (!is.null(out$xi) && lambda_xi > 0)
        valor <- valor + lambda_xi * torch_mean((out$xi - torch_mean(out$xi))^2)
      valor$backward(); opt$step()
    }
    if (!com_validacao) next
    rede$eval()
    pv <- with_no_grad(funcao_perda(rede(xva, tau_treino, idva_t), yva_t))$item()
    hist <- c(hist, pv)
    if (pv < melhor - 1e-5) {
      melhor <- pv; sem_melhora <- 0; melhor_epoca <- ep
      melhor_estado <- lapply(rede$state_dict(), function(p) p$clone())
    } else if ((sem_melhora <- sem_melhora + 1) >= paciencia) break
  }
  if (com_validacao) rede$load_state_dict(melhor_estado)
  rede$eval()
  list(rede = rede, hist_val = hist, melhor_val = melhor,
       epocas = if (com_validacao) length(hist) else epocas,
       melhor_epoca = if (com_validacao) melhor_epoca else epocas, ym = ym, ys = ys)
}

# Conjunto de redes (sementes diferentes); previsão = média dos quantis
# (Vincentização, recomendada por Schulz & Lerch 2022)
treinar_conjunto <- function(..., n_redes = 10, seed0 = 100) {
  lapply(seq_len(n_redes), function(k) treinar_rede(..., seed = seed0 + k))
}
prever_conjunto <- function(conj, X, id = NULL) {
  x <- torch_tensor(X, dtype = torch_float())
  i <- if (is.null(id)) NULL else torch_tensor(as.integer(id), dtype = torch_long())
  q <- 0
  for (r in conj) q <- q + with_no_grad(as_array(r$rede(x, TAU_GRADE, i)$q)) * r$ys + r$ym
  q <- q / length(conj)
  t(apply(q, 1, sort))
}
xi_conjunto <- function(conj, X, id = NULL) {
  x <- torch_tensor(X, dtype = torch_float())
  i <- if (is.null(id)) NULL else torch_tensor(as.integer(id), dtype = torch_long())
  sapply(conj, function(r) with_no_grad(as_array(r$rede(x, TAU_GRADE, i)$xi)))
}

# --- Rede com atenção entre estações (exp6) ------------------------------------------------
# Cada amostra é um DIA; as estações com caso naquele dia são os "tokens" (preditores +
# embedding da estação). Camadas de autoatenção (transformer pre-LN, Vaswani et al. 2017)
# deixam cada estação usar o que o GEFS e as observações dizem nas outras. O logit da
# atenção entre i e j perde beta_h * dist(i, j) (dist em 100 km; beta_h >= 0 aprendido por
# cabeça), um viés espacial no espírito das redes em grafo de Feik et al. (2024) e do
# ALiBi (Press et al. 2022). A cabeça distribucional é rede_distribucional, aplicada token
# a token, então drn/bqn/bgp e a perda são exatamente as mesmas da MLP.

# Linhas (casos) -> arrays dia x estação. id = índice da estação = posição do token.
# ordem: posição de cada linha na sequência de tokens presentes (dia, depois estação).
dias_tokens <- function(X, id, dia, n_est, y = NULL) {
  u <- sort(unique(dia)); di <- match(dia, u); n <- nrow(X); p <- ncol(X)
  A <- array(0, c(length(u), n_est, p))
  A[cbind(rep(di, p), rep(id, p), rep(seq_len(p), each = n))] <- X
  M <- matrix(FALSE, length(u), n_est); M[cbind(di, id)] <- TRUE
  Y <- NULL
  if (!is.null(y)) { Y <- matrix(NA_real_, length(u), n_est); Y[cbind(di, id)] <- y }
  list(X = A, M = M, Y = Y, ordem = order(order(di, id)))
}

# As projeções de saída (o, ffn_saida) começam em ZERO: o bloco começa como identidade, e
# a rede inteira começa igual à MLP (tronco + cabeça). A atenção só entra se reduzir a perda.
bloco_atencao <- nn_module(
  "bloco_atencao",
  initialize = function(dim, cabecas, dropout) {
    self$h <- cabecas; self$dh <- dim %/% cabecas
    self$ln1 <- nn_layer_norm(dim); self$ln2 <- nn_layer_norm(dim)
    self$qkv <- nn_linear(dim, 3 * dim); self$o <- nn_linear(dim, dim)
    self$beta <- nn_parameter(torch_zeros(cabecas))       # softplus(0) = 0,69 por 100 km
    saida <- nn_linear(2 * dim, dim)          # local: registrar duas vezes duplicaria o passo do otimizador
    self$ffn <- nn_sequential(nn_linear(dim, 2 * dim), nn_elu(), nn_dropout(dropout), saida)
    with_no_grad(for (l in list(self$o, saida)) { nn_init_zeros_(l$weight); nn_init_zeros_(l$bias) })
  },
  forward = function(h, ausente, dist) {                  # h: B x S x dim; ausente: B x S
    B <- h$size(1); S <- h$size(2)
    z <- self$qkv(self$ln1(h))$view(c(B, S, 3, self$h, self$dh))
    q <- z[, , 1, , ]$permute(c(1, 3, 2, 4)); k <- z[, , 2, , ]$permute(c(1, 3, 2, 4))
    v <- z[, , 3, , ]$permute(c(1, 3, 2, 4))                                 # B x H x S x dh
    a <- torch_matmul(q, k$transpose(3, 4)) / sqrt(self$dh)
    a <- a - nnf_softplus(self$beta)$view(c(1, self$h, 1, 1)) * dist$view(c(1, 1, S, S))
    a <- nnf_softmax(a$masked_fill(ausente$view(c(B, 1, 1, S)), -1e9), dim = 4)
    o <- torch_matmul(a, v)$permute(c(1, 3, 2, 4))$reshape(c(B, S, self$h * self$dh))
    h <- h + self$o(o)
    h + self$ffn(self$ln2(h))
  }
)

# Tronco = o corpo da MLP (n_in + embedding -> ocultas, ELU); a atenção corrige a
# representação de cada estação; cabeça = só a camada de saída de rede_distribucional.
rede_atencao <- nn_module(
  "rede_atencao",
  initialize = function(n_in, cabeca, n_estacoes, dist, ocultas = c(128, 64), camadas = 2, cabecas = 4,
                        dim_emb = 4, dropout = 0, ...) {
    self$emb <- nn_embedding(n_estacoes + 1, dim_emb)
    tronco <- list(); d <- n_in + dim_emb
    for (h in ocultas) { tronco <- c(tronco, nn_linear(d, h), nn_elu(), nn_dropout(dropout)); d <- h }
    self$tronco <- do.call(nn_sequential, tronco)
    self$dist <- nn_buffer(torch_tensor(dist / 100, dtype = torch_float()))
    self$blocos <- nn_module_list(lapply(seq_len(camadas), function(k) bloco_atencao(d, cabecas, dropout)))
    self$cabeca <- rede_distribucional(d, cabeca, ocultas = integer(0), dropout = 0, ...)
  },
  # x: B x S x n_in, presente: B x S (bool) -> quantis dos tokens presentes (dia, estação)
  forward = function(x, presente, tau) {
    B <- x$size(1); S <- x$size(2)
    ids <- torch_arange(1, S, dtype = torch_long())$unsqueeze(1)$expand(c(B, S))
    h <- self$tronco(torch_cat(list(x, self$emb(ids)), dim = 3))
    ausente <- presente$logical_not()
    for (k in seq_along(self$blocos)) h <- self$blocos[[k]](h, ausente, self$dist)
    self$cabeca(h[presente], tau)
  },
  beta = function() sapply(seq_along(self$blocos), function(k) as_array(nnf_softplus(self$blocos[[k]]$beta)))
)

# Treina UMA rede com atenção; T e Tva vêm de dias_tokens (com y). Tva = NULL: treina
# exatamente `epocas` épocas. lote = nº de dias por passo.
treinar_atencao <- function(T, Tva, cabeca, n_estacoes, dist, seed = 1, epocas = 300, paciencia = 10,
                            lote = 16, lr = 1e-3, wd = 0, tau_treino = TAU_TREINO, ...) {
  torch_manual_seed(seed); set.seed(seed)
  v <- .larguras(tau_treino); v <- v / sum(v) * length(tau_treino)
  ym <- mean(T$Y, na.rm = TRUE); ys <- stats::sd(as.vector(T$Y), na.rm = TRUE)
  tens <- function(T) list(x = torch_tensor(T$X, dtype = torch_float()), m = torch_tensor(T$M, dtype = torch_bool()),
                           y = torch_tensor(ifelse(T$M, (T$Y - ym) / ys, 0), dtype = torch_float()))
  tr <- tens(T)
  rede <- rede_atencao(dim(T$X)[3], cabeca, n_estacoes, dist, ...)
  opt <- otimizador(rede, lr, wd)
  perda_em <- function(D, idx) {
    m <- D$m[idx, , drop = FALSE]
    out <- rede(D$x[idx, , , drop = FALSE], m, tau_treino)
    list(p = .perda_qw(out$q, D$y[idx, , drop = FALSE][m], tau_treino, v), n = m$sum()$item())
  }
  if (!is.null(Tva)) va <- tens(Tva)
  melhor <- Inf; melhor_estado <- NULL; melhor_epoca <- NA_integer_; sem_melhora <- 0; hist <- numeric()
  n <- dim(T$X)[1]
  for (ep in seq_len(epocas)) {
    rede$train()
    for (idx in split(sample.int(n), ceiling(seq_len(n) / lote))) {
      opt$zero_grad(); perda_em(tr, idx)$p$backward(); opt$step()
    }
    if (is.null(Tva)) next
    rede$eval()
    pv <- with_no_grad({
      ps <- lapply(split(seq_len(dim(Tva$X)[1]), ceiling(seq_len(dim(Tva$X)[1]) / 128)), function(i) {
        r <- perda_em(va, i); c(r$p$item() * r$n, r$n)
      })
      s <- Reduce(`+`, ps); s[1] / s[2]
    })
    hist <- c(hist, pv)
    if (pv < melhor - 1e-5) {
      melhor <- pv; sem_melhora <- 0; melhor_epoca <- ep
      melhor_estado <- lapply(rede$state_dict(), function(p) p$clone())
    } else if ((sem_melhora <- sem_melhora + 1) >= paciencia) break
  }
  if (!is.null(Tva)) rede$load_state_dict(melhor_estado)
  rede$eval()
  list(rede = rede, hist_val = hist, melhor_val = melhor,
       melhor_epoca = if (is.null(Tva)) epocas else melhor_epoca, ym = ym, ys = ys, beta = rede$beta())
}

# Conjunto (média dos quantis); devolve as linhas na ordem original de dias_tokens
prever_atencao <- function(conj, T, dias_por_passo = 64) {
  x <- torch_tensor(T$X, dtype = torch_float()); m <- torch_tensor(T$M, dtype = torch_bool())
  partes <- split(seq_len(dim(T$X)[1]), ceiling(seq_len(dim(T$X)[1]) / dias_por_passo))
  q <- 0
  for (r in conj) q <- q + do.call(rbind, lapply(partes, function(i)
    with_no_grad(as_array(r$rede(x[i, , , drop = FALSE], m[i, , drop = FALSE], TAU_GRADE)$q)))) * r$ys + r$ym
  q <- t(apply(q / length(conj), 1, sort))
  q[T$ordem, , drop = FALSE]
}
