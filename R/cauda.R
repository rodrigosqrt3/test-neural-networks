# Funções comuns às análises que combinam os quantis guardados (scripts 37 e 38): leitura dos
# resultados, stack de Caruana, reconstrução da grade de 999 níveis, ξ da validação, recalibração,
# P(Y > t) vetorizada e comparação com bootstrap + Diebold–Mariano.
# Precisa de R/avaliacao.R, R/experimento.R e R/multi.R (TAU_FINO) carregados antes.

ler <- function(p, modelos, padrao = ".") {
  fs <- list.files(p, pattern = "^[^.].*\\.rds$", recursive = TRUE, full.names = TRUE)
  fs <- fs[grepl(sprintf("^(%s)_f[123]\\.rds$", paste(modelos, collapse = "|")), basename(fs)) & grepl(padrao, fs)]
  lapply(fs, function(f) { r <- readRDS(f); r$arq_q <- file.path(dirname(f), sprintf(".quantis_%s_%s.rds", r$config, r$dobra)); r })
}
# pesos de Caruana et al. (2004): seleção gulosa com reposição; `escore` é minimizado na validação
caruana <- function(qs, y, passos = 25, escore = function(q, y) crps_quantis(q, y, TAU_FINO)) {
  cont <- setNames(integer(length(qs)), names(qs)); soma <- 0
  for (p in seq_len(passos)) {
    cv <- vapply(qs, function(q) mean(escore((soma + q) / p, y)), numeric(1))
    j <- which.min(cv); cont[j] <- cont[j] + 1L; soma <- soma + qs[[j]]
  }
  cont / passos
}
# Vincentização: média ponderada dos quantis (`qual` = "q_va" ou "q_te")
vincentizar <- function(q, w, qual) Reduce(`+`, Map(function(x, wi) wi * x[[qual]], q[names(w)], w))
lim_de <- function(d) list(q90 = d$lim_q90, q95 = d$lim_q95, q99 = d$lim_q99, abs35 = d$lim_abs35)

# --- correção leve por cima de um stack (exp14 parte B, exp15) ---------------------------------------
# viés harmônico por estação, ajustado nas linhas `tr` e avaliado em todas
harm <- function(v, est, doy, tr) {
  A <- cbind(1, sin(2 * pi * doy / 365.25), cos(2 * pi * doy / 365.25))
  glob <- stats::lm.fit(A[tr, ], v[tr])$coefficients
  out <- drop(A %*% glob)
  for (s in unique(est)) {
    i <- which(est == s); j <- intersect(i, which(tr))
    if (length(j) >= 60) out[i] <- drop(A[i, , drop = FALSE] %*% stats::lm.fit(A[j, ], v[j])$coefficients)
  }
  out
}
# q'(τ) = m + Xβ + exp(Sγ) (q(τ) − m), por CRPS mínimo nas linhas tr (X[, 2] = viés, começa em 1);
# devolve os quantis corrigidos de todas as linhas e os parâmetros
corrigir <- function(q, y, X, S, tr) {
  m <- rowMeans(q); dq <- q - m
  prev <- function(p, i) m[i] + drop(X[i, , drop = FALSE] %*% p[seq_len(ncol(X))]) +
    exp(drop(S[i, , drop = FALSE] %*% p[-seq_len(ncol(X))])) * dq[i, , drop = FALSE]
  it <- which(tr)
  perda <- function(p) mean(crps_quantis(prev(p, it), y[it], TAU_FINO))
  p0 <- c(0, 1, rep(0, ncol(X) - 2), rep(0, ncol(S)))
  fit <- stats::optim(p0, perda, method = "BFGS", control = list(maxit = 200, reltol = 1e-7))
  list(q = prev(fit$par, seq_len(nrow(q))), par = setNames(fit$par, c(colnames(X), paste0("s_", colnames(S)))),
       conv = fit$convergence)
}

# --- quantis em níveis arbitrários a partir da grade de 99 ------------------------------------------
# Interior: linear entre os níveis 0,01..0,99. Abaixo de 0,01: a extrapolação normal dos scripts 30/33.
# Acima de 0,99: "normal" (= scripts 30/33) ou "gpd": Q(τ) = Q99 + A[(0,01/(1−τ))^ξ − 1], com
# A = (Q99 − Q98)/(1 − 2^−ξ) (a cauda GPD que passa por Q98 e Q99) e ξ fixo por modelo e dobra.
quantis_em <- function(q, niveis, cauda = c("normal", "gpd"), xi = 0, piso = -Inf) {
  cauda <- match.arg(cauda)
  n <- nrow(q); niveis <- matrix(niveis, n, if (is.matrix(niveis)) ncol(niveis) else length(niveis), byrow = !is.matrix(niveis))
  pos <- niveis * 100; lo <- pmin(pmax(floor(pos), 1), 98); w <- pos - lo
  linha <- rep(seq_len(n), ncol(niveis))
  out <- matrix(q[cbind(linha, as.vector(lo))] * (1 - as.vector(w)) + q[cbind(linha, as.vector(lo) + 1)] * as.vector(w), n)
  inf <- niveis < 0.01
  if (any(inf)) {
    b <- (q[, 2] - q[, 1]) / (stats::qnorm(0.02) - stats::qnorm(0.01))
    out[inf] <- (q[, 1] + b * (stats::qnorm(niveis) - stats::qnorm(0.01)))[inf]
  }
  sup <- niveis > 0.99
  if (any(sup)) {
    d1 <- pmax(q[, 99] - q[, 98], 0)
    v <- if (cauda == "normal") {
      q[, 99] + d1 / (stats::qnorm(0.99) - stats::qnorm(0.98)) * (stats::qnorm(niveis) - stats::qnorm(0.99))
    } else if (abs(xi) < 1e-6) {
      q[, 99] + d1 / log(2) * log(0.01 / (1 - niveis))
    } else {
      q[, 99] + d1 / (1 - 2^(-xi)) * ((0.01 / (1 - niveis))^xi - 1)
    }
    out[sup] <- v[sup]
  }
  out <- t(apply(out, 1, sort))
  pmax(out, piso)
}
# ξ de Pickands por caso (níveis 0,96 / 0,98 / 0,99, razão 2 nas probabilidades de excedência),
# resumido pela mediana na VALIDAÇÃO (casos com espaçamentos > 0,1 mm), limitado a [−0,2; 0,6]
xi_validacao <- function(q) {
  d1 <- q[, 99] - q[, 98]; d2 <- q[, 98] - q[, 96]
  ok <- d1 > 0.1 & d2 > 0.1
  if (sum(ok) < 100) return(0)
  min(max(stats::median(log2(d1[ok] / d2[ok])), -0.2), 0.6)
}
# recalibração: níveis de origem G^-1(τ) para cada τ da grade final
niveis_recal <- function(q_va, y_va) {
  u <- pit_quantis(q_va, y_va, TAU_FINO)
  pmin(pmax(stats::quantile(u, TAU_GRADE, type = 8, names = FALSE), 1e-4), 1 - 1e-4)
}
# P(Y > t) de uma grade (n x 999) para t escalar, como prob_excedencia, vetorizado
p_exc <- function(G, t) {
  k <- rowSums(G <= t); n <- ncol(G); p <- numeric(nrow(G))
  p[k == 0] <- 1 - TAU_GRADE[1] / 2; p[k == n] <- (1 - TAU_GRADE[n]) / 2
  m <- which(k > 0 & k < n)
  if (length(m)) {
    a <- G[cbind(m, k[m])]; b <- G[cbind(m, k[m] + 1)]
    f <- ifelse(b > a, (t - a) / (b - a), 0)
    p[m] <- 1 - (TAU_GRADE[k[m]] + f * (TAU_GRADE[2] - TAU_GRADE[1]))
  }
  p
}

comparar <- function(teste, A, B, col, h, estacoes = NULL, rotulo = NULL) {
  a <- teste[modelo == A & horizonte == h]; b <- teste[modelo == B & horizonte == h]
  if (!is.null(estacoes)) { a <- a[estacao %in% estacoes]; b <- b[estacao %in% estacoes] }
  if (!nrow(a) || !nrow(b)) return(NULL)
  chave <- c("horizonte", "dobra", "estacao", "data")
  p <- merge(a[, c(chave, col), with = FALSE], b[, c(chave, col), with = FALSE], by = chave, suffixes = c("_a", "_b"))
  if (nrow(p) != nrow(a) || nrow(p) != nrow(b)) stop(sprintf("%s e %s não foram avaliados nas mesmas linhas (h%d)", A, B, h))
  setnames(p, paste0(col, c("_a", "_b")), c("a", "b"))
  dia <- p[, .(a = sum(a), b = sum(b), d = mean(a - b)), by = data][order(data)]
  if (sum(dia$b) == 0) return(NULL)
  ic <- bootstrap_blocos(dia$a, dia$b, tamanho_bloco = 30, R = 1000)
  data.table(horizonte = h, grupo = if (is.null(rotulo)) "todas" else rotulo, A = A, B = B, escore = col,
             skill = 100 * ic[[1]], ic_inf = 100 * ic[[2]], ic_sup = 100 * ic[[3]], p_dm = diebold_mariano(dia$d)[["p_valor"]])
}
md <- function(dt, dig = 4) {
  dt <- copy(dt)
  for (v in names(dt)) if (is.double(dt[[v]])) set(dt, j = v, value = formatC(dt[[v]], digits = dig, format = "fg"))
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
ic_txt <- function(d) d[, .(horizonte, grupo, comparacao = paste(A, "vs", B), escore, skill = sprintf("%.2f", skill),
                            IC95 = sprintf("[%.1f, %.1f]", ic_inf, ic_sup), p_DM = sprintf("%.3f", p_dm))]
