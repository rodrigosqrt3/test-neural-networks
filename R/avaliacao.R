# Avaliação de previsões probabilísticas representadas por grades de quantis.
#
# Toda previsão (de qualquer modelo) é uma matriz n x K de quantis nos níveis
# TAU_GRADE. Assim CRPS, twCRPS, Brier e PIT são calculados do mesmo jeito para
# todos os modelos.
#
# CRPS(F, y) = 2 * integral_0^1 QS_tau(F^-1(tau), y) dtau   (Gneiting & Ranjan 2011)
# twCRPS com peso 1{z >= t} = CRPS das versões censuradas max(X, t), max(y, t)
# (Allen et al. 2023), cujos quantis são max(Q(tau), t).

TAU_GRADE <- (1:999) / 1000

# Escore quantílico (pinball) elemento a elemento: q é n x K, y tem tamanho n
.pinball <- function(q, y, tau = TAU_GRADE) {
  d <- q - y                                    # recicla y por linha
  sweep((d > 0) * 1, 2, tau) * d                # (1{y < q} - tau) * (q - y)
}

crps_quantis <- function(q, y, tau = TAU_GRADE) {
  2 * rowMeans(.pinball(q, y, tau))
}

# t pode ser escalar ou vetor (limiar sazonal, um por caso)
twcrps_quantis <- function(q, y, t, tau = TAU_GRADE) {
  t <- rep_len(t, length(y))
  crps_quantis(pmax(q, t), pmax(y, t), tau)
}

# P(Y > t) a partir da grade de quantis (interpolação linear em tau)
prob_excedencia <- function(q, t, tau = TAU_GRADE) {
  t <- rep_len(t, nrow(q))
  vapply(seq_len(nrow(q)), function(i) {
    qi <- q[i, ]
    if (t[i] <= qi[1]) return(1 - tau[1] / 2)
    if (t[i] >= qi[length(qi)]) return((1 - tau[length(tau)]) / 2)
    1 - stats::approx(qi, tau, xout = t[i], ties = "ordered")$y
  }, numeric(1))
}

brier <- function(p, y, t) mean((p - (y > t))^2)

# PIT aleatorizado a partir da grade (para histogramas / calibração)
pit_quantis <- function(q, y, tau = TAU_GRADE) {
  vapply(seq_len(nrow(q)), function(i) {
    qi <- q[i, ]
    if (y[i] <= qi[1]) return(stats::runif(1, 0, tau[1]))
    if (y[i] >= qi[length(qi)]) return(stats::runif(1, tau[length(tau)], 1))
    stats::approx(qi, tau, xout = y[i], ties = "ordered")$y
  }, numeric(1))
}

# Limiar sazonal: quantil p da observação de TREINO numa janela de +-janela dias
# em torno de cada dia do ano. Retorna uma função doy -> limiar.
limiar_sazonal <- function(datas_treino, y_treino, p, janela = 15) {
  doy_t <- as.integer(format(datas_treino, "%j"))
  lim <- vapply(1:366, function(d) {
    dist <- pmin(abs(doy_t - d), 366 - abs(doy_t - d))
    stats::quantile(y_treino[dist <= janela], p, na.rm = TRUE, names = FALSE)
  }, numeric(1))
  function(datas) lim[as.integer(format(datas, "%j"))]
}

# Diebold-Mariano com variância HAC (Newey-West) para diferenças de escore
# d = escore_A - escore_B (negativo => A melhor)
diebold_mariano <- function(d, lag = NULL) {
  d <- d[is.finite(d)]
  n <- length(d)
  if (is.null(lag)) lag <- floor(4 * (n / 100)^(2 / 9))
  dc <- d - mean(d)
  gam <- vapply(0:lag, function(k) sum(dc[(k + 1):n] * dc[1:(n - k)]) / n, numeric(1))
  v <- gam[1] + 2 * sum((1 - (1:lag) / (lag + 1)) * gam[-1])
  est <- mean(d) / sqrt(v / n)
  c(media = mean(d), estatistica = est, p_valor = 2 * stats::pnorm(-abs(est)))
}

# IC por bootstrap em blocos (moving block) para a média de um escore ou
# para o skill relativo 1 - mean(a)/mean(b)
bootstrap_blocos <- function(a, b = NULL, tamanho_bloco = 30, R = 1000, seed = 1) {
  set.seed(seed)
  n <- length(a)
  nb <- ceiling(n / tamanho_bloco)
  stat <- function(idx) if (is.null(b)) mean(a[idx]) else 1 - mean(a[idx]) / mean(b[idx])
  reps <- replicate(R, {
    ini <- sample.int(n - tamanho_bloco + 1, nb, replace = TRUE)
    idx <- as.vector(outer(0:(tamanho_bloco - 1), ini, `+`))[1:n]
    stat(idx)
  })
  c(estimativa = stat(seq_len(n)), stats::quantile(reps, c(0.025, 0.975)))
}

# Tabela-resumo: um data.frame por modelo com CRPS, twCRPS e Brier em vários limiares
resumo_escores <- function(lista_quantis, y, limiares) {
  do.call(rbind, lapply(names(lista_quantis), function(m) {
    q <- lista_quantis[[m]]
    linha <- data.frame(modelo = m, CRPS = mean(crps_quantis(q, y)))
    for (nm in names(limiares)) {
      t <- limiares[[nm]]
      linha[[paste0("twCRPS_", nm)]] <- mean(twcrps_quantis(q, y, t))
      linha[[paste0("Brier_", nm)]] <- brier(prob_excedencia(q, t), y, t)
    }
    linha
  }))
}
