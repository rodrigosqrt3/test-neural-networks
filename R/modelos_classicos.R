# Baselines probabilísticos. Todos devolvem uma matriz n x length(TAU_GRADE) de quantis.

# Climatologia: quantis empíricos das observações de TREINO numa janela de
# +-janela dias em torno do dia do ano.
clim_ajustar <- function(datas, y, janela = 15) {
  doy <- as.integer(format(datas, "%j"))
  tab <- t(vapply(1:366, function(d) {
    dist <- pmin(abs(doy - d), 366 - abs(doy - d))
    stats::quantile(y[dist <= janela], TAU_GRADE, na.rm = TRUE, names = FALSE, type = 8)
  }, numeric(length(TAU_GRADE))))
  structure(list(tab = tab), class = "clim")
}
predict.clim <- function(object, datas, ...) object$tab[as.integer(format(datas, "%j")), , drop = FALSE]

# Ensemble bruto: os m membros viram uma distribuição empírica
ens_quantis <- function(membros) {
  m <- ncol(membros)
  s <- t(apply(membros, 1, sort))
  s[, pmin(m, ceiling(TAU_GRADE * m)), drop = FALSE]
}

# EMOS / NGR (Gneiting et al. 2005) ajustada por CRPS mínimo:
#   Y ~ N(mu, sigma),  mu = X b,  log sigma = Z g
emos_ajustar <- function(X, Z, y) {
  X <- cbind(1, X); Z <- cbind(1, Z)
  perda <- function(par) {
    b <- par[seq_len(ncol(X))]; g <- par[-seq_len(ncol(X))]
    mean(scoringRules::crps_norm(y, X %*% b, exp(Z %*% g)))
  }
  ini <- c(stats::coef(stats::lm(y ~ X - 1)), log(stats::sd(y)), rep(0, ncol(Z) - 1))
  fit <- stats::optim(ini, perda, method = "BFGS", control = list(maxit = 2000))
  structure(list(par = fit$par, px = ncol(X), conv = fit$convergence), class = "emos")
}
predict.emos <- function(object, X, Z, ...) {
  X <- cbind(1, X); Z <- cbind(1, Z)
  b <- object$par[seq_len(object$px)]; g <- object$par[-seq_len(object$px)]
  mu <- drop(X %*% b); s <- drop(exp(Z %*% g))
  outer(s, stats::qnorm(TAU_GRADE)) + mu
}

# Quantile regression forest (Meinshausen 2006) via ranger
qrf_ajustar <- function(X, y, seed = 1, ...) {
  ranger::ranger(x = as.data.frame(X), y = y, quantreg = TRUE, num.trees = 1000,
                 min.node.size = 10, seed = seed, num.threads = 6, ...)
}
qrf_quantis <- function(fit, X) {
  q <- stats::predict(fit, as.data.frame(X), type = "quantiles", quantiles = TAU_GRADE)$predictions
  t(apply(q, 1, sort))
}
