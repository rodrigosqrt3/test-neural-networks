# Baseline de série temporal pura: ARIMA na Tmax da A801, sem GEFS.
#
# Uso: Rscript analises/09_baseline_arima.R
# Saída: config "arima" em resultados/exp4/h<H> (split original) e em
#        resultados/exp4_dobra<k>/h<H>, no mesmo formato dos outros modelos.
#
# Modelo (ajustado só no treino de cada dobra):
#   Tmax(t) = m(t) + s(t) * z(t),  m e s^2 = harmônicos anuais (K = 3) por MQO,
#   z(t) ~ ARMA escolhido por auto.arima.
# Na rodada 00 UTC de init (21h local da véspera) conhecemos a Tmax completa até
# init-2 e a parcial de init-1 (obs_vespera, a mesma que a rede usa); o dia
# válido init+H fica H+1 passos à frente. Previsão normal com a variância do ARMA.
source("R/avaliacao.R"); source("R/experimento.R")
suppressMessages(library(forecast))

CORTES <- list(`1` = as.Date(c("2006-12-31", "2008-12-31", "2011-12-31")),
               `2` = as.Date(c("2009-12-31", "2011-12-31", "2015-12-31")))
diaria <- readRDS("dados/processados/a801_tmax_diaria.rds")[order(dia)]
stopifnot(all(diff(diaria$dia) == 1))
harm <- function(datas, K = 3) {
  w <- 2 * pi * as.integer(format(datas, "%j")) / 365.25
  do.call(cbind, lapply(1:K, function(k) cbind(sin(k * w), cos(k * w))))
}
HIST <- 400   # dias de histórico usados no filtro de Kalman a cada rodada

for (dobra in 0:2) for (H in c(1, 3, 7)) {
  dados <- carregar_horizonte(H)
  if (dobra > 0) {
    todos <- rbind(dados$tr, dados$va, dados$te); k <- CORTES[[as.character(dobra)]]
    dados$tr <- todos[valid_date <= k[1]]
    dados$va <- todos[valid_date > k[1] + 7 & valid_date <= k[2]]
    dados$te <- todos[valid_date > k[2] + 7 & valid_date <= k[3]]
  }
  tr <- dados$tr[order(init_date)]; va <- dados$va[order(init_date)]; te <- dados$te[order(init_date)]
  saida <- file.path("resultados", if (dobra > 0) paste0("exp4_dobra", dobra) else "exp4", paste0("h", H))
  dir.create(saida, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(saida, "arima.rds"))) next
  t0 <- Sys.time()

  # sazonalidade da média e da variância: só dias do treino
  Xh <- harm(diaria$dia); emtr <- diaria$dia <= max(tr$valid_date) & !is.na(diaria$tmax)
  fm <- lm(diaria$tmax[emtr] ~ Xh[emtr, ])
  m <- drop(cbind(1, Xh) %*% coef(fm))
  fv <- lm(residuals(fm)^2 ~ Xh[emtr, ])
  s <- sqrt(pmax(drop(cbind(1, Xh) %*% coef(fv)), 0.25))
  z <- (diaria$tmax - m) / s

  fit <- auto.arima(z[diaria$dia <= max(tr$valid_date)], d = 0, seasonal = FALSE,
                    max.p = 5, max.q = 3, stepwise = FALSE, approximation = FALSE)

  prever <- function(d) {
    pos <- match(d$init_date - 1, diaria$dia)       # índice do dia da véspera
    out <- t(vapply(seq_len(nrow(d)), function(i) {
      j <- pos[i]
      zh <- z[max(1, j - HIST):j]
      zh[length(zh)] <- (d$obs_vespera[i] - m[j]) / s[j]   # véspera: só a parcial conhecida
      f <- forecast(Arima(zh, model = fit), h = H + 1)
      alvo <- j + H + 1
      c(m[alvo] + s[alvo] * f$mean[H + 1],
        s[alvo] * (f$upper[H + 1, 1] - f$mean[H + 1]) / qnorm(0.5 + f$level[1] / 200))
    }, numeric(2)))
    stopifnot(all(diaria$dia[pos + H + 1] == d$valid_date))
    outer(out[, 2], qnorm(TAU_GRADE)) + out[, 1]
  }
  q_va <- prever(va); q_te <- prever(te)

  LIM <- list(validacao = limiares_para(tr, va$valid_date), teste = limiares_para(tr, te$valid_date))
  res <- list(config = "arima", horizonte = H, modelo = as.character(fit),
              validacao = cbind(data.table(data = va$valid_date), escores_por_caso(q_va, va$tmax, LIM$validacao)),
              teste = cbind(data.table(data = te$valid_date), escores_por_caso(q_te, te$tmax, LIM$teste)),
              segundos = as.numeric(Sys.time() - t0, units = "secs"))
  saveRDS(list(q_va = q_va, q_te = q_te, data_va = va$valid_date, data_te = te$valid_date),
          file.path(saida, ".grades_arima.rds"))
  saveRDS(res, file.path(saida, "arima.rds"))
  cat(sprintf("dobra %d h%d  %-28s CRPS val %.4f teste %.4f | tw95 %.4f | %.0f s\n", dobra, H,
              as.character(fit), mean(res$validacao$crps), mean(res$teste$crps),
              mean(res$teste$tw_q95), res$segundos))
}
