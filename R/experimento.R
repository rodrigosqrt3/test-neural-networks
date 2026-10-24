# Funções comuns aos experimentos de pós-processamento da Tmax (A801 x GEFSv12):
# amostra comum, limiares, escores por caso e salvamento com checkpoint.
# Mantêm exatamente as mesmas definições de analises/02_experimento_tmax.R,
# para que os resultados de todos os experimentos sejam comparáveis.

library(data.table)

F_GEFS <- c("ens_media", "ens_dp", "ens_min", "ens_max", "doy_sen", "doy_cos")
F_OBS  <- c(F_GEFS, "obs_vespera")
F_LAGS <- c(F_OBS, "obs_lag2", "obs_lag3", "erro_ontem")
NIVEIS_Q <- c(0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99)

# Dados de um horizonte, nos mesmos dias usados pelo experimento 1
carregar_horizonte <- function(H, arquivo = "dados/processados/tabela_tmax_gefs_a801.rds") {
  tab <- readRDS(arquivo)
  membros <- attr(tab, "membros")
  d <- tab[horizon == H & !is.na(bloco) & !is.na(tmax)]
  d <- d[complete.cases(d[, c(F_LAGS, membros), with = FALSE])]
  d[, ano_frac := as.numeric(format(valid_date, "%Y")) +
        (as.numeric(format(valid_date, "%j")) - 1) / 365.25]
  list(tr = d[bloco == "treino"], va = d[bloco == "validacao"], te = d[bloco == "teste"],
       membros = membros)
}

limiares_para <- function(tr, datas) {
  out <- lapply(c(q90 = 0.90, q95 = 0.95, q99 = 0.99),
                function(p) limiar_sazonal(tr$valid_date, tr$tmax, p)(datas))
  out$abs35 <- rep(35, length(datas))
  out
}

escores_por_caso <- function(q, y, lim) {
  out <- data.table(y = y, crps = crps_quantis(q, y), pit = pit_quantis(q, y))
  for (nm in names(lim)) {
    out[[paste0("tw_", nm)]] <- twcrps_quantis(q, y, lim[[nm]])
    out[[paste0("p_", nm)]] <- prob_excedencia(q, lim[[nm]])
    out[[paste0("lim_", nm)]] <- lim[[nm]]
  }
  idx <- vapply(NIVEIS_Q, function(t) which.min(abs(TAU_GRADE - t)), integer(1))
  for (k in seq_along(NIVEIS_Q)) out[[sprintf("q%03d", round(NIVEIS_Q[k] * 100))]] <- q[, idx[k]]
  out
}

# Padronização com média/desvio de `ref`
padronizador <- function(ref, feats) {
  m <- sapply(ref[, ..feats], mean); s <- sapply(ref[, ..feats], sd)
  s[s == 0] <- 1
  function(x) sweep(sweep(as.matrix(x[, ..feats]), 2, m), 2, s, "/")
}

# Correção online do viés: desloca a previsão do caso i pela média dos
# resíduos (obs - mediana prevista) dos casos já verificados até a rodada,
# com peso 0.5^(idade / meia_vida). Um dia local D só é conhecido por
# completo depois das 03 UTC de D+1, então numa rodada 00 UTC de `init` só
# entram casos com valid_date <= init - 2.
correcao_online <- function(init, valid, y, mediana, meia_vida, defasagem = 2) {
  r <- y - mediana
  id <- as.numeric(init); vd <- as.numeric(valid)
  vapply(seq_along(y), function(i) {
    j <- which(vd <= id[i] - defasagem & vd > id[i] - 5 * meia_vida & !is.na(r))
    if (!length(j)) return(0)
    w <- 0.5^((id[i] - vd[j]) / meia_vida)
    sum(w * r[j]) / sum(w)
  }, numeric(1))
}
