# Chuva diária observada nas 43 estações do INMET (fase 2 do PLANO.md), na MESMA janela
# do "dia previsto" do GEFS (script 11): de 00 UTC do dia D a 00 UTC de D+1 (21h local
# da véspera às 21h local). Cada registro das HH:00 acumula a hora anterior, então entram
# os registros com utc em (D 00:00, D+1 00:00]. O dia só vale com >= 22 das 24 horas.
#
# Uso: Rscript analises/23_preparar_chuva.R          (~5 min; lê os 20 anos do INMET)
# Saída: dados/processados/multi_obs_chuva.rds  (estacao, dia, chuva, n_horas)
source("R/dados.R")
t0 <- Sys.time()
est <- readRDS("dados/processados/multi_estacoes.rds")$estacao
horario <- rbindlist(lapply(2000:2019, function(a) {
  d <- as.data.table(readRDS(sprintf("dados/brutos/inmet_rs/%d.rds", a)))
  d[station_code %in% est]
}))
out <- rbindlist(lapply(split(horario, by = "station_code"), function(x) {
  h <- qc_inmet_completo(inmet_para_horario(x))
  h <- unique(h, by = "utc")
  h[, dia := as.Date(utc - 1)]                       # POSIXct em UTC: (D 00:00, D+1 00:00] -> D
  h[, .(chuva = if (sum(!is.na(chuva)) >= 22) sum(chuva, na.rm = TRUE) else NA_real_,
        n_horas = sum(!is.na(chuva))), by = dia][, estacao := x$station_code[1]]
}))
out <- out[!is.na(chuva)]
saveRDS(out[, .(estacao, dia, chuva, n_horas)], "dados/processados/multi_obs_chuva.rds")
cat(sprintf("%d estações, %d dias com chuva medida; %.1f%% dias secos (< 0,2 mm); p99 = %.1f mm; máx = %.1f mm | %.1f min\n",
            uniqueN(out$estacao), nrow(out), 100 * mean(out$chuva < 0.2), quantile(out$chuva, 0.99), max(out$chuva),
            as.numeric(Sys.time() - t0, units = "mins")))
