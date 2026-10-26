# Congela as variáveis da estação A801 no instante de cada rodada (00 UTC):
# temperatura, orvalho, umidade, pressão e tendência de 24 h, vento (u, v),
# chuva e radiação das últimas 24 h, mínima das últimas 24 h.
# Saída: dados/processados/a801_estacao_init.rds (uma linha por init_date).
source("R/dados.R")

horario <- qc_inmet_completo(ler_inmet_completo("dados/brutos/inmet_a801"))
tab <- readRDS("dados/processados/tabela_tmax_gefs_a801.rds")
inits <- sort(unique(tab$init_date))
est <- variaveis_estacao(horario, inits)

cat("Fração de faltantes por variável e ano:\n")
cols <- grep("^est_", names(est), value = TRUE)
print(est[, lapply(.SD, function(x) round(mean(is.na(x)), 2)), by = .(ano = format(init_date, "%Y")),
          .SDcols = cols])
print(summary(est[, ..cols]))

saveRDS(est, "dados/processados/a801_estacao_init.rds")
