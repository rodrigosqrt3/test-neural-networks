# Prepara e CONGELA a tabela de modelagem (GEFSv12 x INMET A801).
# Rodar com o diretório de trabalho na raiz do projeto.
source("R/dados.R")

horario <- qc_inmet(ler_inmet("dados/brutos/inmet_a801"))
diaria  <- tmax_diaria(horario)
gefs    <- ler_gefs("dados/brutos/gefs_v12_tmax_poa_2000_2019.csv")

cat("Casos GEFS descartados (campo anômalo):\n"); print(attr(gefs, "descartados"))
cat(sprintf("A801: %s a %s | dias com Tmax válida: %d de %d\n",
            min(diaria$dia), max(diaria$dia), sum(!is.na(diaria$tmax)), nrow(diaria)))

tab <- atribuir_blocos(montar_tabela(gefs, diaria, horario))
print(tab[, .(n = .N, com_obs = sum(!is.na(tmax) & !is.na(obs_vespera))), by = .(horizon, bloco)])

# diagnóstico rápido: viés e correlação do GEFS bruto vs observação
print(tab[!is.na(tmax), .(vies = mean(ens_media - tmax), rmse = sqrt(mean((ens_media - tmax)^2)),
                          cor = cor(ens_media, tmax), dp_ens = mean(ens_dp)), by = horizon])

dir.create("dados/processados", showWarnings = FALSE)
saveRDS(diaria, "dados/processados/a801_tmax_diaria.rds")
saveRDS(tab, "dados/processados/tabela_tmax_gefs_a801.rds")
