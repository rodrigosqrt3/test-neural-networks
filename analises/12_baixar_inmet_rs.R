# Baixa os dados horários das estações automáticas do INMET no RS com o rmet
# (CRAN) e congela em dados/brutos/inmet_rs/<ano>.rds. Pode ser retomado.
# Uso: Rscript analises/12_baixar_inmet_rs.R [ano_ini] [ano_fim]
suppressMessages(library(rmet))
args <- commandArgs(TRUE)
anos <- if (length(args) >= 2) as.integer(args[1]):as.integer(args[2]) else 2000:2019

est <- inmet_stations()
rs <- est[est$state == "RS" & substr(est$station_code, 1, 1) == "A", ]   # automáticas
dir.create("dados/brutos/inmet_rs", recursive = TRUE, showWarnings = FALSE)
saveRDS(rs, "dados/brutos/inmet_rs/estacoes.rds")
cat(nrow(rs), "estações automáticas no RS\n")

for (a in anos) {
  arq <- sprintf("dados/brutos/inmet_rs/%d.rds", a)
  if (file.exists(arq)) next
  d <- inmet_get(years = a, stations = rs$station_code, tz = "UTC", quiet = TRUE)
  saveRDS(d, arq)
  cat(a, ":", nrow(d), "linhas,", length(unique(d$station_code)), "estações\n")
}
