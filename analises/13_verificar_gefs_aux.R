# Verifica os arquivos do GEFS baixados. Ignora arquivos com menos de 2 min.
#   1. o arquivo abre e tem H1/H3/H7  -> se não, APAGA (o baixador refaz);
#   2. a Tmax do dia previsto na célula de Porto Alegre (máximo das 4 janelas de
#      6 h) é comparada com dados/brutos/gefs_v12_tmax_poa_2000_2019.csv;
#   3. cada variável está numa faixa física plausível.
# 2 e 3 só são REGISTRADOS (resultados/tabelas/gefs_aux_verificacao.csv), sem
# apagar: em ~4% dos dias (inverno) o CSV antigo dá Tmax 0,02–0,8 °C acima (ele
# usou outra janela para o "dia"), e há campos corrompidos na própria NOAA
# (Tmax > 50 °C; o ler_gefs() já os descartava), que baixar de novo não conserta.
# Uso: Rscript analises/13_verificar_gefs_aux.R
suppressMessages(library(data.table))
FAIXAS <- list(tmax2m = c(-15, 50), t2m = c(-15, 50), chuva = c(0, 400), q2m = c(0, 0.035),
               u10 = c(-40, 40), v10 = c(-40, 40), pnmm = c(94000, 105000), nuvens = c(0, 100),
               radsw = c(0, 1400), cape = c(0, 10000), pwat = c(0, 100), fsens = c(-600, 1200),
               flat = c(-600, 1500), solo = c(0, 1), t850 = c(-25, 40), q850 = c(0, 0.03),
               u850 = c(-60, 60), v850 = c(-60, 60), z500 = c(4500, 6200))

csv <- fread("dados/brutos/gefs_v12_tmax_poa_2000_2019.csv")
csv[, `:=`(chave = paste(member, format(as.Date(init_date), "%Y%m%d"), horizon), tc = forecast_tmax - 273.15)]
setkey(csv, chave)

fs <- list.files("dados/brutos/gefs_aux", pattern = "^[0-9]{8}\\.rds$", recursive = TRUE, full.names = TRUE)
fs <- fs[difftime(Sys.time(), file.mtime(fs), units = "mins") > 2]
problemas <- list(); quebrados <- character(); n_cmp <- 0; dif <- numeric()
for (f in fs) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(x) || !all(c("H1", "H3", "H7") %in% names(x)) || !length(x$H1)) { quebrados <- c(quebrados, f); next }
  motivo <- tryCatch({
    i <- which.min(abs(x$grade$lat + 30.05)); j <- which.min(abs(x$grade$lon + 51.17))
    m <- basename(dirname(dirname(f))); dia <- sub(".rds", "", basename(f), fixed = TRUE)
    for (H in c(1, 3, 7)) {
      cmp <- x[[paste0("H", H)]]
      ref <- csv[.(paste(m, dia, H)), tc]
      nossa <- max(sapply(cmp[grep("^tmax2m", names(cmp))], function(z) z[i, j]))
      if (length(ref) == 1 && !is.na(ref) && !is.na(nossa)) {
        n_cmp <- n_cmp + 1; dif <- c(dif, nossa - ref)
        if (abs(nossa - ref) > 1) stop(sprintf("Tmax POA H%d difere do CSV (%.2f vs %.2f)", H, nossa, ref))
      }
      for (nm in names(cmp)) {
        v <- sub("_h[0-9]+$", "", nm); z <- cmp[[nm]]
        if (all(is.na(z))) next                                   # campo ausente na NOAA
        r <- range(z, na.rm = TRUE)
        if (r[1] < FAIXAS[[v]][1] || r[2] > FAIXAS[[v]][2])
          stop(sprintf("%s H%d fora da faixa [%.4g, %.4g]", nm, H, r[1], r[2]))
      }
    }
    NA_character_
  }, error = function(e) conditionMessage(e))
  if (!is.na(motivo)) problemas[[f]] <- motivo
}
cat(length(fs), "arquivos verificados |", length(quebrados), "quebrados (apagados)\n")
cat(sprintf("Tmax POA vs CSV: %d comparações | idênticas (até 0,01 °C) %.1f%% | |dif| > 0,5 °C: %d\n",
            n_cmp, 100 * mean(abs(dif) <= 0.01), sum(abs(dif) > 0.5)))
cat(length(problemas), "arquivos com valor suspeito (registrados, não apagados)\n")
if (length(problemas)) {
  dir.create("resultados/tabelas", recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(arquivo = names(problemas), motivo = unlist(problemas)),
         "resultados/tabelas/gefs_aux_verificacao.csv")
  print(head(unlist(problemas), 10))
}
if (length(quebrados)) { file.remove(quebrados); cat("apagados; rode o baixador de novo para refazê-los\n") }
