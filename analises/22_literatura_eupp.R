# Compara os nossos modelos na EUPPBench (tarefa R2F) com os escores publicados.
# Uso: Rscript analises/22_literatura_eupp.R [pasta]   (padrão resultados/eupp; depois do script 18)
# Saída: resultados/<pasta>_literatura.md
#
# Feik, Lerch & Stühmer (2024, arXiv:2407.11050), Tabela 1, tarefa R2F (treino nos reforecasts,
# teste nas previsões de 2017–18), CRPS médio de T2m nas 122 estações, 2 casas decimais.
# O CRPS do ensemble bruto é a conferência de que o nosso conjunto de teste é o mesmo deles.
library(data.table)
args <- commandArgs(TRUE)
pasta <- if (length(args)) args[[1]] else "resultados/eupp"
nome <- basename(pasta)
FEIK <- data.table(horizonte = c(24, 72, 120),
                   ENS = c(1.12, 1.18, 1.38), DRN = c(0.61, 0.79, 1.11), SMRY = c(0.62, 0.79, 1.10),
                   DS = c(0.61, 0.78, 1.14), GAT = c(0.60, 0.78, 1.09))
res <- fread(sprintf("resultados/tabelas/%s_resumo.csv", nome))[, .(horizonte, modelo, CRPS)]
ens <- rbindlist(lapply(c(24, 72, 120), function(h) {
  f <- file.path(pasta, paste0("h", h), ".ensemble_bruto.rds")
  if (file.exists(f)) data.table(horizonte = h, modelo = "ensemble bruto", CRPS = mean(readRDS(f)$crps_ens))
}))
nossos <- dcast(rbind(ens, res[modelo %in% c("emos_loc", "emos_bst", "qrf", "gbm", "drn_mais", "stack")]),
                modelo ~ horizonte, value.var = "CRPS")
lit <- melt(FEIK, id.vars = "horizonte", variable.name = "modelo", value.name = "CRPS")
lit <- dcast(lit[, modelo := paste0(modelo, " (Feik et al.)")], modelo ~ horizonte, value.var = "CRPS")
tab <- rbind(nossos, lit, fill = TRUE)
fmt <- function(dt) {
  dt <- copy(dt)
  for (v in names(dt)[-1]) set(dt, j = v, value = sprintf("%.3f", dt[[v]]))
  c(paste0("| modelo | ", paste(paste0(names(dt)[-1], " h"), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
# skill do nosso melhor contra a DRN e a GAT publicadas
melhor <- res[modelo %in% c("drn_mais", "stack"), .SD[which.min(CRPS)], by = horizonte]
comp <- merge(melhor, FEIK[, .(horizonte, DRN, GAT, ENS)], by = "horizonte")
comp[, `:=`(vs_DRN = sprintf("%.1f%%", 100 * (1 - CRPS / DRN)), vs_GAT = sprintf("%.1f%%", 100 * (1 - CRPS / GAT)))]
ensv <- merge(ens, FEIK[, .(horizonte, ENS)], by = "horizonte")
L <- c("# EUPPBench (R2F): os nossos modelos e a literatura", "",
       "CRPS médio de T2m (°C), previsões de 2017–2018. Nós: 118 estações (2 sem observação nos Zarr e 2 com",
       "menos de 365 casos de treino ficam de fora); Feik et al. (2024, Tabela 1): 122 estações, 2 casas decimais.",
       "Se o CRPS do ensemble bruto bate com o deles, os conjuntos de teste são equivalentes.", "",
       fmt(tab), "",
       "## Conferência do ensemble bruto", "",
       paste0("| passo | nosso | Feik et al. |"), "|---|---|---|",
       sprintf("| %d h | %.3f | %.2f |", ensv$horizonte, ensv$CRPS, ensv$ENS), "",
       "## O nosso melhor modelo contra as redes publicadas", "",
       "| passo | modelo | CRPS | vs DRN | vs GAT (rede em grafo) |", "|---|---|---|---|---|",
       sprintf("| %d h | %s | %.3f | %s | %s |", comp$horizonte, comp$modelo, comp$CRPS, comp$vs_DRN, comp$vs_GAT), "",
       "Cuidados: (1) diferenças abaixo de ~0,01 não se resolvem com 2 casas decimais; (2) Feik et al. usam os",
       "mesmos preditores de Höhlein et al. (2024), e os nossos são parecidos mas não idênticos; (3) não há IC para",
       "os números publicados, então a comparação é pontual.")
writeLines(L, sprintf("resultados/%s_literatura.md", nome))
cat(L, sep = "\n")
