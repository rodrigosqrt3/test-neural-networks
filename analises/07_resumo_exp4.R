# Resumo do experimento 4: piso de ruído entre sementes e ganho da informação nova.
# Uso: Rscript analises/07_resumo_exp4.R [pasta_exp4]   (padrão resultados/exp4)
# Saída: resultados/exp4_sementes.md
source("R/avaliacao.R")
library(data.table)

args <- commandArgs(TRUE)
p4 <- if (length(args)) args[[1]] else "resultados/exp4"
ler <- function(pasta, h, nome) {
  f <- file.path(pasta, paste0("h", h), paste0(nome, ".rds"))
  if (file.exists(f)) readRDS(f)$teste else NULL
}
linhas <- character(); add <- function(...) linhas <<- c(linhas, sprintf(...))
add("# Exp. 4: sementes e informação nova (teste 2016–2019)\n")
add("CRPS e tw95 médios no teste. Conjunto 'base' = exp3/base (seed0 100) + x4_base_s2..s5.\n")

for (h in c(1, 3, 7)) {
  base <- c(list(ler("resultados/exp3", h, "base")), lapply(2:5, function(s) ler(p4, h, sprintf("x4_base_s%d", s))))
  tudo <- lapply(1:5, function(s) ler(p4, h, sprintf("x4_tudo_s%d", s)))
  base <- Filter(Negate(is.null), base); tudo <- Filter(Negate(is.null), tudo)
  if (!length(base) || !length(tudo)) next
  cb <- sapply(base, function(d) mean(d$crps)); ct <- sapply(tudo, function(d) mean(d$crps))
  wb <- sapply(base, function(d) mean(d$tw_q95)); wt <- sapply(tudo, function(d) mean(d$tw_q95))
  add("\n## Horizonte %d\n", h)
  add("| conjunto | n sementes | CRPS média | CRPS dp entre sementes | min–max | tw95 média | tw95 dp |")
  add("|---|---|---|---|---|---|---|")
  add("| base | %d | %.4f | %.4f | %.4f–%.4f | %.4f | %.5f |", length(cb), mean(cb), sd(cb), min(cb), max(cb), mean(wb), sd(wb))
  add("| tudo | %d | %.4f | %.4f | %.4f–%.4f | %.4f | %.5f |", length(ct), mean(ct), sd(ct), min(ct), max(ct), mean(wt), sd(wt))
  # comparação pareada por dia, com a média dos escores entre sementes
  mb <- Reduce(`+`, lapply(base, `[[`, "crps")) / length(base)
  mt <- Reduce(`+`, lapply(tudo, `[[`, "crps")) / length(tudo)
  ic <- bootstrap_blocos(mt, mb, tamanho_bloco = 30, R = 1000)
  dm <- diebold_mariano(mt - mb)
  add("\nSkill CRPS de tudo vs base (médias entre sementes): %.2f%% [%.1f, %.1f], p DM = %.3f\n",
      100 * ic[[1]], 100 * ic[[2]], 100 * ic[[3]], dm[["p_valor"]])
  # cada bloco de informação contra a média das sementes base
  add("| config | CRPS | diferença vs base (desvios-padrão entre sementes) |")
  add("|---|---|---|")
  for (nm in c("x4_lags", "x4_traj", "x4_lagens", "x4_estacao", "x4_tudo_largo", "x4_emos", "x4_qrf")) {
    d <- ler(p4, h, nm); if (is.null(d)) next
    add("| %s | %.4f | %+.1f |", nm, mean(d$crps), (mean(d$crps) - mean(cb)) / sd(cb))
  }
}
writeLines(linhas, "resultados/exp4_sementes.md")
cat(linhas, sep = "\n")
