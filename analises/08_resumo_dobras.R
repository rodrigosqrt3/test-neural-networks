# Resumo da origem móvel: junta os testes das 3 dobras (2009–2011, 2012–2015,
# 2016–2019) e compara os modelos principais em 11 anos de teste.
# Uso: Rscript analises/08_resumo_dobras.R
# Saída: resultados/exp4_dobras.md
#
# Dobra 0 = split original: emos vem do exp. 1, base (semente 1) do exp. 3.
# Com duas sementes, o escore por caso é a média entre elas.
source("R/avaliacao.R")
library(data.table)

SUF <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1]] else ""   # "_smoke" para teste
arquivo <- function(dobra, h, nome) {
  if (dobra > 0) return(file.path(sprintf("resultados/exp4_dobra%d%s/h%d", dobra, SUF, h), paste0(nome, ".rds")))
  pasta <- switch(nome, emos = "exp1", x4_base_s1 = "exp3", paste0("exp4", SUF))
  if (nome == "x4_base_s1") nome <- "base"
  file.path(sprintf("resultados/%s/h%d", pasta, h), paste0(nome, ".rds"))
}
ler <- function(dobra, h, nomes) {
  fs <- vapply(nomes, function(n) arquivo(dobra, h, n), character(1))
  fs <- fs[file.exists(fs)]
  if (!length(fs)) return(NULL)
  ds <- lapply(fs, function(f) readRDS(f)$teste)
  d <- copy(ds[[1]])
  for (col in c("crps", "tw_q95", "tw_q99", "p_q95")) set(d, j = col, value = Reduce(`+`, lapply(ds, `[[`, col)) / length(ds))
  set(d, j = "dobra", value = dobra); set(d, j = "n_sementes", value = length(ds))
  d
}
MODELOS <- list(emos = "emos", arima = "arima", x4_emos = "x4_emos", x4_qrf = "x4_qrf",
                base = c("x4_base_s1", "x4_base_s2"), estacao = "x4_estacao",
                tudo = c("x4_tudo_s1", "x4_tudo_s2"), lagens = "x4_lagens",
                # 10_combinacoes.R: base = média dos QUANTIS das 2 sementes; c_* = base + outro
                base_2sem = "base_2sem", c_base_lagens = "c_base_lagens", c_base_tudo = "c_base_tudo",
                c_base_qrf = "c_base_qrf", c_base_emos = "c_base_emos", c_base_arima = "c_base_arima")

linhas <- character(); add <- function(...) linhas <<- c(linhas, sprintf(...))
add("# Origem móvel: 3 dobras, 11 anos de teste (2009–2019)\n")
add("Limiares sazonais e padronizações refeitos com o treino de cada dobra. Skill vs EMOS (definição do exp. 1) com bootstrap em blocos de 30 dias sobre os casos das 3 dobras juntos.\n")

for (h in c(1, 3, 7)) {
  dd <- lapply(MODELOS, function(m) rbindlist(lapply(0:2, function(k) ler(k, h, m))))
  dd <- Filter(function(d) nrow(d) > 0, dd)
  if (!"emos" %in% names(dd)) next
  ref <- dd$emos
  add("\n## Horizonte %d\n", h)
  add("| modelo | dobras | n | CRPS | tw95 | tw99 | exc95 obs | exc95 prev | skill CRPS %% [IC] | p DM | skill tw95 %% [IC] |")
  add("|---|---|---|---|---|---|---|---|---|---|---|")
  for (nm in names(dd)) {
    d <- dd[[nm]]
    chave <- paste(d$dobra, d$data); r <- ref[match(chave, paste(ref$dobra, ref$data))]
    ok <- !is.na(r$crps)
    s1 <- bootstrap_blocos(d$crps[ok], r$crps[ok]); s2 <- bootstrap_blocos(d$tw_q95[ok], r$tw_q95[ok])
    p <- if (nm == "emos") NA else diebold_mariano(d$crps[ok] - r$crps[ok])[["p_valor"]]
    add("| %s | %s | %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.2f [%.1f, %.1f] | %.3f | %.2f [%.1f, %.1f] |",
        nm, paste(sort(unique(d$dobra)), collapse = ","), nrow(d), mean(d$crps), mean(d$tw_q95), mean(d$tw_q99),
        mean(d$y > d$lim_q95), mean(d$p_q95),
        100 * s1[[1]], 100 * s1[[2]], 100 * s1[[3]], p, 100 * s2[[1]], 100 * s2[[2]], 100 * s2[[3]])
  }
  # consistência entre períodos: skill CRPS vs EMOS em cada dobra
  add("\nSkill CRPS (%%) vs EMOS por período de teste:\n")
  add("| modelo | 2009–2011 | 2012–2015 | 2016–2019 |")
  add("|---|---|---|---|")
  for (nm in setdiff(names(dd), "emos")) {
    v <- vapply(1:3, function(i) {
      k <- c(1, 2, 0)[i]
      a <- dd[[nm]][dobra == k]; b <- ref[dobra == k]
      if (!nrow(a) || !nrow(b)) return(NA_real_)
      b <- b[match(a$data, b$data)]
      100 * (1 - mean(a$crps) / mean(b$crps, na.rm = TRUE))
    }, numeric(1))
    add("| %s | %s |", nm, paste(ifelse(is.na(v), "—", sprintf("%.2f", v)), collapse = " | "))
  }
}
writeLines(linhas, sprintf("resultados/exp4_dobras%s.md", SUF))
cat(linhas, sep = "\n")

