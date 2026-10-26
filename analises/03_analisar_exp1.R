# Analisa os resultados de 02_experimento_tmax.R.
# Uso: Rscript analises/03_analisar_exp1.R [pasta]   (padrão resultados/exp1)
# Saídas: resultados/tabelas/exp1_*.csv, resultados/figuras/exp1_*.png e
#         resultados/exp1_resumo.md

source("R/avaliacao.R")
library(data.table)

# Uso: Rscript analises/03_analisar_exp1.R [pasta1 pasta2 ...] [--saida=nome]
args <- commandArgs(TRUE)
nome_saida <- sub("^--saida=", "", grep("^--saida=", args, value = TRUE))
pastas <- grep("^--", args, value = TRUE, invert = TRUE)
if (!length(pastas)) pastas <- "resultados/exp1"
pasta <- paste(pastas, collapse = " + ")
if (!length(nome_saida)) nome_saida <- "exp1"
sufixo <- if (any(grepl("smoke", pastas))) "_smoke" else ""
dir.create("resultados/tabelas", recursive = TRUE, showWarnings = FALSE)
dir.create("resultados/figuras", recursive = TRUE, showWarnings = FALSE)

arquivos <- unlist(lapply(pastas, list.files, pattern = "^[^.].*\\.rds$", recursive = TRUE,
                          full.names = TRUE))
res <- lapply(arquivos, readRDS)
cat(length(res), "configurações encontradas em", pasta, "\n")
na_amostra <- vapply(res, function(r) isTRUE(r$validacao_na_amostra), logical(1))
names(na_amostra) <- vapply(res, function(r) r$config, character(1))

qs <- function(q, y, tau) ((y < q) - tau) * (q - y)

# --- tabela de escores médios ---------------------------------------------------
resumir <- function(r, bloco) {
  d <- r[[bloco]]
  data.table(
    horizonte = r$horizonte, config = r$config, bloco = bloco, n = nrow(d),
    CRPS = mean(d$crps), tw90 = mean(d$tw_q90), tw95 = mean(d$tw_q95), tw99 = mean(d$tw_q99),
    tw35 = mean(d$tw_abs35),
    Brier95 = mean((d$p_q95 - (d$y > d$lim_q95))^2),
    Brier99 = mean((d$p_q99 - (d$y > d$lim_q99))^2),
    QS95 = mean(qs(d$q095, d$y, 0.95)), QS99 = mean(qs(d$q099, d$y, 0.99)),
    # calibração: cobertura do intervalo central de 80% e frequência observada
    # vs prevista de exceder o q95 sazonal (razão 1 = calibrado na cauda)
    cobertura80 = mean(d$pit > 0.1 & d$pit < 0.9),
    pit_inf5 = mean(d$pit < 0.05), pit_sup5 = mean(d$pit > 0.95),   # ideal: 0.05 cada
    exc95_obs = mean(d$y > d$lim_q95), exc95_prev = mean(d$p_q95),
    segundos = if (is.null(r$segundos)) NA_real_ else r$segundos
  )
}
tab <- rbindlist(c(lapply(res, resumir, "validacao"), lapply(res, resumir, "teste")))

# --- skill vs EMOS com IC por bootstrap em blocos e Diebold-Mariano --------------
compara <- function(r, ref, bloco, coluna) {
  a <- r[[bloco]][[coluna]]; b <- ref[[bloco]][[coluna]]
  stopifnot(identical(r[[bloco]]$data, ref[[bloco]]$data))
  ic <- bootstrap_blocos(a, b, tamanho_bloco = 30, R = 1000)
  dm <- diebold_mariano(a - b)
  data.table(skill = ic[[1]], ic_inf = ic[[2]], ic_sup = ic[[3]], p_dm = dm[["p_valor"]])
}
comp <- rbindlist(lapply(res, function(r) {
  ref <- Filter(function(x) x$config == "emos" && x$horizonte == r$horizonte, res)
  if (!length(ref) || r$config == "emos") return(NULL)
  rbindlist(lapply(c("crps", "tw_q95", "tw_q99"), function(col) {
    cbind(data.table(horizonte = r$horizonte, config = r$config, escore = col),
          compara(r, ref[[1]], "teste", col))
  }))
}))

tabela_csv <- function(nome) sprintf("resultados/tabelas/%s_%s%s.csv", nome_saida, nome, sufixo)
fwrite(tab, tabela_csv("escores"))
fwrite(comp, tabela_csv("skill_vs_emos"))

# --- seleção pela validação (sem olhar o teste) ----------------------------------
selecao <- rbindlist(lapply(split(tab, tab$horizonte), function(t) {
  # modelos que treinaram com a validação não podem ser escolhidos por ela
  v <- t[bloco == "validacao" & !config %in% names(na_amostra)[na_amostra]]
  rbindlist(lapply(c("CRPS", "tw95"), function(criterio) {
    melhor <- v$config[which.min(v[[criterio]])]
    cbind(data.table(criterio_validacao = criterio, escolhido = melhor),
          t[bloco == "teste" & config == melhor, .(horizonte, CRPS, tw95, tw99, Brier95)])
  }))
}))

# --- figuras ----------------------------------------------------------------------
fig <- function(nome, expr, w = 1600, h = 900) {
  png(sprintf("resultados/figuras/%s_%s%s.png", nome_saida, nome, sufixo), width = w, height = h, res = 150)
  on.exit(dev.off())
  expr
}

# Fronteira CRPS x twCRPS95 para o peso na cauda (gama)
fig("pareto_gama", {
  hs <- sort(unique(tab$horizonte))
  par(mfrow = c(1, length(hs)), mar = c(4.5, 4.5, 3, 1))
  for (h in hs) {
    t <- tab[horizonte == h & bloco == "teste" & !config %in% c("climatologia", "ensemble_bruto")]
    plot(t$CRPS, t$tw95, pch = 19, col = "grey70", xlab = "CRPS (teste)",
         ylab = "twCRPS q95 (teste)", main = paste("Horizonte", h, "dia(s)"))
    for (fam in c("bqn", "bgp")) {
      cfg <- c(paste0(fam, "_obs"), paste0(fam, "_obs_g", c(1, 4, 9, 19)))
      s <- t[match(cfg, config)][!is.na(CRPS)]
      lines(s$CRPS, s$tw95, type = "b", pch = 19, col = if (fam == "bqn") "#2471A3" else "#C0392B")
    }
    e <- t[config == "emos"]
    points(e$CRPS, e$tw95, pch = 17, cex = 1.5)
    legend("topleft", c("BQN, gama 0-19", "BGP, gama 0-19", "EMOS", "outros"),
           col = c("#2471A3", "#C0392B", "black", "grey70"), pch = c(19, 19, 17, 19), bty = "n")
  }
})

# Histogramas PIT (calibração) para modelos-chave, teste, horizonte 1
fig("pit_h1", {
  chave <- c("ensemble_bruto", "emos", "qrf_obs", "drn_obs", "bqn_obs", "bgp_obs")
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  for (cfg in chave) {
    r <- Filter(function(x) x$config == cfg && x$horizonte == min(tab$horizonte), res)
    if (!length(r)) next
    hist(r[[1]]$teste$pit, breaks = seq(0, 1, 0.05), main = cfg, xlab = "PIT",
         col = "grey80", border = "white", freq = FALSE)
    abline(h = 1, lty = 2)
  }
})

# Diagrama de confiabilidade para exceder o q95 sazonal
fig("confiabilidade_q95", {
  chave <- c("emos", "qrf_obs", "drn_obs", "bgp_obs", "bgp_obs_g9")
  cores <- c("black", "#7F8C8D", "#2471A3", "#C0392B", "#E67E22")
  hs <- sort(unique(tab$horizonte))
  par(mfrow = c(1, length(hs)), mar = c(4.5, 4.5, 3, 1))
  for (h in hs) {
    plot(0:1, 0:1, type = "n", xlab = "Probabilidade prevista", ylab = "Frequência observada",
         main = paste("Exceder q95 — horizonte", h))
    abline(0, 1, lty = 2)
    for (k in seq_along(chave)) {
      r <- Filter(function(x) x$config == chave[k] && x$horizonte == h, res)
      if (!length(r)) next
      d <- r[[1]]$teste
      bins <- cut(d$p_q95, c(0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1), include.lowest = TRUE)
      a <- d[, .(p = mean(p_q95), o = mean(y > lim_q95), n = .N), by = bins][n >= 10][order(p)]
      lines(a$p, a$o, type = "b", pch = 19, col = cores[k])
    }
    legend("topleft", chave, col = cores, pch = 19, bty = "n", cex = 0.8)
  }
})

# Distribuição de xi estimado (BGP)
xi <- rbindlist(lapply(res, function(r) {
  if (is.null(r$xi_teste)) return(NULL)
  data.table(horizonte = r$horizonte, config = r$config, xi_media = mean(r$xi_teste),
             xi_p05 = quantile(r$xi_teste, 0.05), xi_p95 = quantile(r$xi_teste, 0.95))
}))
fwrite(xi, tabela_csv("xi"))

# --- relatório em markdown ----------------------------------------------------------
md <- function(dt, digitos = 4) {
  dt <- copy(dt)
  for (j in names(dt)) if (is.numeric(dt[[j]])) set(dt, j = j, value = round(dt[[j]], digitos))
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
teste <- tab[bloco == "teste", .(horizonte, config, CRPS, tw95, tw99, Brier95, QS99,
                                 cobertura80, pit_inf5, pit_sup5, exc95_obs,
                                 exc95_prev)][order(horizonte, CRPS)]
skill <- comp[escore %in% c("crps", "tw_q95")][, .(horizonte, config, escore,
  skill = round(100 * skill, 2), ic = sprintf("[%.1f, %.1f]", 100 * ic_inf, 100 * ic_sup),
  p_dm = signif(p_dm, 2))][order(horizonte, escore, -skill)]
linhas <- c(
  sprintf("# Resultados (%s) — Tmax A801, pós-processamento do GEFSv12", nome_saida), "",
  "Colunas: pit_inf5 e pit_sup5 = fração de observações abaixo do quantil 5% e acima do 95% previstos (ideal 0,05); exc95 = frequência observada e prevista de exceder o q95 sazonal do treino.", "",
  sprintf("Gerado em %s a partir de `%s` (%d configurações).", format(Sys.time()), pasta, length(res)), "",
  "## Escores no teste (2016–2019), menor é melhor", "", md(teste), "",
  "## Skill (%) em relação ao EMOS, IC 95% por bootstrap em blocos de 30 dias, p de Diebold–Mariano", "",
  md(skill, 3), "",
  "## Modelo escolhido pela validação (2012–2015) e seu desempenho no teste", "", md(selecao), "",
  "## xi estimado nas redes BGP (teste)", "", if (nrow(xi)) md(xi, 3) else "(nenhuma)", "",
  sprintf("Figuras: `resultados/figuras/%s_*.png`.", nome_saida)
)
relatorio <- sprintf("resultados/%s_resumo%s.md", nome_saida, sufixo)
writeLines(linhas, relatorio)
cat("Relatório:", relatorio, "\n")
print(selecao)
