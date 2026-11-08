# Analisa o experimento 5 (multi-estação).
# Uso: Rscript analises/16_analisar_exp5.R [pasta]     (padrão resultados/exp5)
# Saídas: resultados/exp5_resumo.md e resultados/tabelas/exp5_*.csv (sufixo _smoke no teste)
#
# Inferência: estações no mesmo dia têm erros correlacionados, então a unidade de
# reamostragem é o DIA, não o caso. Os escores são somados por dia (todas as
# estações e dobras) e o skill 1 - S_modelo / S_ref tem IC por bootstrap em blocos
# de 30 dias; o p de Diebold–Mariano usa a série diária das diferenças médias
# (variância HAC). Antes de comparar, verifica-se que os dois modelos foram
# avaliados exatamente nas mesmas linhas.
source("R/avaliacao.R")
library(data.table)

args <- commandArgs(TRUE)
pasta <- if (length(args)) args[[1]] else "resultados/exp5"
suf <- if (grepl("smoke", pasta)) "_smoke" else ""
REF <- "emos_loc"
ORDEM <- c("clim", "emos_loc", "emos_loc_obs", "qrf", "drn", "bgp")

fs <- list.files(pasta, pattern = "^[^.].*\\.rds$", recursive = TRUE, full.names = TRUE)
res <- lapply(fs, readRDS)
teste <- rbindlist(lapply(res, function(r) cbind(r$teste, modelo = r$config, dobra = r$dobra, horizonte = r$horizonte)))
cat(length(res), "resultados em", pasta, "\n")

qs <- function(q, y, tau) ((y < q) - tau) * (q - y)
resumo <- teste[, .(n = .N, CRPS = mean(crps), tw90 = mean(tw_q90), tw95 = mean(tw_q95), tw99 = mean(tw_q99),
                    Brier95 = mean((p_q95 - (y > lim_q95))^2), Brier99 = mean((p_q99 - (y > lim_q99))^2),
                    QS99 = mean(qs(q099, y, 0.99)), cob80 = mean(pit > 0.1 & pit < 0.9),
                    pit_inf5 = mean(pit < 0.05), pit_sup5 = mean(pit > 0.95),
                    exc95_obs = mean(y > lim_q95), exc95_prev = mean(p_q95)),
                by = .(horizonte, modelo)]

# skill de A em relação a B num escore, com reamostragem por dia
comparar <- function(A, B, col, filtro = NULL) {
  a <- teste[modelo == A]; b <- teste[modelo == B]
  if (!is.null(filtro)) { a <- a[eval(filtro)]; b <- b[eval(filtro)] }
  chave <- c("horizonte", "dobra", "estacao", "data")
  p <- merge(a[, c(chave, col), with = FALSE], b[, c(chave, col), with = FALSE], by = chave, suffixes = c("_a", "_b"))
  if (nrow(p) != nrow(a) || nrow(p) != nrow(b)) stop(sprintf("%s e %s não foram avaliados nas mesmas linhas", A, B))
  setnames(p, paste0(col, c("_a", "_b")), c("a", "b"))
  dia <- p[, .(a = sum(a), b = sum(b), d = mean(a - b)), by = data][order(data)]
  ic <- bootstrap_blocos(dia$a, dia$b, tamanho_bloco = 30, R = 1000)
  data.table(skill = 100 * ic[[1]], ic_inf = 100 * ic[[2]], ic_sup = 100 * ic[[3]],
             p_dm = diebold_mariano(dia$d)[["p_valor"]], dias = nrow(dia))
}

modelos <- intersect(ORDEM, unique(teste$modelo))
hs <- sort(unique(teste$horizonte))
skill <- rbindlist(lapply(hs, function(h) rbindlist(lapply(setdiff(modelos, REF), function(m)
  rbindlist(lapply(c("crps", "tw_q95", "tw_q99"), function(col)
    cbind(horizonte = h, modelo = m, escore = col, comparar(m, REF, col, bquote(horizonte == .(h))))))))))
# a nossa rede contra as alternativas mais fortes
duelos <- rbindlist(lapply(hs, function(h) rbindlist(lapply(c("drn", "qrf", "emos_loc_obs"), function(b) {
  if (!all(c("bgp", b) %in% modelos)) return(NULL)
  rbindlist(lapply(c("crps", "tw_q95"), function(col)
    cbind(horizonte = h, A = "bgp", B = b, escore = col, comparar("bgp", b, col, bquote(horizonte == .(h))))))
}))))
por_dobra <- teste[, .(crps = mean(crps)), by = .(horizonte, dobra, modelo)]
por_dobra <- merge(por_dobra, por_dobra[modelo == REF, .(horizonte, dobra, ref = crps)], by = c("horizonte", "dobra"))
por_dobra[, skill := 100 * (1 - crps / ref)]
por_estacao <- teste[, .(n = .N, crps = mean(crps)), by = .(horizonte, estacao, modelo)]
por_estacao <- merge(por_estacao, por_estacao[modelo == REF, .(horizonte, estacao, ref = crps)], by = c("horizonte", "estacao"))
por_estacao[, skill := 100 * (1 - crps / ref)]

dir.create("resultados/tabelas", recursive = TRUE, showWarnings = FALSE)
fwrite(resumo, sprintf("resultados/tabelas/exp5_escores%s.csv", suf))
fwrite(skill, sprintf("resultados/tabelas/exp5_skill%s.csv", suf))
fwrite(duelos, sprintf("resultados/tabelas/exp5_duelos%s.csv", suf))
fwrite(por_estacao, sprintf("resultados/tabelas/exp5_por_estacao%s.csv", suf))

# --- relatório -------------------------------------------------------------------------------
md <- function(dt, dig = 4) {
  dt <- copy(dt)
  for (v in names(dt)) if (is.double(dt[[v]])) set(dt, j = v, value = formatC(dt[[v]], digits = dig, format = "fg"))
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
ic_txt <- function(d) d[, .(horizonte, modelo = if ("modelo" %in% names(d)) modelo else paste(A, "vs", B), escore,
                            skill = sprintf("%.2f", skill), IC95 = sprintf("[%.1f, %.1f]", ic_inf, ic_sup),
                            p_DM = sprintf("%.3f", p_dm))]
ordenar <- function(d) d[order(horizonte, match(modelo, ORDEM))]
L <- c(sprintf("# Experimento 5 — Tmax, 43 estações do RS x GEFSv12 (gerado em %s)", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
       "Testes de 2009–2011, 2012–2015 e 2016–2019 juntos (cada dobra só com os modelos ajustados no seu treino).",
       sprintf("Referência dos skills: `%s`. IC por bootstrap em blocos de 30 dias sobre a série DIÁRIA; p de Diebold–Mariano (HAC).", REF), "",
       "## Escores (menor é melhor; cob80, pit_*: ideal 0,80 / 0,05; exc95: obs ≈ prev)", "",
       md(ordenar(resumo)), "",
       "## Skill (%) em relação ao EMOS local", "", md(ic_txt(skill[order(horizonte, escore, -skill)])), "",
       "## BGP (nossa rede) contra as alternativas mais fortes", "", md(ic_txt(duelos)), "",
       "## Skill CRPS (%) por dobra", "",
       md(dcast(por_dobra[modelo != REF], horizonte + modelo ~ dobra, value.var = "skill")[order(horizonte, match(modelo, ORDEM))], 3), "",
       "## Por estação: skill CRPS (%) vs EMOS local (mediana, mín., máx., estações com ganho)", "",
       md(por_estacao[modelo != REF, .(estacoes = .N, mediana = median(skill), min = min(skill), max = max(skill),
                                        com_ganho = sum(skill > 0)), by = .(horizonte, modelo)][order(horizonte, match(modelo, ORDEM))], 3))
writeLines(L, sprintf("resultados/exp5_resumo%s.md", suf))
cat("relatório:", sprintf("resultados/exp5_resumo%s.md", suf), "\n")
