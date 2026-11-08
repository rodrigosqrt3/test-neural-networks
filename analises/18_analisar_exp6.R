# Analisa o experimento 6 e o compara com o exp5 (mesmas linhas de teste). Serve também
# para a EUPPBench (pasta resultados/eupp, sem exp5 para comparar).
# Uso: Rscript analises/18_analisar_exp6.R [pasta]     (padrão resultados/exp6)
# Saídas: resultados/<pasta>_resumo.md e resultados/tabelas/<pasta>_*.csv
#
# 1. Combinação "stack": seleção de conjunto (Caruana et al. 2004), gulosa e com
#    reposição, 25 passos: a cada passo entra o componente que mais reduz o CRPS da
#    média dos quantis (Vincentização) na VALIDAÇÃO. A validação vem da etapa 1 (modelos
#    ajustados só no treino), então os pesos não veem o teste. Os quantis guardados têm
#    99 níveis; a grade de 999 é reconstruída por interpolação (cauda extrapolada no
#    espaço normal), e os componentes são reavaliados pela MESMA reconstrução ("_99"),
#    para que a comparação com o stack não dependa dela.
# 2. Inferência como no script 16: unidade de reamostragem = DIA; IC por bootstrap em
#    blocos de 30 dias; Diebold–Mariano na série diária; mesmas linhas conferidas.
source("R/avaliacao.R"); source("R/experimento.R"); source("R/multi.R")   # TAU_FINO

args <- commandArgs(TRUE)
pasta <- if (length(args)) args[[1]] else "resultados/exp6"
SMOKE <- grepl("smoke", pasta)
nome <- basename(pasta)
pasta5 <- if (!startsWith(nome, "exp6")) "" else if (SMOKE) "resultados/exp5_smoke" else "resultados/exp5"
REF <- "emos_loc"
COMPONENTES <- c("emos_bst", "gbm", "qrf", "drn_mais", "att")
ESCADA <- c("drn_exp5", "drn_tv", "drn_mais", "att", "stack")
ORDEM <- c("clim", "emos_loc_exp5", "emos_loc", "emos_bst", "qrf", "gbm", "drn_exp5", "bgp_exp5",
           "drn_tv", "drn_mais", "att", "stack")
DUELO_A <- c("stack", "att", "drn_mais"); DUELO_B <- c("gbm", "qrf", "emos_bst")
if (startsWith(nome, "exp7")) {                 # chuva: referência = climatologia (CRPSS habitual)
  REF <- "clim"
  COMPONENTES <- c("qrf", "drn_c3", "bqn_c3", "bgp_c3", "bgp", "bgp_xi")
  ESCADA <- c("drn_c3", "bqn_c3", "bgp_c3", "bgp", "bgp_xi", "stack")
  ORDEM <- c("clim", "qrf", "drn_c3", "bqn_c3", "bgp_c3", "bgp", "bgp_xi", "stack")
  DUELO_A <- c("stack", "bgp_xi", "bgp", "bgp_c3", "bqn_c3", "drn_c3"); DUELO_B <- c("qrf", "bqn_c3", "drn_c3")
}

ler <- function(p, filtro = ".*") {
  fs <- list.files(p, pattern = "^[^.].*\\.rds$", recursive = TRUE, full.names = TRUE)
  fs <- fs[grepl(filtro, basename(fs)) & !grepl("/_", fs)]     # subpastas "_*" são backups
  lapply(fs, readRDS)
}
res <- ler(pasta)
res5 <- ler(pasta5, "^(drn|bgp|emos_loc)_f")
for (k in seq_along(res5)) res5[[k]]$config <- paste0(res5[[k]]$config, "_exp5")
para_tab <- function(rs) rbindlist(lapply(rs, function(r) cbind(r$teste, modelo = r$config, dobra = r$dobra, horizonte = r$horizonte)))
teste <- rbind(para_tab(res), para_tab(res5), fill = TRUE)
cat(length(res), "resultados em", pasta, "+", length(res5), "do exp5\n")

# --- 1. stack -------------------------------------------------------------------------------
# 99 quantis (TAU_FINO) -> grade TAU_GRADE de 999
grade_de_99 <- function(q) {
  pos <- TAU_GRADE * 100; lo <- pmin(pmax(floor(pos), 1), 98); w <- pos - lo
  out <- sweep(q[, lo, drop = FALSE], 2, 1 - w, `*`) + sweep(q[, lo + 1, drop = FALSE], 2, w, `*`)
  z <- stats::qnorm(TAU_GRADE)
  inf <- TAU_GRADE < 0.01; sup <- TAU_GRADE > 0.99
  dz1 <- stats::qnorm(0.02) - stats::qnorm(0.01); dz2 <- stats::qnorm(0.99) - stats::qnorm(0.98)
  out[, inf] <- q[, 1] + outer((q[, 2] - q[, 1]) / dz1, z[inf] - stats::qnorm(0.01))
  out[, sup] <- q[, 99] + outer((q[, 99] - q[, 98]) / dz2, z[sup] - stats::qnorm(0.99))
  t(apply(out, 1, sort))
}
caruana <- function(qs, y, passos = 25) {
  cont <- setNames(integer(length(qs)), names(qs)); soma <- 0
  for (p in seq_len(passos)) {
    cv <- vapply(qs, function(q) mean(crps_quantis((soma + q) / p, y, TAU_FINO)), numeric(1))
    j <- which.min(cv); cont[j] <- cont[j] + 1L; soma <- soma + qs[[j]]
  }
  cont / passos
}
lim_de <- function(d) list(q90 = d$lim_q90, q95 = d$lim_q95, q99 = d$lim_q99, abs35 = d$lim_abs35)

pesos <- list(); extra <- list()
for (h in sort(unique(teste$horizonte))) for (dob in sort(unique(teste$dobra))) {
  disp <- intersect(COMPONENTES, teste[horizonte == h & dobra == dob, unique(modelo)])
  if (length(disp) < 2) next
  arq <- function(m) file.path(pasta, paste0("h", h), sprintf(".quantis_%s_%s.rds", m, dob))
  q <- setNames(lapply(disp, function(m) readRDS(arq(m))), disp)
  r0 <- Filter(function(r) r$config == disp[1] && r$horizonte == h && r$dobra == dob, res)[[1]]
  w <- caruana(lapply(q, `[[`, "q_va"), r0$validacao$y)
  pesos[[length(pesos) + 1]] <- data.table(horizonte = h, dobra = dob, t(w))
  base <- r0$teste[, .(estacao, data)]
  qte <- Reduce(`+`, Map(function(x, wi) wi * x$q_te, q, w))
  extra[[length(extra) + 1]] <- cbind(base, escores_por_caso(grade_de_99(qte), r0$teste$y, lim_de(r0$teste)),
                                      modelo = "stack", dobra = dob, horizonte = h)
  for (m in disp) extra[[length(extra) + 1]] <- cbind(base, escores_por_caso(grade_de_99(q[[m]]$q_te), r0$teste$y, lim_de(r0$teste)),
                                                      modelo = paste0(m, "_99"), dobra = dob, horizonte = h)
  rm(q); gc()
}
pesos <- rbindlist(pesos, fill = TRUE)
teste <- rbind(teste, rbindlist(extra), fill = TRUE)

# --- 2. escores e comparações ----------------------------------------------------------------
qs <- function(q, y, tau) ((y < q) - tau) * (q - y)
resumo <- teste[, .(n = .N, CRPS = mean(crps), tw90 = mean(tw_q90), tw95 = mean(tw_q95), tw99 = mean(tw_q99),
                    tw35 = mean(tw_abs35), Brier95 = mean((p_q95 - (y > lim_q95))^2),
                    QS99 = mean(qs(q099, y, 0.99)), cob80 = mean(pit > 0.1 & pit < 0.9),
                    exc95_obs = mean(y > lim_q95), exc95_prev = mean(p_q95)),
                by = .(horizonte, modelo)]

comparar <- function(A, B, col, h) {
  a <- teste[modelo == A & horizonte == h]; b <- teste[modelo == B & horizonte == h]
  chave <- c("horizonte", "dobra", "estacao", "data")
  p <- merge(a[, c(chave, col), with = FALSE], b[, c(chave, col), with = FALSE], by = chave, suffixes = c("_a", "_b"))
  if (nrow(p) != nrow(a) || nrow(p) != nrow(b)) stop(sprintf("%s e %s não foram avaliados nas mesmas linhas (h%d)", A, B, h))
  setnames(p, paste0(col, c("_a", "_b")), c("a", "b"))
  dia <- p[, .(a = sum(a), b = sum(b), d = mean(a - b)), by = data][order(data)]
  if (sum(dia$b) == 0)                           # limiar nunca atingido (ex.: 35 °C às 00 UTC na Europa)
    return(data.table(skill = NA_real_, ic_inf = NA_real_, ic_sup = NA_real_, p_dm = NA_real_, dias = nrow(dia)))
  ic <- bootstrap_blocos(dia$a, dia$b, tamanho_bloco = 30, R = 1000)
  data.table(skill = 100 * ic[[1]], ic_inf = 100 * ic[[2]], ic_sup = 100 * ic[[3]],
             p_dm = diebold_mariano(dia$d)[["p_valor"]], dias = nrow(dia))
}
hs <- sort(unique(teste$horizonte))
completos <- function(h) {                       # modelos com as mesmas dobras que a referência
  dref <- teste[modelo == REF & horizonte == h, sort(unique(dobra))]
  teste[horizonte == h, .(ok = identical(sort(unique(dobra)), dref)), by = modelo][ok == TRUE, modelo]
}
principais <- function(h) intersect(ORDEM, completos(h))
skill <- rbindlist(lapply(hs, function(h) rbindlist(lapply(setdiff(principais(h), REF), function(m)
  rbindlist(lapply(c("crps", "tw_q95", "tw_q99", "tw_abs35"), function(col)
    cbind(horizonte = h, modelo = m, escore = col, comparar(m, REF, col, h))))))))
# a escada: cada degrau contra o anterior
escada <- rbindlist(lapply(hs, function(h) {
  e <- intersect(ESCADA, completos(h))
  if (length(e) < 2) return(NULL)
  rbindlist(lapply(2:length(e), function(k) rbindlist(lapply(c("crps", "tw_q95"), function(col)
    cbind(horizonte = h, A = e[k], B = e[k - 1], escore = col, comparar(e[k], e[k - 1], col, h))))))
}))
# o melhor modelo único e o stack contra os baselines fortes
duelos <- rbindlist(lapply(hs, function(h) {
  ok <- completos(h)
  rbindlist(lapply(intersect(DUELO_A, ok), function(A)
    rbindlist(lapply(setdiff(intersect(DUELO_B, ok), A), function(B)
      rbindlist(lapply(c("crps", "tw_q95", "tw_q99", "tw_abs35"), function(col)
        cbind(horizonte = h, A = A, B = B, escore = col, comparar(A, B, col, h))))))))
}))
# reconstrução da grade de 99 quantis: quanto muda o CRPS de cada componente
recon <- merge(resumo[modelo %in% COMPONENTES, .(horizonte, modelo, CRPS)],
               resumo[grepl("_99$", modelo), .(horizonte, modelo = sub("_99$", "", modelo), CRPS_99 = CRPS)],
               by = c("horizonte", "modelo"))
recon[, dif_pct := 100 * (CRPS_99 / CRPS - 1)]

por_dobra <- teste[modelo %in% ORDEM, .(crps = mean(crps)), by = .(horizonte, dobra, modelo)]
por_dobra <- merge(por_dobra, por_dobra[modelo == REF, .(horizonte, dobra, ref = crps)], by = c("horizonte", "dobra"))
por_dobra[, skill := 100 * (1 - crps / ref)]
por_estacao <- teste[modelo %in% ORDEM, .(crps = mean(crps)), by = .(horizonte, estacao, modelo)]
por_estacao <- merge(por_estacao, por_estacao[modelo == REF, .(horizonte, estacao, ref = crps)], by = c("horizonte", "estacao"))
por_estacao[, skill := 100 * (1 - crps / ref)]
# em graus: a mediana como previsão pontual e o intervalo de 80%
graus <- teste[modelo %in% ORDEM, .(erro_abs = mean(abs(q050 - y)), ate1C = 100 * mean(abs(q050 - y) <= 1),
                                     ate2C = 100 * mean(abs(q050 - y) <= 2), larg80 = mean(q090 - q010),
                                     cob80 = 100 * mean(y >= q010 & y <= q090)), by = .(horizonte, modelo)]
# o que a validação escolheu
escolhas <- rbindlist(lapply(res, function(r) {
  i <- r$info
  data.table(horizonte = r$horizonte, dobra = r$dobra, modelo = r$config,
             receita = if (is.null(i$receita)) NA_character_ else i$receita,
             epocas = if (is.null(i$melhor_epoca)) NA_character_ else paste(i$melhor_epoca, collapse = " "),
             xi = if (is.null(i$xi)) NA_real_ else i$xi,
             arvores = if (is.null(i$arvores_media)) NA_character_ else paste0(i$arvores_media, "/", i$arvores_var),
             iter_bst = if (is.null(i$iteracoes)) NA_real_ else stats::median(i$iteracoes),
             minutos = i$segundos / 60)
}))

dir.create("resultados/tabelas", recursive = TRUE, showWarnings = FALSE)
for (nm in c("resumo", "skill", "escada", "duelos", "por_estacao", "graus", "pesos", "escolhas"))
  fwrite(get(nm), sprintf("resultados/tabelas/%s_%s.csv", nome, nm))

# --- 3. relatório -----------------------------------------------------------------------------
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
cabecalho <- if (startsWith(nome, "exp7")) c(
  sprintf("# Experimento 7 — chuva diária, cauda GPD (gerado em %s)", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
  "Chuva diária (00–00 UTC) nas estações do INMET no RS x GEFSv12; testes de 2009–2011, 2012–2015 e 2016–2019 juntos.",
  "Redes censuradas em 0, em duas etapas. Referência dos skills: `clim`. `tw35`/`tw_abs35` = acima de 35 mm.") else
  if (startsWith(nome, "eupp")) c(
  sprintf("# EUPPBench, tarefa R2F — os modelos do exp6 (gerado em %s)", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
  "T2m em 122 estações europeias; treino nos reforecasts do ECMWF (1997–2016), teste nas previsões de 2017–2018.",
  "`horizonte` = passo em horas. Modelos em duas etapas (escolhas na validação 2013–16, final em todos os reforecasts).",
  "Referência dos skills: `emos_loc`. Comparação com a literatura: resultados/eupp_literatura.md.") else c(
  sprintf("# Experimento 6 — modelo de referência da fase 1 (gerado em %s)", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
  "Tmax, estações do INMET no RS x GEFSv12; testes de 2009–2011, 2012–2015 e 2016–2019 juntos.",
  "Todos os modelos do exp6 em duas etapas (escolhas na validação, modelo final em treino + validação).",
  "`*_exp5`: os modelos do exp5 (só treino), nas mesmas linhas. Referência dos skills: `emos_loc` do exp6.")
L <- c(cabecalho,
       "IC por bootstrap em blocos de 30 dias sobre a série DIÁRIA; p de Diebold–Mariano (HAC).", "",
       "## Escores (menor é melhor; cob80 ideal 0,80; exc95: obs ≈ prev)", "",
       md(ordenar(resumo[modelo %in% ORDEM])), "",
       "## Escada: cada degrau contra o anterior (skill de A sobre B)", "", md(ic_txt(escada)), "",
       "## Skill (%) em relação ao EMOS local", "", md(ic_txt(skill[order(horizonte, escore, -skill)])), "",
       "## O melhor contra os baselines fortes (CRPS)", "", md(ic_txt(duelos)), "",
       "## Em graus (mediana como previsão pontual; intervalo de 80%)", "", md(ordenar(graus), 3), "",
       "## Skill CRPS (%) por dobra", "",
       md(dcast(por_dobra[modelo != REF], horizonte + modelo ~ dobra, value.var = "skill")[order(horizonte, match(modelo, ORDEM))], 3), "",
       "## Por estação: skill CRPS (%) vs EMOS local", "",
       md(por_estacao[modelo != REF, .(estacoes = .N, mediana = median(skill), min = min(skill), max = max(skill),
                                        com_ganho = sum(skill > 0)), by = .(horizonte, modelo)][order(horizonte, match(modelo, ORDEM))], 3), "",
       "## Pesos do stack (escolhidos na validação)", "", md(pesos, 2), "",
       "## Reconstrução de 99 quantis: variação do CRPS de cada componente (%)", "", md(recon, 3), "",
       "## O que a validação escolheu (receita, épocas por semente, árvores média/variância, iterações do EMOS-bst)", "",
       md(escolhas[order(horizonte, dobra, match(modelo, ORDEM))], 3))
writeLines(L, sprintf("resultados/%s_resumo.md", nome))
cat("relatório:", sprintf("resultados/%s_resumo.md", nome), "\n")
