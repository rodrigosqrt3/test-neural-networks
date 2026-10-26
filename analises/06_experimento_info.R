# Experimento 4: a rede melhora com INFORMAÇÃO NOVA? E qual é o piso de ruído?
#
# Uso:  Rscript analises/06_experimento_info.R <horizonte>
#       variáveis de ambiente: N_REDES (padrão 10), TORCH_THREADS (padrão 4), SMOKE=1
#       requer dados/processados/a801_estacao_init.rds (analises/05_preparar_estacao.R)
#
# Motivação: no exp. 3 as variações de treino (pesos, grade, online, ...) ficaram
# a ~0.1–0.5% umas das outras, provavelmente dentro do ruído entre sementes. O
# GEFS que usamos é só a Tmax de 5 membros num ponto; o ganho deve vir de mais
# informação, não de outra perda. Modelo-base em tudo: BGP, preditores "obs",
# dropout 0, treino 2000–2011, early stopping em 2012–2015 (o "base" do exp. 3).
#
#   x4_base_s2..s5  base com outras sementes (seed0 200..500): piso de ruído
#   x4_lags         + Tmax de 2 e 3 dias antes e erro de ontem (com dropout 0)
#   x4_traj         + média/dp do ensemble da MESMA rodada nos outros prazos
#   x4_lagens       + rodadas mais antigas para o mesmo dia (ensemble defasado:
#                   h1 <- h3 de init-2 e h7 de init-6; h3 <- h7 de init-4; h7 não tem)
#   x4_estacao      + estado da estação às 00 UTC (T, orvalho, UR, pressão e
#                   tendência 24 h, vento u/v, chuva, radiação e mínima de 24 h)
#   x4_tudo_s1..s5  tudo junto, 5 sementes
#   x4_tudo_largo   tudo, camadas 128-64
#   x4_emos, x4_qrf baselines com a MESMA informação de x4_tudo (regra 4)
#
# Mesma amostra (mesmos dias) dos exp. 1 e 3: faltas nos novos preditores são
# imputadas com a média do treino + indicador. Resultados em resultados/exp4/.

source("R/avaliacao.R"); source("R/modelos_classicos.R"); source("R/modelos_nn.R")
source("R/experimento.R")

args <- commandArgs(TRUE)
H <- as.integer(if (length(args)) args[[1]] else 1)
# Dobra temporal (origem móvel). 0 = split original (teste 2016–2019).
# 1: treino ..2006, validação 2007–2008, teste 2009–2011
# 2: treino ..2009, validação 2010–2011, teste 2012–2015
# (os primeiros anos da A801 têm muitas falhas: com treino só até 2004 sobravam
# ~900 casos, e a dobra mediria "pouco dado", não "informação nova")
# Nas dobras 1 e 2 roda só o conjunto principal de modelos (ver o fim do script);
# com a dobra 0 dá 11 anos de teste (2009–2019) para os mesmos modelos.
DOBRA <- as.integer(if (length(args) >= 2) args[[2]] else 0)
CORTES <- list(`1` = as.Date(c("2006-12-31", "2008-12-31", "2011-12-31")),
               `2` = as.Date(c("2009-12-31", "2011-12-31", "2015-12-31")))
SMOKE <- Sys.getenv("SMOKE") == "1"
N_REDES <- if (SMOKE) 1L else as.integer(Sys.getenv("N_REDES", "10"))
EPOCAS <- if (SMOKE) 15L else 400L
torch_set_num_threads(as.integer(Sys.getenv("TORCH_THREADS", "4")))

pasta <- paste0(if (DOBRA > 0) paste0("exp4_dobra", DOBRA) else "exp4", if (SMOKE) "_smoke" else "")
saida <- file.path("resultados", pasta, paste0("h", H))
dir.create(saida, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) {
  linha <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(linha, "\n"); cat(linha, "\n", file = file.path(saida, "log.txt"), append = TRUE)
}

# --- preditores novos ----------------------------------------------------------------
tab <- readRDS("dados/processados/tabela_tmax_gefs_a801.rds")
est <- readRDS("dados/processados/a801_estacao_init.rds")
F_EST <- setdiff(names(est), "init_date")

# mesma rodada, outros prazos
outros <- setdiff(c(1, 3, 7), H)
F_TRAJ <- character()
for (k in outros) {
  x <- tab[horizon == k, .(init_date, m = ens_media, s = ens_dp)]
  setnames(x, c("m", "s"), paste0(c("traj_media_h", "traj_dp_h"), k))
  F_TRAJ <- c(F_TRAJ, paste0(c("traj_media_h", "traj_dp_h"), k))
  assign(paste0("traj_", k), x)
}
# rodadas anteriores para o mesmo dia válido
defas <- switch(as.character(H), "1" = c(3, 7), "3" = 7, "7" = numeric())
F_LAGENS <- character()
for (k in defas) {
  x <- tab[horizon == k, .(valid_date, m = ens_media, s = ens_dp)]
  setnames(x, c("m", "s"), paste0(c("lagens_media_h", "lagens_dp_h"), k))
  F_LAGENS <- c(F_LAGENS, paste0(c("lagens_media_h", "lagens_dp_h"), k))
  assign(paste0("lagens_", k), x)
}

anexar <- function(d) {
  d <- merge(d, est, by = "init_date", all.x = TRUE, sort = FALSE)
  for (k in outros) d <- merge(d, get(paste0("traj_", k)), by = "init_date", all.x = TRUE, sort = FALSE)
  for (k in defas) d <- merge(d, get(paste0("lagens_", k)), by = "valid_date", all.x = TRUE, sort = FALSE)
  d[order(init_date)]
}

dados <- carregar_horizonte(H)
if (DOBRA > 0) {
  todos <- rbind(dados$tr, dados$va, dados$te)
  k <- CORTES[[as.character(DOBRA)]]; gap <- 7
  dados$tr <- todos[valid_date <= k[1]]
  dados$va <- todos[valid_date > k[1] + gap & valid_date <= k[2]]
  dados$te <- todos[valid_date > k[2] + gap & valid_date <= k[3]]
}
tr <- anexar(dados$tr); va <- anexar(dados$va); te <- anexar(dados$te)
stopifnot(identical(tr$valid_date, dados$tr[order(init_date)]$valid_date),
          identical(va$valid_date, dados$va[order(init_date)]$valid_date),
          identical(te$valid_date, dados$te[order(init_date)]$valid_date))

# imputação com a média do treino + indicadores de falta (por grupo)
NOVOS <- c(F_EST, F_TRAJ, F_LAGENS)
medias <- sapply(tr[, ..NOVOS], mean, na.rm = TRUE)
imputar <- function(d) {
  f_est <- as.numeric(!complete.cases(d[, F_EST, with = FALSE]))
  f_gefs <- as.numeric(!complete.cases(d[, c(F_TRAJ, F_LAGENS), with = FALSE]))
  d[, `:=`(falta_est = f_est, falta_gefs = f_gefs)]
  for (v in NOVOS) set(d, which(is.na(d[[v]])), v, medias[[v]])
  d
}
tr <- imputar(tr); va <- imputar(va); te <- imputar(te)
trva <- rbind(tr, va)
log_msg(sprintf("faltas no treino: estação %.1f%%, gefs defasado/trajetória %.1f%%",
                100 * mean(tr$falta_est), 100 * mean(tr$falta_gefs)))

F_X_EST <- c(F_OBS, F_EST, "falta_est")
F_X_TRAJ <- c(F_OBS, F_TRAJ, "falta_gefs")
F_X_LAGENS <- c(F_OBS, F_LAGENS, "falta_gefs")
F_TUDO <- unique(c(F_LAGS, F_EST, F_TRAJ, F_LAGENS, "falta_est", "falta_gefs"))

LIM <- list(validacao = limiares_para(tr, va$valid_date), teste = limiares_para(tr, te$valid_date))
log_msg(sprintf("horizonte %d, dobra %d | treino %d, validacao %d, teste %d | %d preditores em 'tudo' | N_REDES=%d%s",
                H, DOBRA, nrow(tr), nrow(va), nrow(te), length(F_TUDO), N_REDES, if (SMOKE) " (SMOKE)" else ""))

feito <- function(nome) file.exists(file.path(saida, paste0(nome, ".rds")))
salvar <- function(nome, q_va, q_te, extra = list()) {
  # grades de quantis inteiras (para combinar modelos depois sem retreinar);
  # o ponto no início do nome faz 03_analisar_exp1.R ignorá-las
  saveRDS(list(q_va = q_va, q_te = q_te, data_va = va$valid_date, data_te = te$valid_date),
          file.path(saida, paste0(".grades_", nome, ".rds")))
  res <- c(list(
    config = nome, horizonte = H,
    validacao = cbind(data.table(data = va$valid_date), escores_por_caso(q_va, va$tmax, LIM$validacao)),
    teste = cbind(data.table(data = te$valid_date), escores_por_caso(q_te, te$tmax, LIM$teste))
  ), extra)
  saveRDS(res, file.path(saida, paste0(nome, ".rds")))
  d <- res$teste
  log_msg(sprintf("  %-16s CRPS val %.4f teste %.4f | tw95 %.4f | exc95 obs %.3f prev %.3f",
                  nome, mean(res$validacao$crps), mean(d$crps), mean(d$tw_q95),
                  mean(d$y > d$lim_q95), mean(d$p_q95)))
}

treinar <- function(nome, feats, ...) {
  if (feito(nome)) return(invisible())
  padr <- padronizador(tr, feats)
  t0 <- Sys.time()
  conj <- treinar_conjunto(Xtr = padr(tr), ytr = tr$tmax, Xva = padr(va), yva = va$tmax,
                           n_redes = N_REDES, epocas = EPOCAS, cabeca = "bgp", dropout = 0, ...)
  salvar(nome, prever_conjunto(conj, padr(va)), prever_conjunto(conj, padr(te)),
         list(segundos = as.numeric(Sys.time() - t0, units = "secs"),
              xi_media = mean(xi_conjunto(conj, padr(te))),
              melhor_epoca = sapply(conj, `[[`, "melhor_epoca"), preditores = feats))
}

# --- baselines com a mesma informação (rápidos: rodam primeiro) ------------------------
if (!feito("x4_emos")) {
  t0 <- Sys.time()
  xe <- c("ens_media", "obs_vespera", "doy_sen", "doy_cos", F_EST, F_TRAJ, F_LAGENS)
  ze <- function(x) cbind(log(x$ens_dp + 0.05), x$doy_sen, x$doy_cos)
  pe <- function(a, b) {
    p <- padronizador(a, xe)
    predict(emos_ajustar(p(a), ze(a), a$tmax), p(b), ze(b))
  }
  salvar("x4_emos", pe(tr, va), pe(trva, te), list(segundos = as.numeric(Sys.time() - t0, units = "secs")))
}
if (!feito("x4_qrf")) {
  t0 <- Sys.time()
  pq <- function(a, b) qrf_quantis(qrf_ajustar(as.matrix(a[, ..F_TUDO]), a$tmax), as.matrix(b[, ..F_TUDO]))
  salvar("x4_qrf", pq(tr, va), pq(trva, te), list(segundos = as.numeric(Sys.time() - t0, units = "secs")))
}

# --- dobras 1 e 2: só o conjunto principal --------------------------------------------
if (DOBRA > 0) {
  # EMOS com a mesma definição do exp. 1 (referência dos skills)
  if (!feito("emos")) {
    xe0 <- c("ens_media", "obs_vespera", "doy_sen", "doy_cos")
    ze0 <- function(x) cbind(log(x$ens_dp + 0.05), x$doy_sen, x$doy_cos)
    pe0 <- function(a, b) predict(emos_ajustar(as.matrix(a[, ..xe0]), ze0(a), a$tmax),
                                  as.matrix(b[, ..xe0]), ze0(b))
    salvar("emos", pe0(tr, va), pe0(trva, te))
  }
  treinar("x4_base_s1", F_OBS, seed0 = 100)
  treinar("x4_tudo_s1", F_TUDO, seed0 = 100)
  treinar("x4_estacao", F_X_EST)
  treinar("x4_base_s2", F_OBS, seed0 = 200)
  treinar("x4_tudo_s2", F_TUDO, seed0 = 200)
  # acrescentado depois: o ensemble defasado foi o que ajudou no h1 do split original
  if (length(F_LAGENS)) treinar("x4_lagens", F_X_LAGENS)
  log_msg("horizonte ", H, ", dobra ", DOBRA, " concluído")
  quit(save = "no")
}

# --- redes: primeiro o que responde à pergunta principal ---------------------------------
treinar("x4_tudo_s1", F_TUDO, seed0 = 100)
treinar("x4_estacao", F_X_EST)
treinar("x4_base_s2", F_OBS, seed0 = 200)
treinar("x4_tudo_s2", F_TUDO, seed0 = 200)
treinar("x4_traj", F_X_TRAJ)
if (length(F_LAGENS)) treinar("x4_lagens", F_X_LAGENS)
treinar("x4_lags", F_LAGS)
treinar("x4_base_s3", F_OBS, seed0 = 300)
treinar("x4_tudo_s3", F_TUDO, seed0 = 300)
treinar("x4_tudo_largo", F_TUDO, ocultas = c(128, 64))
treinar("x4_base_s4", F_OBS, seed0 = 400)
treinar("x4_tudo_s4", F_TUDO, seed0 = 400)
treinar("x4_base_s5", F_OBS, seed0 = 500)
treinar("x4_tudo_s5", F_TUDO, seed0 = 500)

log_msg("horizonte ", H, " concluído")
