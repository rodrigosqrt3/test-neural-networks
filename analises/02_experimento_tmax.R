# Experimento 1: pós-processamento da Tmax do GEFSv12 para a A801 (Porto Alegre).
#
# Uso:  Rscript analises/02_experimento_tmax.R <horizonte>     (1, 3 ou 7)
#       variáveis de ambiente: N_REDES (padrão 10), TORCH_THREADS (padrão 4),
#                              SMOKE=1 (teste rápido: 1 rede, poucas épocas)
#
# Hipóteses testadas (todas avaliadas nos MESMOS dias de validação e teste):
#   H1  redes distribucionais (DRN, BQN, BGP) vs baselines (climatologia,
#       ensemble bruto, EMOS, QRF)
#   H2  cauda GPD explícita (BGP) vs quantis Bernstein (BQN)
#   H3  perda ponderada na cauda (qwCRPS, gama) melhora a cauda? a que custo
#       no CRPS? (fronteira de Pareto)
#   H4  xi global vs xi(x) com encolhimento lambda
#   H5  onde começa a cauda (tau0)
#   H6  quais preditores ajudam (só GEFS; + Tmax da véspera; + defasagens e
#       erro de ontem)
#
# Cada configuração é salva em resultados/exp1/h<h>/<config>.rds ao terminar;
# rodar de novo pula o que já existe (retoma após interrupção). Salvamos
# escores POR CASO (não as grades de quantis inteiras), para permitir
# bootstrap e testes de Diebold-Mariano depois, em 03_analisar_exp1.R.

source("R/avaliacao.R"); source("R/modelos_classicos.R"); source("R/modelos_nn.R")
library(data.table)

args <- commandArgs(TRUE)
H <- as.integer(if (length(args)) args[[1]] else 1)
BATERIA <- if (length(args) >= 2) args[[2]] else "exp1"   # "exp1" ou "exp2"
SMOKE <- Sys.getenv("SMOKE") == "1"
N_REDES <- if (SMOKE) 1L else as.integer(Sys.getenv("N_REDES", "10"))
EPOCAS <- if (SMOKE) 15L else 400L
torch_set_num_threads(as.integer(Sys.getenv("TORCH_THREADS", "4")))

saida <- file.path("resultados", if (SMOKE) "exp1_smoke" else "exp1", paste0("h", H))
dir.create(saida, recursive = TRUE, showWarnings = FALSE)
arquivo_log <- file.path(saida, "log.txt")
log_msg <- function(...) {
  linha <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(linha, "\n")
  cat(linha, "\n", file = arquivo_log, append = TRUE)
}

# --- dados ---------------------------------------------------------------
F_GEFS <- c("ens_media", "ens_dp", "ens_min", "ens_max", "doy_sen", "doy_cos")
F_OBS  <- c(F_GEFS, "obs_vespera")
F_LAGS <- c(F_OBS, "obs_lag2", "obs_lag3", "erro_ontem")
F_MEM  <- c(paste0("membro_", 1:5), "doy_sen", "doy_cos", "obs_vespera")  # membros ordenados
F_POOL <- c(F_OBS, "h_1", "h_3", "h_7")                                   # todos os prazos juntos
CONJUNTOS <- list(gefs = F_GEFS, obs = F_OBS, lags = F_LAGS, mem = F_MEM, pool = F_POOL)

tab <- readRDS("dados/processados/tabela_tmax_gefs_a801.rds")
membros <- attr(tab, "membros")
ordenados <- t(apply(as.matrix(tab[, ..membros]), 1, sort, na.last = TRUE))
for (k in 1:5) tab[[paste0("membro_", k)]] <- ordenados[, k]
for (h in c(1, 3, 7)) tab[[paste0("h_", h)]] <- as.numeric(tab$horizon == h)

# amostra comum a todas as configurações de um horizonte
amostra <- function(h) {
  x <- tab[horizon == h & !is.na(bloco) & !is.na(tmax)]
  x[complete.cases(x[, c(F_LAGS, membros), with = FALSE])]
}
d <- amostra(H)
tr <- d[bloco == "treino"]; va <- d[bloco == "validacao"]; te <- d[bloco == "teste"]
trva <- rbind(tr, va)
log_msg(sprintf("horizonte %d | treino %d, validacao %d, teste %d | N_REDES=%d%s",
                H, nrow(tr), nrow(va), nrow(te), N_REDES, if (SMOKE) " (SMOKE)" else ""))

# limiares sazonais do TREINO (+ um absoluto de impacto)
limiares_para <- function(datas) {
  out <- lapply(c(q90 = 0.90, q95 = 0.95, q99 = 0.99),
                function(p) limiar_sazonal(tr$valid_date, tr$tmax, p)(datas))
  out$abs35 <- rep(35, length(datas))
  out
}
LIM <- list(validacao = limiares_para(va$valid_date), teste = limiares_para(te$valid_date))
NIVEIS_Q <- c(0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99)

# Escores por caso a partir da grade de quantis
escores_por_caso <- function(q, y, lim) {
  out <- data.table(y = y, crps = crps_quantis(q, y), pit = pit_quantis(q, y))
  for (nm in names(lim)) {
    out[[paste0("tw_", nm)]] <- twcrps_quantis(q, y, lim[[nm]])
    out[[paste0("p_", nm)]] <- prob_excedencia(q, lim[[nm]])
    out[[paste0("lim_", nm)]] <- lim[[nm]]
  }
  idx <- vapply(NIVEIS_Q, function(t) which.min(abs(TAU_GRADE - t)), integer(1))
  for (k in seq_along(NIVEIS_Q)) out[[sprintf("q%03d", round(NIVEIS_Q[k] * 100))]] <- q[, idx[k]]
  out
}

# padronização com média e desvio do TREINO
padronizador <- function(feats) {
  m <- sapply(tr[, ..feats], mean); s <- sapply(tr[, ..feats], sd)
  function(x) sweep(sweep(as.matrix(x[, ..feats]), 2, m), 2, s, "/")
}

salvar <- function(nome, prev_va, prev_te, extra = list()) {
  res <- c(list(
    config = nome, horizonte = H,
    validacao = cbind(data.table(data = va$valid_date), escores_por_caso(prev_va, va$tmax, LIM$validacao)),
    teste = cbind(data.table(data = te$valid_date), escores_por_caso(prev_te, te$tmax, LIM$teste))
  ), extra)
  saveRDS(res, file.path(saida, paste0(nome, ".rds")))
  log_msg(sprintf("  %-22s CRPS val %.4f | teste %.4f | twCRPS95 teste %.4f",
                  nome, mean(res$validacao$crps), mean(res$teste$crps), mean(res$teste$tw_q95)))
}
feito <- function(nome) file.exists(file.path(saida, paste0(nome, ".rds")))

# --- baselines: ajuste em treino (-> validação) e em treino+val (-> teste) ---
roda_baseline <- function(nome, ajusta_prevê) {
  if (feito(nome)) return(invisible())
  t0 <- Sys.time()
  prev_va <- ajusta_prevê(tr, va)
  prev_te <- ajusta_prevê(trva, te)
  salvar(nome, prev_va, prev_te, list(segundos = as.numeric(Sys.time() - t0, units = "secs")))
}

roda_baseline("climatologia", function(a, b) predict(clim_ajustar(a$valid_date, a$tmax), b$valid_date))
roda_baseline("ensemble_bruto", function(a, b) ens_quantis(as.matrix(b[, ..membros])))
xe <- c("ens_media", "obs_vespera", "doy_sen", "doy_cos")
ze <- function(x) cbind(log(x$ens_dp + 0.05), x$doy_sen, x$doy_cos)
roda_baseline("emos", function(a, b) {
  fit <- emos_ajustar(as.matrix(a[, ..xe]), ze(a), a$tmax)
  predict(fit, as.matrix(b[, ..xe]), ze(b))
})
for (cj in c("obs", "lags")) {
  roda_baseline(paste0("qrf_", cj), function(a, b) {
    f <- CONJUNTOS[[cj]]
    qrf_quantis(qrf_ajustar(as.matrix(a[, ..f]), a$tmax), as.matrix(b[, ..f]))
  })
}

# --- redes ----------------------------------------------------------------
configs <- list(
  # H1/H2: redes básicas (preditores "obs")
  drn_obs          = list(cabeca = "drn", feats = "obs"),
  bqn_obs          = list(cabeca = "bqn", feats = "obs"),
  bgp_obs          = list(cabeca = "bgp", feats = "obs"),
  # H6: preditores
  drn_gefs         = list(cabeca = "drn", feats = "gefs"),
  bgp_gefs         = list(cabeca = "bgp", feats = "gefs"),
  drn_lags         = list(cabeca = "drn", feats = "lags"),
  bgp_lags         = list(cabeca = "bgp", feats = "lags"),
  # H3: peso na cauda (gama) — fronteira CRPS x twCRPS
  bqn_obs_g1       = list(cabeca = "bqn", feats = "obs", gama = 1),
  bqn_obs_g4       = list(cabeca = "bqn", feats = "obs", gama = 4),
  bqn_obs_g9       = list(cabeca = "bqn", feats = "obs", gama = 9),
  bqn_obs_g19      = list(cabeca = "bqn", feats = "obs", gama = 19),
  bgp_obs_g1       = list(cabeca = "bgp", feats = "obs", gama = 1),
  bgp_obs_g4       = list(cabeca = "bgp", feats = "obs", gama = 4),
  bgp_obs_g9       = list(cabeca = "bgp", feats = "obs", gama = 9),
  bgp_obs_g19      = list(cabeca = "bgp", feats = "obs", gama = 19),
  # H4: xi dependente de x, com encolhimento
  bgp_obs_xix_l01  = list(cabeca = "bgp", feats = "obs", xi_modo = "x", lambda_xi = 0.1),
  bgp_obs_xix_l1   = list(cabeca = "bgp", feats = "obs", xi_modo = "x", lambda_xi = 1),
  bgp_obs_xix_l10  = list(cabeca = "bgp", feats = "obs", xi_modo = "x", lambda_xi = 10),
  # H5: início da cauda
  bgp_obs_tau85    = list(cabeca = "bgp", feats = "obs", tau0 = 0.85),
  bgp_obs_tau95    = list(cabeca = "bgp", feats = "obs", tau0 = 0.95),
  # combinações
  bgp_lags_g4      = list(cabeca = "bgp", feats = "lags", gama = 4),
  bgp_lags_xix_l1  = list(cabeca = "bgp", feats = "lags", xi_modo = "x", lambda_xi = 1)
)

# Segunda bateria (Rscript ... <h> exp2), resultados nas mesmas pastas:
#   H7  perda híbrida: corpo por escore quantílico + cauda por verossimilhança
#       GPD, limiar e cauda estimados juntos num passo
#   H8  membros do ensemble (ordenados) como entrada
#   H9  um modelo para os três prazos (3x mais dados)
#   H10 sensibilidade à arquitetura
configs2 <- list(
  bgp_obs_nll        = list(cabeca = "bgp", feats = "obs", perda = "hibrida"),
  bgp_obs_nll_l025   = list(cabeca = "bgp", feats = "obs", perda = "hibrida", lambda_nll = 0.25),
  bgp_obs_nll_xix    = list(cabeca = "bgp", feats = "obs", perda = "hibrida",
                            xi_modo = "x", lambda_xi = 1),
  drn_mem            = list(cabeca = "drn", feats = "mem"),
  bgp_mem            = list(cabeca = "bgp", feats = "mem"),
  drn_pool           = list(cabeca = "drn", feats = "pool"),
  bgp_pool           = list(cabeca = "bgp", feats = "pool"),
  bgp_pool_g4        = list(cabeca = "bgp", feats = "pool", gama = 4),
  bgp_obs_largo      = list(cabeca = "bgp", feats = "obs", ocultas = c(128, 64)),
  bgp_obs_fundo      = list(cabeca = "bgp", feats = "obs", ocultas = c(64, 64, 32)),
  bgp_obs_drop0      = list(cabeca = "bgp", feats = "obs", dropout = 0),
  bgp_obs_drop3      = list(cabeca = "bgp", feats = "obs", dropout = 0.3),
  drn_obs_largo      = list(cabeca = "drn", feats = "obs", ocultas = c(128, 64))
)
if (BATERIA == "exp2") configs <- configs2

for (nome in names(configs)) {
  if (feito(nome)) { log_msg("  (já existe) ", nome); next }
  cfg <- configs[[nome]]
  # "pool": treina com treino/validação de TODOS os prazos; avalia neste prazo
  if (cfg$feats == "pool") {
    todos <- rbindlist(lapply(c(1, 3, 7), amostra))
    tr_fit <- todos[bloco == "treino"]; va_fit <- todos[bloco == "validacao"]
  } else {
    tr_fit <- tr; va_fit <- va
  }
  f <- CONJUNTOS[[cfg$feats]]
  m <- sapply(tr_fit[, ..f], mean); s <- sapply(tr_fit[, ..f], sd)
  s[s == 0] <- 1
  padr <- function(x) sweep(sweep(as.matrix(x[, ..f]), 2, m), 2, s, "/")
  cfg$feats <- NULL
  t0 <- Sys.time()
  conj <- do.call(treinar_conjunto, c(
    list(Xtr = padr(tr_fit), ytr = tr_fit$tmax, Xva = padr(va_fit), yva = va_fit$tmax,
         n_redes = N_REDES, epocas = EPOCAS), cfg))
  extra <- list(
    segundos = as.numeric(Sys.time() - t0, units = "secs"),
    epocas = sapply(conj, `[[`, "epocas"),
    melhor_val = sapply(conj, `[[`, "melhor_val")
  )
  if (cfg$cabeca == "bgp") extra$xi_teste <- rowMeans(as.matrix(xi_conjunto(conj, padr(te))))
  salvar(nome, prever_conjunto(conj, padr(va)), prever_conjunto(conj, padr(te)), extra)
}
log_msg("horizonte ", H, " concluído")
