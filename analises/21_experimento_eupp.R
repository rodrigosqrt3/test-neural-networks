# Experimento EUPPBench (critério 2 da fase 1): os modelos do exp6 na tarefa R2F de
# Feik et al. (2024): T2m em 122 estações europeias; treino nos reforecasts do ECMWF
# (11 membros, 1997–2016), teste nas previsões operacionais de 2017–2018 (51 membros).
#
# Uso:  Rscript analises/21_experimento_eupp.R <passo>     (24, 72 ou 120 horas)
#   variáveis de ambiente: N_REDES (5), TORCH_THREADS, QRF_THREADS, GBM_THREADS, QRF_ARVORES (200), SMOKE=1
#   SMOKE=1: 6 estações, 1 rede, 3 épocas, poucas árvores (só testa o código)
# Requer dados/processados/euppbench/eupp_{ref,fc}_h<passo>.csv (analises/20_converter_euppbench.py).
#
# Divisão: treino = reforecasts com data válida até 2012; validação = 2013–2016 (7 dias de
# intervalo), usada para escolher épocas/iterações/receitas e os pesos do stack; o modelo
# final é reajustado em todos os reforecasts (duas etapas, como no exp6) e prevê 2017–18.
# Uma "dobra" só (eupp). A atenção ficou de fora: no exp6 a validação escolheu 0 camadas
# em 7 das 9 dobras.
#
# Saída: resultados/eupp/h<passo>/<modelo>_eupp.rds e .quantis_*, no formato do exp6
# (analisados por analises/18_analisar_exp6.R resultados/eupp e por 22_literatura_eupp.R).

source("R/avaliacao.R"); source("R/modelos_classicos.R"); source("R/modelos_nn.R")
source("R/experimento.R"); source("R/multi.R")

args <- commandArgs(TRUE)
PASSO <- as.integer(if (length(args)) args[[1]] else 24)
SMOKE <- Sys.getenv("SMOKE") == "1"
N_REDES <- if (SMOKE) 1L else as.integer(Sys.getenv("N_REDES", "5"))
EPOCAS <- if (SMOKE) 3L else 300L
QRF_ARVORES <- if (SMOKE) 50L else as.integer(Sys.getenv("QRF_ARVORES", "200"))
QRF_THREADS <- as.integer(Sys.getenv("QRF_THREADS", "2"))
GBM_THREADS <- as.integer(Sys.getenv("GBM_THREADS", "2"))
torch_set_num_threads(as.integer(Sys.getenv("TORCH_THREADS", "2")))
# as receitas do exp6; com ~400 mil casos de treino, 1 semente por candidata basta para escolher
REDES <- list(
  base  = list(ocultas = c(128, 64), dropout = 0, dim_emb = 4, lote = 256, lr = 1e-3, paciencia = 20, n_busca = 1),
  drop  = list(ocultas = c(128, 64), dropout = 0.1, dim_emb = 4, lote = 256, lr = 1e-3, paciencia = 20, wd = 1e-4, n_busca = 1),
  larga = list(ocultas = c(256, 128, 64), dropout = 0.1, dim_emb = 8, lote = 256, lr = 5e-4, paciencia = 20, wd = 1e-4,
               n_busca = 1))

saida <- file.path("resultados", if (SMOKE) "eupp_smoke" else "eupp", paste0("h", PASSO))
dir.create(saida, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) {
  linha <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(linha, "\n"); cat(linha, "\n", file = file.path(saida, "log.txt"), append = TRUE)
}

# --- dados no formato da pipeline (tmax = alvo, ens_media/ens_dp = T2m do ensemble) ---------
ler <- function(tipo) {
  d <- fread(sprintf("dados/processados/euppbench/eupp_%s_h%03d.csv", tipo, PASSO))
  d[, `:=`(estacao = as.character(station_id), valid_date = as.Date(substr(as.character(valid), 1, 10)),
           tmax = obs, ens_media = t2m_media, ens_dp = t2m_dp, fonte = tipo)]
}
d <- rbind(ler("ref"), ler("fc"))
d[!is.finite(ens_dp), ens_dp := stats::median(d$ens_dp, na.rm = TRUE)]   # usado pelo EMOS, fora da imputação
todas <- d[, .(lat = lat[1], lon = lon[1]), by = estacao][order(estacao)]
if (SMOKE) { todas <- todas[unique(round(seq(1, .N, length.out = min(6, .N))))]; d <- d[estacao %in% todas$estacao] }
d[, id := match(estacao, todas$estacao)]
doy <- as.integer(format(d$valid_date, "%j"))
d[, `:=`(doy_sen = sin(2 * pi * doy / 365.25), doy_cos = cos(2 * pi * doy / 365.25),
         alt_dif = altitude - model_altitude)]
for (v in c("tp6_media", "tp6_dp", "cp6_media", "cp6_dp", "cape_media", "cape_dp", "sd_media", "sd_dp"))
  set(d, j = v, value = log1p(pmax(d[[v]], 0)))
setorder(d, valid_date, estacao)
FEATS <- c(setdiff(grep("_(media|dp)$", names(d), value = TRUE), c("ens_media", "ens_dp")), "t2m_min", "t2m_max",
           "lat", "lon", "altitude", "model_altitude", "model_orography", "alt_dif", "doy_sen", "doy_cos")
dist <- distancias_km(todas)

sp <- list(tr = d[fonte == "ref" & valid_date <= as.Date("2012-12-31")],
           va = d[fonte == "ref" & valid_date > as.Date("2012-12-31") + 7],
           te = d[fonte == "fc"])
# como no exp5/6 (MIN_TREINO): a estação entra se tiver >= 365 casos de treino. Ficam de fora
# 6438 e 7072 (85 e 136 casos até 2012), cuja climatologia nem cobre o ano todo.
ok <- sp$tr[, .N, by = estacao][N >= MIN_TREINO, estacao]
for (b in c("tr", "va", "te")) sp[[b]] <- sp[[b]][estacao %in% ok]
sp <- imputar(adicionar_clim(sp), c(FEATS, F_CLIM), list())
sp$estacoes <- sort(ok)
sp$f6 <- sp$feats
lim <- list(va = limiares_estacao(sp$tr, sp$va), te = limiares_estacao(sp$tr, sp$te))
log_msg(sprintf("passo %d h | %d estações | treino %d, validação %d, teste %d | %d preditores | N_REDES=%d%s",
                PASSO, length(sp$estacoes), nrow(sp$tr), nrow(sp$va), nrow(sp$te), length(sp$f6), N_REDES,
                if (SMOKE) " (SMOKE)" else ""))
log_msg(sprintf("CRPS do ensemble bruto no teste: %.4f", mean(sp$te$crps_ens)))

MODELOS <- list(
  clim     = function(sp) duas_etapas(modelo_clim, sp),
  emos_loc = function(sp) duas_etapas(function(s) modelo_emos_loc(s, c("ens_media", "doy_sen", "doy_cos")), sp),
  emos_bst = function(sp) modelo_emos_bst(sp, sp$f6, m_max = if (SMOKE) 50 else 4000),   # no exp6, 2000 ainda encostou no teto
  gbm      = function(sp) modelo_gbm(sp, sp$f6, GBM_THREADS, n_max = if (SMOKE) 100 else 3000),
  qrf      = function(sp) duas_etapas(function(s) modelo_qrf(s, s$f6, QRF_ARVORES, QRF_THREADS), sp),
  drn_mais = function(sp) modelo_rede2(sp, sp$f6, "drn", nrow(todas), N_REDES, EPOCAS, REDES)
)

for (m in names(MODELOS)) {
  arq <- file.path(saida, sprintf("%s_eupp.rds", m))
  if (file.exists(arq)) next
  t0 <- Sys.time()
  mod <- MODELOS[[m]](sp)
  va <- avaliar_em_blocos(mod$va, sp$va, lim$va)
  te <- avaliar_em_blocos(mod$te, sp$te, lim$te)
  seg <- as.numeric(Sys.time() - t0, units = "secs")
  saveRDS(list(q_va = va$quantis, q_te = te$quantis), file.path(saida, sprintf(".quantis_%s_eupp.rds", m)))
  saveRDS(list(config = m, dobra = "eupp", horizonte = PASSO, validacao = va$escores, teste = te$escores,
               info = c(mod$info, list(segundos = seg, n_treino = nrow(sp$tr), n_validacao = nrow(sp$va),
                                       estacoes = sp$estacoes))),
          arq)
  log_msg(sprintf("  %-9s | CRPS val %.4f teste %.4f | tw95 teste %.4f | %.0f s", m,
                  mean(va$escores$crps), mean(te$escores$crps), mean(te$escores$tw_q95), seg))
}
# o ensemble bruto, para conferir com a literatura (mesmas linhas de teste)
saveRDS(sp$te[, .(estacao, data = valid_date, crps_ens)], file.path(saida, ".ensemble_bruto.rds"))
log_msg("passo ", PASSO, " concluído")
