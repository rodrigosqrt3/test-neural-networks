# Experimento 5: pós-processamento multi-estação da Tmax (43 estações do INMET no
# RS x GEFSv12), horizontes 1, 3 e 7, três dobras temporais (teste 2009–2011,
# 2012–2015, 2016–2019). Fase 1.3 do PLANO.md.
#
# Uso:  Rscript analises/15_experimento_multi.R <horizonte>
#   variáveis de ambiente: N_REDES (10), TORCH_THREADS (2), QRF_THREADS (2), SMOKE=1
#   SMOKE=1: 4 estações antigas, 1 rede, 3 épocas, QRF com 50 árvores (só testa o código)
# Requer dados/processados/multi_*.rds (analises/14_preparar_multiestacao.R).
#
# Modelos (todos ajustados só no treino de cada dobra; mesmas linhas avaliadas):
#   clim           climatologia da estação (+-15 dias)
#   emos_loc       EMOS local (Gneiting et al. 2005): mu ~ média do ensemble + sazonalidade,
#                  log sigma ~ log dp do ensemble; um ajuste por estação. REFERÊNCIA dos skills
#   emos_loc_obs   EMOS local + Tmax da véspera, erro de ontem e ensemble defasado
#   qrf            quantile regression forest global com todos os preditores
#   drn            rede distribucional normal (Rasp & Lerch 2018) com embedding de estação
#   bgp            nossa rede: corpo Bernstein + cauda GPD, embedding de estação
# As redes usam TODOS os preditores: ensemble (Tmax, T2m, chuva), ~23 campos do membro
# de controle, ensemble defasado, observações na hora da rodada, estação, sazonalidade.
#
# Saída: resultados/exp5/h<H>/<modelo>_<dobra>.rds (escores por caso de validação e
# teste) e .quantis_<modelo>_<dobra>.rds (99 quantis por caso, para combinações).
# Hiperparâmetros das redes fixados ANTES de ver o teste (herdados dos exp. 1–4,
# camadas maiores pelo volume de dados); nada é escolhido olhando o teste.

source("R/avaliacao.R"); source("R/modelos_classicos.R"); source("R/modelos_nn.R")
source("R/experimento.R"); source("R/multi.R")

args <- commandArgs(TRUE)
H <- as.integer(if (length(args)) args[[1]] else 1)
SMOKE <- Sys.getenv("SMOKE") == "1"
N_REDES <- if (SMOKE) 1L else as.integer(Sys.getenv("N_REDES", "10"))
EPOCAS <- if (SMOKE) 3L else 300L
QRF_ARVORES <- if (SMOKE) 50L else 500L
QRF_THREADS <- as.integer(Sys.getenv("QRF_THREADS", "2"))
torch_set_num_threads(as.integer(Sys.getenv("TORCH_THREADS", "2")))
REDE <- list(ocultas = c(128, 64), dropout = 0, dim_emb = 4, lote = 256, lr = 1e-3, paciencia = 20)

saida <- file.path("resultados", if (SMOKE) "exp5_smoke" else "exp5", paste0("h", H))
dir.create(saida, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) {
  linha <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(linha, "\n"); cat(linha, "\n", file = file.path(saida, "log.txt"), append = TRUE)
}

tab <- montar_tabela_multi(H, estacoes = if (SMOKE) c("A801", "A802", "A803", "A805"))
FEATS <- c(F_ENS, F_CTRL, tab$f_defas, F_OBS_INIT, F_EST, F_FIXOS)
GRUPOS_FALTA <- list(falta_obs = F_OBS_INIT, falta_est = F_EST, falta_ctrl = F_CTRL, falta_defas = tab$f_defas)
X_EMOS <- c("ens_media", "doy_sen", "doy_cos")
X_EMOS_OBS <- c(X_EMOS, "obs_vespera", "erro_ontem", grep("^defas_media", tab$f_defas, value = TRUE))
log_msg(sprintf("horizonte %d | %d casos, %d estações, %d preditores | N_REDES=%d%s",
                H, nrow(tab$d), uniqueN(tab$d$estacao), length(FEATS), N_REDES, if (SMOKE) " (SMOKE)" else ""))

MODELOS <- list(
  clim         = function(sp) modelo_clim(sp),
  emos_loc     = function(sp) modelo_emos_loc(sp, X_EMOS),
  emos_loc_obs = function(sp) modelo_emos_loc(sp, X_EMOS_OBS),
  qrf          = function(sp) modelo_qrf(sp, sp$feats, QRF_ARVORES, QRF_THREADS),
  bgp          = function(sp) do.call(modelo_rede, c(list(sp, sp$feats, "bgp", tab$n_estacoes, N_REDES, EPOCAS), REDE)),
  drn          = function(sp) do.call(modelo_rede, c(list(sp, sp$feats, "drn", tab$n_estacoes, N_REDES, EPOCAS), REDE))
)

for (dobra in names(DOBRAS)) {
  feitos <- file.exists(file.path(saida, sprintf("%s_%s.rds", names(MODELOS), dobra)))
  if (all(feitos)) next
  sp <- dividir_dobra(tab$d, dobra)
  if (!nrow(sp$te) || !nrow(sp$va)) { log_msg(dobra, ": sem casos de validação/teste, pulando"); next }
  estacoes <- sp$estacoes
  sp <- imputar(sp, FEATS, GRUPOS_FALTA)
  lim <- list(va = limiares_estacao(sp$tr, sp$va), te = limiares_estacao(sp$tr, sp$te))
  log_msg(sprintf("%s | %d estações | treino %d, validação %d, teste %d casos", dobra,
                  uniqueN(sp$tr$estacao), nrow(sp$tr), nrow(sp$va), nrow(sp$te)))
  for (m in names(MODELOS)) {
    arq <- file.path(saida, sprintf("%s_%s.rds", m, dobra))
    if (file.exists(arq)) next
    t0 <- Sys.time()
    mod <- MODELOS[[m]](sp)
    va <- avaliar_em_blocos(mod$va, sp$va, lim$va)
    te <- avaliar_em_blocos(mod$te, sp$te, lim$te)
    seg <- as.numeric(Sys.time() - t0, units = "secs")
    saveRDS(list(q_va = va$quantis, q_te = te$quantis),
            file.path(saida, sprintf(".quantis_%s_%s.rds", m, dobra)))
    saveRDS(list(config = m, dobra = dobra, horizonte = H, validacao = va$escores, teste = te$escores,
                 info = c(mod$info, list(segundos = seg, n_treino = nrow(sp$tr), estacoes = estacoes))),
            arq)                                    # por último: sua existência marca "feito"
    log_msg(sprintf("  %-13s %s | CRPS val %.4f teste %.4f | tw95 teste %.4f | %.0f s", m, dobra,
                    mean(va$escores$crps), mean(te$escores$crps), mean(te$escores$tw_q95), seg))
  }
}
log_msg("horizonte ", H, " concluído")
