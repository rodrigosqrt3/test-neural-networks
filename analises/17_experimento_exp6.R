# Experimento 6: o modelo de referência da fase 1 (PLANO.md 1.3). Mesmas estações,
# dobras e linhas de teste do exp5, para comparar direto com ele.
#
# Uso:  Rscript analises/17_experimento_exp6.R <horizonte>
#   variáveis de ambiente: N_REDES (5), TORCH_THREADS (3), QRF_THREADS (3), GBM_THREADS (3), SMOKE=1
#   SMOKE=1: 4 estações, 1 rede, 3 épocas, 50 árvores (só testa o código)
#
# O que muda em relação ao exp5 (cada item é um degrau, medido separadamente):
#   1. Duas etapas (R/multi.R): o que é escolhido na validação (épocas, iterações) fica
#      fixo e o modelo é reajustado em treino + validação antes de prever o teste. Vale
#      para TODOS os modelos, então a comparação continua com a mesma informação.
#   2. Preditores novos: Tmax do controle por janela de 6 h e climatologia da estação
#      (média, dp e anomalia do ensemble; só do treino). Os erros recentes do GEFS
#      (F_ERROS) foram testados no piloto na validação da f1 e pioraram (-0,4 a -1,1%):
#      ficam fora.
#   3. Baselines da literatura com TODOS os preditores: EMOS com boosting (Messner et al.
#      2017), QRF (Taillardat et al. 2016) e boosting de árvores heterocedástico (xgboost).
#   4. Rede com atenção entre estações (cada dia é um conjunto de estações; viés de
#      distância aprendido; nº de camadas escolhido na validação).
#
# Piloto (validação da f1, 2007–08, anterior a todos os testes; h1, 1 semente; CRPS):
#   DRN do exp5 0,901 | + janelas 0,884 | + climatologia 0,887 | + janelas e clim 0,887 |
#   + erros recentes 0,902–0,911 | atenção 0 / 1 / 2 camadas 0,882 / 0,899 / 0,920
#   (com 5 estações não há vizinhas o bastante; nas dobras seguintes a validação decide).
#
# Modelos:
#   clim, emos_loc   como no exp5 (emos_loc é a REFERÊNCIA dos skills), em duas etapas
#   emos_bst         EMOS local com boosting, todos os preditores
#   qrf              QRF global, todos os preditores
#   gbm              xgboost global: média + variância (resíduos fora da amostra)
#   drn_tv           a DRN do exp5 (mesmos preditores e receita), só que em duas etapas
#   drn_mais         DRN + preditores novos, receita escolhida na validação
#   att              rede com atenção entre estações, preditores novos
# A combinação (seleção de conjunto de Caruana et al. 2004, pesos na validação) é feita
# depois, em analises/18_analisar_exp6.R.
#
# Saída: resultados/exp6/h<H>/<modelo>_<dobra>.rds e .quantis_<modelo>_<dobra>.rds, no
# formato do exp5. q_va vem da etapa 1 (ajuste só no treino), q_te do modelo final.

source("R/avaliacao.R"); source("R/modelos_classicos.R"); source("R/modelos_nn.R")
source("R/experimento.R"); source("R/multi.R")

args <- commandArgs(TRUE)
H <- as.integer(if (length(args)) args[[1]] else 1)
SMOKE <- Sys.getenv("SMOKE") == "1"
N_REDES <- if (SMOKE) 1L else as.integer(Sys.getenv("N_REDES", "5"))
EPOCAS <- if (SMOKE) 3L else 300L
QRF_ARVORES <- if (SMOKE) 50L else 500L
QRF_THREADS <- as.integer(Sys.getenv("QRF_THREADS", "3"))
GBM_THREADS <- as.integer(Sys.getenv("GBM_THREADS", "3"))
torch_set_num_threads(as.integer(Sys.getenv("TORCH_THREADS", "3")))
# Receitas das redes. drn_tv usa só a do exp5; nas outras, a receita é escolhida na
# validação de cada dobra (R/multi.R: .escolher). Na atenção, o nº de camadas é uma das
# escolhas: com 0 camadas a rede é a MLP treinada em lotes de dias, então a atenção só
# entra onde a validação mostra que ajuda.
REDE5 <- list(exp5 = list(ocultas = c(128, 64), dropout = 0, dim_emb = 4, lote = 256, lr = 1e-3, paciencia = 20))
REDES6 <- list(
  base  = list(ocultas = c(128, 64), dropout = 0, dim_emb = 4, lote = 256, lr = 1e-3, paciencia = 20),
  drop  = list(ocultas = c(128, 64), dropout = 0.1, dim_emb = 4, lote = 256, lr = 1e-3, paciencia = 20, wd = 1e-4),
  larga = list(ocultas = c(256, 128, 64), dropout = 0.1, dim_emb = 8, lote = 256, lr = 5e-4, paciencia = 20, wd = 1e-4))
ATENCOES <- list(
  c0 = list(ocultas = c(128, 64), camadas = 0, dim_emb = 4, lote = 32, lr = 1e-3, paciencia = 20, n_busca = 1),
  c1 = list(ocultas = c(128, 64), camadas = 1, cabecas = 2, dim_emb = 4, lote = 32, lr = 3e-4, wd = 1e-3,
            paciencia = 20, n_busca = 1),
  c2 = list(ocultas = c(128, 64), camadas = 2, cabecas = 4, dim_emb = 4, dropout = 0.1, lote = 32, lr = 3e-4,
            wd = 1e-3, paciencia = 20, n_busca = 1))

saida <- file.path("resultados", if (SMOKE) "exp6_smoke" else "exp6", paste0("h", H))
dir.create(saida, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) {
  linha <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(linha, "\n"); cat(linha, "\n", file = file.path(saida, "log.txt"), append = TRUE)
}

tab <- montar_tabela_multi(H, estacoes = if (SMOKE) c("A801", "A802", "A803", "A805"))
FEATS5 <- c(F_ENS, F_CTRL, tab$f_defas, F_OBS_INIT, F_EST, F_FIXOS)
FEATS6 <- c(FEATS5, F_JANELAS, F_CLIM)
GRUPOS5 <- list(falta_obs = F_OBS_INIT, falta_est = F_EST, falta_ctrl = F_CTRL, falta_defas = tab$f_defas)
GRUPOS6 <- GRUPOS5
X_EMOS <- c("ens_media", "doy_sen", "doy_cos")
log_msg(sprintf("horizonte %d | %d casos, %d estações, %d preditores (exp5: %d) | N_REDES=%d%s",
                H, nrow(tab$d), uniqueN(tab$d$estacao), length(FEATS6), length(FEATS5), N_REDES,
                if (SMOKE) " (SMOKE)" else ""))

rede <- function(sp, f, cabeca, receitas) modelo_rede2(sp, f, cabeca, tab$n_estacoes, N_REDES, EPOCAS, receitas)
atencao <- function(sp, cabeca) modelo_atencao(sp, sp$f6, cabeca, tab$n_estacoes, tab$dist, N_REDES, EPOCAS, ATENCOES)
MODELOS <- list(
  clim     = function(sp) duas_etapas(modelo_clim, sp),
  emos_loc = function(sp) duas_etapas(function(s) modelo_emos_loc(s, X_EMOS), sp),
  emos_bst = function(sp) modelo_emos_bst(sp, sp$f6, m_max = if (SMOKE) 50 else 2000),
  gbm      = function(sp) modelo_gbm(sp, sp$f6, GBM_THREADS, n_max = if (SMOKE) 100 else 3000),
  qrf      = function(sp) duas_etapas(function(s) modelo_qrf(s, s$f6, QRF_ARVORES, QRF_THREADS), sp),
  drn_tv   = function(sp) rede(sp, sp$f5, "drn", REDE5),
  drn_mais = function(sp) rede(sp, sp$f6, "drn", REDES6),
  att      = function(sp) atencao(sp, "drn")
)

for (dobra in names(DOBRAS)) {
  feitos <- file.exists(file.path(saida, sprintf("%s_%s.rds", names(MODELOS), dobra)))
  if (all(feitos)) next
  sp <- dividir_dobra(tab$d, dobra)
  if (!nrow(sp$te) || !nrow(sp$va)) { log_msg(dobra, ": sem casos de validação/teste, pulando"); next }
  estacoes <- sp$estacoes
  sp <- imputar(adicionar_clim(sp), FEATS6, GRUPOS6)
  sp$f6 <- sp$feats
  sp$f5 <- intersect(c(FEATS5, names(GRUPOS5)), sp$feats)     # exatamente o conjunto do exp5
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
                 info = c(mod$info, list(segundos = seg, n_treino = nrow(sp$tr), n_validacao = nrow(sp$va),
                                         estacoes = estacoes))),
            arq)                                    # por último: sua existência marca "feito"
    log_msg(sprintf("  %-9s %s | CRPS val %.4f teste %.4f | tw95 teste %.4f | %.0f s", m, dobra,
                    mean(va$escores$crps), mean(te$escores$crps), mean(te$escores$tw_q95), seg))
  }
}
log_msg("horizonte ", H, " concluído")
