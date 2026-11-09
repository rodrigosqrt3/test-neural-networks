# Experimento 7 (fase 2): chuva diária nas estações do INMET no RS x GEFSv12, h1/h3/h7,
# mesmas dobras do exp5/6. A pergunta é se uma cauda GPD na rede (corpo Bernstein + GPD,
# "bgp") melhora a previsão de chuva FORTE em relação às cabeças sem cauda explícita.
#
# Uso:  Rscript analises/24_experimento_chuva.R <horizonte>
#   variáveis de ambiente: N_REDES (5), TORCH_THREADS, QRF_THREADS, SMOKE=1
# Requer dados/processados/multi_obs_chuva.rds (analises/23_preparar_chuva.R).
#
# A distribuição da chuva tem massa em zero. Todas as redes são CENSURADAS em 0: os quantis
# são max(Q(tau), 0) dentro da perda (CRPS), o que dá P(Y = 0) = P(Q latente <= 0) e mantém
# o escore próprio (como no EMOS censurado de Scheuerer 2014 / Scheuerer & Hamill 2015).
#
# Modelos (todos em duas etapas, como no exp6; mesmas linhas avaliadas):
#   clim     climatologia da estação (+-15 dias); REFERÊNCIA dos skills (CRPSS habitual em chuva)
#   qrf      quantile regression forest (Taillardat et al. 2019 para chuva), todos os preditores
#   drn_c3   normal censurada na raiz cúbica (a forma do EMOS de chuva clássico), como rede
#   bqn_c3   Bernstein censurado na raiz cúbica (Bremnes 2020)
#   bgp_c3   Bernstein + cauda GPD, na raiz cúbica: isola o efeito da cauda na mesma escala
#   bgp      Bernstein + cauda GPD na escala ORIGINAL (mm): a cauda GPD de fato da chuva
#   bgp_xi   o mesmo, com xi(x) dependendo dos preditores e encolhido para a média
# O stack vem depois, no script 18.
#
# Avaliação: as mesmas colunas do exp6. Limiares: quantis sazonais q90/q95/q99 da chuva da
# estação (com os dias secos) e o absoluto de 35 mm (coluna "abs35").

source("R/avaliacao.R"); source("R/modelos_classicos.R"); source("R/modelos_nn.R")
source("R/experimento.R"); source("R/multi.R")

args <- commandArgs(TRUE)
H <- as.integer(if (length(args)) args[[1]] else 1)
SMOKE <- Sys.getenv("SMOKE") == "1"
N_REDES <- if (SMOKE) 1L else as.integer(Sys.getenv("N_REDES", "5"))
EPOCAS <- if (SMOKE) 3L else 300L
QRF_ARVORES <- if (SMOKE) 50L else 500L
QRF_THREADS <- as.integer(Sys.getenv("QRF_THREADS", "2"))
torch_set_num_threads(as.integer(Sys.getenv("TORCH_THREADS", "2")))
# clip: corte do gradiente (norma 5). Na rodada de 25/10 o bgp_xi deu NaN em h1 e h7, f3; a
# causa era o ramo descartado do torch_where na cauda GPD com xi ~ 0 (0/0 no gradiente),
# corrigida em R/modelos_nn.R. A saída da rede não muda; só o gradiente deixa de ser NaN.
REDE <- list(base = list(ocultas = c(128, 64), dropout = 0, dim_emb = 4, lote = 256, lr = 1e-3, paciencia = 20, clip = 5))
C3 <- list(f = function(y) pmax(y, 0)^(1 / 3), inv = function(t) pmax(t, 0)^3)

saida <- file.path("resultados", if (SMOKE) "exp7_smoke" else "exp7", paste0("h", H))
dir.create(saida, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) {
  linha <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(linha, "\n"); cat(linha, "\n", file = file.path(saida, "log.txt"), append = TRUE)
}

# a tabela do exp5/6 com o alvo trocado: a coluna `tmax` passa a ser a CHUVA observada
# (o resto da pipeline usa esse nome para o alvo). Ficam os dias com Tmax e chuva medidas.
tab <- montar_tabela_multi(H, estacoes = if (SMOKE) c("A801", "A802", "A803", "A805"))
obs <- readRDS("dados/processados/multi_obs_chuva.rds")
d <- merge(tab$d, obs[, .(estacao, valid_date = dia, chuva)], by = c("estacao", "valid_date"))
d[, tmax_obs := tmax][, tmax := chuva][, chuva := NULL]
F_JANELAS_CHUVA <- c("chuva_h06", "chuva_h12", "chuva_h18", "chuva_h24")
for (v in F_JANELAS_CHUVA) set(d, j = v, value = log1p(pmax(d[[v]], 0)))
setorder(d, valid_date, estacao)

# climatologia da chuva da estação no treino da dobra: P(chuva >= 1 mm) e média, +-15 dias
adicionar_clim_chuva <- function(sp) {
  tabs <- lapply(split(sp$tr, by = "estacao"), function(a) {
    doy <- as.integer(format(a$valid_date, "%j"))
    t(vapply(1:366, function(k) {
      y <- a$tmax[pmin(abs(doy - k), 366 - abs(doy - k)) <= 15]
      c(mean(y >= 1), log1p(mean(y)))
    }, numeric(2)))
  })
  aplicar <- function(x) {
    x <- copy(x); doy <- as.integer(format(x$valid_date, "%j"))
    p <- m <- rep(NA_real_, nrow(x))
    for (e in names(tabs)) { i <- which(x$estacao == e); p[i] <- tabs[[e]][doy[i], 1]; m[i] <- tabs[[e]][doy[i], 2] }
    x[, `:=`(clim_p1mm = p, clim_media = m)]
  }
  sp$tr <- aplicar(sp$tr); sp$va <- aplicar(sp$va); sp$te <- aplicar(sp$te)
  sp
}
FEATS <- c(F_ENS, F_CTRL, F_JANELAS_CHUVA, tab$f_defas, F_OBS_INIT, F_EST, F_FIXOS, "clim_p1mm", "clim_media")
GRUPOS <- list(falta_obs = F_OBS_INIT, falta_est = F_EST, falta_ctrl = F_CTRL, falta_defas = tab$f_defas)
log_msg(sprintf("horizonte %d | %d casos, %d estações, %d preditores | %.0f%% dias com >= 1 mm | N_REDES=%d%s",
                H, nrow(d), uniqueN(d$estacao), length(FEATS), 100 * mean(d$tmax >= 1), N_REDES,
                if (SMOKE) " (SMOKE)" else ""))

rede <- function(sp, cabeca, transf, receita = REDE) modelo_rede2(sp, sp$feats, cabeca, tab$n_estacoes, N_REDES, EPOCAS,
                                                                 receita, transf = transf, piso = 0)
# xi(x): a forma da cauda depende das condições, encolhida para a média do lote com
# peso lambda_xi (R/modelos_nn.R; análogo à penalização do ERF). Candidata à fase 2.
REDE_XI <- list(xi_x = c(REDE$base, list(xi_modo = "x", lambda_xi = 0.1)))
MODELOS <- list(
  clim   = function(sp) duas_etapas(modelo_clim, sp),
  qrf    = function(sp) duas_etapas(function(s) modelo_qrf(s, s$feats, QRF_ARVORES, QRF_THREADS), sp),
  drn_c3 = function(sp) rede(sp, "drn", C3),
  bqn_c3 = function(sp) rede(sp, "bqn", C3),
  bgp_c3 = function(sp) rede(sp, "bgp", C3),
  bgp    = function(sp) rede(sp, "bgp", NULL),
  bgp_xi = function(sp) rede(sp, "bgp", NULL, REDE_XI)
)

for (dobra in names(DOBRAS)) {
  feitos <- file.exists(file.path(saida, sprintf("%s_%s.rds", names(MODELOS), dobra)))
  if (all(feitos)) next
  sp <- dividir_dobra(d, dobra)
  if (!nrow(sp$te) || !nrow(sp$va)) { log_msg(dobra, ": sem casos de validação/teste, pulando"); next }
  estacoes <- sp$estacoes
  sp <- imputar(adicionar_clim_chuva(sp), FEATS, GRUPOS)
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
    saveRDS(list(q_va = va$quantis, q_te = te$quantis), file.path(saida, sprintf(".quantis_%s_%s.rds", m, dobra)))
    saveRDS(list(config = m, dobra = dobra, horizonte = H, validacao = va$escores, teste = te$escores,
                 info = c(mod$info, list(segundos = seg, n_treino = nrow(sp$tr), n_validacao = nrow(sp$va),
                                         estacoes = estacoes))),
            arq)
    log_msg(sprintf("  %-7s %s | CRPS val %.4f teste %.4f | tw95 teste %.4f | tw35 teste %.4f | %.0f s", m, dobra,
                    mean(va$escores$crps), mean(te$escores$crps), mean(te$escores$tw_q95),
                    mean(te$escores$tw_abs35), seg))
  }
}
log_msg("horizonte ", H, " concluído")
