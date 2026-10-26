# Experimento 3: deriva do viés do GEFS e calibração da cauda inferior.
#
# Uso:  Rscript analises/04_experimento_deriva.R <horizonte>
#       variáveis de ambiente: N_REDES (padrão 10), TORCH_THREADS (padrão 4), SMOKE=1
#
# Motivação: o viés do GEFS em relação à A801 passou de ~+0.1 °C (2000–2010)
# para ~-0.4 °C (2011–2019), e todos os modelos do experimento 1 subestimam as
# excedências do q95 no teste. Estratégias comparadas (modelo-base: BGP com os
# preditores "obs" e dropout 0, o melhor do experimento 1):
#
#   base           treino 2000–2011, early stopping em 2012–2015 (referência)
#   trva           retreina em 2000–2015 com o nº de épocas escolhido no base
#   peso4, peso8   2000–2015 com peso 0.5^(idade/meia-vida), meia-vida 4 ou 8 anos
#   tendencia      base + ano como preditor
#   *_online30/90  correção online do viés com os erros dos dias já observados
#                  (meia-vida de 30 ou 90 dias), aplicada depois da previsão
#   grade, grade_g20  grade de treino fina também na cauda INFERIOR (e grau 20),
#                  para corrigir a distorção do PIT perto de 0 vista no exp. 1
#   combo          peso4 + grade + online90
#
# Atenção: modelos treinados com 2012–2015 têm escores de validação dentro da
# amostra; só o teste (2016–2019) é comparável. Resultados em resultados/exp3/.

source("R/avaliacao.R"); source("R/modelos_classicos.R"); source("R/modelos_nn.R")
source("R/experimento.R")

args <- commandArgs(TRUE)
H <- as.integer(if (length(args)) args[[1]] else 1)
SMOKE <- Sys.getenv("SMOKE") == "1"
N_REDES <- if (SMOKE) 1L else as.integer(Sys.getenv("N_REDES", "10"))
EPOCAS <- if (SMOKE) 15L else 400L
torch_set_num_threads(as.integer(Sys.getenv("TORCH_THREADS", "4")))

saida <- file.path("resultados", if (SMOKE) "exp3_smoke" else "exp3", paste0("h", H))
dir.create(saida, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) {
  linha <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(linha, "\n"); cat(linha, "\n", file = file.path(saida, "log.txt"), append = TRUE)
}

dados <- carregar_horizonte(H)
tr <- dados$tr; va <- dados$va; te <- dados$te
trva <- rbind(tr, va)
LIM <- list(validacao = limiares_para(tr, va$valid_date), teste = limiares_para(tr, te$valid_date))
log_msg(sprintf("horizonte %d | treino %d, validacao %d, teste %d | N_REDES=%d%s",
                H, nrow(tr), nrow(va), nrow(te), N_REDES, if (SMOKE) " (SMOKE)" else ""))

feito <- function(nome) file.exists(file.path(saida, paste0(nome, ".rds")))
salvar <- function(nome, q_va, q_te, extra = list()) {
  res <- c(list(
    config = nome, horizonte = H,
    validacao = cbind(data.table(data = va$valid_date), escores_por_caso(q_va, va$tmax, LIM$validacao)),
    teste = cbind(data.table(data = te$valid_date), escores_por_caso(q_te, te$tmax, LIM$teste))
  ), extra)
  saveRDS(res, file.path(saida, paste0(nome, ".rds")))
  d <- res$teste
  log_msg(sprintf("  %-22s CRPS teste %.4f | tw95 %.4f | exc95 obs %.3f prev %.3f | PIT<0.05 %.3f",
                  nome, mean(d$crps), mean(d$tw_q95), mean(d$y > d$lim_q95), mean(d$p_q95),
                  mean(d$pit < 0.05)))
}

# Aplica a correção online a um par de grades (validação + teste em sequência)
salvar_online <- function(nome, q_va, q_te, meia_vida, extra = list()) {
  if (feito(nome)) return(invisible())
  todos <- rbind(va, te)
  mediana <- c(q_va[, 500], q_te[, 500])            # TAU_GRADE[500] = 0.5
  c_todos <- correcao_online(todos$init_date, todos$valid_date, todos$tmax, mediana, meia_vida)
  n_va <- nrow(va)
  salvar(nome, q_va + c_todos[seq_len(n_va)], q_te + c_todos[-seq_len(n_va)],
         c(extra, list(correcao_teste = c_todos[-seq_len(n_va)])))
}

# Grades de quantis guardadas em memória para as correções online
grades <- list()
cache <- function(nome) file.path(saida, paste0(".grades_", nome, ".rds"))

treinar <- function(nome, ajuste, feats, pesos = NULL, epocas = EPOCAS,
                    usar_validacao = TRUE, ...) {
  if (file.exists(cache(nome))) {
    grades[[nome]] <<- readRDS(cache(nome))
    return(invisible())
  }
  padr <- padronizador(ajuste, feats)
  t0 <- Sys.time()
  conj <- treinar_conjunto(
    Xtr = padr(ajuste), ytr = ajuste$tmax,
    Xva = if (usar_validacao) padr(va) else NULL, yva = if (usar_validacao) va$tmax else NULL,
    pesos = pesos, n_redes = N_REDES, epocas = epocas, cabeca = "bgp", dropout = 0, ...
  )
  g <- list(q_va = prever_conjunto(conj, padr(va)), q_te = prever_conjunto(conj, padr(te)),
            melhor_epoca = sapply(conj, `[[`, "melhor_epoca"),
            xi = mean(xi_conjunto(conj, padr(te))),
            segundos = as.numeric(Sys.time() - t0, units = "secs"))
  saveRDS(g, cache(nome))
  grades[[nome]] <<- g
  salvar(nome, g$q_va, g$q_te, list(segundos = g$segundos, xi_media = g$xi,
                                    validacao_na_amostra = !usar_validacao))
}

pesos_recencia <- function(x, meia_vida_anos) {
  idade <- as.numeric(max(x$valid_date) - x$valid_date) / 365.25
  0.5^(idade / meia_vida_anos)
}

TAU_SIMETRICA <- c((1:100 - 0.5) / 1000, (10:89 + 0.5) / 100, 0.9 + (1:100 - 0.5) / 1000)

# --- modelos treinados ---------------------------------------------------------
treinar("base", tr, F_OBS)
ep <- max(10L, as.integer(stats::median(grades$base$melhor_epoca)))
log_msg("  épocas para os modelos sem validação: ", ep)
treinar("trva", trva, F_OBS, epocas = ep, usar_validacao = FALSE)
treinar("peso4", trva, F_OBS, pesos = pesos_recencia(trva, 4), epocas = ep,
        usar_validacao = FALSE)
treinar("peso8", trva, F_OBS, pesos = pesos_recencia(trva, 8), epocas = ep,
        usar_validacao = FALSE)
treinar("tendencia", tr, c(F_OBS, "ano_frac"))
treinar("grade", tr, F_OBS, tau_treino = TAU_SIMETRICA)
treinar("grade_g20", tr, F_OBS, tau_treino = TAU_SIMETRICA, grau = 20)
treinar("peso4_grade", trva, F_OBS, pesos = pesos_recencia(trva, 4), epocas = ep,
        usar_validacao = FALSE, tau_treino = TAU_SIMETRICA)

# --- correções online (sem treino adicional) -----------------------------------------
for (m in c("base", "trva", "peso4")) {
  for (mv in c(30, 90)) {
    salvar_online(sprintf("%s_online%d", m, mv), grades[[m]]$q_va, grades[[m]]$q_te, mv,
                  list(validacao_na_amostra = m != "base"))
  }
}
salvar_online("combo", grades$peso4_grade$q_va, grades$peso4_grade$q_te, 90,
              list(validacao_na_amostra = TRUE))

# EMOS com e sem correção online, como referência clássica
xe <- c("ens_media", "obs_vespera", "doy_sen", "doy_cos")
ze <- function(x) cbind(log(x$ens_dp + 0.05), x$doy_sen, x$doy_cos)
emos_va <- predict(emos_ajustar(as.matrix(tr[, ..xe]), ze(tr), tr$tmax), as.matrix(va[, ..xe]), ze(va))
emos_te <- predict(emos_ajustar(as.matrix(trva[, ..xe]), ze(trva), trva$tmax), as.matrix(te[, ..xe]), ze(te))
if (!feito("emos_ref")) salvar("emos_ref", emos_va, emos_te)
salvar_online("emos_online90", emos_va, emos_te, 90)

log_msg("horizonte ", H, " concluído")
