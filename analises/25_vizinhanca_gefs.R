# H6 do exp8: o GEFS EM VOLTA da estação, não só no ponto. Para cada rodada (membro de
# controle), horizonte e estação, num quadrado de 7x7 pontos de grade (0,25°, ~±0,75°)
# centrado no ponto de grade mais próximo da estação, e para cada campo-chave:
#   viz_<campo>_m   média do quadrado
#   viz_<campo>_s   desvio padrão (heterogeneidade: frente, borda de nuvem, costa)
#   viz_<campo>_dx  leste menos oeste (média da coluna de borda), gradiente zonal
#   viz_<campo>_dy  norte menos sul, gradiente meridional (advecção quente/fria)
# Resumos da vizinhança no lugar de uma CNN sobre o recorte: com ~10^5 casos, resumos
# bem escolhidos costumam render o mesmo com muito menos risco de sobreajuste.
#
# Uso: Rscript analises/25_vizinhanca_gefs.R       (~5-15 min; lê as 7.305 rodadas do controle)
#      TESTE=1 -> só 2010 (verificação rápida, saída *_teste.rds)
#      FONTE=oper -> GEFS operacional 2020-09-24 a 2025-12-31 (exp14), saída *_oper.rds
# Saída: dados/processados/multi_vizinhanca.rds (estacao, init_date, horizon, viz_*)
library(data.table)
TESTE <- Sys.getenv("TESTE") == "1"
OPER <- Sys.getenv("FONTE") == "oper"
CAMPOS <- c("tmax2m_h18", "t2m_h18", "t850_h18", "q2m_h18", "nuvens_h18", "radsw_h18", "pnmm_h18", "z500_h18")
R <- 3                                                   # raio: quadrado de 2R+1 = 7 pontos
t0 <- Sys.time()

est <- readRDS("dados/processados/multi_estacoes.rds")
anos <- if (TESTE) 2010 else if (OPER) 2020:as.integer(Sys.getenv("ANO_FIM", "2025")) else 2000:2019
dias <- seq(if (OPER) as.Date("2020-09-24") else as.Date(sprintf("%d-01-01", min(anos))),
            as.Date(sprintf("%d-12-31", max(anos))), by = "day")
arq <- function(d) sprintf("dados/brutos/%s/c00/%s/%s.rds", if (OPER) "gefs_oper" else "gefs_aux",
                           format(d, "%Y"), format(d, "%Y%m%d"))
g0 <- readRDS(arq(dias[file.exists(arq(dias))][1]))
lon <- g0$grade$lon; lat <- g0$grade$lat              # matrizes de campo: [lat, lon]
ii <- vapply(est$lat, function(y) which.min(abs(lat - y)), 1L)
jj <- vapply(est$lon, function(x) which.min(abs(lon - x)), 1L)
ii <- pmin(pmax(ii, R + 1), length(lat) - R); jj <- pmin(pmax(jj, R + 1), length(lon) - R)
desl <- -R:R
# índices lineares do quadrado (n_est x 49) e das bordas leste/oeste/norte/sul (n_est x 7)
lin <- function(i, j) (j - 1) * length(lat) + i
Q <- t(mapply(function(i, j) as.vector(outer(i + desl, j + desl, lin)), ii, jj))
L <- t(mapply(function(i, j) lin(i + desl, j + R), ii, jj)); O <- t(mapply(function(i, j) lin(i + desl, j - R), ii, jj))
N <- t(mapply(function(i, j) lin(ifelse(lat[1] > lat[2], i - R, i + R), j + desl), ii, jj))
S <- t(mapply(function(i, j) lin(ifelse(lat[1] > lat[2], i + R, i - R), j + desl), ii, jj))
resumo <- function(M) {
  q <- matrix(M[Q], nrow(Q)); m <- rowMeans(q)
  cbind(m = m, s = sqrt(pmax(rowMeans(q^2) - m^2, 0)),
        dx = rowMeans(matrix(M[L], nrow(L))) - rowMeans(matrix(M[O], nrow(O))),
        dy = rowMeans(matrix(M[N], nrow(N))) - rowMeans(matrix(M[S], nrow(S))))
}

out <- vector("list", length(dias))
for (k in seq_along(dias)) {
  f <- arq(dias[k]); if (!file.exists(f)) next
  x <- readRDS(f)
  out[[k]] <- rbindlist(lapply(c(1, 3, 7), function(H) {
    campos <- x[[paste0("H", H)]]
    b <- data.table(estacao = est$estacao, init_date = dias[k], horizon = H)
    for (cp in CAMPOS) {
      M <- campos[[cp]]
      r <- if (is.null(M)) matrix(NA_real_, nrow(est), 4) else resumo(M)
      nm <- paste0("viz_", sub("_h18$", "", cp), "_", c("m", "s", "dx", "dy"))
      for (c in 1:4) set(b, j = nm[c], value = r[, c])
    }
    b
  }))
  if (k %% 1000 == 0) cat(format(Sys.time(), "%H:%M:%S"), k, "de", length(dias), "\n")
}
out <- rbindlist(out)
saveRDS(out, sprintf("dados/processados/multi_vizinhanca%s.rds", if (TESTE) "_teste" else if (OPER) "_oper" else ""))
cat(sprintf("pronto: %d linhas x %d colunas | NA em %.2f%% das células | %.1f min\n", nrow(out), ncol(out),
            100 * mean(is.na(as.matrix(out[, -(1:3)]))), as.numeric(Sys.time() - t0, units = "mins")))
