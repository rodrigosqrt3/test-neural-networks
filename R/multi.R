# Experimento multi-estação (fase 1.3 do PLANO.md): tabela de modelagem,
# dobras temporais, imputação, limiares por estação, avaliação em blocos e os
# modelos (climatologia, EMOS local, QRF, redes com embedding de estação).
#
# Princípios (README.md, regras 1–5):
#   - Todo modelo é ajustado SÓ no bloco de treino da dobra. A validação serve
#     apenas para o early stopping das redes; o teste é tocado uma vez. Assim
#     todos os modelos usam exatamente a mesma informação.
#   - Imputação, padronização e limiares sazonais vêm só do treino.
#   - Uma estação entra numa dobra se tiver >= MIN_TREINO casos de treino; as
#     linhas de validação/teste são as MESMAS para todos os modelos.
# Depende de R/avaliacao.R, R/modelos_classicos.R, R/modelos_nn.R e R/experimento.R.

library(data.table)

# treino até [1] | validação até [2] | teste até [3]; GAP_DIAS de intervalo entre blocos
DOBRAS <- list(f1 = as.Date(c("2006-12-31", "2008-12-31", "2011-12-31")),
               f2 = as.Date(c("2009-12-31", "2011-12-31", "2015-12-31")),
               f3 = as.Date(c("2011-12-31", "2015-12-31", "2019-12-31")))
GAP_DIAS <- 7
MIN_TREINO <- 365
TAU_FINO <- (1:99) / 100                 # quantis guardados por caso (para combinações)

# --- preditores ------------------------------------------------------------------------
F_ENS <- c("ens_media", "ens_dp", "ens_min", "ens_max", "t2m_ens_media", "t2m_ens_dp",
           "chuva_ens_media", "chuva_ens_dp", "chuva_p1mm")
F_CTRL <- c("ctrl_tmax_dia", "ctrl_chuva_dia", "t2m_h12", "t2m_h18", "q2m_h18", "u10_h18", "v10_h18",
            "pnmm_h06", "pnmm_h18", "nuvens_h12", "nuvens_h18", "radsw_h12", "radsw_h18", "cape_h18",
            "pwat_h18", "fsens_h18", "flat_h18", "solo_h12", "t850_h18", "q850_h18", "u850_h18",
            "v850_h18", "z500_h18")
F_EST <- c("est_t00", "est_td00", "est_ur00", "est_p00", "est_u", "est_v", "est_chuva24",
           "est_rad24", "est_tmin24", "est_dp24")
F_OBS_INIT <- c("obs_vespera", "obs_lag2", "obs_lag3", "erro_ontem")
F_FIXOS <- c("lat", "lon", "altitude", "doy_sen", "doy_cos")
LOG1P <- c("chuva_ens_media", "chuva_ens_dp", "ctrl_chuva_dia", "cape_h18")   # caudas longas

# exp6: Tmax do controle por janela de 6 h, erros recentes do GEFS (na estação, nas
# vizinhas e no estado) e climatologia da estação (esta só dentro da dobra: adicionar_clim)
F_JANELAS <- c("tmax2m_h06", "tmax2m_h12", "tmax2m_h18", "tmax2m_h24")
F_ERROS <- c("erro_h_7d", "erro_h_30d", "eabs_h_30d", "erro_ontem_viz", "erro_ontem_reg")
F_CLIM <- c("clim_media", "clim_dp", "ens_anom")

# rodadas anteriores que preveem o mesmo dia (ensemble defasado)
defasagens <- function(H) switch(as.character(H), "1" = c(3, 7), "3" = 7, numeric())

# distâncias entre estações em km (haversine), na ordem de `est`
distancias_km <- function(est) {
  la <- est$lat * pi / 180; lo <- est$lon * pi / 180
  a <- outer(la, la, function(x, y) sin((x - y) / 2)^2) +
    outer(la, la, function(x, y) cos(x) * cos(y)) * outer(lo, lo, function(x, y) sin((x - y) / 2)^2)
  2 * 6371 * asin(sqrt(pmin(a, 1)))
}

# média móvel de n valores terminando em cada posição, ignorando NA (NA se < min_n válidos)
.media_movel <- function(x, n, min_n) {
  ok <- !is.na(x); s <- cumsum(ifelse(ok, x, 0)); k <- cumsum(ok)
  atras <- function(v) c(rep(0, n), v)[seq_along(v)]
  m <- (s - atras(s)) / (k - atras(k))
  m[(k - atras(k)) < min_n] <- NA
  m
}

# Erros recentes do GEFS conhecidos na rodada 00 UTC de init. Um dia local D só se
# fecha às 03 UTC de D+1, então entram os dias D <= init - 2 (como em correcao_online).
#   erro_h_7d, erro_h_30d  média de (obs - média do ensemble no horizonte H) em
#                          [init-8, init-2] e [init-31, init-2]: o viés recente
#   eabs_h_30d             média do |erro| nos mesmos 30 dias: a escala recente
#   erro_ontem_viz, _reg   erro_ontem médio das 4 estações mais próximas / de todas
#                          (as observações das outras estações também chegam às 00 UTC)
.erros_recentes <- function(d, g, od, oi, est, H) {
  e <- merge(g[horizon == H, .(estacao, dia = valid_date, prev = ens_media)], od[, .(estacao, dia, tmax)],
             by = c("estacao", "dia"))[!is.na(tmax) & !is.na(prev)]
  dias <- seq(min(d$init_date) - 40, max(d$init_date), by = "day")
  rol <- e[, {
    x <- rep(NA_real_, length(dias)); i <- match(dia, dias); ok <- !is.na(i)
    x[i[ok]] <- (tmax - prev)[ok]
    list(init_date = dias + 2, erro_h_7d = .media_movel(x, 7, 4), erro_h_30d = .media_movel(x, 30, 15),
         eabs_h_30d = .media_movel(abs(x), 30, 15))
  }, by = estacao]
  d <- merge(d, rol, by = c("estacao", "init_date"), all.x = TRUE)

  p1 <- g[horizon == 1, .(estacao, init_date = init_date + 2, prev_vespera = ens_media)]
  eo <- merge(oi[, .(estacao, init_date, obs_vespera)], p1, by = c("estacao", "init_date"))
  eo <- eo[!is.na(obs_vespera) & !is.na(prev_vespera), .(estacao, init_date, e = obs_vespera - prev_vespera)]
  W <- dcast(eo, init_date ~ estacao, value.var = "e")
  M <- as.matrix(W[, -1]); D <- distancias_km(est[match(colnames(M), estacao)])
  media <- function(A) { r <- rowMeans(A, na.rm = TRUE); r[!is.finite(r)] <- NA; r }
  reg <- media(M)
  viz <- rbindlist(lapply(seq_len(ncol(M)), function(j) {
    k <- order(D[j, ])[2:5]                                 # 4 mais próximas (a 1ª é ela mesma)
    data.table(estacao = colnames(M)[j], init_date = W$init_date,
               erro_ontem_viz = media(M[, k, drop = FALSE]), erro_ontem_reg = reg)
  }))
  merge(d, viz, by = c("estacao", "init_date"), all.x = TRUE)
}

# --- tabela de modelagem ------------------------------------------------------------------
# Uma linha por (estação, rodada) com alvo observado e ensemble disponível.
# extras = TRUE acrescenta F_ERROS (exp6); as linhas são as mesmas nos dois casos.
montar_tabela_multi <- function(H, estacoes = NULL, extras = FALSE) {
  est <- readRDS("dados/processados/multi_estacoes.rds")
  g <- readRDS("dados/processados/multi_gefs.rds")
  od <- readRDS("dados/processados/multi_obs_dia.rds")
  oi <- readRDS("dados/processados/multi_obs_init.rds")
  if (!is.null(estacoes)) {
    est <- est[estacao %in% estacoes]; g <- g[estacao %in% estacoes]
    od <- od[estacao %in% estacoes]; oi <- oi[estacao %in% estacoes]
  }
  d <- g[horizon == H]
  d <- merge(d, od[, .(estacao, valid_date = dia, tmax)], by = c("estacao", "valid_date"))
  d <- d[!is.na(tmax) & !is.na(ens_media)]
  f_defas <- character()
  for (k in defasagens(H)) {
    nm <- paste0(c("defas_media_h", "defas_dp_h"), k)
    x <- g[horizon == k, .(estacao, valid_date, ens_media, ens_dp)]
    setnames(x, c("ens_media", "ens_dp"), nm)
    d <- merge(d, x, by = c("estacao", "valid_date"), all.x = TRUE)
    f_defas <- c(f_defas, nm)
  }
  d <- merge(d, oi, by = c("estacao", "init_date"), all.x = TRUE)
  # erro de ontem: Tmax parcial da véspera menos a previsão de 1 dia (rodada init-2)
  # para esse mesmo dia; conhecido na hora da rodada
  p1 <- g[horizon == 1, .(estacao, init_date = init_date + 2, prev_vespera = ens_media)]
  d <- merge(d, p1, by = c("estacao", "init_date"), all.x = TRUE)
  d[, erro_ontem := obs_vespera - prev_vespera][, prev_vespera := NULL]
  if (extras) d <- .erros_recentes(d, g, od, oi, est, H)
  d <- merge(d, est[, .(estacao, lat, lon, altitude)], by = "estacao")
  doy <- as.integer(format(d$valid_date, "%j"))
  d[, `:=`(doy_sen = sin(2 * pi * doy / 365.25), doy_cos = cos(2 * pi * doy / 365.25))]
  for (v in LOG1P) set(d, j = v, value = log1p(pmax(d[[v]], 0)))
  todas <- readRDS("dados/processados/multi_estacoes.rds")[order(estacao)]
  d[, id := match(estacao, todas$estacao)]            # índice fixo do embedding
  setorder(d, valid_date, estacao)
  list(d = d, f_defas = f_defas, n_estacoes = nrow(todas),
       dist = distancias_km(todas))                     # linhas/colunas na ordem de id
}

# Climatologia da estação (média e dp da Tmax em +-15 dias do dia do ano), estimada só
# no treino da dobra, e a anomalia da média do ensemble em relação a ela (exp6)
adicionar_clim <- function(sp) {
  tabs <- lapply(split(sp$tr, by = "estacao"), function(a) {
    doy <- as.integer(format(a$valid_date, "%j"))
    t(vapply(1:366, function(k) {
      y <- a$tmax[pmin(abs(doy - k), 366 - abs(doy - k)) <= 15]
      c(mean(y), stats::sd(y))
    }, numeric(2)))
  })
  aplicar <- function(x) {
    x <- copy(x); doy <- as.integer(format(x$valid_date, "%j"))
    m <- s <- rep(NA_real_, nrow(x))
    for (e in names(tabs)) {
      i <- which(x$estacao == e); m[i] <- tabs[[e]][doy[i], 1]; s[i] <- tabs[[e]][doy[i], 2]
    }
    x[, `:=`(clim_media = m, clim_dp = s, ens_anom = ens_media - m)]
  }
  sp$tr <- aplicar(sp$tr); sp$va <- aplicar(sp$va); sp$te <- aplicar(sp$te)
  sp
}

dividir_dobra <- function(d, dobra) {
  k <- DOBRAS[[dobra]]
  tr <- d[valid_date <= k[1]]
  ok <- tr[, .N, by = estacao][N >= MIN_TREINO, estacao]
  list(tr = tr[estacao %in% ok],
       va = d[valid_date > k[1] + GAP_DIAS & valid_date <= k[2] & estacao %in% ok],
       te = d[valid_date > k[2] + GAP_DIAS & valid_date <= k[3] & estacao %in% ok],
       estacoes = sort(ok))
}

# Faltas: média do treino + um indicador por grupo de preditores.
# Devolve os três blocos imputados e a lista final de preditores.
imputar <- function(sp, feats, grupos) {
  medias <- vapply(feats, function(v) mean(sp$tr[[v]], na.rm = TRUE), numeric(1))
  medias[!is.finite(medias)] <- 0
  grupos <- Filter(function(g) length(intersect(g, feats)) > 0, grupos)
  aplicar <- function(x) {
    x <- copy(x)
    for (g in names(grupos)) {
      cols <- intersect(grupos[[g]], feats)
      set(x, j = g, value = as.numeric(!complete.cases(x[, cols, with = FALSE])))
    }
    for (v in feats) set(x, i = which(is.na(x[[v]])), j = v, value = medias[[v]])
    x
  }
  list(tr = aplicar(sp$tr), va = aplicar(sp$va), te = aplicar(sp$te), feats = c(feats, names(grupos)))
}

# Limiares sazonais q90/q95/q99 POR ESTAÇÃO, estimados no treino (+ 35 °C absoluto)
limiares_estacao <- function(tr, x) {
  ps <- c(q90 = 0.90, q95 = 0.95, q99 = 0.99)
  out <- lapply(ps, function(p) rep(NA_real_, nrow(x)))
  for (s in unique(x$estacao)) {
    a <- tr[estacao == s]; i <- which(x$estacao == s)
    for (nm in names(ps)) out[[nm]][i] <- limiar_sazonal(a$valid_date, a$tmax, ps[[nm]])(x$valid_date[i])
  }
  out$abs35 <- rep(35, nrow(x))
  out
}

# Avalia uma previsão em blocos de linhas (a grade de 999 quantis de todos os
# casos não cabe na memória). prever(i) devolve a grade das linhas i de x.
avaliar_em_blocos <- function(prever, x, lim, bloco = 5000L) {
  partes <- lapply(split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / bloco)), function(i) {
    q <- prever(i)
    stopifnot(nrow(q) == length(i), ncol(q) == length(TAU_GRADE), !anyNA(q))
    e <- escores_por_caso(q, x$tmax[i], lapply(lim, `[`, i))
    list(e = cbind(data.table(estacao = x$estacao[i], data = x$valid_date[i]), e),
         q = round(q[, round(TAU_FINO * 1000)], 2))       # TAU_GRADE[k] = k / 1000
  })
  list(escores = rbindlist(lapply(partes, `[[`, "e")), quantis = do.call(rbind, lapply(partes, `[[`, "q")))
}

# --- modelos -------------------------------------------------------------------------------
# Cada um recebe os blocos imputados (sp) e devolve list(va = prever, te = prever, info)
# com prever(i) -> grade de quantis das linhas i.

por_estacao <- function(fits, x, f) function(i) {
  out <- matrix(NA_real_, length(i), length(TAU_GRADE))
  xi <- x[i]
  for (s in unique(xi$estacao)) { j <- which(xi$estacao == s); out[j, ] <- f(fits[[s]], xi[j]) }
  out
}

# Climatologia da estação: quantis empíricos do treino em +-15 dias do dia do ano
modelo_clim <- function(sp) {
  fits <- lapply(split(sp$tr, by = "estacao"), function(a) clim_ajustar(a$valid_date, a$tmax))
  f <- function(fit, x) predict(fit, x$valid_date)
  list(va = por_estacao(fits, sp$va, f), te = por_estacao(fits, sp$te, f), info = list())
}

# EMOS local (Gneiting et al. 2005), um ajuste por estação por CRPS mínimo:
#   mu = X b (X = xcols), log sigma = g0 + g1 log(dp do ensemble)
modelo_emos_loc <- function(sp, xcols) {
  Z <- function(x) cbind(log(x$ens_dp + 0.05))
  X <- function(x) as.matrix(x[, xcols, with = FALSE])
  fits <- lapply(split(sp$tr, by = "estacao"), function(a) emos_ajustar(X(a), Z(a), a$tmax))
  f <- function(fit, x) predict(fit, X(x), Z(x))
  conv <- vapply(fits, `[[`, 0L, "conv")
  list(va = por_estacao(fits, sp$va, f), te = por_estacao(fits, sp$te, f),
       info = list(nao_convergiu = names(conv)[conv != 0]))
}

# Quantile regression forest global (Meinshausen 2006; Taillardat et al. 2016):
# todas as estações juntas, identificadas por lat/lon/altitude
modelo_qrf <- function(sp, feats, num.trees, threads) {
  fit <- ranger::ranger(x = as.data.frame(sp$tr[, feats, with = FALSE]), y = sp$tr$tmax, quantreg = TRUE,
                        num.trees = num.trees, min.node.size = 10, seed = 1, num.threads = threads)
  f <- function(x) function(i) qrf_quantis(fit, as.matrix(x[i, feats, with = FALSE]))
  list(va = f(sp$va), te = f(sp$te), info = list(oob_mse = fit$prediction.error))
}

# Rede distribucional com embedding de estação, conjunto de n_redes (Vincentização)
modelo_rede <- function(sp, feats, cabeca, n_estacoes, n_redes, epocas, seed0 = 100, ...) {
  padr <- padronizador(sp$tr, feats)
  conj <- treinar_conjunto(Xtr = padr(sp$tr), ytr = sp$tr$tmax, Xva = padr(sp$va), yva = sp$va$tmax,
                           id_tr = sp$tr$id, id_va = sp$va$id, n_estacoes = n_estacoes,
                           cabeca = cabeca, n_redes = n_redes, epocas = epocas, seed0 = seed0, ...)
  f <- function(x) { X <- padr(x); function(i) prever_conjunto(conj, X[i, , drop = FALSE], x$id[i]) }
  info <- list(melhor_epoca = vapply(conj, `[[`, 0, "melhor_epoca"),
               melhor_val = vapply(conj, `[[`, 0, "melhor_val"))
  if (cabeca == "bgp") {
    k <- seq_len(min(2000L, nrow(sp$te)))
    info$xi <- mean(xi_conjunto(conj, padr(sp$te)[k, , drop = FALSE], sp$te$id[k]))
  }
  list(va = f(sp$va), te = f(sp$te), info = info)
}

# --- exp6: ajuste em duas etapas ------------------------------------------------------------
# (1) No treino, com a validação para escolher o que precisa ser escolhido (épocas,
#     nº de iterações/árvores) -> prevê a VALIDAÇÃO (usada só para os pesos da combinação).
# (2) De novo em treino + validação, com essas escolhas fixas -> prevê o TESTE.
# O teste continua intocado; o modelo final só ganha os anos mais próximos do teste.
# Imputação, padronização da etapa 1 e limiares continuam vindo só do treino.

com_tv <- function(sp) { sp$tr <- rbind(sp$tr, sp$va); sp }

# para modelos sem nada a escolher na validação
duas_etapas <- function(f, sp) {
  a <- f(sp); b <- f(com_tv(sp))
  list(va = a$va, te = b$te, info = list(etapa1 = a$info, final = b$info))
}

.crps_norm <- function(y, mu, s) {
  z <- (y - mu) / s
  s * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
}

# EMOS com boosting (Messner, Mayr & Zeileis 2017; o "EMOS-loc-bst" de Rasp & Lerch 2018):
# mu e log sigma lineares em TODOS os preditores, por boosting componente a componente do
# CRPS normal. A cada iteração, para cada parâmetro, o preditor que melhor ajusta (MQO) o
# gradiente negativo anda nu nessa direção; os demais coeficientes ficam em zero, o que
# faz a seleção de preditores. Parte do EMOS simples (mu ~ ens_media, sigma constante).
emos_bst_ajustar <- function(X, y, m, nu = 0.1, Xva = NULL, yva = NULL) {
  usar <- apply(X, 2, stats::sd) > 1e-8                   # constantes na estação (lat, lon, ...)
  cen <- colMeans(X[, usar, drop = FALSE]); esc <- apply(X[, usar, drop = FALSE], 2, stats::sd)
  Z <- function(A) cbind(1, sweep(sweep(A[, usar, drop = FALSE], 2, cen), 2, esc, "/"))
  Xs <- Z(X); ss <- colSums(Xs^2)
  base <- intersect(c("ens_media", "t2m_media"), colnames(X)[usar])[1]   # média do ensemble do alvo
  if (is.na(base)) stop("emos_bst: falta a média do ensemble (ens_media ou t2m_media) entre os preditores")
  j0 <- 1 + match(base, colnames(X)[usar])
  b0 <- stats::coef(stats::lm.fit(Xs[, c(1, j0)], y))
  b <- numeric(ncol(Xs)); b[c(1, j0)] <- b0
  g <- numeric(ncol(Xs)); g[1] <- log(stats::sd(y - Xs %*% b))
  em <- drop(Xs %*% b); es <- drop(Xs %*% g)
  if (!is.null(Xva)) { Xv <- Z(Xva); vm <- drop(Xv %*% b); vs <- drop(Xv %*% g); crps_va <- numeric(m) }
  melhor <- function(u) { c <- drop(crossprod(Xs, u)) / ss; j <- which.max(c^2 * ss); c(j, nu * c[j]) }
  for (k in seq_len(m)) {
    s <- exp(es); z <- (y - em) / s
    pm <- melhor(2 * stats::pnorm(z) - 1)                          # -dCRPS/dmu
    ps <- melhor(-s * (2 * stats::dnorm(z) - 1 / sqrt(pi)))        # -dCRPS/dlog sigma
    b[pm[1]] <- b[pm[1]] + pm[2]; g[ps[1]] <- g[ps[1]] + ps[2]
    em <- em + pm[2] * Xs[, pm[1]]; es <- es + ps[2] * Xs[, ps[1]]
    if (!is.null(Xva)) {
      vm <- vm + pm[2] * Xv[, pm[1]]; vs <- vs + ps[2] * Xv[, ps[1]]
      crps_va[k] <- mean(.crps_norm(yva, vm, exp(vs)))
    }
  }
  structure(list(b = b, g = g, Z = Z, crps_va = if (!is.null(Xva)) crps_va), class = "emos_bst")
}
predict.emos_bst <- function(object, X, ...) {
  A <- object$Z(X)
  outer(drop(exp(A %*% object$g)), stats::qnorm(TAU_GRADE)) + drop(A %*% object$b)
}

# nº de iterações escolhido no CRPS de validação de cada estação. No exp6 (m_max = 1000,
# nu = 0,1) a escolha bateu no teto em quase todas as estações: o padrão agora anda mais
# rápido (nu = 0,3) e vai mais longe (m_max = 2000), ~6x o caminho anterior.
modelo_emos_bst <- function(sp, feats, m_max = 2000, nu = 0.3) {
  X <- function(x) as.matrix(x[, feats, with = FALSE])
  tv <- com_tv(sp)$tr
  m <- fits <- finais <- list()
  for (s in unique(sp$tr$estacao)) {
    a <- sp$tr[estacao == s]; v <- sp$va[estacao == s]
    if (nrow(v)) {
      cv <- emos_bst_ajustar(X(a), a$tmax, m_max, nu, Xva = X(v), yva = v$tmax)$crps_va
      m[[s]] <- which.min(cv)
    } else m[[s]] <- m_max %/% 2
    fits[[s]] <- emos_bst_ajustar(X(a), a$tmax, m[[s]], nu)
    finais[[s]] <- emos_bst_ajustar(X(tv[estacao == s]), tv[estacao == s]$tmax, m[[s]], nu)
  }
  f <- function(fit, x) predict(fit, X(x))
  list(va = por_estacao(fits, sp$va, f), te = por_estacao(finais, sp$te, f),
       info = list(iteracoes = unlist(m)))
}

# Boosting de árvores (xgboost) global, normal heterocedástica: árvores para a média
# (erro quadrático; nº de árvores pela validação) e outras para a variância, ajustadas
# nos quadrados dos resíduos FORA DA AMOSTRA do treino (4 blocos de tempo), com perda
# gama e ligação log, i.e. E[r^2 | x] = sigma^2(x).
GBM_PAR <- list(eta = 0.05, max_depth = 6, subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 20)

.xgb <- function(X, y, objetivo, n, threads, Xva = NULL, yva = NULL) {
  par <- c(GBM_PAR, list(objective = objetivo, nthread = threads))
  dtr <- xgboost::xgb.DMatrix(X, label = y)
  if (is.null(Xva)) return(xgboost::xgb.train(par, dtr, nrounds = n, verbose = 0))
  xgboost::xgb.train(par, dtr, nrounds = n, watchlist = list(va = xgboost::xgb.DMatrix(Xva, label = yva)),
                     early_stopping_rounds = 50, verbose = 0)
}
.xgb_prever <- function(fit, X, n) predict(fit, X, iterationrange = c(1, n + 1))

# resíduos fora da amostra: cada bloco de tempo previsto pelo modelo da média ajustado nos outros
.residuos_oof <- function(X, y, datas, n, threads, blocos = 4) {
  b <- cut(rank(as.numeric(datas), ties.method = "first"), blocos, labels = FALSE)
  r <- numeric(length(y))
  for (k in seq_len(blocos)) {
    i <- b == k
    r[i] <- y[i] - .xgb_prever(.xgb(X[!i, , drop = FALSE], y[!i], "reg:squarederror", n, threads), X[i, , drop = FALSE], n)
  }
  r
}

modelo_gbm <- function(sp, feats, threads, n_max = 3000) {
  X <- function(x) as.matrix(x[, feats, with = FALSE])
  Xtr <- X(sp$tr); Xva <- X(sp$va)
  # etapa 1: escolhas na validação
  media <- .xgb(Xtr, sp$tr$tmax, "reg:squarederror", n_max, threads, Xva, sp$va$tmax)
  n_m <- media$best_iteration
  mu_va <- .xgb_prever(media, Xva, n_m)
  r2 <- .residuos_oof(Xtr, sp$tr$tmax, sp$tr$valid_date, n_m, threads)^2 + 1e-3
  vari <- .xgb(Xtr, r2, "reg:gamma", n_max, threads, Xva, (sp$va$tmax - mu_va)^2 + 1e-3)
  n_v <- vari$best_iteration
  s_va <- sqrt(.xgb_prever(vari, Xva, n_v))
  # etapa 2: treino + validação com n_m e n_v fixos
  tv <- com_tv(sp)$tr; Xtv <- X(tv)
  media2 <- .xgb(Xtv, tv$tmax, "reg:squarederror", n_m, threads)
  r2 <- .residuos_oof(Xtv, tv$tmax, tv$valid_date, n_m, threads)^2 + 1e-3
  vari2 <- .xgb(Xtv, r2, "reg:gamma", n_v, threads)
  grade <- function(mu, s) outer(s, stats::qnorm(TAU_GRADE)) + mu
  te <- function(i) { Xi <- X(sp$te[i]); grade(.xgb_prever(media2, Xi, n_m), sqrt(.xgb_prever(vari2, Xi, n_v))) }
  list(va = function(i) grade(mu_va[i], s_va[i]), te = te, info = list(arvores_media = n_m, arvores_var = n_v))
}

# Receita (hiperparâmetros) escolhida na validação DESTA dobra: cada candidata treina
# n_busca sementes só no treino (a receita pode trazer o seu n_busca) e fica a de menor
# perda de validação média (a perda é o CRPS do alvo padronizado com média/dp do treino,
# comparável entre candidatas).
.escolher <- function(receitas, treinar, n_busca = 2) {
  nb <- vapply(receitas, function(r) if (is.null(r$n_busca)) n_busca else r$n_busca, 0)
  receitas <- lapply(receitas, function(r) { r$n_busca <- NULL; r })
  if (length(receitas) == 1) return(list(receita = receitas[[1]], nome = names(receitas)[1], val = NA))
  val <- vapply(seq_along(receitas), function(j)
    mean(vapply(seq_len(nb[j]), function(k) treinar(receitas[[j]], 900 + k)$melhor_val, 0)), 0)
  names(val) <- names(receitas)
  list(receita = receitas[[which.min(val)]], nome = names(receitas)[which.min(val)], val = val)
}

# Rede (MLP) em duas etapas: cada semente treina com early stopping na validação e
# depois é treinada de novo, do zero, em treino + validação pelo mesmo nº de épocas.
# receitas: lista nomeada de listas de argumentos (treinar_rede / rede_distribucional).
modelo_rede2 <- function(sp, feats, cabeca, n_estacoes, n_redes, epocas, receitas, seed0 = 100) {
  padr <- padronizador(sp$tr, feats); Xtr <- padr(sp$tr); Xva <- padr(sp$va)
  treinar <- function(r, seed, X = Xtr, y = sp$tr$tmax, id = sp$tr$id, Xv = Xva, ep = epocas)
    do.call(treinar_rede, c(list(Xtr = X, ytr = y, Xva = Xv, yva = sp$va$tmax, id_tr = id, id_va = sp$va$id,
                                 n_estacoes = n_estacoes, cabeca = cabeca, epocas = ep, seed = seed), r))
  esc <- .escolher(receitas, treinar)
  conj <- lapply(seq_len(n_redes), function(k) treinar(esc$receita, seed0 + k))
  tv <- com_tv(sp)$tr; padr2 <- padronizador(tv, feats); Xtv <- padr2(tv)
  final <- lapply(seq_along(conj), function(k)
    treinar(esc$receita, seed0 + k, Xtv, tv$tmax, tv$id, NULL, max(3L, conj[[k]]$melhor_epoca)))
  f <- function(cj, p, x) { X <- p(x); function(i) prever_conjunto(cj, X[i, , drop = FALSE], x$id[i]) }
  info <- list(receita = esc$nome, val_receitas = esc$val,
               melhor_epoca = vapply(conj, `[[`, 0, "melhor_epoca"),
               melhor_val = vapply(conj, `[[`, 0, "melhor_val"))
  if (cabeca == "bgp") {
    k <- seq_len(min(2000L, nrow(sp$te)))
    info$xi <- mean(xi_conjunto(final, padr2(sp$te)[k, , drop = FALSE], sp$te$id[k]))
  }
  list(va = f(conj, padr, sp$va), te = f(final, padr2, sp$te), info = info)
}

# Rede com atenção entre estações (modelos_nn.R), também em duas etapas. A previsão das
# linhas i usa como contexto TODAS as linhas do mesmo bloco nos dias dessas linhas.
modelo_atencao <- function(sp, feats, cabeca, n_estacoes, dist, n_redes, epocas, receitas, seed0 = 200) {
  padr <- padronizador(sp$tr, feats)
  dias <- function(x, p) dias_tokens(p(x), x$id, x$valid_date, n_estacoes, x$tmax)
  Ttr <- dias(sp$tr, padr); Tva <- dias(sp$va, padr)
  treinar <- function(r, seed, T = Ttr, Tv = Tva, ep = epocas)
    do.call(treinar_atencao, c(list(T, Tv, cabeca, n_estacoes, dist, seed = seed, epocas = ep), r))
  esc <- .escolher(receitas, treinar)
  conj <- lapply(seq_len(n_redes), function(k) treinar(esc$receita, seed0 + k))
  tv <- com_tv(sp)$tr; padr2 <- padronizador(tv, feats)
  Ttv <- dias(tv, padr2); rm(Ttr, Tva)
  final <- lapply(seq_along(conj), function(k)
    treinar(esc$receita, seed0 + k, Ttv, NULL, max(3L, conj[[k]]$melhor_epoca)))
  f <- function(cj, p, x) {
    X <- p(x)
    function(i) {
      j <- which(x$valid_date %in% unique(x$valid_date[i]))     # contexto: o dia inteiro
      q <- prever_atencao(cj, dias_tokens(X[j, , drop = FALSE], x$id[j], x$valid_date[j], n_estacoes))
      q[match(i, j), , drop = FALSE]
    }
  }
  list(va = f(conj, padr, sp$va), te = f(final, padr2, sp$te),
       info = list(receita = esc$nome, val_receitas = esc$val,
                   melhor_epoca = vapply(conj, `[[`, 0, "melhor_epoca"),
                   melhor_val = vapply(conj, `[[`, 0, "melhor_val"),
                   beta_dist = lapply(final, `[[`, "beta")))
}
