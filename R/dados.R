# Leitura e preparação dos dados: INMET A801 (observação) e GEFSv12 (previsão).

library(data.table)

# --- INMET: arquivos anuais horários (UTC) -> data.table horário ---------------
ler_inmet <- function(pasta) {
  arqs <- list.files(pasta, pattern = "\\.CSV$", full.names = TRUE, ignore.case = TRUE)
  rbindlist(lapply(arqs, function(a) {
    d <- fread(a, skip = 8, sep = ";", dec = ",", encoding = "Latin-1",
               na.strings = c("-9999", "-9999,0", ""), header = TRUE, fill = TRUE)
    nm <- names(d)
    col <- function(padrao) nm[grepl(padrao, nm, ignore.case = TRUE)][1]
    data <- gsub("/", "-", d[[1]])
    hora <- substr(gsub("[^0-9]", "", d[[2]]), 1, 2)          # "00:00" ou "0000 UTC"
    num <- function(x) as.numeric(gsub(",", ".", as.character(x)))
    data.table(
      utc  = as.POSIXct(paste(data, hora), format = "%Y-%m-%d %H", tz = "UTC"),
      tmax_h = num(d[[col("TEMPERATURA M.XIMA NA HORA")]]),
      temp_h = num(d[[col("BULBO SECO")]])
    )
  }))[!is.na(utc)][order(utc)]
}

# Controle de qualidade simples: faixa física e coerência com a temperatura horária
qc_inmet <- function(h) {
  h[tmax_h < -5 | tmax_h > 45, tmax_h := NA]
  h[!is.na(temp_h) & !is.na(tmax_h) & (tmax_h < temp_h - 1 | tmax_h > temp_h + 8), tmax_h := NA]
  h[]
}

# Tmax diária no DIA LOCAL (America/Sao_Paulo, com horário de verão até 2019).
# Cada registro horário é o máximo da hora anterior, então o registro das HH:00
# cobre (HH-1, HH]; atribuímos ao dia local do início do intervalo.
tmax_diaria <- function(h, min_horas = 20) {
  h <- copy(h)[, dia := as.Date(format(utc - 1, tz = "America/Sao_Paulo"))]
  h[, .(tmax = if (sum(!is.na(tmax_h)) >= min_horas) max(tmax_h, na.rm = TRUE) else NA_real_,
        n_horas = sum(!is.na(tmax_h))), by = dia]
}

# Tmax "conhecida" do dia local anterior à inicialização 00 UTC de init_date:
# só registros com utc <= init_date 00:00 UTC (sem olhar o futuro).
tmax_parcial_vespera <- function(h, init_dates, min_horas = 15) {
  h <- copy(h)[, dia := as.Date(format(utc - 1, tz = "America/Sao_Paulo"))]
  setkey(h, dia)
  vapply(init_dates, function(d0) {
    corte <- as.POSIXct(paste(d0, "00:00"), tz = "UTC")
    x <- h[.(d0 - 1)][utc <= corte, tmax_h]
    if (sum(!is.na(x)) >= min_horas) max(x, na.rm = TRUE) else NA_real_
  }, numeric(1))
}

# --- INMET completo: todas as variáveis horárias usadas como preditores ---------
# (padrões sem acento: o cabeçalho vem em Latin-1)
ler_inmet_completo <- function(pasta) {
  arqs <- list.files(pasta, pattern = "\\.CSV$", full.names = TRUE, ignore.case = TRUE)
  rbindlist(lapply(arqs, function(a) {
    d <- fread(a, skip = 8, sep = ";", dec = ",", encoding = "Latin-1",
               na.strings = c("-9999", "-9999,0", ""), header = TRUE, fill = TRUE)
    nm <- names(d)
    col <- function(padrao) {
      k <- nm[grepl(padrao, nm, ignore.case = TRUE)][1]
      if (is.na(k)) return(rep(NA_real_, nrow(d)))
      as.numeric(gsub(",", ".", as.character(d[[k]])))
    }
    data <- gsub("/", "-", d[[1]])
    hora <- substr(gsub("[^0-9]", "", d[[2]]), 1, 2)
    data.table(
      utc    = as.POSIXct(paste(data, hora), format = "%Y-%m-%d %H", tz = "UTC"),
      chuva  = col("PRECIPITA"),
      pres   = col("PRESSAO ATMOSFERICA AO NIVEL"),
      rad    = col("RADIACAO GLOBAL"),
      temp   = col("BULBO SECO"),
      orv    = col("PONTO DE ORVALHO"),
      ur     = col("UMIDADE RELATIVA DO AR"),
      vdir   = col("VENTO, DIRE"),
      vvel   = col("VENTO, VELOCIDADE")
    )
  }))[!is.na(utc)][order(utc)]
}

qc_inmet_completo <- function(h) {
  h[!(chuva >= 0 & chuva <= 150), chuva := NA]
  h[!(pres >= 950 & pres <= 1050), pres := NA]
  h[!(rad >= 0 & rad <= 5000), rad := NA]
  h[!(temp >= -5 & temp <= 45), temp := NA]
  h[!(orv >= -15 & orv <= 35) | orv > temp + 0.5, orv := NA]
  h[!(ur >= 1 & ur <= 100), ur := NA]
  h[!(vdir >= 0 & vdir <= 360), vdir := NA]
  h[!(vvel >= 0 & vvel <= 40), vvel := NA]
  h[]
}

# Estado da atmosfera na estação no instante da rodada (00 UTC de init_date,
# 21h local da véspera): só registros com utc <= init 00 UTC.
variaveis_estacao <- function(h, init_dates) {
  grade <- data.table(utc = seq(min(h$utc), max(h$utc), by = "hour"))
  h <- h[grade, on = "utc"]
  idx <- match(as.POSIXct(paste(init_dates, "00:00"), tz = "UTC"), h$utc)
  ultimo <- function(x, i, janela = 2) {        # valor em i, ou o mais recente até i-janela
    if (is.na(i)) return(NA_real_)
    v <- x[max(1, i - janela):i]; v <- v[!is.na(v)]
    if (length(v)) v[length(v)] else NA_real_
  }
  soma24 <- function(x, i, min_n) {
    if (is.na(i) || i < 24) return(NA_real_)
    v <- x[(i - 23):i]
    if (sum(!is.na(v)) >= min_n) sum(v, na.rm = TRUE) else NA_real_
  }
  u <- -h$vvel * sin(h$vdir * pi / 180); v <- -h$vvel * cos(h$vdir * pi / 180)
  media3 <- function(x, i) if (is.na(i) || i < 3) NA_real_ else mean(x[(i - 2):i], na.rm = TRUE)
  # radiação: horas noturnas costumam vir como NA -> contam como 0
  rad0 <- fifelse(is.na(h$rad) & !is.na(h$temp), 0, h$rad)
  out <- data.table(
    init_date = init_dates,
    est_t00   = vapply(idx, function(i) ultimo(h$temp, i), numeric(1)),
    est_td00  = vapply(idx, function(i) ultimo(h$orv, i), numeric(1)),
    est_ur00  = vapply(idx, function(i) ultimo(h$ur, i), numeric(1)),
    est_p00   = vapply(idx, function(i) ultimo(h$pres, i), numeric(1)),
    est_p24   = vapply(idx, function(i) if (is.na(i) || i <= 24) NA_real_ else ultimo(h$pres, i - 24), numeric(1)),
    est_u     = vapply(idx, function(i) media3(u, i), numeric(1)),
    est_v     = vapply(idx, function(i) media3(v, i), numeric(1)),
    est_chuva24 = vapply(idx, function(i) soma24(h$chuva, i, 20), numeric(1)),
    est_rad24 = vapply(idx, function(i) soma24(rad0, i, 20), numeric(1)),
    est_tmin24 = vapply(idx, function(i) {
      if (is.na(i) || i < 24) return(NA_real_)
      x <- h$temp[(i - 23):i]; if (sum(!is.na(x)) >= 18) min(x, na.rm = TRUE) else NA_real_
    }, numeric(1))
  )
  out[is.nan(est_u), est_u := NA][is.nan(est_v), est_v := NA]
  out[, est_dp24 := est_p00 - est_p24][, est_p24 := NULL]
  out[, est_chuva24 := log1p(est_chuva24)]
  out[]
}

# --- Multi-estação (INMET via rmet + GEFS na caixa do RS) ----------------------

# Dados horários do rmet -> nomes usados pelas funções acima (qc_inmet,
# tmax_diaria, tmax_parcial_vespera, qc_inmet_completo, variaveis_estacao)
inmet_para_horario <- function(d) {
  d <- as.data.table(d)
  d[, .(utc = datetime, tmax_h = temp_max_prev_c, temp_h = temp_dry_c,
        chuva = precip_mm, pres = pressure_station_hpa, rad = radiation_kjm2,
        temp = temp_dry_c, orv = temp_dew_c, ur = humid_rel_pct,
        vdir = wind_dir_deg, vvel = wind_speed_ms)][order(utc)]
}

# Faixas físicas dos campos do GEFS (fora disso: campo corrompido na NOAA -> NA)
FAIXAS_GEFS <- list(tmax2m = c(-15, 50), t2m = c(-15, 50), chuva = c(0, 400), q2m = c(0, 0.035),
                    u10 = c(-40, 40), v10 = c(-40, 40), pnmm = c(94000, 105000), nuvens = c(0, 100),
                    radsw = c(0, 1400), cape = c(0, 10000), pwat = c(0, 100), fsens = c(-600, 1200),
                    flat = c(-600, 1500), solo = c(0, 1), t850 = c(-25, 40), q850 = c(0, 0.03),
                    u850 = c(-60, 60), v850 = c(-60, 60), z500 = c(4500, 6200))

# Pesos de interpolação bilinear: matriz n_estações x n_células da grade
# (células na ordem de as.vector(m), m = matriz lat x lon do arquivo do GEFS)
pesos_bilineares <- function(lon_grade, lat_grade, lon, lat) {
  nl <- length(lat_grade); W <- matrix(0, length(lon), nl * length(lon_grade))
  for (s in seq_along(lon)) {
    j0 <- findInterval(lon[s], lon_grade)                      # lon crescente
    i0 <- findInterval(-lat[s], -lat_grade)                    # lat decrescente
    if (j0 < 1 || j0 >= length(lon_grade) || i0 < 1 || i0 >= nl) stop("estação fora da grade")
    fx <- (lon[s] - lon_grade[j0]) / (lon_grade[j0 + 1] - lon_grade[j0])
    fy <- (lat[s] - lat_grade[i0]) / (lat_grade[i0 + 1] - lat_grade[i0])
    for (a in 0:1) for (b in 0:1)
      W[s, (i0 + a) + (j0 + b - 1) * nl] <- (if (a) fy else 1 - fy) * (if (b) fx else 1 - fx)
  }
  W
}

# Interpola todos os campos de uma lista de matrizes (um H de um arquivo) nas
# estações; células NA (ex.: solo no mar) saem da média e os pesos são renormalizados
interpolar_campos <- function(campos, W) {
  M <- vapply(campos, as.vector, numeric(length(campos[[1]])))
  ok <- !is.na(M); M[!ok] <- 0
  v <- (W %*% M) / (W %*% ok)
  v[!is.finite(v)] <- NA
  for (nm in colnames(v)) {
    f <- FAIXAS_GEFS[[sub("_h[0-9]+$", "", nm)]]
    v[, nm][v[, nm] < f[1] | v[, nm] > f[2]] <- NA
  }
  v
}

# --- GEFS: formato longo (init, valid, member, horizon, forecast_tmax em K) -----
ler_gefs <- function(arquivo) {
  g <- fread(arquivo)
  g[, `:=`(init_date = as.Date(init_date), valid_date = as.Date(valid_date),
           tc = forecast_tmax - 273.15)]
  # caso anômalo conhecido (campo GRIB corrompido): descarta o ensemble inteiro do caso
  ruins <- unique(g[tc > 50 | tc < -10, .(init_date, horizon)])
  if (nrow(ruins)) g <- g[!ruins, on = .(init_date, horizon)]
  attr(g, "descartados") <- ruins
  g
}

# Tabela de modelagem: uma linha por (init_date, horizon)
montar_tabela <- function(g, diaria, horario) {
  w <- dcast(g, init_date + valid_date + horizon ~ member, value.var = "tc")
  mem <- setdiff(names(w), c("init_date", "valid_date", "horizon"))
  M <- as.matrix(w[, ..mem])
  w[, `:=`(ens_media = rowMeans(M), ens_dp = apply(M, 1, sd),
           ens_min = apply(M, 1, min), ens_max = apply(M, 1, max))]
  w <- merge(w, diaria[, .(valid_date = dia, tmax)], by = "valid_date", all.x = TRUE)
  inits <- sort(unique(w$init_date))
  lag <- data.table(init_date = inits, obs_vespera = tmax_parcial_vespera(horario, inits))
  # Tmax de dias locais completos antes da rodada (00 UTC = 21h local da véspera)
  lag[, obs_lag2 := diaria$tmax[match(init_date - 2, diaria$dia)]]
  lag[, obs_lag3 := diaria$tmax[match(init_date - 3, diaria$dia)]]
  # "erro de ontem": Tmax observada na véspera menos a previsão de 1 dia
  # (rodada de init_date - 2) para esse mesmo dia. Conhecido na hora da rodada.
  prev_ontem <- w[horizon == 1, .(init_date = init_date + 2, prev_vespera = ens_media)]
  lag <- merge(lag, prev_ontem, by = "init_date", all.x = TRUE)
  lag[, erro_ontem := obs_vespera - prev_vespera][, prev_vespera := NULL]
  w <- merge(w, lag, by = "init_date", all.x = TRUE)
  doy <- as.integer(format(w$valid_date, "%j"))
  w[, `:=`(doy_sen = sin(2 * pi * doy / 365.25), doy_cos = cos(2 * pi * doy / 365.25))]
  attr(w, "membros") <- mem
  w[order(horizon, init_date)]
}

# Blocos temporais pela data VÁLIDA, com intervalo (gap) de 7 dias entre blocos
atribuir_blocos <- function(w, fim_treino = as.Date("2011-12-31"),
                            fim_val = as.Date("2015-12-31"), gap = 7) {
  w[, bloco := fcase(
    valid_date <= fim_treino, "treino",
    valid_date > fim_treino + gap & valid_date <= fim_val, "validacao",
    valid_date > fim_val + gap, "teste",
    default = NA_character_)]
  w[]
}
