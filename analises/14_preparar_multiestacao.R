# Congela a base multi-estação (fase 1.2 do PLANO.md): 43 estações automáticas
# do INMET no RS x GEFSv12 (caixa do RS), 2000–2019, horizontes 1, 3 e 7.
#
# Uso: Rscript analises/14_preparar_multiestacao.R
#      TESTE=1 -> só A801 e A803 em 2010 (verificação rápida, saída *_teste.rds)
#      FONTE=oper -> GEFS OPERACIONAL 2020-09-24 a 2025-12-31 (script 39) e INMET 2020–25, para o exp14;
#                    saídas multi_obs_dia_oper, multi_obs_init_oper e multi_gefs_oper (mesmas colunas)
#
# Saídas em dados/processados/:
#   multi_estacoes.rds  código, nome, lat, lon, altitude
#   multi_obs_dia.rds   Tmax do DIA LOCAL por estação (mesma definição da A801)
#   multi_obs_init.rds  o que se sabe na rodada 00 UTC: Tmax parcial da véspera,
#                       Tmax de 2 e 3 dias antes, estado da estação (variaveis_estacao)
#   multi_gefs.rds      por estação, rodada e H: os ~29 campos do membro de controle
#                       interpolados (bilinear), Tmax e chuva do dia previsto, e
#                       média/dp/mín/máx do ensemble de 5 membros (Tmax, T2m, chuva)
source("R/dados.R")
TESTE <- Sys.getenv("TESTE") == "1"
OPER <- Sys.getenv("FONTE") == "oper"
suf <- if (TESTE) "_teste" else if (OPER) "_oper" else ""
t0 <- Sys.time()
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

est <- as.data.table(readRDS("dados/brutos/inmet_rs/estacoes.rds"))[
  , .(estacao = station_code, nome = station_name, lat = latitude, lon = longitude, altitude = elevation)]
anos <- if (OPER) 2020:2025 else 2000:2019
if (TESTE) { est <- est[estacao %in% c("A801", "A803")]; anos <- 2010 }

# --- 1. observações ---------------------------------------------------------------
horario <- rbindlist(lapply(anos, function(a) {
  d <- as.data.table(readRDS(sprintf("dados/brutos/inmet_rs/%d.rds", a)))
  d[station_code %in% est$estacao]
}))
obs_dia <- list(); obs_init <- list()
for (s in est$estacao) {
  h <- inmet_para_horario(horario[station_code == s])
  if (!nrow(h)) next
  diaria <- tmax_diaria(qc_inmet(copy(h)))
  inits <- seq(min(diaria$dia) + 3, max(diaria$dia) + 1, by = "day")
  oi <- data.table(estacao = s, init_date = inits,
                   obs_vespera = tmax_parcial_vespera(qc_inmet(copy(h)), inits))
  oi[, obs_lag2 := diaria$tmax[match(init_date - 2, diaria$dia)]]
  oi[, obs_lag3 := diaria$tmax[match(init_date - 3, diaria$dia)]]
  # preditores da estação: QC do rmet + faixas próprias (desde 26/10; ver script 35)
  hp <- inmet_para_horario(qc_rmet(horario[station_code == s]))
  oi <- merge(oi, variaveis_estacao(qc_inmet_completo(hp), inits), by = "init_date")
  obs_dia[[s]] <- diaria[, .(estacao = s, dia, tmax, n_horas)]
  obs_init[[s]] <- oi
  msg(s, ":", sum(!is.na(diaria$tmax)), "dias com Tmax")
}
obs_dia <- rbindlist(obs_dia); obs_init <- rbindlist(obs_init, use.names = TRUE)

# --- 2. GEFS nas estações ------------------------------------------------------------
dias <- seq(if (OPER) as.Date("2020-09-24") else as.Date(sprintf("%d-01-01", min(anos))),
            as.Date(sprintf("%d-12-31", max(anos))), by = "day")
arq <- function(m, d) sprintf("dados/brutos/%s/%s/%s/%s.rds", if (OPER) "gefs_oper" else "gefs_aux", m,
                              format(d, "%Y"), format(d, "%Y%m%d"))
g0 <- readRDS(arq("c00", dias[1]))
W <- pesos_bilineares(g0$grade$lon, g0$grade$lat, est$lon, est$lat)

# por membro: Tmax e chuva do dia previsto (máx./soma das 4 janelas) e T2m às 18 UTC;
# para o controle, todos os campos
ler_membro <- function(m, todos = FALSE) {
  out <- vector("list", length(dias))
  for (k in seq_along(dias)) {
    f <- arq(m, dias[k]); if (!file.exists(f)) next
    x <- readRDS(f)
    out[[k]] <- rbindlist(lapply(c(1, 3, 7), function(H) {
      v <- interpolar_campos(x[[paste0("H", H)]], W)
      tx <- v[, grep("^tmax2m_", colnames(v)), drop = FALSE]
      ch <- v[, grep("^chuva_", colnames(v)), drop = FALSE]
      base <- data.table(estacao = est$estacao, init_date = dias[k], horizon = H,
                         tmax_dia = apply(tx, 1, max), chuva_dia = rowSums(ch), t2m_18 = v[, "t2m_h18"])
      if (todos) cbind(base, as.data.table(v)) else base
    }))
    if (k %% 1000 == 0) msg(m, ":", k, "de", length(dias))
  }
  rbindlist(out)
}
ctrl <- ler_membro("c00", todos = TRUE)
membros <- rbindlist(c(list(ctrl[, .(estacao, init_date, horizon, tmax_dia, chuva_dia, t2m_18, membro = "c00")]),
                       lapply(c("p01", "p02", "p03", "p04"), function(m) ler_membro(m)[, membro := m])))
ens <- membros[, .(ens_media = mean(tmax_dia), ens_dp = sd(tmax_dia), ens_min = min(tmax_dia),
                   ens_max = max(tmax_dia), n_membros = sum(!is.na(tmax_dia)),
                   t2m_ens_media = mean(t2m_18), t2m_ens_dp = sd(t2m_18),
                   chuva_ens_media = mean(chuva_dia), chuva_ens_dp = sd(chuva_dia),
                   chuva_p1mm = mean(chuva_dia > 1)),
               by = .(estacao, init_date, horizon)]
gefs <- merge(ctrl, ens, by = c("estacao", "init_date", "horizon"))
gefs[, valid_date := init_date + horizon]
setnames(gefs, c("tmax_dia", "chuva_dia"), c("ctrl_tmax_dia", "ctrl_chuva_dia"))

# conferência: ensemble de Tmax na A801 vs a tabela antiga (mesma fonte)
if (!OPER && "A801" %in% est$estacao) {
  velha <- readRDS("dados/processados/tabela_tmax_gefs_a801.rds")[, .(init_date, horizon, ens_media_velha = ens_media)]
  cmp <- merge(gefs[estacao == "A801", .(init_date, horizon, ens_media)], velha, by = c("init_date", "horizon"))
  msg(sprintf("A801: ens_media interpolada vs ponto da tabela antiga: dif média %.3f °C, cor %.4f (n = %d)",
              mean(cmp$ens_media - cmp$ens_media_velha, na.rm = TRUE),
              cor(cmp$ens_media, cmp$ens_media_velha, use = "complete.obs"), nrow(cmp)))
}

dir.create("dados/processados", showWarnings = FALSE)
if (!OPER) saveRDS(est, sprintf("dados/processados/multi_estacoes%s.rds", suf))
saveRDS(obs_dia, sprintf("dados/processados/multi_obs_dia%s.rds", suf))
saveRDS(obs_init, sprintf("dados/processados/multi_obs_init%s.rds", suf))
saveRDS(gefs, sprintf("dados/processados/multi_gefs%s.rds", suf))
msg(sprintf("pronto: %d estações | gefs %d linhas x %d colunas | %.1f min", nrow(est), nrow(gefs), ncol(gefs),
            as.numeric(Sys.time() - t0, units = "mins")))
