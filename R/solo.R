# Memória da superfície conhecida na rodada 00 UTC de init (exp13). Hipótese física: solo seco ->
# menos evapotranspiração -> mais calor sensível -> Tmax mais alta, e o GEFS conhece mal o estado do
# solo (o RS fica na bacia do Prata, um dos "hotspots" de acoplamento solo–atmosfera).
# Depende de R/multi.R (.media_movel, distancias_km).
#
#   H19  F_CHUVA_ANT  chuva antecedente observada nas estações do INMET. Como nos erros recentes, um dia
#        local D só se fecha às 03 UTC de D+1, então entram os dias D <= init - 2:
#          chuva_ant_7d/30d/90d  log1p da soma em [init-8, init-2] etc. (média dos dias válidos x n,
#                                NA se < 70% válidos)
#          dias_sem_chuva        dias desde o último dia com >= 1 mm (até 90)
#          chuva_ant_30d_viz/_reg  chuva_ant_30d médio das 4 estações mais próximas / de todas
#   H20  F_SOLO  umidade do solo de reanálise (MERRA-2, NASA POWER; script 40) defasada DEFAS_SOLO dias,
#        para respeitar a latência de um produto operacional:
#          gwet_sup, gwet_raiz   superfície (0–5 cm) e zona de raízes em init - DEFAS_SOLO
#          gwet_raiz_d30         variação da zona de raízes nos 30 dias anteriores

F_CHUVA_ANT <- c("chuva_ant_7d", "chuva_ant_30d", "chuva_ant_90d", "dias_sem_chuva", "chuva_ant_30d_viz", "chuva_ant_30d_reg")
F_SOLO <- c("gwet_sup", "gwet_raiz", "gwet_raiz_d30")
DEFAS_SOLO <- 7

# uma linha por (estação, init_date) para todas as datas de `inits`
memoria_superficie <- function(inits, obs_chuva = readRDS("dados/processados/multi_obs_chuva.rds"),
                               solo = readRDS("dados/processados/multi_solo.rds"),
                               est = readRDS("dados/processados/multi_estacoes.rds")) {
  dias <- seq(min(inits) - 100, max(inits), by = "day")
  ch <- obs_chuva[!is.na(chuva) & dia %in% dias]
  rol <- ch[, {
    x <- rep(NA_real_, length(dias)); x[match(dia, dias)] <- chuva
    soma <- function(n) log1p(.media_movel(x, n, ceiling(0.7 * n)) * n)
    molhado <- ifelse(!is.na(x) & x >= 1, seq_along(x), NA_integer_)
    ultimo <- cummax(ifelse(is.na(molhado), 0L, molhado))
    list(init_date = dias + 2, chuva_ant_7d = soma(7), chuva_ant_30d = soma(30), chuva_ant_90d = soma(90),
         dias_sem_chuva = ifelse(ultimo > 0, pmin(seq_along(x) - ultimo, 90), NA_real_))
  }, by = estacao]
  rol <- rol[init_date %in% inits]

  # vizinhas e estado (as observações das outras estações também já chegaram)
  W <- dcast(rol, init_date ~ estacao, value.var = "chuva_ant_30d")
  M <- as.matrix(W[, -1]); D <- distancias_km(est[match(colnames(M), estacao)])
  media <- function(A) { r <- rowMeans(A, na.rm = TRUE); r[!is.finite(r)] <- NA; r }
  reg <- media(M)
  viz <- rbindlist(lapply(seq_len(ncol(M)), function(j)
    data.table(estacao = colnames(M)[j], init_date = W$init_date,
               chuva_ant_30d_viz = media(M[, order(D[j, ])[2:5], drop = FALSE]), chuva_ant_30d_reg = reg)))
  rol <- merge(rol, viz, by = c("estacao", "init_date"), all = TRUE)

  s <- solo[, .(estacao, dia, gwet_sup, gwet_raiz)]
  s30 <- s[, .(estacao, dia = dia + 30, raiz30 = gwet_raiz)]
  s <- merge(s, s30, by = c("estacao", "dia"), all.x = TRUE)[, gwet_raiz_d30 := gwet_raiz - raiz30][, raiz30 := NULL]
  s[, init_date := dia + DEFAS_SOLO][, dia := NULL]
  out <- merge(CJ(estacao = est$estacao, init_date = inits), rol, by = c("estacao", "init_date"), all.x = TRUE)
  merge(out, s, by = c("estacao", "init_date"), all.x = TRUE)
}
