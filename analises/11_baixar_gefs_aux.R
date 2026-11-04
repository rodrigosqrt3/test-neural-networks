# Baixa variáveis do GEFSv12 reforecast (AWS, bucket público noaa-gefs-retrospective)
# recortadas numa caixa em torno do RS, só as mensagens GRIB necessárias, por
# intervalo de bytes (arquivos .idx). Não faz parte das análises (regra 6): gera
# dados/brutos/gefs_aux/<membro>/<ano>/<AAAAMMDD>.rds e pode ser retomado.
#
# Uso: Rscript analises/11_baixar_gefs_aux.R <ano_ini> <ano_fim> [membros] [conjunto]
#      membros: "c00" (padrão) ou "p01,p02,p03,p04"
#      conjunto: "completo" (padrão, ~25 campos por dia) ou "chave" (Tmax, T2m, chuva)
#
# Para cada rodada 00 UTC e horizonte H (dias 1, 3, 7) o "dia previsto" é o
# intervalo [24H, 24H+24] h de previsão (21h local da véspera às 21h local do
# dia), coberto por 4 janelas de 6 h nos campos de máximo/acumulado/média.
suppressMessages({ library(curl); library(terra) })
terraOptions(progress = 0)

args <- commandArgs(TRUE)
ANOS <- if (length(args) >= 2) as.integer(args[1]):as.integer(args[2]) else 2010L
MEMBROS <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else "c00"
CONJUNTO <- if (length(args) >= 4) args[4] else "completo"
HORIZONTES <- c(1, 3, 7)
CAIXA <- ext(-58.625, -48.375, -34.625, -26.375)          # RS + margem (bordas de célula de 0,25°)
BASE <- "https://noaa-gefs-retrospective.s3.amazonaws.com/GEFSv12/reforecast"
SAIDA <- "dados/brutos/gefs_aux"

# Especificação: arquivo, nível (como no .idx), tipo e horas (relativas a 24H)
#   tipo "inst": campo instantâneo em 24H+h; "jan": janela de 6 h que termina em 24H+h
janelas <- c(6, 12, 18, 24)
ESPEC <- list(
  tmax2m  = list("tmax_2m",     "2 m above ground",       "jan",  janelas, "max"),
  chuva   = list("apcp_sfc",    "surface",                "jan",  janelas, "acc"),
  t2m     = list("tmp_2m",      "2 m above ground",       "inst", c(12, 18)),
  q2m     = list("spfh_2m",     "2 m above ground",       "inst", 18),
  u10     = list("ugrd_hgt",    "10 m above ground",      "inst", 18),
  v10     = list("vgrd_hgt",    "10 m above ground",      "inst", 18),
  pnmm    = list("pres_msl",    "mean sea level",         "inst", c(6, 18)),
  nuvens  = list("tcdc_eatm",   "entire atmosphere",      "jan",  c(12, 18), "ave"),
  radsw   = list("dswrf_sfc",   "surface",                "jan",  c(12, 18), "ave"),
  cape    = list("cape_sfc",    "surface",                "inst", 18),
  pwat    = list("pwat_eatm",   "entire atmosphere (considered as a single layer)", "inst", 18),
  fsens   = list("shtfl_sfc",   "surface",                "jan",  18, "ave"),
  flat    = list("lhtfl_sfc",   "surface",                "jan",  18, "ave"),
  solo    = list("soilw_bgrnd", "0-0.1 m below ground",   "inst", 12),
  t850    = list("tmp_pres",    "850 mb",                 "inst", 18),
  q850    = list("spfh_pres",   "850 mb",                 "inst", 18),
  u850    = list("ugrd_pres",   "850 mb",                 "inst", 18),
  v850    = list("vgrd_pres",   "850 mb",                 "inst", 18),
  z500    = list("hgt_pres_abv700mb", "500 mb",           "inst", 18)   # níveis < 700 hPa ficam nesse arquivo
)
if (CONJUNTO == "chave") ESPEC <- ESPEC[c("tmax2m", "chuva", "t2m")]

descritor <- function(e, horas) {
  if (e[[3]] == "inst") sprintf("%s:%d hour fcst", e[[2]], horas)
  else sprintf("%s:%d-%d hour %s fcst", e[[2]], horas - 6, horas, e[[5]])
}

pool <- new_pool(host_con = 12, total_con = 24)
buscar <- function(urls, ranges = NULL) {         # baixa em paralelo; devolve lista de raw
  out <- vector("list", length(urls))
  for (k in seq_along(urls)) local({
    k <- k; h <- new_handle(failonerror = TRUE)
    if (!is.null(ranges)) handle_setheaders(h, Range = ranges[k])
    # out[k] <- list(...) e não out[[k]] <- NULL: atribuir NULL com [[ ]] APAGA o
    # elemento e desloca os seguintes (as mensagens trocariam de lugar)
    curl_fetch_multi(urls[k], handle = h, pool = pool,
                     done = function(r) out[k] <<- list(if (r$status_code %in% c(200, 206)) r$content else NULL),
                     fail = function(e) out[k] <<- list(NULL))
  })
  multi_run(pool = pool)
  stopifnot(length(out) == length(urls))
  out
}
# erros de rede: tenta de novo; erro de especificação (mensagem não encontrada) não adianta repetir
com_tentativas <- function(f, n = 4) {
  for (i in seq_len(n)) {
    r <- try(f(), silent = TRUE)
    if (!inherits(r, "try-error")) return(r)
    if (grepl("não encontrada", r)) break
    Sys.sleep(5 * i)
  }
  stop(r)
}

baixar_rodada <- function(data, membro) {
  ano <- format(data, "%Y"); ymd <- format(data, "%Y%m%d")
  arq <- file.path(SAIDA, membro, ano, paste0(ymd, ".rds"))
  if (file.exists(arq)) return("já existe")
  pasta <- sprintf("%s/%s/%s00/%s/Days:1-10/", BASE, ano, ymd, membro)
  arquivos <- unique(vapply(ESPEC, `[[`, "", 1))
  urls_grib <- sprintf("%s%s_%s00_%s.grib2", pasta, arquivos, ymd, membro)
  idxs <- buscar(paste0(urls_grib, ".idx"))
  if (any(vapply(idxs, is.null, NA))) stop("idx ausente")
  idxs <- setNames(lapply(idxs, function(x) strsplit(rawToChar(x), "\n")[[1]]), arquivos)
  # localiza cada mensagem pedida: (variável, H, rótulo) -> intervalo de bytes
  pedidos <- list()
  for (nm in names(ESPEC)) for (H in HORIZONTES) {
    e <- ESPEC[[nm]]; linhas <- idxs[[e[[1]]]]
    ini <- as.numeric(sapply(strsplit(linhas, ":"), `[`, 2))
    for (h in e[[4]]) {
      d <- descritor(e, 24 * H + h)
      i <- which(endsWith(sub(":ENS=.*$", "", linhas), d))
      # alguns arquivos da NOAA têm cada mensagem gravada duas vezes: usa a primeira;
      # se a mensagem não existe mesmo, o campo fica NA (não perde o dia inteiro)
      if (!length(i)) {
        pedidos[[length(pedidos) + 1]] <- list(var = nm, H = H, h = h, url = NA_character_, range = NA_character_)
        next
      }
      i <- i[1]
      fim <- if (i < length(linhas)) ini[i + 1] - 1 else NA
      pedidos[[length(pedidos) + 1]] <- list(var = nm, H = H, h = h, url = urls_grib[match(e[[1]], arquivos)],
                                             range = sprintf("bytes=%.0f-%s", ini[i], if (is.na(fim)) "" else sprintf("%.0f", fim)))
    }
  }
  origem <- vapply(pedidos, `[[`, "", "url")
  existe <- !is.na(origem)
  msgs <- vector("list", length(pedidos))
  msgs[existe] <- buscar(origem[existe], vapply(pedidos[existe], `[[`, "", "range"))
  # pedaço que falhou: se o GRIB não existe no bucket (404, mesmo com .idx), o
  # campo fica NA; outra falha (rede) faz a rodada inteira ser tentada de novo
  falhou <- which(existe & vapply(msgs, is.null, NA))
  for (u in unique(origem[falhou])) {
    st <- curl_fetch_memory(u, handle = new_handle(nobody = TRUE))$status_code
    if (st != 404) stop("falha em mensagens")
    existe[origem == u] <- FALSE
  }
  origem[!existe] <- NA
  # um GRIB por arquivo de origem (as grades diferem: níveis acima de 700 hPa vêm
  # em 0,5°), lido de uma vez, recortado e levado à grade de 0,25° da caixa
  camadas <- vector("list", length(pedidos)); molde <- NULL
  for (u in unique(origem[existe])) {
    k <- which(origem == u)
    tmp <- tempfile(fileext = ".grib2")
    con <- file(tmp, "wb"); for (m in msgs[k]) writeBin(m, con); close(con)
    rr <- crop(rast(tmp), CAIXA, snap = "out")
    stopifnot(nlyr(rr) == length(k))
    if (is.null(molde) && all(abs(res(rr) - 0.25) < 1e-6)) molde <- crop(rr[[1]], CAIXA)
    camadas[k] <- as.list(rr)
    unlink(tmp)
  }
  if (is.null(molde)) stop("nenhum campo em 0,25°")
  r <- rast(lapply(camadas, function(x) {
    if (is.null(x)) return(init(molde, NA))                  # campo ausente no arquivo da NOAA
    if (compareGeom(x, molde, stopOnError = FALSE)) x else resample(x, molde, method = "bilinear")
  }))
  v <- values(r, mat = TRUE)                                  # células x camadas
  rotulo <- vapply(pedidos, function(p) sprintf("%s_h%02d", p$var, p$h), "")
  out <- list(init = data, membro = membro, grade = list(lon = xFromCol(r), lat = yFromRow(r)))
  for (H in HORIZONTES) {
    k <- which(vapply(pedidos, `[[`, 0, "H") == H)
    out[[paste0("H", H)]] <- setNames(lapply(k, function(j) matrix(v[, j], nrow = nrow(r), byrow = TRUE)), rotulo[k])
  }
  dir.create(dirname(arq), recursive = TRUE, showWarnings = FALSE)
  # grava num temporário e renomeia: uma interrupção nunca deixa arquivo pela metade
  tmp_rds <- paste0(arq, ".tmp")
  saveRDS(out, tmp_rds)
  file.rename(tmp_rds, arq)
  sprintf("%d mensagens, %.1f MB%s", sum(existe), sum(lengths(msgs)) / 2^20,
          if (any(!existe)) paste0(" | ausentes na NOAA: ", paste(rotulo[!existe], collapse = ",")) else "")
}

log_arq <- file.path(SAIDA, sprintf("log_%d_%d_%s.txt", min(ANOS), max(ANOS), paste(MEMBROS, collapse = "")))
dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)
for (membro in MEMBROS) for (ano in ANOS) {
  datas <- seq(as.Date(sprintf("%d-01-01", ano)), as.Date(sprintf("%d-12-31", ano)), by = "day")
  if (Sys.getenv("PILOTO") != "") datas <- datas[seq_len(as.integer(Sys.getenv("PILOTO")))]
  for (d in as.list(datas)) {
    t0 <- Sys.time()
    msg <- tryCatch(com_tentativas(function() baixar_rodada(d, membro)),
                    error = function(e) paste("ERRO:", conditionMessage(e)))
    linha <- sprintf("%s %s %s | %s | %.1f s", format(Sys.time(), "%H:%M:%S"), membro, d, msg,
                     as.numeric(Sys.time() - t0, units = "secs"))
    cat(linha, "\n"); cat(linha, "\n", file = log_arq, append = TRUE)
  }
}
