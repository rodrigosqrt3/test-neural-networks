# Combinação de previsões (sem treinar nada): média ponderada dos quantis
# (Vincentização) de dois modelos, q = w * q_A + (1 - w) * q_B, com w em
# {0, 0.1, ..., 1} escolhido pelo CRPS na VALIDAÇÃO de cada dobra e horizonte.
#
# Uso: Rscript analises/10_combinacoes.R     (depois de 06 e 09)
# Saída: configs "base_2sem" e "c_*" nas mesmas pastas de exp4 / exp4_dobra<k>,
#        lidas por 08_resumo_dobras.R.
#
# A rede usou a validação para o early stopping, então o CRPS de validação dela é
# um pouco otimista: isso puxa w para a rede, ou seja, contra a combinação.
source("R/avaliacao.R"); source("R/experimento.R")

pasta <- function(dobra, H) file.path("resultados", if (dobra > 0) paste0("exp4_dobra", dobra) else "exp4", paste0("h", H))
arq_grade <- function(dobra, H, nome) {
  if (dobra == 0 && nome == "x4_base_s1") return(sprintf("resultados/exp3/h%d/.grades_base.rds", H))
  file.path(pasta(dobra, H), paste0(".grades_", nome, ".rds"))
}
# média dos quantis entre sementes (as que existirem)
grade <- function(dobra, H, nomes) {
  fs <- vapply(nomes, function(n) arq_grade(dobra, H, n), character(1)); fs <- fs[file.exists(fs)]
  if (!length(fs)) return(NULL)
  gs <- lapply(fs, readRDS)
  list(q_va = Reduce(`+`, lapply(gs, `[[`, "q_va")) / length(gs),
       q_te = Reduce(`+`, lapply(gs, `[[`, "q_te")) / length(gs), n = length(gs))
}
COMPONENTES <- list(base = c("x4_base_s1", "x4_base_s2"), tudo = c("x4_tudo_s1", "x4_tudo_s2"),
                    lagens = "x4_lagens", qrf = "x4_qrf", emos = "x4_emos", arima = "arima")
PARES <- list(c_base_qrf = c("base", "qrf"), c_base_emos = c("base", "emos"),
              c_base_arima = c("base", "arima"), c_base_lagens = c("base", "lagens"),
              c_base_tudo = c("base", "tudo"))
W <- seq(0, 1, by = 0.1)

lim_de <- function(d) list(q90 = d$lim_q90, q95 = d$lim_q95, q99 = d$lim_q99, abs35 = d$lim_abs35)
salvar <- function(dobra, H, nome, q_va, q_te, ref, extra = list()) {
  res <- c(list(config = nome, horizonte = H,
                validacao = cbind(data.table(data = ref$validacao$data),
                                  escores_por_caso(q_va, ref$validacao$y, lim_de(ref$validacao))),
                teste = cbind(data.table(data = ref$teste$data),
                              escores_por_caso(q_te, ref$teste$y, lim_de(ref$teste)))), extra)
  saveRDS(res, file.path(pasta(dobra, H), paste0(nome, ".rds")))
  res
}

linhas <- list()
for (dobra in 0:2) for (H in c(1, 3, 7)) {
  ref <- readRDS(file.path(pasta(dobra, H), "arima.rds"))    # y, datas e limiares da dobra
  comp <- lapply(COMPONENTES, function(n) grade(dobra, H, n))
  comp <- Filter(Negate(is.null), comp)
  stopifnot(nrow(comp$base$q_te) == nrow(ref$teste))
  # base só com a média das 2 sementes (para separar "mais redes" de "outro modelo")
  b <- salvar(dobra, H, "base_2sem", comp$base$q_va, comp$base$q_te, ref, list(n_sementes = comp$base$n))
  for (nm in names(PARES)) {
    p <- PARES[[nm]]
    if (!all(p %in% names(comp))) next
    A <- comp[[p[1]]]; B <- comp[[p[2]]]
    cv <- vapply(W, function(w) mean(crps_quantis(w * A$q_va + (1 - w) * B$q_va, ref$validacao$y)), numeric(1))
    w <- W[which.min(cv)]
    r <- salvar(dobra, H, nm, w * A$q_va + (1 - w) * B$q_va, w * A$q_te + (1 - w) * B$q_te, ref, list(peso = w))
    linhas[[length(linhas) + 1]] <- data.table(dobra, H, combinacao = nm, peso_base = w,
                                               crps_teste = mean(r$teste$crps), crps_base_2sem = mean(b$teste$crps))
  }
}
tab <- rbindlist(linhas)
tab[, ganho_pct := round(100 * (1 - crps_teste / crps_base_2sem), 2)]
print(tab)
fwrite(tab, "resultados/tabelas/exp4_combinacoes_pesos.csv")
