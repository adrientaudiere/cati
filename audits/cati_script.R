# ============================================================
# External audit script for cati v0.99.5
# Generated: 2026-03-30
# Self-contained and re-runnable.
# ============================================================

devtools::load_all("/home/adrien/cati", quiet = TRUE)

audit_results <- data.frame(
  fn = character(),
  test = character(),
  status = character(),
  message = character(),
  stringsAsFactors = FALSE
)

audit_call <- function(fn_name, expr, test_label) {
  result <- tryCatch(
    withCallingHandlers(
      {
        force(expr)
        list(status = "OK", message = "")
      },
      warning = function(w) {
        invokeRestart("muffleWarning")
        list(status = "WARNING", message = conditionMessage(w))
      }
    ),
    error = function(e) list(status = "ERROR", message = conditionMessage(e))
  )
  data.frame(
    fn = fn_name,
    test = test_label,
    status = result$status,
    message = result$message,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# TEST DATA
# ============================================================
set.seed(42)

N <- 120 # individuals
Nsp <- 6 # species
Ncom <- 4 # communities

sp_vec <- factor(rep(paste0("sp", 1:Nsp), each = N / Nsp))
plot_vec <- factor(rep(paste0("site", 1:Ncom), times = N / Ncom))

# Full traits matrix with NAs
traits_na <- matrix(
  c(rnorm(N, 50, 10), rnorm(N, 2, 0.5), rnorm(N, 100, 20), runif(N, 0, 1)),
  nrow = N,
  ncol = 4,
  dimnames = list(NULL, c("Tr1", "Tr2", "Tr3", "Tr4"))
)
traits_na[sample(N * 4, 20)] <- NA # sprinkle NAs

# Complete traits matrix (no NAs)
traits_complete <- traits_na
traits_complete[is.na(traits_complete)] <- 0

# Single trait (as matrix)
traits_1col <- traits_na[, 1, drop = FALSE]

# Community matrix: species × sites
comm_sp <- table(sp_vec, plot_vec)
class(comm_sp) <- "matrix"
comm_sp <- comm_sp + 0L # integer matrix

# Species-mean traits (for RaoRel / dist-based functions)
traits_sp_mean <- apply(traits_complete, 2, function(x) tapply(x, sp_vec, mean))

# pairwise distance matrix (sp × sp), halved squared Euclidean for Rao
mat_dist <- (as.matrix(dist(traits_sp_mean))^2) / 2

# factors for partvar
genus_vec <- factor(sub("sp[0-9]", "genX", as.character(sp_vec)))
factors_mat <- cbind(
  site = as.factor(as.character(plot_vec)),
  species = sp_vec,
  genus = genus_vec
)

# Small dataset for slow/heavy functions
N_small <- 30
sp_small <- factor(rep(paste0("sp", 1:3), each = N_small / 3))
pl_small <- factor(rep(paste0("s", 1:3), times = N_small / 3))
tr_small <- matrix(
  rnorm(N_small * 2),
  nrow = N_small,
  ncol = 2,
  dimnames = list(NULL, c("A", "B"))
)

# ============================================================
# NND METRICS: CVNND, MNND, MinNND, SDNND, SDND, MND
# ============================================================

audit_results <- rbind(
  audit_results,
  audit_call(
    "CVNND",
    CVNND(traits_na[, 1], na.rm = TRUE),
    "single trait with NAs"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "CVNND",
    CVNND(traits_na[, 1], div_range = TRUE, na.rm = TRUE),
    "div_range = TRUE"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "CVNND",
    CVNND(traits_complete, na.rm = TRUE),
    "multi-trait matrix"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "CVNND",
    CVNND(traits_complete, scale.tr = FALSE, na.rm = TRUE),
    "scale.tr = FALSE"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "MNND",
    MNND(traits_na[, 1], na.rm = TRUE),
    "single trait with NAs"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call("MNND", MNND(traits_complete, na.rm = TRUE), "multi-trait")
)

audit_results <- rbind(
  audit_results,
  audit_call("MinNND", MinNND(traits_na[, 1], na.rm = TRUE), "single trait")
)

audit_results <- rbind(
  audit_results,
  audit_call("SDNND", SDNND(traits_na[, 1], na.rm = TRUE), "single trait")
)

audit_results <- rbind(
  audit_results,
  audit_call("SDND", SDND(traits_na[, 1], na.rm = TRUE), "single trait")
)

audit_results <- rbind(
  audit_results,
  audit_call("MND", MND(traits_na[, 1], na.rm = TRUE), "single trait")
)

# ============================================================
# SumBL / MinMaxMST
# ============================================================

tr_nomiss <- na.omit(traits_na[1:60, 1, drop = FALSE])

audit_results <- rbind(
  audit_results,
  audit_call("SumBL", SumBL(tr_nomiss, gower.dist = TRUE), "gower.dist = TRUE")
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "SumBL",
    SumBL(tr_nomiss, gower.dist = FALSE),
    "gower.dist = FALSE"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "SumBL",
    SumBL(tr_nomiss, method.hclust = "complete"),
    "method.hclust = complete"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "MinMaxMST",
    MinMaxMST(tr_nomiss, gower.dist = TRUE),
    "gower.dist = TRUE"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "MinMaxMST",
    MinMaxMST(tr_nomiss, gower.dist = FALSE),
    "gower.dist = FALSE"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "MinMaxMST",
    MinMaxMST(tr_nomiss, gower.dist = FALSE, scale.tr = FALSE),
    "scale.tr = FALSE"
  )
)

# ============================================================
# Tstats
# ============================================================

res_tstats <- suppressWarnings(Tstats(
  traits_na,
  ind.plot = plot_vec,
  sp = sp_vec,
  nperm = 9,
  printprogress = FALSE
))

audit_results <- rbind(
  audit_results,
  audit_call(
    "Tstats",
    Tstats(
      traits_na,
      ind.plot = plot_vec,
      sp = sp_vec,
      nperm = 9,
      printprogress = FALSE
    ),
    "normal traits with NAs"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "Tstats",
    Tstats(
      traits_na,
      ind.plot = plot_vec,
      sp = sp_vec,
      nperm = NULL,
      printprogress = FALSE
    ),
    "nperm = NULL"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "Tstats",
    Tstats(
      traits_na,
      ind.plot = plot_vec,
      sp = sp_vec,
      nperm = 9,
      SE = 2,
      printprogress = FALSE
    ),
    "SE = 2"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "Tstats",
    Tstats(
      traits_na,
      ind.plot = plot_vec,
      sp = sp_vec,
      nperm = 9,
      independantTraits = FALSE,
      printprogress = FALSE
    ),
    "independantTraits = FALSE"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call("print.Tstats", print(res_tstats), "print method")
)

audit_results <- rbind(
  audit_results,
  audit_call("summary.Tstats", summary(res_tstats), "summary method")
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "sum_Tstats",
    sum_Tstats(res_tstats, type = "binary"),
    "type = binary"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "sum_Tstats",
    sum_Tstats(res_tstats, type = "percent"),
    "type = percent"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call("sum_Tstats", sum_Tstats(res_tstats, type = "site"), "type = site")
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "sum_Tstats",
    sum_Tstats(res_tstats, type = "p.value"),
    "type = p.value"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call("sum_Tstats", sum_Tstats(res_tstats, type = "all"), "type = all")
)

audit_results <- rbind(
  audit_results,
  audit_call("ses.Tstats", ses.Tstats(res_tstats), "normal")
)

pdf(tempfile())
audit_results <- rbind(
  audit_results,
  audit_call("plot.Tstats", plot(res_tstats, type = "normal"), "type = normal")
)
audit_results <- rbind(
  audit_results,
  audit_call("plot.Tstats", plot(res_tstats, type = "simple"), "type = simple")
)
audit_results <- rbind(
  audit_results,
  audit_call(
    "plot.Tstats",
    plot(res_tstats, type = "simple_range"),
    "type = simple_range"
  )
)
audit_results <- rbind(
  audit_results,
  audit_call(
    "plot.Tstats",
    plot(res_tstats, type = "barplot"),
    "type = barplot"
  )
)
audit_results <- rbind(
  audit_results,
  audit_call("barplot.Tstats", barplot(res_tstats), "normal")
)
dev.off()

# ============================================================
# ses / as.listofindex / ses.listofindex / plot.listofindex
# ============================================================

audit_results <- rbind(
  audit_results,
  audit_call(
    "ses",
    ses(res_tstats$Tstats$T_IP.IC, res_tstats$Tstats$T_IP.IC_nm),
    "matrix obs"
  )
)

loi <- as.listofindex(list(res_tstats))
audit_results <- rbind(
  audit_results,
  audit_call(
    "as.listofindex",
    as.listofindex(list(res_tstats)),
    "single Tstats"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call("ses.listofindex", ses.listofindex(loi), "from Tstats")
)

pdf(tempfile())
audit_results <- rbind(
  audit_results,
  audit_call("plot.listofindex", plot(loi), "normal")
)
dev.off()

audit_results <- rbind(
  audit_results,
  audit_call("Pval", Pval(res_tstats), "from Tstats")
)

# ============================================================
# ComIndex
# ============================================================

funct_basic <- c("mean(x, na.rm = TRUE)", "sd(x, na.rm = TRUE)")

res_ci_reg <- suppressWarnings(
  ComIndex(
    traits = traits_na,
    index = funct_basic,
    sp = sp_vec,
    nullmodels = "regional.ind",
    ind.plot = plot_vec,
    nperm = 9,
    printprogress = FALSE
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "ComIndex",
    ComIndex(
      traits = traits_na,
      index = funct_basic,
      sp = sp_vec,
      nullmodels = "regional.ind",
      ind.plot = plot_vec,
      nperm = 9,
      printprogress = FALSE
    ),
    "regional.ind"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "ComIndex",
    ComIndex(
      traits = traits_na,
      index = "mean(x, na.rm = TRUE)",
      sp = sp_vec,
      nullmodels = "local",
      ind.plot = plot_vec,
      nperm = 9,
      printprogress = FALSE
    ),
    "local"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "ComIndex",
    ComIndex(
      traits = traits_na,
      index = "mean(x, na.rm = TRUE)",
      sp = sp_vec,
      nullmodels = "regional.pop",
      ind.plot = plot_vec,
      nperm = 9,
      printprogress = FALSE
    ),
    "regional.pop"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "ComIndex",
    ComIndex(
      traits = traits_na,
      index = "mean(x, na.rm = TRUE)",
      sp = sp_vec,
      nullmodels = "regional.pop.prab",
      ind.plot = plot_vec,
      nperm = 9,
      printprogress = FALSE
    ),
    "regional.pop.prab"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "ComIndex",
    ComIndex(
      traits = traits_na,
      index = "mean(x, na.rm = TRUE)",
      sp = sp_vec,
      nullmodels = "regional.ind",
      ind.plot = plot_vec,
      nperm = NULL,
      printprogress = FALSE
    ),
    "nperm = NULL"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call("print.ComIndex", print(res_ci_reg), "print method")
)

audit_results <- rbind(
  audit_results,
  audit_call("summary.ComIndex", summary(res_ci_reg), "summary method")
)

audit_results <- rbind(
  audit_results,
  audit_call("Pval", Pval(res_ci_reg), "from ComIndex")
)

pdf(tempfile())
audit_results <- rbind(
  audit_results,
  audit_call("plot.ComIndex", plot(res_ci_reg), "type = normal")
)
dev.off()

# ============================================================
# ComIndexMulti
# NOTE: ComIndexMulti crashes with >1 index — documented as bug
# ============================================================

audit_results <- rbind(
  audit_results,
  data.frame(
    fn = "ComIndexMulti",
    test = "multiple indices",
    status = "ERROR",
    message = "Known bug: crashes with >1 index ('list' cannot be coerced to 'double')",
    stringsAsFactors = FALSE
  )
)

res_cim <- tryCatch(
  suppressWarnings(
    ComIndexMulti(
      traits = traits_na,
      index = "mean(x, na.rm = TRUE)",
      by.factor = plot_vec,
      sp = sp_vec,
      nullmodels = "regional.ind",
      ind.plot = plot_vec,
      nperm = 9,
      printprogress = FALSE
    )
  ),
  error = function(e) NULL
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "ComIndexMulti",
    ComIndexMulti(
      traits = traits_na,
      index = "mean(x, na.rm = TRUE)",
      by.factor = plot_vec,
      sp = sp_vec,
      nullmodels = "regional.ind",
      ind.plot = plot_vec,
      nperm = 9,
      printprogress = FALSE
    ),
    "single index with by.factor"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call("print.ComIndexMulti", print(res_cim), "print method")
)

audit_results <- rbind(
  audit_results,
  audit_call("summary.ComIndexMulti", summary(res_cim), "summary method")
)

pdf(tempfile())
audit_results <- rbind(
  audit_results,
  audit_call("plot.ComIndexMulti", plot(res_cim), "type = normal")
)
dev.off()

# ============================================================
# IndexByGroups / samplingSubsetData / AbToInd
# ============================================================

audit_results <- rbind(
  audit_results,
  audit_call(
    "IndexByGroups",
    IndexByGroups(c("mean(x)", "sd(x)"), "pop"),
    "two metrics"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "samplingSubsetData",
    samplingSubsetData(
      d = traits_na,
      sampUnit = sp_vec,
      nperm = 3,
      prop = c(50, 100)
    ),
    "type = proportion"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "samplingSubsetData",
    samplingSubsetData(
      d = traits_na,
      sampUnit = sp_vec,
      nperm = 2,
      type = "count",
      Size = c(5, 10)
    ),
    "type = count"
  )
)

# AbToInd: needs a species×sites abundance matrix and species-level traits
comm_t <- t(comm_sp) # sites × species
tr_sp <- matrix(
  rnorm(Nsp * 4),
  nrow = Nsp,
  ncol = 4,
  dimnames = list(paste0("sp", 1:Nsp), c("T1", "T2", "T3", "T4"))
)
audit_results <- rbind(
  audit_results,
  audit_call("AbToInd", AbToInd(traits = tr_sp, com = comm_t), "count type")
)

# ============================================================
# RaoRel
# ============================================================

audit_results <- rbind(
  audit_results,
  audit_call(
    "RaoRel",
    RaoRel(
      sample = comm_sp,
      dfunc = mat_dist,
      dphyl = NULL,
      weight = FALSE,
      Jost = FALSE,
      structure = NULL
    ),
    "weight = FALSE"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "RaoRel",
    RaoRel(
      sample = comm_sp,
      dfunc = mat_dist,
      dphyl = NULL,
      weight = TRUE,
      Jost = FALSE,
      structure = NULL
    ),
    "weight = TRUE"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "RaoRel",
    RaoRel(
      sample = comm_sp,
      dfunc = mat_dist,
      dphyl = NULL,
      weight = FALSE,
      Jost = TRUE,
      structure = NULL
    ),
    "Jost = TRUE"
  )
)

# ============================================================
# Fred (requires no NAs)
# ============================================================

audit_results <- rbind(
  audit_results,
  audit_call(
    "Fred",
    Fred(traits = traits_complete[, 1:3], ind.plot = plot_vec),
    "3 traits no NAs"
  )
)

# ============================================================
# partvar / barPartvar / piePartvar
# ============================================================

res_pv <- suppressWarnings(
  partvar(traits = traits_na, factors = factors_mat, printprogress = FALSE)
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "partvar",
    partvar(traits = traits_na, factors = factors_mat, printprogress = FALSE),
    "3 nested factors"
  )
)

pdf(tempfile())
audit_results <- rbind(
  audit_results,
  audit_call("barPartvar", barPartvar(res_pv), "normal")
)
audit_results <- rbind(
  audit_results,
  audit_call("piePartvar", piePartvar(res_pv), "normal")
)
dev.off()

# ============================================================
# decompCTRE / traitflex.anova / plot.traitflex / print.traitflex
# ============================================================

res_dc <- suppressWarnings(
  decompCTRE(
    traits = traits_na,
    sp = sp_vec,
    ind.plot = plot_vec,
    printprogress = FALSE
  )
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "decompCTRE",
    decompCTRE(
      traits = traits_na,
      sp = sp_vec,
      ind.plot = plot_vec,
      printprogress = FALSE
    ),
    "normal"
  )
)

pdf(tempfile())
audit_results <- rbind(
  audit_results,
  audit_call("barplot.decompCTRE", barplot(res_dc), "resume = TRUE")
)
audit_results <- rbind(
  audit_results,
  audit_call(
    "barplot.decompCTRE",
    barplot(res_dc, resume = FALSE),
    "resume = FALSE"
  )
)
dev.off()

tf <- res_dc[[1]]
audit_results <- rbind(
  audit_results,
  audit_call("print.traitflex", print(tf), "normal")
)

pdf(tempfile())
audit_results <- rbind(
  audit_results,
  audit_call("plot.traitflex", plot(tf), "normal")
)
audit_results <- rbind(
  audit_results,
  audit_call(
    "plot.traitflex",
    plot(tf, use.percentage = FALSE),
    "use.percentage = FALSE"
  )
)
dev.off()

# traitflex.anova directly
specif_avg <- tapply(traits_complete[, 1], plot_vec, mean)
const_avg <- rep(mean(traits_complete[, 1]), length(specif_avg))
audit_results <- rbind(
  audit_results,
  audit_call(
    "traitflex.anova",
    traitflex.anova(~1, specif.avg = specif_avg, const.avg = const_avg),
    "intercept-only formula"
  )
)

# ============================================================
# RandCom
# ============================================================

audit_results <- rbind(
  audit_results,
  audit_call("RandCom", RandCom(), "defaults")
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "RandCom",
    RandCom(Ncom = 5, Nsp = 10, Nind.com = 50, Filter = "Internal"),
    "Internal filter"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call("RandCom", RandCom(Filter = "External"), "External filter")
)

audit_results <- rbind(
  audit_results,
  audit_call("RandCom", RandCom(Filter = "Both"), "Both filters")
)

# ============================================================
# Plotting helpers: plotRandtest, plotSpVar, plotDistri,
#                   plotSpPop, plotSESvar, plotCorTstats
# ============================================================

pdf(tempfile())

audit_results <- rbind(
  audit_results,
  audit_call("plotRandtest", plotRandtest(res_tstats), "from Tstats")
)

audit_results <- rbind(
  audit_results,
  audit_call(
    "plotSpVar",
    plotSpVar(traits_na, sp_vec, plot_vec),
    "traits with NAs"
  )
)

audit_results <- rbind(
  audit_results,
  audit_call("plotDistri", plotDistri(traits_na, plot_vec), "by plot")
)

audit_results <- rbind(
  audit_results,
  audit_call("plotSpPop", plotSpPop(traits_na, sp_vec, plot_vec), "normal")
)

audit_results <- rbind(
  audit_results,
  audit_call("plotSESvar", plotSESvar(res_tstats), "from Tstats")
)

audit_results <- rbind(
  audit_results,
  audit_call("plotCorTstats", plotCorTstats(res_tstats), "from Tstats")
)

dev.off()

# ============================================================
# RESULTS
# ============================================================

write.csv(
  audit_results,
  "/home/adrien/cati/audits/cati_results.csv",
  row.names = FALSE
)
cat("\nAudit complete:", nrow(audit_results), "tests run\n")
cat("OK:     ", sum(audit_results$status == "OK"), "\n")
cat("WARNING:", sum(audit_results$status == "WARNING"), "\n")
cat("ERROR:  ", sum(audit_results$status == "ERROR"), "\n")
