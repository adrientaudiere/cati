data(finch.ind)

test_that("Tstats returns correct class and structure", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  expect_s3_class(res, "Tstats")
  expect_true(all(
    c("Tstats", "variances", "traits", "ind.plot", "sp", "call") %in% names(res)
  ))
  expect_true(all(res$Tstats$T_IP.IC >= 0, na.rm = TRUE))
  expect_true(all(res$Tstats$T_IC.IR >= 0, na.rm = TRUE))
  expect_true(all(res$Tstats$T_PC.PR >= 0, na.rm = TRUE))
})

test_that("Tstats with nperm = NULL returns only observed values", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = NULL,
    printprogress = FALSE
  )
  expect_s3_class(res, "Tstats")
  expect_null(res$Tstats$T_IP.IC_nm)
  expect_null(res$Tstats$T_IC.IR_nm)
  expect_null(res$Tstats$T_PC.PR_nm)
})

test_that("Tstats with independantTraits = FALSE works", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE,
    independantTraits = FALSE
  )
  expect_s3_class(res, "Tstats")
})

test_that("sum_Tstats works with all type options", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  expect_no_error(sum_Tstats(res, type = "binary"))
  expect_no_error(sum_Tstats(res, type = "percent"))
  expect_no_error(sum_Tstats(res, type = "site"))
  expect_no_error(sum_Tstats(res, type = "p.value"))
  expect_no_error(sum_Tstats(res, type = "all"))
})

test_that("ses.Tstats returns a list with ses/ses.inf/ses.sup per index", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  ses_res <- ses.Tstats(res)
  expect_true(is.list(ses_res))
  expect_true(all(c("ses", "ses.inf", "ses.sup") %in% names(ses_res[[1]])))
})

test_that("Tstats print and summary S3 methods do not error", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  expect_no_error(print(res))
  expect_no_error(summary(res))
})

test_that("Tstats with SE > 0 works", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE,
    SE = 5
  )
  expect_s3_class(res, "Tstats")
})

test_that("null models keep each plot's own individuals or pool (GH-12, GH-18)", {
  set.seed(1)
  plots <- c("pa", "pb", "pc")
  d <- expand.grid(
    ind = 1:5,
    sp = c("x", "y"),
    plot = plots,
    stringsAsFactors = FALSE
  )
  d$trait <- rnorm(nrow(d), mean = match(d$plot, plots) * 50)
  d <- d[sample(nrow(d)), ] # rows not sorted by plot
  res <- Tstats(
    as.matrix(d$trait),
    ind.plot = factor(d$plot),
    sp = factor(d$sp),
    nperm = 5,
    printprogress = FALSE
  )
  # local null model only permutes within plots: community variance is unchanged
  for (i in 1:5) {
    expect_equal(res$variances$var_IC_nm1[i, 1, ], res$variances$var_IC[, 1])
  }

  pools <- lapply(c(1, 100), \(s) matrix(rnorm(400, 0, s), ncol = 1))
  res2 <- Tstats(
    as.matrix(d$trait[d$plot != "pc"]),
    ind.plot = factor(d$plot[d$plot != "pc"]),
    sp = factor(d$sp[d$plot != "pc"]),
    reg.pool = pools,
    nperm = 20,
    printprogress = FALSE
  )
  m <- colMeans(res2$variances$var_IC_nm2[, 1, ])
  expect_lt(m[1], 10)
  expect_gt(m[2], 1000)
})

test_that("T_IP.IC uses only the populations of its own plot (GH-11)", {
  set.seed(1)
  d <- expand.grid(
    ind = 1:5,
    sp = c("a", "b"),
    plot = c("p1", "p2", "p10", "A(1)"),
    stringsAsFactors = FALSE
  )
  d$trait <- rnorm(
    nrow(d),
    mean = ifelse(d$sp == "a", 0, 20),
    sd = ifelse(d$plot == "p10", 10, 1)
  )
  res <- Tstats(
    as.matrix(d$trait),
    ind.plot = factor(d$plot),
    sp = factor(d$sp),
    nperm = NULL,
    printprogress = FALSE
  )
  by_hand <- sapply(levels(factor(d$plot)), function(p) {
    x <- d[d$plot == p, ]
    mean(tapply(x$trait, x$sp, var)) / var(x$trait)
  })
  expect_equal(unname(res$Tstats$T_IP.IC[, 1]), unname(by_hand))
})

test_that("sum_Tstats rows T_IC.IR.inf and T_PC.PR.sup report their own statistic (GH-16)", {
  set.seed(1)
  plots <- paste0("c", 1:5)
  trs <- c("A", "B")
  obs <- function() matrix(0, 5, 2, dimnames = list(plots, trs))
  null <- function() {
    a <- array(rnorm(199 * 2 * 5), c(199, 2, 5))
    colnames(a) <- trs
    a
  }
  x <- list(
    Tstats = list(
      T_IP.IC = obs(),
      T_IC.IR = obs(),
      T_PC.PR = obs(),
      T_IP.IC_nm = null(),
      T_IC.IR_nm = null(),
      T_PC.PR_nm = null()
    ),
    nb.ind_plot = table(factor(plots)),
    namestraits = trs,
    call = quote(Tstats())
  )
  class(x) <- "Tstats"
  x$Tstats$T_IP.IC["c1", "A"] <- -10
  x$Tstats$T_IC.IR[, "A"] <- c(10, 10, -10, 10, 10)

  expect_equal(unname(sum_Tstats(x, type = "site")["T_IC.IR.inf", "A"]), "c3")
  expect_equal(
    unname(sum_Tstats(x, type = "binary")[
      c("T_IC.IR.sup", "T_PC.PR.sup"),
      "A"
    ]),
    c(TRUE, FALSE)
  )
  all <- sum_Tstats(x, type = "all")
  expect_equal(unname(all$sites["T_IC.IR.inf", "A"]), "c3")
  expect_false(unname(all$binary["T_PC.PR.sup", "A"]))
})

test_that("perm_pval handles ties, NA null values and untestable cells", {
  expect_equal(perm_pval(1, rep(1, 9), "inf"), 1)
  expect_equal(perm_pval(1, rep(1, 9), "sup"), 1)
  expect_true(is.na(perm_pval(NA, 1:9, "sup")))
  expect_true(is.na(perm_pval(1, rep(NA, 9), "sup")))
  expect_equal(perm_pval(10, c(1:4, rep(NA, 5)), "sup"), 1 / 5)
  expect_equal(perm_pval(10, c(1:4, rep(NA, 5)), "inf"), 1)
})

test_that("sum_Tstats and Pval agree on tail names", {
  set.seed(1)
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  pv <- sum_Tstats(res, type = "p.value")
  P <- Pval(res)
  expect_equal(
    unname(pv[paste("T_IC.IR.inf", rownames(res$Tstats$T_IC.IR)), 1]),
    unname(P[[2]][[1]]$pval.inf[, 1])
  )
})
