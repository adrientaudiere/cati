data(finch.ind)

test_that("ses returns a list with ses, ses.inf, ses.sup", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  ses_res <- ses(res$Tstats$T_IP.IC, res$Tstats$T_IP.IC_nm)
  expect_named(ses_res, c("ses", "ses.inf", "ses.sup"))
  expect_true(is.numeric(ses_res$ses))
  expect_true(is.numeric(ses_res$ses.inf))
  expect_true(is.numeric(ses_res$ses.sup))
})

test_that("ses confidence bounds bracket the SES", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  ses_res <- ses(res$Tstats$T_IP.IC, res$Tstats$T_IP.IC_nm)
  expect_true(all(ses_res$ses.inf <= ses_res$ses.sup, na.rm = TRUE))
})

test_that("as.listofindex converts a Tstats object to a list", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  loi <- as.listofindex(list(res))
  expect_true(is.list(loi))
  expect_true(length(loi) > 0)
})

test_that("ses.listofindex returns ses/ses.inf/ses.sup for each index", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  loi <- as.listofindex(list(res))
  ses_loi <- ses.listofindex(loi)
  expect_true(is.list(ses_loi))
  expect_true(all(c("ses", "ses.inf", "ses.sup") %in% names(ses_loi[[1]])))
})

test_that("ses and Pval use the null layout when sizes collide (GH-20)", {
  set.seed(3)
  ind.plot <- factor(rep(paste0("P", 1:3), each = 12))
  sp <- factor(rep(c("a", "b", "c"), length.out = 36))
  traits <- matrix(
    rnorm(108),
    ncol = 3,
    dimnames = list(NULL, paste0("tr", 1:3))
  )

  # 3 plots x 3 traits: Tstats null is permutations x traits x plots
  ts <- Tstats(traits, ind.plot, sp, nperm = 20, printprogress = FALSE)
  nm <- ts$Tstats$T_IC.IR_nm
  expected <- (ts$Tstats$T_IC.IR - apply(nm, c(3, 2), mean)) /
    apply(nm, c(3, 2), sd)
  expect_no_warning(s <- ses.Tstats(ts))
  expect_equal(unname(unclass(s[[2]]$ses)), unname(expected))
  expect_no_warning(s_bare <- ses(ts$Tstats$T_IC.IR, nm, null.dims = c(3, 2)))
  expect_equal(unname(unclass(s_bare$ses)), unname(expected))

  p_inf <- sapply(1:3, \(t) {
    sapply(1:3, \(g) (1 + sum(nm[, t, g] <= ts$Tstats$T_IC.IR[g, t])) / 21)
  })
  P <- Pval(ts)[[2]]
  expect_length(P, 3)
  expect_equal(unname(sapply(P, \(x) as.vector(x$pval.inf))), p_inf)

  # 3 plots x 2 traits, nperm = 3: ComIndex null is traits x plots x permutations
  ci <- ComIndex(
    traits[, 1:2],
    index = "tapply(x, ind.plot, mean)",
    nullmodels = "regional.ind",
    ind.plot = ind.plot,
    sp = sp,
    nperm = 3,
    printprogress = FALSE
  )
  nm <- ci$null[[1]]
  expected <- (ci$obs[[1]] - apply(nm, c(2, 1), mean)) / apply(nm, c(2, 1), sd)
  expect_equal(
    unname(unclass(ses.listofindex(as.listofindex(ci))[[1]]$ses)),
    unname(expected)
  )
})

test_that("ses accepts a transposed obs with a null from as.listofindex", {
  res <- Tstats(
    traits.finch,
    ind.plot = ind.plot.finch,
    sp = sp.finch,
    nperm = 9,
    printprogress = FALSE
  )
  loi <- as.listofindex(list(res))
  s <- ses(loi[[1]], loi[[2]])$ses
  expect_equal(ses(t(loi[[1]]), loi[[2]])$ses, t(s))
})
