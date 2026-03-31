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
