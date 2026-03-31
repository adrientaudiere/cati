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
