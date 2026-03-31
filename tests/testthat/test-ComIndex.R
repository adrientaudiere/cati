data(finch.ind)

test_that("ComIndex with regional.ind returns correct class and structure", {
  res <- ComIndex(
    traits = traits.finch,
    index = c("mean(x, na.rm = TRUE)", "sd(x, na.rm = TRUE)"),
    sp = sp.finch,
    nullmodels = "regional.ind",
    ind.plot = ind.plot.finch,
    nperm = 9,
    printprogress = FALSE
  )
  expect_s3_class(res, "ComIndex")
  expect_true(all(
    c("obs", "null", "list.index", "traits", "ind.plot", "sp") %in% names(res)
  ))
  expect_named(res$obs, c("mean(x, na.rm = TRUE)", "sd(x, na.rm = TRUE)"))
})

test_that("ComIndex with local nullmodel works", {
  res <- ComIndex(
    traits = traits.finch,
    index = "mean(x, na.rm = TRUE)",
    sp = sp.finch,
    nullmodels = "local",
    ind.plot = ind.plot.finch,
    nperm = 9,
    printprogress = FALSE
  )
  expect_s3_class(res, "ComIndex")
})

test_that("ComIndex with regional.pop nullmodel works", {
  res <- ComIndex(
    traits = traits.finch,
    index = "mean(x, na.rm = TRUE)",
    sp = sp.finch,
    nullmodels = "regional.pop",
    ind.plot = ind.plot.finch,
    nperm = 9,
    printprogress = FALSE
  )
  expect_s3_class(res, "ComIndex")
})

test_that("ComIndex with regional.pop.prab nullmodel works", {
  res <- ComIndex(
    traits = traits.finch,
    index = "mean(x, na.rm = TRUE)",
    sp = sp.finch,
    nullmodels = "regional.pop.prab",
    ind.plot = ind.plot.finch,
    nperm = 9,
    printprogress = FALSE
  )
  expect_s3_class(res, "ComIndex")
})

test_that("ComIndex with nperm = NULL returns no null model", {
  res <- ComIndex(
    traits = traits.finch,
    index = "mean(x, na.rm = TRUE)",
    sp = sp.finch,
    nullmodels = "regional.ind",
    ind.plot = ind.plot.finch,
    nperm = NULL,
    printprogress = FALSE
  )
  expect_s3_class(res, "ComIndex")
  expect_null(res$null)
})

test_that("ComIndex print and summary S3 methods do not error", {
  res <- ComIndex(
    traits = traits.finch,
    index = "mean(x, na.rm = TRUE)",
    sp = sp.finch,
    nullmodels = "regional.ind",
    ind.plot = ind.plot.finch,
    nperm = 9,
    printprogress = FALSE
  )
  expect_no_error(print(res))
  expect_no_error(summary(res))
})
