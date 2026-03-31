data(finch.ind)

test_that("IndexByGroups wraps each metric in tapply", {
  result <- IndexByGroups(c("mean(x)", "sd(x)"), "pop")
  expect_type(result, "character")
  expect_length(result, 2)
  expect_true(all(grepl("tapply", result)))
  expect_true(all(grepl("pop", result)))
})

test_that("IndexByGroups works with a single metric", {
  result <- IndexByGroups("mean(x, na.rm = TRUE)", "site")
  expect_length(result, 1)
  expect_match(result, "tapply")
})

test_that("samplingSubsetData returns list with nperm elements", {
  res <- samplingSubsetData(
    d = traits.finch,
    sampUnit = sp.finch,
    nperm = 3,
    prop = c(50, 100)
  )
  expect_true(is.list(res))
  expect_length(res, 3)
})

test_that("samplingSubsetData inner lists match number of proportions", {
  res <- samplingSubsetData(
    d = traits.finch,
    sampUnit = sp.finch,
    nperm = 2,
    prop = c(25, 50, 75)
  )
  expect_length(res[[1]], 3)
})

test_that("samplingSubsetData with type = 'count' works", {
  res <- samplingSubsetData(
    d = traits.finch,
    sampUnit = sp.finch,
    nperm = 2,
    type = "count",
    Size = c(5, 10)
  )
  expect_true(is.list(res))
})
