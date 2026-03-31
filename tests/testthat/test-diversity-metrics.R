data(finch.ind)

# ── NND metrics ──────────────────────────────────────────────────────────────

test_that("CVNND works on a single trait vector", {
  val <- CVNND(traits.finch[, 1], na.rm = TRUE)
  expect_length(val, 1)
  expect_true(is.numeric(val) && val >= 0)
})

test_that("CVNND with div_range = TRUE returns numeric", {
  val <- CVNND(traits.finch[, 1], div_range = TRUE, na.rm = TRUE)
  expect_length(val, 1)
  expect_true(is.numeric(val))
})

test_that("CVNND works on a multi-trait matrix", {
  val <- CVNND(traits.finch, na.rm = TRUE)
  expect_length(val, 1)
  expect_true(is.numeric(val) && val >= 0)
})

test_that("CVNND with scale.tr = FALSE works", {
  val <- CVNND(traits.finch, scale.tr = FALSE, na.rm = TRUE)
  expect_length(val, 1)
  expect_true(is.numeric(val))
})

test_that("MNND returns a non-negative numeric scalar", {
  val <- MNND(traits.finch[, 1], na.rm = TRUE)
  expect_length(val, 1)
  expect_true(is.numeric(val) && val >= 0)
})

test_that("MinNND returns a non-negative numeric scalar", {
  val <- MinNND(traits.finch[, 1], na.rm = TRUE)
  expect_length(val, 1)
  expect_true(is.numeric(val) && val >= 0)
})

test_that("SDNND returns a non-negative numeric scalar", {
  val <- SDNND(traits.finch[, 1], na.rm = TRUE)
  expect_length(val, 1)
  expect_true(is.numeric(val) && val >= 0)
})

test_that("SDND returns a numeric scalar", {
  val <- SDND(traits.finch[, 1], na.rm = TRUE)
  expect_length(val, 1)
  expect_true(is.numeric(val))
})

test_that("MND returns a numeric scalar", {
  val <- MND(traits.finch[, 1], na.rm = TRUE)
  expect_length(val, 1)
  expect_true(is.numeric(val))
})

# ── SumBL ─────────────────────────────────────────────────────────────────────

test_that("SumBL with gower.dist = TRUE returns non-negative numeric", {
  tr <- na.omit(traits.finch[1:100, 1, drop = FALSE])
  val <- SumBL(tr)
  expect_length(val, 1)
  expect_true(is.numeric(val) && val >= 0)
})

test_that("SumBL with gower.dist = FALSE returns non-negative numeric", {
  tr <- na.omit(traits.finch[1:100, 1, drop = FALSE])
  val <- SumBL(tr, gower.dist = FALSE)
  expect_length(val, 1)
  expect_true(is.numeric(val) && val >= 0)
})

test_that("SumBL with different hclust methods works", {
  tr <- na.omit(traits.finch[1:100, 1, drop = FALSE])
  for (m in c("average", "complete", "single")) {
    expect_no_error(SumBL(tr, gower.dist = TRUE, method.hclust = m))
  }
})

# ── MinMaxMST ─────────────────────────────────────────────────────────────────

test_that("MinMaxMST with gower.dist = TRUE returns value in (0, 1]", {
  tr <- na.omit(traits.finch[1:50, 1, drop = FALSE])
  val <- MinMaxMST(tr)
  expect_length(val, 1)
  expect_true(is.numeric(val) && val > 0 && val <= 1)
})

test_that("MinMaxMST with gower.dist = FALSE works", {
  tr <- na.omit(traits.finch[1:50, 1, drop = FALSE])
  val <- MinMaxMST(tr, gower.dist = FALSE)
  expect_length(val, 1)
  expect_true(is.numeric(val))
})

test_that("MinMaxMST with scale.tr = FALSE works", {
  tr <- na.omit(traits.finch[1:50, 1, drop = FALSE])
  val <- MinMaxMST(tr, gower.dist = FALSE, scale.tr = FALSE)
  expect_length(val, 1)
  expect_true(is.numeric(val))
})
