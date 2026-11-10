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

test_that("regional.pop null values of a trait come from that trait only", {
  set.seed(1)
  tr <- cbind(t1 = runif(60), t2 = runif(60, 100, 101))
  plot <- factor(rep(c("P1", "P2", "P3"), each = 20))
  sp <- factor(rep(c("a", "b", "c", "d"), 15))
  res <- ComIndex(
    traits = tr,
    index = "mean(x)",
    sp = sp,
    nullmodels = "regional.pop",
    ind.plot = plot,
    nperm = 9,
    printprogress = FALSE
  )
  expect_true(all(res$null[[1]][1, , ] <= 1))
  expect_true(all(res$null[[1]][2, , ] >= 100))
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

test_that("ComIndexMulti null values of one trait use only that trait (GH-13)", {
  set.seed(1)
  tr <- cbind(t1 = runif(60), t2 = runif(60, 100, 101))
  plot <- factor(rep(c("P1", "P2", "P3"), each = 20))
  sp <- factor(rep(c("a", "b", "c", "d"), 15))
  for (nm in c("local", "regional.ind")) {
    res <- ComIndexMulti(
      tr,
      index = c("mean(x[, 1])", "mean(x[, 2])"),
      by.factor = plot,
      sp = sp,
      ind.plot = plot,
      nullmodels = c(nm, nm),
      nperm = 9,
      printprogress = FALSE
    )
    expect_true(all(res$null[[1]] >= 0 & res$null[[1]] <= 1))
    expect_true(all(res$null[[2]] >= 100 & res$null[[2]] <= 101))
  }
})

test_that("ComIndex and ComIndexMulti do not depend on row order (GH-15)", {
  set.seed(1)
  plot <- factor(rep(c("P1", "P2", "P3"), each = 20))
  sp <- factor(rep(c("a", "b", "c", "d"), 15))
  tr <- cbind(
    t1 = c(a = 10, b = 20, c = 30, d = 40)[as.character(sp)] + rnorm(60),
    t2 = rnorm(60)
  )
  ord <- sample(60)
  run <- function(i, nm, index = "tapply(x, sp, mean)") {
    ComIndex(
      tr[i, ],
      index = index,
      sp = sp[i],
      ind.plot = plot[i],
      nullmodels = nm,
      nperm = 2,
      printprogress = FALSE
    )
  }
  expect_equal(run(ord, "local")$obs, run(1:60, "local")$obs)
  # regional.pop works on population means, x is one value per population
  expect_equal(
    run(ord, "regional.pop", "mean(x)")$obs,
    run(1:60, "regional.pop", "mean(x)")$obs
  )
  res <- run(ord, "local")
  expect_equal(as.character(res$sp), as.character(sp[ord][order(plot[ord])]))

  pop <- factor(paste(sp, plot, sep = "_"))
  multi <- function(i) {
    ComIndexMulti(
      tr[i, ],
      index = "mean(x[, 1])",
      by.factor = pop[i],
      sp = sp[i],
      ind.plot = plot[i],
      nullmodels = "local",
      nperm = 2,
      printprogress = FALSE
    )
  }
  expect_equal(multi(ord)$obs, multi(1:60)$obs)
})

test_that("ComIndexMulti works with a single trait and with regional.pop.prab", {
  set.seed(1)
  tr <- cbind(t1 = runif(60), t2 = runif(60, 100, 101))
  plot <- factor(rep(c("P1", "P2", "P3"), each = 20))
  sp <- factor(rep(c("a", "b", "c", "d"), 15))

  r1 <- ComIndexMulti(
    tr[, 1, drop = FALSE],
    index = "mean(x[, 1])",
    by.factor = plot,
    sp = sp,
    ind.plot = plot,
    nullmodels = "local",
    nperm = 5,
    printprogress = FALSE
  )
  expect_equal(r1$obs[[1]], as.vector(tapply(tr[, 1], plot, mean)))
  expect_equal(r1$null[[1]][, 1], r1$obs[[1]])
  for (nm in c("regional.ind", "regional.pop", "regional.pop.prab")) {
    expect_no_error(ComIndexMulti(
      tr[, 1, drop = FALSE],
      index = "mean(x[, 1])",
      by.factor = plot,
      sp = sp,
      ind.plot = plot,
      nullmodels = nm,
      nperm = 5,
      printprogress = FALSE
    ))
  }

  r2 <- ComIndexMulti(
    tr,
    index = "mean(x[, 1])",
    by.factor = plot,
    sp = sp,
    ind.plot = plot,
    nullmodels = "regional.pop.prab",
    nperm = 20,
    printprogress = FALSE
  )
  pop_mean <- tapply(tr[, 1], paste(sp, plot), mean)
  expect_equal(
    r2$obs[[1]],
    as.vector(tapply(pop_mean, sub(".* ", "", names(pop_mean)), mean))
  )
  expect_true(all(
    r2$null[[1]] >= min(pop_mean) & r2$null[[1]] <= max(pop_mean)
  ))
})

test_that("ComIndex works with a single trait", {
  set.seed(1)
  tr <- cbind(t1 = runif(60))
  plot <- factor(rep(c("P1", "P2", "P3"), each = 20))
  sp <- factor(rep(c("a", "b", "c", "d"), 15))
  r <- ComIndex(
    tr,
    index = "tapply(x, ind.plot, mean)",
    nullmodels = "local",
    ind.plot = plot,
    sp = sp,
    nperm = 5,
    printprogress = FALSE
  )
  expect_equal(as.vector(r$obs[[1]]), as.vector(tapply(tr[, 1], plot, mean)))
})

test_that("regional.pop and regional.pop.prab group population values with ind.plot and sp", {
  set.seed(1)
  tr <- cbind(t1 = runif(60), t2 = runif(60))
  plot <- factor(rep(c("P1", "P2", "P3"), each = 20))
  sp <- factor(rep(c("a", "b", "c", "d"), 15))
  pop_mean <- apply(tr, 2, \(x) tapply(x, paste(sp, plot), mean))
  pop_plot <- sub(".* ", "", rownames(pop_mean))
  expected <- apply(pop_mean, 2, \(x) tapply(x, pop_plot, mean))

  for (nm in c("regional.pop", "regional.pop.prab")) {
    r <- ComIndex(
      tr,
      index = c("tapply(x, ind.plot, mean)", "tapply(x, sp, length)"),
      nullmodels = c(nm, nm),
      ind.plot = plot,
      sp = sp,
      nperm = 5,
      printprogress = FALSE
    )
    expect_equal(unname(r$obs[[1]]), unname(expected))
    expect_equal(dim(r$null[[1]]), c(2, 3, 5))
    expect_true(all(r$obs[[2]] == 3))
  }
})
