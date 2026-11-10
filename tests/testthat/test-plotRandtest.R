test_that("plotRandtest tests the right null values for Tstats and ComIndex", {
  set.seed(2)
  ind.plot <- factor(rep(paste0("P", 1:3), each = 12))
  sp <- factor(rep(c("a", "b", "c"), length.out = 36))
  traits <- matrix(rnorm(72), ncol = 2, dimnames = list(NULL, c("t1", "t2")))

  calls <- list()
  local_mocked_bindings(
    as.randtest = function(sim, obs, alter) {
      calls[[length(calls) + 1]] <<- list(sim = sim, obs = obs)
      ade4::as.randtest(sim, obs, alter)
    }
  )
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  ci <- ComIndex(
    traits,
    index = "tapply(x, ind.plot, mean)",
    nullmodels = "regional.ind",
    ind.plot = ind.plot,
    sp = sp,
    nperm = 9,
    printprogress = FALSE
  )
  nm <- ci$null[[1]]
  for (input in list(ci, as.listofindex(ci))) {
    calls <- list()
    plotRandtest(input)
    # trait 1: one call per community, then the mean over communities
    expect_equal(calls[[2]]$sim, nm[1, 2, ])
    expect_equal(calls[[2]]$obs, ci$obs[[1]][2, 1])
    expect_equal(calls[[4]]$sim, colMeans(nm[1, , ]))
    expect_equal(calls[[4]]$obs, mean(ci$obs[[1]][, 1]))
  }

  ts <- Tstats(traits, ind.plot, sp, nperm = 9, printprogress = FALSE)
  calls <- list()
  plotRandtest(ts)
  nm <- ts$Tstats$T_IP.IC_nm
  expect_equal(as.vector(calls[[3]]$sim), nm[, 1, 3])
  expect_equal(calls[[4]]$sim, rowMeans(nm[, 1, ]))
})
