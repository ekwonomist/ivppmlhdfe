test_that("one fixed effect: the engine is the exact weighted within transform", {
  set.seed(3)
  N <- 300
  g <- ivppmlhdfe:::.dense(sample(1:17, N, TRUE))
  w <- runif(N, 0.2, 3)
  V <- cbind(rnorm(N), runif(N), rexp(N))
  S <- ivppmlhdfe:::.fe_setup(list(g))
  ivppmlhdfe:::.fe_weights(S, w)
  E <- ivppmlhdfe:::.fe_partial(S, V, 1e-12)$E
  expected <- V - apply(V, 2, function(v) ave(w * v, g, FUN = sum) / ave(w, g, FUN = sum))
  expect_lt(max(abs(E - expected)), 1e-12)
})

test_that("three fixed effects: Schur-PCG residuals equal the dense dummy regression", {
  d <- sim_panel(n = 6, T = 4, seed = 5, cls = "C")
  fes <- list(ivppmlhdfe:::.interact(list(d$exp, d$imp)),
              ivppmlhdfe:::.interact(list(d$exp, d$year)),
              ivppmlhdfe:::.interact(list(d$imp, d$year)))
  w <- d$w * exp(rnorm(nrow(d)))
  V <- cbind(d$x1, d$x2, d$z, log(d$y + 1))
  S <- ivppmlhdfe:::.fe_setup(fes)
  ivppmlhdfe:::.fe_weights(S, w)
  E <- ivppmlhdfe:::.fe_partial(S, V, 1e-14)$E
  R <- dummy_resid(V, fes, w)
  expect_lt(max(abs(E - R)) / max(abs(V)), 1e-9)
})

test_that("unbalanced groups (rowsum path) agree with the dummy regression", {
  d <- sim_panel(n = 7, T = 5, seed = 6, cls = "B")
  d <- d[-c(3, 17, 40, 41, 99), ]
  fes <- list(ivppmlhdfe:::.interact(list(d$exp, d$year)), ivppmlhdfe:::.interact(list(d$imp, d$year)))
  w <- runif(nrow(d), 0.1, 5)
  V <- cbind(d$x1, d$z)
  S <- ivppmlhdfe:::.fe_setup(fes)
  ivppmlhdfe:::.fe_weights(S, w)
  E <- ivppmlhdfe:::.fe_partial(S, V, 1e-14)$E
  expect_lt(max(abs(E - dummy_resid(V, fes, w))) / max(abs(V)), 1e-9)
})

test_that("warm-start identity: M_w2 (M_w1 v) = M_w2 v", {
  d <- sim_panel(n = 6, T = 5, seed = 7, cls = "C")
  fes <- list(ivppmlhdfe:::.interact(list(d$exp, d$imp)),
              ivppmlhdfe:::.interact(list(d$exp, d$year)),
              ivppmlhdfe:::.interact(list(d$imp, d$year)))
  V <- cbind(d$x1, d$z)
  S <- ivppmlhdfe:::.fe_setup(fes)
  ivppmlhdfe:::.fe_weights(S, runif(nrow(d), 0.5, 2))
  E1 <- ivppmlhdfe:::.fe_partial(S, V, 1e-14)$E
  ivppmlhdfe:::.fe_weights(S, runif(nrow(d), 0.5, 2))
  a <- ivppmlhdfe:::.fe_partial(S, E1, 1e-14)$E
  b <- ivppmlhdfe:::.fe_partial(S, V, 1e-14)$E
  expect_lt(max(abs(a - b)) / max(abs(V)), 1e-10)
})
