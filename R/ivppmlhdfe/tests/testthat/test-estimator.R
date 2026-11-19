test_that("PPML without IV equals a direct Poisson fit on dummies (coef and robust SE)", {
  d <- sim_panel(n = 10, T = 6, seed = 2)
  r <- ivppmlhdfe(y ~ x1 + x2 | id + year, d, tol = 1e-10, verbose = -1)
  g <- stats::glm(y ~ x1 + x2 + factor(id) + factor(year), family = stats::quasipoisson(), data = d,
                  control = stats::glm.control(epsilon = 1e-13, maxit = 100))
  expect_equal(unname(coef(r)), unname(coef(g)[c("x1", "x2")]), tolerance = 1e-8)
  Xf <- stats::model.matrix(g)
  mu <- stats::fitted(g)
  N <- nrow(d)
  Hi <- solve(crossprod(Xf, Xf * mu))
  sc <- Xf * (d$y - mu)
  V <- N / (N - 1) * Hi %*% crossprod(sc) %*% Hi
  expect_equal(unname(r$se), unname(sqrt(diag(V))[c("x1", "x2")]), tolerance = 1e-6)
  # no fixed effects at all
  r0 <- ivppmlhdfe(y ~ x1 + x2, d, tol = 1e-10, verbose = -1)
  g0 <- stats::glm(y ~ x1 + x2, family = stats::quasipoisson(), data = d,
                   control = stats::glm.control(epsilon = 1e-13, maxit = 100))
  expect_equal(unname(coef(r0)), unname(coef(g0)[-1]), tolerance = 1e-8)
  expect_equal(r0$constant, unname(coef(g0)[1]), tolerance = 1e-8)
})

test_that("just-identified IV equals 2SLS on dummies at the converged weights", {
  d <- sim_panel(n = 12, T = 6, seed = 4)
  r <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, d, tol = 1e-10, verbose = -1)
  mu <- predict(r, type = "mu")
  zz <- log(mu) - 1 + d$y / mu
  D <- stats::model.matrix(~ factor(id) + factor(year) - 1, d)
  Xf <- cbind(x2 = d$x2, x1 = d$x1, D)
  Zf <- cbind(x2 = d$x2, z = d$z, D)
  b <- solve(crossprod(Zf, Xf * mu), crossprod(Zf, zz * mu))
  expect_equal(unname(coef(r)), unname(b[1:2]), tolerance = 1e-7)
  # the additive moment E[q (y - mu)] = 0 holds at the solution
  mom <- crossprod(Zf, d$y - mu)
  expect_lt(max(abs(mom)) / sum(d$y), 1e-7)
})

test_that("overidentified IV: moment conditions of the 2SLS fixed point", {
  d <- sim_panel(n = 12, T = 6, seed = 8)
  r <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z + z2, d, tol = 1e-10, verbose = -1)
  mu <- predict(r, type = "mu")
  zz <- log(mu) - 1 + d$y / mu
  D <- stats::model.matrix(~ factor(id) + factor(year) - 1, d)
  Xf <- cbind(d$x2, d$x1, D)
  Zf <- cbind(d$x2, d$z, d$z2, D)
  Pz <- Zf %*% solve(crossprod(Zf, Zf * mu), crossprod(Zf, Xf * mu))
  b <- solve(crossprod(Pz, Xf * mu), crossprod(Pz, zz * mu))
  expect_equal(unname(coef(r)), b[1:2], tolerance = 1e-7)
})

test_that("accel on and off agree (coef 1e-8, SE 1e-6)", {
  d <- sim_panel(n = 9, T = 6, seed = 9, cls = "C")
  f <- y ~ x2 | exp^imp + exp^year + imp^year | x1 ~ z + z2
  a <- ivppmlhdfe(f, d, vcov = ~pair, accel = TRUE, verbose = -1)
  b <- ivppmlhdfe(f, d, vcov = ~pair, accel = FALSE, verbose = -1)
  expect_lt(max(abs(coef(a) - coef(b)) / pmax(abs(coef(b)), 1e-3)), 1e-8)
  expect_lt(max(abs(a$se - b$se) / b$se), 1e-6)
  for (s in c("w", "t", "r", "wt")) {
    c1 <- ivppmlhdfe(f, d, vcov = ~pair, accel = s, verbose = -1)
    expect_lt(max(abs(coef(c1) - coef(b)) / pmax(abs(coef(b)), 1e-3)), 1e-8)
  }
})

test_that("singletons and fixed-effect separation are dropped and counted as in Stata", {
  d <- sim_panel(n = 10, T = 5, seed = 11)
  extra <- d[1, ]; extra$id <- 999                     # a singleton id
  d <- rbind(d, extra)
  d$y[d$id == 3] <- 0                                   # an all-zero id: separated by its FE
  r <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, d, verbose = -1)
  expect_equal(r$num_singletons, 6L)
  expect_equal(r$num_separated, 0L)
  expect_equal(r$N, nrow(d) - 6L)
  expect_equal(sum(r$esample), nrow(d) - 6L)
})

test_that("frequency weights equal replicated data", {
  d <- sim_panel(n = 8, T = 5, seed = 12)
  rr <- d[rep(seq_len(nrow(d)), d$fw), ]
  a <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, d, weights = "fw", weight_type = "fweight", verbose = -1)
  b <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, rr, verbose = -1)
  expect_equal(a$N, sum(d$fw))
  expect_equal(coef(a), coef(b), tolerance = 1e-8)
  expect_equal(a$se, b$se, tolerance = 1e-6)
})

test_that("cluster, two-way cluster, pweights, offset and exposure run; predict is consistent", {
  d <- sim_panel(n = 7, T = 5, seed = 13, cls = "B")
  f <- y ~ x2 | exp^year + imp^year | x1 ~ z
  r1 <- ivppmlhdfe(f, d, vcov = ~pair, weights = "w", verbose = -1)
  r2 <- ivppmlhdfe(f, d, vcov = ~exp + imp, verbose = -1)
  expect_equal(r2$N_clust, c(7L, 7L))
  expect_true(all(is.finite(r1$se)) && all(is.finite(r2$se)))
  d$expo <- exp(d$off)
  r3 <- ivppmlhdfe(f, d, offset = "off", verbose = -1)
  r4 <- ivppmlhdfe(f, d, exposure = "expo", verbose = -1)
  expect_equal(coef(r3), coef(r4), tolerance = 1e-10)
  eta <- predict(r3, type = "eta")
  expect_equal(eta, predict(r3, type = "xb") + predict(r3, type = "d") + d$off, tolerance = 1e-10)
  expect_equal(predict(r3, type = "mu"), exp(eta), tolerance = 1e-12)
  expect_output(print(summary(r3)), "x1")
  expect_equal(nobs(r3), nrow(d))
  expect_equal(dim(vcov(r3)), c(2L, 2L))
})

test_that("collinear regressors and instruments are dropped", {
  d <- sim_panel(n = 8, T = 5, seed = 14)
  d$x2b <- 2 * d$x2
  d$zb <- d$z + 0 * d$x1
  d$zc <- 3 * d$z
  r <- ivppmlhdfe(y ~ x2 + x2b | id + year | x1 ~ z + zc, d, verbose = -1)
  expect_equal(names(coef(r)), c("x2", "x1"))
  expect_equal(r$instruments, "z")
  # a regressor constant within id is absorbed by the id fixed effect
  d$xid <- rnorm(8)[d$id]
  d$z3 <- d$z + d$z2
  t0 <- proc.time()[["elapsed"]]
  r2 <- ivppmlhdfe(y ~ x2 + xid | id + year | x1 ~ z + z2 + z3, d, verbose = -1)
  expect_lt(proc.time()[["elapsed"]] - t0, 10)
  expect_equal(names(coef(r2)), c("x2", "x1"))
  expect_equal(r2$instruments, c("z", "z2"))
})
