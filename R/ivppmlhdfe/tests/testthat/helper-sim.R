# Small simulated gravity-style panels for the tests.
sim_panel <- function(n = 12, T = 6, seed = 1, cls = "A", zeros = 0) {
  set.seed(seed)
  if (cls == "A") {
    d <- expand.grid(year = seq_len(T), id = seq_len(n))
    a <- 0.3 * rnorm(n)[d$id] + 0.2 * rnorm(T)[d$year]
  } else {
    d <- expand.grid(year = seq_len(T), imp = seq_len(n), exp = seq_len(n))
    d <- d[d$exp != d$imp, ]
    d$pair <- (d$exp - 1) * n + d$imp
    ot <- matrix(0.3 * rnorm(n * T), n, T)
    dt <- matrix(0.3 * rnorm(n * T), n, T)
    a <- ot[cbind(d$exp, d$year)] + dt[cbind(d$imp, d$year)]
    if (cls == "C") a <- a + 0.5 * rnorm(n * n)[d$pair]
  }
  N <- nrow(d)
  d$z <- rnorm(N)
  d$z2 <- rnorm(N)
  d$x2 <- 0.5 * rnorm(N)
  e <- rnorm(N)
  d$x1 <- 0.8 * d$z + 0.3 * d$z2 + e + 0.4 * a
  d$y <- exp(0.5 * d$x1 + 0.3 * d$x2 + a - 0.5 + 0.5 * e + sqrt(0.75) * rnorm(N))
  if (zeros > 0) d$y[runif(N) < zeros] <- 0
  d$w <- runif(N, 0.5, 2)
  d$fw <- 1 + (seq_len(N) %% 3)
  d$off <- 0.2 * rnorm(N)
  rownames(d) <- NULL
  d
}

# Weighted least squares residuals of the columns of V on a dense dummy design.
dummy_resid <- function(V, fe_list, w) {
  D <- do.call(cbind, lapply(fe_list, function(g) stats::model.matrix(~ factor(g) - 1)))
  apply(V, 2, function(v) stats::lm.wfit(D, v, w)$residuals)
}
