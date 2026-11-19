# Sample preparation: singletons, separation, collinearity.

# Iterated singleton drop, reghdfe5 style. With sep_fe = TRUE the outcome acts as an
# iweight (ppmlhdfe's trick): groups whose outcome sum is zero are dropped together with
# singleton groups, and both count as "singletons" (reghdfe5's num_singletons; Stata's
# ivppmlhdfe v0.9.4 reports them there). Returns the logical keep mask.
.drop_singletons <- function(fes, y, sep_fe = TRUE, alive = rep(TRUE, length(y))) {
  F <- length(fes)
  idx <- which(alive)
  nmax <- vapply(fes, function(g) max(g), 0L)
  remaining <- F
  i <- 0L
  ndrop <- 0L
  while (remaining > 0L) {
    i <- i + 1L
    g <- ((i - 1L) %% F) + 1L
    if (length(idx) < 2L) stop("ivppmlhdfe: fewer than 2 observations left after dropping singletons")
    codes <- fes[[g]][idx]
    cnt <- tabulate(codes, nmax[g])
    drop <- cnt[codes] == 1L
    if (sep_fe) {
      ys <- numeric(nmax[g])
      s <- rowsum(y[idx], codes, reorder = TRUE)
      ys[as.integer(rownames(s))] <- s
      drop <- drop | (ys[codes] == 0)
    }
    nd <- sum(drop)
    if (nd == 0L) {
      remaining <- remaining - 1L
    } else {
      remaining <- F - 1L
      idx <- idx[!drop]
      ndrop <- ndrop + nd
    }
  }
  keep <- logical(length(y))
  keep[idx] <- TRUE
  list(keep = keep, ndrop = ndrop)
}

# Two-part collinearity test (as the Julia backend v0.9.4, _select_not_collinear):
# (a) absorption: weighted norm of the partialled column > rtol * max raw weighted norm;
# (b) within-block rank: QR with limited pivoting (LINPACK dqrdc2), which keeps the columns
#     in their order and moves a column to the end only when it is (numerically) a linear
#     combination of the columns before it, so the FIRST of two collinear columns is kept,
#     as Stata's invsym-based rule (reghdfe_rmcoll) does.
# Returns the kept column indices (sorted).
.select_not_collinear <- function(A_partial, A_raw, w, rtol = 1e-9, rtol_rank = 1e-7) {
  K <- ncol(A_partial)
  if (K == 0L) return(integer(0))
  sw <- sqrt(pmax(w, 0))
  raw_n <- sqrt(.colSums((A_raw * sw)^2, nrow(A_raw), K))
  max_raw <- max(raw_n)
  if (max_raw == 0) return(integer(0))
  part_n <- sqrt(.colSums((A_partial * sw)^2, nrow(A_partial), K))
  nz <- which(part_n > rtol * max_raw)
  if (length(nz) <= 1L) return(nz)
  # rtol_rank = 1e-7 (R's lm default): the partial-out is iterative, so an exact linear relation
  # between columns survives only up to the inner tolerance (about 1e-9 relative)
  q <- qr(A_partial[, nz, drop = FALSE] * sw, tol = rtol_rank, LAPACK = FALSE)
  sort(nz[q$pivot[seq_len(q$rank)]])
}

# ReLU separation test (Correia, Guimaraes and Zylkin 2019; the method of ppmlhdfe):
# u = 1{y = 0}; repeat: weighted least squares of u on the regressors and the fixed
# effects with weight M = ceil(eps^-1/2) on y > 0 and 1 on y = 0 (least squares with
# equality constraints by weighting), xbd = fitted values, u <- max(xbd, 0) on y = 0;
# stop when no boundary fitted value is negative (separated: xbd > 0), or when no
# boundary residual is negative. Reimplemented from the paper and ppmlhdfe's description.
# Returns the indices (into y) of separated observations.
.relu_separation <- function(y, X, fes, tol = 1e-4, maxiter = 100L, verbose = 0L) {
  boundary <- which(y == 0)
  interior <- which(y > 0)
  if (!length(boundary) || ncol(X) == 0L) return(integer(0))
  N <- length(y)
  regtol <- max(tol^2, 1e-13)
  M <- ceiling(.Machine$double.eps^(-1 / 2))
  w <- rep(M, N)
  w[boundary] <- 1
  S <- .fe_setup(fes)
  .fe_weights(S, w)
  sw <- sqrt(w)
  norm1 <- function(A) if (nrow(A)) max(.colSums(abs(A), nrow(A), ncol(A))) else 0
  nu <- norm1(X[boundary, , drop = FALSE])
  nc <- norm1(X[interior, , drop = FALSE])
  u <- as.numeric(y == 0)
  utilde <- u; u_last <- u
  xtilde <- X
  sep <- integer(0)
  for (iter in seq_len(maxiter)) {
    if (iter > 1L) utilde <- u + utilde - u_last
    u_last <- u
    pr <- .fe_partial(S, cbind(utilde, xtilde), regtol)$E
    utilde <- pr[, 1L]
    xtilde <- pr[, -1L, drop = FALSE]
    q <- qr(xtilde * sw, LAPACK = TRUE)
    b <- qr.coef(q, utilde * sw)
    b[is.na(b)] <- 0
    resid <- as.vector(utilde - xtilde %*% b)
    yy <- sum(u * u)
    epsilon <- yy / M * max(if (nc > 0) nu / nc else 1, 1)^2
    epsilon <- max(epsilon, 1e-8)
    delta <- epsilon + tol
    xbd <- u - resid
    xbd[interior] <- 0
    xb <- xbd[boundary]
    xb[xb >= -0.1 * delta & xb <= delta] <- 0
    xbd[boundary] <- xb
    if (all(xb >= 0)) {
      sep <- boundary[xb > 0]
      break
    }
    resid[abs(resid) < 1e-8] <- 0
    if (min(resid[boundary]) >= 0) {
      rb <- resid[boundary]
      xb[rb > delta] <- 0
      sep <- boundary[xb > 0]
      break
    }
    if (iter == maxiter) {
      if (verbose > -1L) message("(ReLU separation check: maximum number of iterations reached; no observations dropped)")
      return(integer(0))
    }
    u[boundary] <- pmax(xb, 0)
  }
  if (length(sep) && verbose > -1L)
    message(sprintf("(ReLU method dropped %d separated observation%s in %d iterations)",
                    length(sep), if (length(sep) > 1L) "s" else "", iter))
  sep
}
