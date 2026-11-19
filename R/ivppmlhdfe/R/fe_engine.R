# Fixed-effect engine: exact elimination of the FE with the most levels, then
# Jacobi-preconditioned conjugate gradients on the remaining system
#     S = D2' W M1 D2        (never formed)
# for all data columns in lockstep (one pass over the data serves every column).
# Port of the Julia sandbox solver schur.jl (2026-11-19), which follows fereg's
# fedemean (Schmieder 2026, NBER WP 36017, pp. 20-24).
#
# Everything is vectorised: group sums are column sums of a reshaped matrix when the
# groups are balanced (a gather by a precomputed permutation, then .colSums), and
# rowsum() otherwise; group means are gathered back by integer indexing. There is no
# loop over observations.

# ---------------------------------------------------------------------------------
# Grouping helper
# ---------------------------------------------------------------------------------

# codes: integer vector with values 1..n, every level present.
.mkgrp <- function(codes, n = max(codes)) {
  codes <- as.integer(codes)
  cnt <- tabulate(codes, n)
  if (any(cnt == 0L)) stop("internal: group codes are not dense")
  sorted <- !is.unsorted(codes)
  bal <- all(cnt == cnt[1L])
  perm <- if (bal && !sorted) order(codes, method = "radix") else NULL
  list(codes = codes, n = n, cnt = cnt, sorted = sorted, bal = bal, perm = perm,
       size = cnt[1L], N = length(codes))
}

# Group sums of a vector (length n result) or of the columns of a matrix (n x nc).
.gsum <- function(G, X) {
  if (is.null(dim(X))) {
    if (G$bal) {
      xs <- if (G$sorted) X else X[G$perm]
      return(.colSums(xs, G$size, G$n))
    }
    return(as.vector(rowsum(X, G$codes, reorder = TRUE)))
  }
  nc <- ncol(X)
  if (nc == 0L) return(matrix(0, G$n, 0L))
  if (G$bal) {
    xs <- if (G$sorted) X else X[G$perm, , drop = FALSE]
    return(matrix(.colSums(xs, G$size, G$n * nc), G$n, nc))
  }
  r <- rowsum(X, G$codes, reorder = TRUE)
  dimnames(r) <- NULL
  r
}

# Dense integer codes 1..n for a vector (order of first appearance is irrelevant).
.dense <- function(x) {
  u <- unique(x)
  match(x, u)
}

# Dense codes of the interaction of several vectors (a list of atomic vectors).
.interact <- function(cols) {
  if (length(cols) == 1L) return(.dense(cols[[1L]]))
  code <- .dense(cols[[1L]])
  for (j in seq_along(cols)[-1L]) {
    cj <- .dense(cols[[j]])
    key <- (code - 1) * max(cj) + cj          # double: no integer overflow
    code <- .dense(key)
  }
  code
}

# ---------------------------------------------------------------------------------
# Solver setup (structure once per sample) and weight refill (once per IRLS step)
# ---------------------------------------------------------------------------------

# fes: list of dense integer code vectors (one per FE), all of length N.
.fe_setup <- function(fes) {
  F <- length(fes)
  N <- length(fes[[1L]])
  nlev <- vapply(fes, function(g) max(g), 0L)
  big <- which.max(nlev)                       # ties: the first
  g1 <- .mkgrp(fes[[big]], nlev[big])
  rest <- seq_len(F)[-big]
  R <- length(rest)
  lev <- vector("list", R)        # per remaining FE: grouping object (codes 1..m_r)
  levo <- vector("list", R)       # codes offset into 1..m
  rows <- vector("list", R)       # rows of the m-vector that belong to FE r
  keys <- vector("list", R)       # (g1, level) combinations for the Jacobi diagonal
  off <- 0L
  for (j in seq_len(R)) {
    f <- rest[j]
    mr <- nlev[f]
    lev[[j]] <- .mkgrp(fes[[f]], mr)
    levo[[j]] <- fes[[f]] + off
    rows[[j]] <- off + seq_len(mr)
    key <- (g1$codes - 1) * mr + fes[[f]]
    ukey <- unique(key)
    if (length(ukey) == N) {
      # every (g1, level) cell holds one observation: C_{g,l} = w_i, summed by level directly
      keys[[j]] <- list(unique = TRUE)
    } else {
      kidx <- match(key, ukey)
      keys[[j]] <- list(unique = FALSE, grp = .mkgrp(kidx, length(ukey)),
                        g = as.integer((ukey - 1) %/% mr) + 1L,
                        l = .mkgrp(as.integer((ukey - 1) %% mr) + 1L, mr))
    }
    off <- off + mr
  }
  S <- new.env(parent = emptyenv())
  S$N <- N; S$F <- F; S$big <- big; S$g1 <- g1; S$G <- nlev[big]
  S$R <- R; S$lev <- lev; S$levo <- levo; S$rows <- rows; S$keys <- keys; S$m <- off
  S$nlev <- nlev
  S$passes <- 0L; S$restarts <- 0L; S$cg_iters <- 0L
  S
}

.fe_weights <- function(S, w) {
  S$w <- w
  n1 <- .gsum(S$g1, w)
  S$n1inv <- ifelse(n1 > 0, 1 / n1, 0)
  if (S$R > 0L) {
    m <- S$m
    d2 <- numeric(m)
    diagS <- numeric(m)
    n1inv_obs <- NULL
    for (j in seq_len(S$R)) {
      rj <- S$rows[[j]]
      d2[rj] <- .gsum(S$lev[[j]], w)
      k <- S$keys[[j]]
      # diag(S)_l = d2_l - sum_g C_{g,l}^2 / n1_g,  C_{g,l} = sum of w over obs in (g, l)
      if (k$unique) {
        if (is.null(n1inv_obs)) n1inv_obs <- S$n1inv[S$g1$codes]
        diagS[rj] <- -.gsum(S$lev[[j]], w * w * n1inv_obs)
      } else {
        C <- .gsum(k$grp, w)
        diagS[rj] <- -.gsum(k$l, C * C * S$n1inv[k$g])
      }
    }
    diagS <- diagS + d2
    S$pinv <- ifelse(diagS > 1e-12 * d2, 1 / diagS, 0)
    S$d2inv <- ifelse(d2 > 0, 1 / d2, 0)
  }
  invisible(S)
}

# ---------------------------------------------------------------------------------
# Operator and residual passes
# ---------------------------------------------------------------------------------

# Q = S P  (P: m x nc)
.fe_apply <- function(S, P) {
  U <- P[S$levo[[1L]], , drop = FALSE]
  if (S$R > 1L) for (j in 2:S$R) U <- U + P[S$levo[[j]], , drop = FALSE]
  sb <- .gsum(S$g1, U * S$w) * S$n1inv
  E <- U - .gexpand(S$g1, sb)
  WE <- E * S$w
  Q <- matrix(0, S$m, ncol(P))
  for (j in seq_len(S$R)) Q[S$rows[[j]], ] <- .gsum(S$lev[[j]], WE)
  S$passes <- S$passes + 1L
  Q
}

# E = M1 (V - D2 X);  Gr = D2' W E;  en2 = colSums(W E^2)
.fe_resid <- function(S, V, X = NULL, vnorm = FALSE) {
  U <- V
  if (!is.null(X)) for (j in seq_len(S$R)) U <- U - X[S$levo[[j]], , drop = FALSE]
  WU <- U * S$w
  vn2 <- if (vnorm) .colSums(WU * U, nrow(U), ncol(U)) else NULL
  sb <- .gsum(S$g1, WU) * S$n1inv
  E <- U - .gexpand(S$g1, sb)
  WE <- E * S$w
  en2 <- .colSums(WE * E, nrow(E), ncol(E))
  Gr <- NULL
  if (S$R > 0L) {
    Gr <- matrix(0, S$m, ncol(V))
    for (j in seq_len(S$R)) Gr[S$rows[[j]], ] <- .gsum(S$lev[[j]], WE)
  }
  S$passes <- S$passes + 1L
  list(E = E, Gr = Gr, en2 = en2, vn2 = vn2)
}

.cscale <- function(M, s) M * rep(s, each = nrow(M))

# Group means back to observations: sb[codes, ] (rep(each =) when the groups are sorted and balanced)
.gexpand <- function(G, sb) {
  if (G$sorted && G$bal) return(rep(as.vector(sb), each = G$size))
  sb[G$codes, , drop = FALSE]
}

# ---------------------------------------------------------------------------------
# Partial out every column of V (N x nc) under the current weights.
# Stopping rule per column: ||A'e||_2 <= tol * sqrt(F) * ||sqrt(W) e||_2 with A the
# column-normalised design (so ||A'e||^2 = sum_l Gr_l^2 / d2_l), plus fereg's
# certificate: a column that passes on the recursive CG residual must pass again on a
# freshly computed residual, else CG restarts from the current iterate.
# ---------------------------------------------------------------------------------
.fe_partial <- function(S, V, tol, maxiter = 10000L) {
  V <- as.matrix(V)
  storage.mode(V) <- "double"
  nc <- ncol(V)
  if (nc == 0L) return(list(E = V, iters = integer(0), conv = logical(0)))
  if (S$R == 0L) {
    r <- .fe_resid(S, V)
    return(list(E = r$E, iters = rep(0L, nc), conv = rep(TRUE, nc)))
  }
  m <- S$m
  thr2 <- (tol * sqrt(S$F))^2
  d2inv <- S$d2inv; pinv <- S$pinv
  gs <- function(A) .colSums(A * A * d2inv, m, ncol(A))
  X <- matrix(0, m, nc)
  r0 <- .fe_resid(S, V, NULL, vnorm = TRUE)
  B <- r0$Gr
  c0 <- r0$en2
  # absolute floor: a column (numerically) inside the span of the fixed effects stops once its
  # residual is below 1e-11 of the input's weighted norm (the relative rule alone never stops)
  fl2 <- (1e-11)^2 * r0$vn2
  Rr <- B
  P <- matrix(0, m, nc)
  iters <- integer(nc)
  active <- gs(Rr) > thr2 * c0 & c0 > fl2
  # a column that passes at X = 0 is final: its certificate would recompute exactly r0
  conv <- !active
  E <- r0$E
  rz <- numeric(nc)
  total <- 0L
  restarts <- 0L
  repeat {
    # (re)start CG for the active columns from the current X with residual Rr
    a <- which(active)
    if (length(a)) {
      Zp <- Rr[, a, drop = FALSE] * pinv
      P[, a] <- Zp
      rz[a] <- .colSums(Rr[, a, drop = FALSE] * Zp, m, length(a))
      bad <- a[!(rz[a] > 0)]
      if (length(bad)) { active[bad] <- FALSE; P[, bad] <- 0 }
    }
    while (any(active) && total < maxiter) {
      a <- which(active)
      Pa <- P[, a, drop = FALSE]
      Qa <- .fe_apply(S, Pa)
      total <- total + 1L
      pq <- .colSums(Pa * Qa, m, length(a))
      ok <- pq > 0
      if (!all(ok)) {
        bad <- a[!ok]
        active[bad] <- FALSE
        P[, bad] <- 0
        a <- a[ok]; Pa <- Pa[, ok, drop = FALSE]; Qa <- Qa[, ok, drop = FALSE]; pq <- pq[ok]
        if (!length(a)) next
      }
      alpha <- rz[a] / pq
      Xa <- X[, a, drop = FALSE] + .cscale(Pa, alpha)
      Ra <- Rr[, a, drop = FALSE] - .cscale(Qa, alpha)
      X[, a] <- Xa
      Rr[, a] <- Ra
      iters[a] <- iters[a] + 1L
      xb <- .colSums(Xa * (B[, a, drop = FALSE] + Ra), m, length(a))
      e2 <- pmax(c0[a] - xb, 0)
      done <- gs(Ra) <= thr2 * e2 | e2 <= fl2[a]
      if (any(done)) {
        dn <- a[done]
        active[dn] <- FALSE
        P[, dn] <- 0
      }
      keep <- !done
      if (any(keep)) {
        ak <- a[keep]
        Rk <- Ra[, keep, drop = FALSE]
        Zk <- Rk * pinv
        s <- .colSums(Rk * Zk, m, length(ak))
        beta <- s / rz[ak]
        rz[ak] <- s
        P[, ak] <- Zk + .cscale(P[, ak, drop = FALSE], beta)
      }
    }
    # certificate: recompute the residual and the gradient from scratch (columns not yet final)
    nconv <- which(!conv)
    if (!length(nconv)) break
    last <- .fe_resid(S, V[, nconv, drop = FALSE], X[, nconv, drop = FALSE])
    E[, nconv] <- last$E
    pass <- gs(last$Gr) <= thr2 * last$en2 | last$en2 <= fl2[nconv]
    conv[nconv[pass]] <- TRUE
    fail <- nconv[!pass]
    if (!length(fail) || total >= maxiter || restarts >= 50L) break
    Rr[, fail] <- last$Gr[, !pass, drop = FALSE]
    active[fail] <- TRUE
    restarts <- restarts + 1L
  }
  S$restarts <- S$restarts + restarts
  S$cg_iters <- S$cg_iters + total
  list(E = E, iters = iters, conv = conv)
}
