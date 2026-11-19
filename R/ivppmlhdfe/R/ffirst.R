# First-stage statistics (ffirst = TRUE), following the conventions pinned by the Stata track
# (v010_build/stata/SPEC_ffirst.md): the numbers ivreg2/ivreghdfe's `ffirst` reports for the
# linear IV problem at the converged IRLS weights omega = w * mu, with the fixed effects
# partialled out under omega. Reimplemented from the formulas (Kleibergen and Paap 2006;
# Sanderson and Windmeijer 2016; Cragg and Donald 1993), not from GPL code.

# Number of connected components of the bipartite graph whose nodes are the levels of two
# fixed effects and whose edges are the observations (min-label propagation).
.group_min <- function(v, g, n) {
  o <- order(g, v)
  go <- g[o]
  first <- !duplicated(go)
  r <- numeric(n)
  r[go[first]] <- v[o][first]
  r
}
.n_components <- function(a, b) {
  na <- max(a); nb <- max(b)
  key <- unique((a - 1) * nb + b)
  ea <- as.integer((key - 1) %/% nb) + 1L
  eb <- as.integer((key - 1) %% nb) + 1L
  la <- seq_len(na) + 0
  repeat {
    lb <- .group_min(la[ea], eb, nb)
    la2 <- pmin(.group_min(lb[eb], ea, na), la)
    if (identical(la2, la)) break
    la <- la2
  }
  length(unique(la))
}

# Absorbed degrees of freedom, reghdfe's rule (SPEC_ffirst section 4).
.df_absorb <- function(fes, clusters = NULL, has_fe = TRUE) {
  F <- length(fes)
  G <- vapply(fes, max, 0L)
  M <- integer(F)
  nested <- logical(F)
  if (!is.null(clusters) && length(clusters)) {
    for (g in seq_len(F)) {
      for (cl in clusters) {
        if (length(unique((fes[[g]] - 1) * (max(cl) + 1) + cl)) == G[g]) { nested[g] <- TRUE; break }
      }
    }
    M[nested] <- G[nested]
  }
  idx <- which(!nested)
  if (length(idx) > 1L) M[idx[-1L]] <- 1L
  if (any(nested) && length(idx) >= 1L) M[idx] <- 1L
  if (length(idx) > 1L) {
    for (i in seq_len(length(idx) - 1L)) for (j in (i + 1L):length(idx)) {
      g <- idx[i]; h <- idx[j]
      M[h] <- max(M[h], .n_components(fes[[g]], fes[[h]]))
    }
  }
  sum(G - M)
}

.msqrt <- function(A) {
  e <- eigen(A, symmetric = TRUE)
  ev <- pmax(e$values, 0)
  e$vectors %*% (sqrt(ev) * t(e$vectors))
}

.ffirst <- function(Xdm, Zdm, w, K_exog, endog, inst, N, N_eff, vtype, clusters, is_fw,
                    w_user, mu, S, fes, has_fe = TRUE) {
  Kw <- K_exog; Ke <- length(endog); Lz <- length(inst); L <- Kw + Lz
  if (Ke == 0L || Lz == 0L) return(NULL)
  Wx <- Xdm[, seq_len(Kw), drop = FALSE]
  Xe <- Xdm[, Kw + seq_len(Ke), drop = FALSE]
  Z <- Zdm[, Kw + seq_len(Lz), drop = FALSE]
  R <- cbind(Z, Wx)
  om <- w
  psi <- if (is_fw) om^2 / w_user else om^2
  Nn <- N_eff
  cl <- identical(vtype, "cluster") && !is.null(clusters) && length(clusters)
  df_a <- if (has_fe) .df_absorb(fes, if (cl) clusters else NULL) else 0L
  sdof <- if (df_a > 0) df_a else 1
  terms <- list()
  G <- NA_real_
  if (cl) {
    nC <- length(clusters)
    for (mask in seq_len(2L^nC - 1L)) {
      members <- which(bitwAnd(mask, 2L^(seq_len(nC) - 1L)) > 0L)
      cid <- if (length(members) == 1L) clusters[[members]] else .interact(clusters[members])
      terms[[length(terms) + 1L]] <- list(sign = if (length(members) %% 2L == 1L) 1 else -1, cid = cid)
    }
    G <- min(vapply(clusters, function(c) length(unique(c)), 0L))
  }
  meat <- function(Sc) {
    if (!cl) return(crossprod(Sc * psi, Sc))
    WS <- Sc * om
    M <- matrix(0, ncol(Sc), ncol(Sc))
    for (t in terms) M <- M + t$sign * crossprod(rowsum(WS, t$cid, reorder = FALSE))
    M
  }
  wls <- function(Rm, y) {
    Rw <- Rm * om
    RtR <- crossprod(Rw, Rm)
    b <- solve(RtR, crossprod(Rw, y))
    list(b = b, e = y - Rm %*% b, RtR = RtR)
  }
  wald_Z <- function(e, RtR, b) {
    A <- solve(RtR)
    V <- A %*% meat(R * as.vector(e)) %*% A
    bz <- b[seq_len(Lz)]
    sum(bz * solve(V[seq_len(Lz), seq_len(Lz), drop = FALSE], bz))
  }
  to_F <- function(W0, df1) {
    if (cl) return(c(F = W0 / (Nn - 1) * (Nn - L - sdof) * (G - 1) / G / df1, df2 = G - 1))
    df2 <- Nn - L - sdof
    c(F = W0 / Nn * df2 / df1, df2 = df2)
  }
  first <- vector("list", Ke)
  names(first) <- endog
  Pi <- matrix(0, L, Ke)
  for (k in seq_len(Ke)) {
    f <- wls(R, Xe[, k])
    Pi[, k] <- f$b
    Fk <- to_F(wald_Z(f$e, f$RtR, f$b), Lz)
    er <- if (Kw) wls(Wx, Xe[, k])$e else Xe[, k]
    pr2 <- 1 - sum(om * f$e^2) / sum(om * er^2)
    first[[k]] <- c(F = Fk[["F"]], df1 = Lz, df2 = Fk[["df2"]],
                    p = stats::pf(Fk[["F"]], Lz, Fk[["df2"]], lower.tail = FALSE), partialR2 = pr2)
  }
  Xhat <- R %*% Pi
  for (k in seq_len(Ke)) {
    o <- first[[k]]
    if (Ke == 1L) {
      first[[k]] <- c(o, SWF = o[["F"]], SWdf1 = o[["df1"]], SWdf2 = o[["df2"]], SWp = o[["p"]])
      next
    }
    oth <- setdiff(seq_len(Ke), k)
    Sh <- cbind(Xhat[, oth, drop = FALSE], Wx)
    Sa <- cbind(Xe[, oth, drop = FALSE], Wx)
    Shw <- Sh * om
    b1 <- solve(crossprod(Shw, Sh), crossprod(Shw, Xe[, k]))
    e1 <- Xe[, k] - Sa %*% b1
    f2 <- wls(R, e1)
    df1 <- Lz - Ke + 1
    Fs <- to_F(wald_Z(f2$e, f2$RtR, f2$b), df1)
    first[[k]] <- c(o, SWF = Fs[["F"]], SWdf1 = df1, SWdf2 = Fs[["df2"]],
                    SWp = stats::pf(Fs[["F"]], df1, Fs[["df2"]], lower.tail = FALSE))
  }
  # system statistics: W partialled out of the endogenous regressors and the instruments
  if (Kw) {
    Y <- wls(Wx, Xe)$e
    Zb <- wls(Wx, Z)$e
  } else {
    Y <- Xe; Zb <- Z
  }
  Zbw <- Zb * om
  A <- crossprod(Zbw, Zb)
  B <- crossprod(Y * om, Y)
  C <- crossprod(Zbw, Y)
  PiZ <- solve(A, C)
  LA <- t(chol(A))
  LB <- t(chol(B))
  Theta <- t(LA) %*% PiZ %*% solve(t(LB))
  sv <- svd(Theta, nu = Lz, nv = Ke)
  U <- sv$u; V <- sv$v; s <- sv$d
  kk <- Ke - 1L
  ru <- (kk + 1L):Lz; rv <- (kk + 1L):Ke
  u12 <- U[seq_len(kk), ru, drop = FALSE]; u22 <- U[ru, ru, drop = FALSE]
  v12 <- V[seq_len(kk), rv, drop = FALSE]; v22 <- V[rv, rv, drop = FALSE]
  Aperp <- rbind(u12, u22) %*% solve(u22) %*% .msqrt(u22 %*% t(u22))
  Bperp <- .msqrt(v22 %*% t(v22)) %*% solve(t(v22)) %*% t(rbind(v12, v22))
  Tm <- kronecker(Bperp, t(Aperp))
  lam <- Tm %*% as.vector(Theta)
  J <- kronecker(solve(LB), solve(LA))
  rk <- function(res) {
    Sc <- do.call(cbind, lapply(seq_len(Ke), function(k) Zb * res[, k]))
    Vt <- Tm %*% J %*% meat(Sc) %*% t(J) %*% t(Tm)
    sum(lam * solve(Vt, lam))
  }
  idstat <- rk(Y)
  rk_w <- rk(Y - Zb %*% PiZ)
  iddf <- Lz - Ke + 1
  widstat <- if (cl) rk_w / (Nn - 1) * (Nn - L - sdof) * (G - 1) / G / Lz else rk_w / Nn * (Nn - L - sdof) / Lz
  cd <- min(s^2 / (1 - s^2))
  cdf <- cd * (Nn - sdof - L) / Lz
  fm <- do.call(rbind, first)
  rownames(fm) <- endog
  list(first = first, first_matrix = fm,
       idstat = idstat, iddf = iddf, idp = stats::pchisq(idstat, iddf, lower.tail = FALSE),
       cdf = cdf, widstat = widstat, ccorr = s, df_a = df_a, sdof = sdof,
       N = Nn, G = if (cl) G else NA_real_, L = L, Ke = Ke, Lz = Lz)
}

print_ffirst <- function(ff, digits = 4L) {
  cat("\nFirst-stage statistics (excluded instruments):\n")
  fm <- ff$first_matrix
  show <- cbind(F = fm[, "F"], df1 = fm[, "df1"], df2 = fm[, "df2"], `Pr(>F)` = fm[, "p"],
                `partial R2` = fm[, "partialR2"], `SW F` = fm[, "SWF"], `SW Pr(>F)` = fm[, "SWp"])
  print(signif(show, digits))
  cat(sprintf("Underidentification (Kleibergen-Paap rk LM): %.4g, chi2(%d), p = %.4g\n",
              ff$idstat, as.integer(ff$iddf), ff$idp))
  cat(sprintf("Weak identification: Cragg-Donald Wald F = %.4g, Kleibergen-Paap rk Wald F = %.4g\n",
              ff$cdf, ff$widstat))
  invisible(ff)
}
