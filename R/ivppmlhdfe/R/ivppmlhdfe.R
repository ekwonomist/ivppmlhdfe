# ivppmlhdfe: IV-PPML with high-dimensional fixed effects (R port of the Stata command
# ivppmlhdfe v0.9.4 with the v1.0.0 fast path). See ?ivppmlhdfe.

.EPS1 <- .Machine$double.eps          # Mata epsilon(1)
.EPS100 <- 2^-46                       # Mata epsilon(100): spacing of doubles at 100

.parse_accel <- function(accel) {
  if (isTRUE(accel)) return(list(w = TRUE, t = TRUE, r = TRUE, label = "wtr"))
  if (isFALSE(accel)) return(list(w = FALSE, t = FALSE, r = FALSE, label = "none"))
  if (is.character(accel)) {
    a <- tolower(paste(accel, collapse = ""))
    bad <- setdiff(strsplit(gsub("[ ,]", "", a), "")[[1L]], c("w", "t", "r", "s"))
    if (length(bad)) stop("ivppmlhdfe: accel must be TRUE, FALSE or a combination of 'w', 't', 'r'")
    out <- list(w = grepl("w", a), t = grepl("t", a), r = grepl("r", a))
    out$label <- paste0(c("w", "t", "r")[unlist(out[c("w", "t", "r")])], collapse = "")
    if (!nzchar(out$label)) out$label <- "none"
    return(out)
  }
  stop("ivppmlhdfe: accel must be TRUE, FALSE or a combination of 'w', 't', 'r'")
}

.parse_sep <- function(separation) {
  if (is.null(separation) || !length(separation)) return(character(0))
  s <- tolower(separation)
  if (any(s %in% c("none", "no", "off"))) return(character(0))
  if (any(s %in% c("default", "standard", "on", "auto"))) s <- c("fe", "relu")
  if (any(s %in% c("all", "full"))) s <- c("fe", "relu", "mu")
  s[s == "ir"] <- "relu"
  s <- setdiff(s, "simplex")           # subsumed by relu
  bad <- setdiff(s, c("fe", "relu", "mu"))
  if (length(bad)) stop("ivppmlhdfe: unknown separation method(s): ", paste(bad, collapse = ", "),
                        " (valid: fe, relu, mu)")
  unique(s)
}

# Weighted 2SLS on partialled data (no constant).
.tsls <- function(zdm, Xdm, Zdm, w) {
  Zw <- Zdm * w
  Pi <- solve(crossprod(Zw, Zdm), crossprod(Zw, Xdm))
  Xhat <- Zdm %*% Pi
  Xw <- Xhat * w
  b <- solve(crossprod(Xw, Xdm), crossprod(Xw, zdm))
  list(b = as.vector(b), resid = as.vector(zdm - Xdm %*% b), Xhat = Xhat)
}

# fereg's staged inner tolerance (fereg-core/src/ppml/mod.rs:204-247, as the Julia sandbox).
.tol_fereg <- function(iter, eps, inner_tol, tol, start_tol, target_tol) {
  k_next <- iter + 1
  fr <- min(k_next - 1, 7) / 7
  staged <- if (k_next >= 8) target_tol else exp((1 - fr) * log(start_tol) + fr * log(target_tol))
  prog <- if (iter > 1) {
    if (eps <= tol) target_tol else min(max(0.1 * eps, target_tol), start_tol)
  } else start_tol
  max(min(staged, prog, inner_tol), target_tol)
}

# The IRLS-IV loop. Returns the state at convergence.
.irls <- function(y, X, Zi, K_exog, w_user, off_c, S, ws_init, tol, target_tol, maxit, acc,
                  sep_mu, verbose) {
  N <- length(y); K <- ncol(X); n_inst <- ncol(Zi)
  mean_y <- sum(w_user * y) / sum(w_user)
  mu <- 0.5 * (y + mean_y)
  mu <- pmax(mu, pmax(0.05 * y, 1e-3))
  eta <- log(mu) + off_c
  ypos <- y > 0
  y0 <- which(!ypos)
  wylogy <- sum((w_user * y * log(y))[ypos])
  wy <- w_user * y
  start_tol <- 1e-4
  inner_tol <- max(start_tol, tol)
  alt_tol <- start_tol
  converged <- FALSE
  deviance <- NA_real_
  eps <- NA_real_
  ok <- 0L
  nsh <- 0L
  prev_sh <- FALSE
  b <- numeric(K)
  old_eta <- eta
  sep_mask <- logical(N)
  num_sep_mu <- 0L
  D <- NULL
  z_old <- NULL
  tr_tol <- numeric(0); tr_cg <- integer(0)
  t_part <- 0
  ex <- seq_len(K_exog)
  xi <- 1L + seq_len(K)
  zi <- 1L + K + seq_len(n_inst)
  if (verbose > 0L) cat(sprintf("IRLS-IV iterations (N = %d, K = %d, L = %d)\n", N, K, K_exog + n_inst))
  iter <- 0L
  for (iter in seq_len(maxit)) {
    if (any(!is.finite(mu)) || any(!is.finite(eta)))
      stop(sprintf("ivppmlhdfe: mu has infinite or missing values on iteration %d (runaway divergence); ", iter),
           "try separation = \"all\" or a simpler fixed-effect structure")
    if (iter > 10L && K > 0L && (max(abs(b)) > 1e6 || max(abs(eta)) > 30))
      stop(sprintf("ivppmlhdfe: IRLS diverged at iteration %d (max|b| = %.2e, max|eta| = %.2e)",
                   iter, max(abs(b)), max(abs(eta))))
    z <- eta - off_c - 1 + y / mu
    if (sep_mu && length(y0)) z[y0] <- eta[y0] - off_c[y0] - 1
    irls_w <- pmax(w_user * mu, 1e-20)
    tp0 <- proc.time()[["elapsed"]]
    .fe_weights(S, irls_w)
    if (acc$w && !is.null(D)) {
      D[, 1L] <- D[, 1L] + (z - z_old)
    } else if (acc$w && !is.null(ws_init)) {
      D <- cbind(z, ws_init)
      ws_init <- NULL
    } else {
      D <- cbind(z, X, Zi)
    }
    if (acc$w) z_old <- z
    cg0 <- S$cg_iters
    D <- .fe_partial(S, D, inner_tol)$E
    t_part <- t_part + (proc.time()[["elapsed"]] - tp0)
    tr_tol <- c(tr_tol, inner_tol); tr_cg <- c(tr_cg, S$cg_iters - cg0)
    zdm <- D[, 1L]
    Xdm <- D[, xi, drop = FALSE]
    Zdm <- cbind(Xdm[, ex, drop = FALSE], D[, zi, drop = FALSE])
    b_old <- b
    fit <- tryCatch(.tsls(zdm, Xdm, Zdm, irls_w), error = function(e)
      stop("ivppmlhdfe: 2SLS solve failed at iteration ", iter,
           " (rank-deficient design: check for collinear regressors or instruments)", call. = FALSE))
    b <- fit$b
    resid <- fit$resid
    beta_change <- max(abs(b - b_old))
    if (!prev_sh) old_eta <- eta
    eta <- z - resid + off_c
    if (sep_mu && iter > 1L && length(y0)) {
      log_septol <- log(1e-6)
      adj <- if (any(ypos)) log_septol + min(min(eta[ypos]) + 5, 0) else log_septol
      sep_mask <- sep_mask | (eta <= adj & !ypos)
      num_sep_mu <- sum(sep_mask)
    }
    mu <- exp(eta)
    if (num_sep_mu > 0L) mu[sep_mask] <- 0
    mu <- pmax(mu, .EPS100)
    old_dev <- deviance
    dh <- sum(w_user * (mu - y)) + wylogy - sum(wy * eta)       # wy = 0 where y = 0
    if (2 * dh / N < .EPS1) dh <- 0
    if (abs(dh) < .EPS1) dh <- 0
    deviance <- max(2 * dh, 0)
    is_sh <- FALSE
    if (iter > 1L) {
      delta <- old_dev - deviance
      if (deviance < 0.1 * delta) delta <- deviance
      eps <- abs(delta) / max(min(deviance, old_dev), 0.1)
      if (eps < tol) {
        if (inner_tol <= 1.1 * target_tol || S$F <= 1L) {
          ok <- ok + 1L
          converged <- ok >= 1L
        }
      } else if (delta < 0 && nsh < 2L) {
        eta <- 0.9 * old_eta + 0.1 * eta
        if (nsh > 0L) eta <- pmax(eta, -10)
        mu <- pmax(exp(eta), .EPS100)
        is_sh <- TRUE
        ok <- 0L
        nsh <- nsh + 1L
      } else {
        ok <- 0L
        nsh <- 0L
      }
    }
    if (verbose > 0L) {
      cat(sprintf("Iter %3d:  dev = %-11.5e", iter, deviance))
      if (iter > 1L) cat(sprintf("  eps = %-9.4e", eps))
      cat(sprintf("  db = %-9.4e  tol = %5.0e%s%s\n", beta_change, inner_tol,
                  if (is_sh) "  H" else "", if (ok > 0L) "  O" else ""))
    }
    prev_sh <- is_sh
    if (is_sh) {
      deviance <- old_dev
      next
    }
    if (converged) break
    if (acc$t) {
      inner_tol <- .tol_fereg(iter, if (iter > 1L) eps else Inf, inner_tol, tol, start_tol, target_tol)
    } else if (iter > 1L && eps < inner_tol) {
      inner_tol <- max(min(0.1 * inner_tol, alt_tol), target_tol)
      alt_tol <- 10^-ceiling(log10(1 / max(0.1 * eps, .EPS1)))
    }
  }
  if (!converged)
    stop(sprintf("ivppmlhdfe: failed to converge in %d iterations (eps = %.4e); ", maxit, eps),
         "try a larger maxit or a looser tol")
  if (verbose > 0L) cat(sprintf("Converged in %d iterations (tol = %.4e)\n", iter, tol))
  list(b = b, resid = resid, D = D, zdm = zdm, Xdm = Xdm, Zdm = Zdm, z_last = z, irls_w = irls_w,
       eta = eta, mu = mu, iterations = iter, deviance = deviance, eps = eps, inner_tol = inner_tol,
       num_sep_mu = num_sep_mu, sep_mask = sep_mask, tol_trace = tr_tol, cg_trace = tr_cg,
       t_partial = t_part)
}

#' @export
ivppmlhdfe <- function(fml, data, vcov = "robust", weights = NULL,
                       weight_type = c("pweight", "fweight"), offset = NULL, exposure = NULL,
                       accel = TRUE, ffirst = FALSE, tol = 1e-8, itol = NULL, maxit = 10000L,
                       separation = c("fe", "relu"), keep_singletons = FALSE, verbose = 0L) {
  cl <- match.call()
  t_start <- proc.time()[["elapsed"]]
  data <- as.data.frame(data)
  weight_type <- match.arg(weight_type)
  acc <- .parse_accel(accel)
  sep <- .parse_sep(separation)
  pf <- .parse_fml(fml)
  vc <- .parse_vcov(vcov)
  verbose <- as.integer(verbose)
  if (!is.numeric(tol) || tol <= 0 || tol >= 1) stop("ivppmlhdfe: tol must be in (0, 1)")
  target_tol <- if (is.null(itol) || itol <= 0) max(1e-12, 0.1 * tol) else itol
  note <- function(...) if (verbose >= 0L) message(...)
  N0 <- nrow(data)

  # ---- 1. Columns ----
  yname <- paste(deparse(pf$y), collapse = "")
  y_all <- as.numeric(eval(pf$y, data, environment(fml)))
  if (length(y_all) != N0) stop("ivppmlhdfe: the dependent variable must have one value per row")
  Xe <- .part_matrix(pf$exog, data)
  Xn <- .part_matrix(pf$endog, data)
  Zi <- .part_matrix(pf$inst, data)
  exog <- colnames(Xe); endog <- colnames(Xn); inst <- colnames(Zi)
  if (length(exog) + length(endog) == 0L)
    stop("ivppmlhdfe: at least one regressor is required (an empty model is not supported)")
  if (length(inst) < length(endog)) stop("ivppmlhdfe: equation not identified (fewer instruments than endogenous regressors)")
  if (length(inst) && length(intersect(inst, c(exog, endog))))
    stop("ivppmlhdfe: an instrument cannot equal a regressor: ",
         paste(intersect(inst, c(exog, endog)), collapse = ", "))
  fe_vars <- unique(unlist(pf$fe))
  miss_fe <- setdiff(fe_vars, names(data))
  if (length(miss_fe)) stop("ivppmlhdfe: fixed-effect variable(s) not found: ", paste(miss_fe, collapse = ", "))
  cl_vars <- unique(unlist(vc$vars))
  miss_cl <- setdiff(cl_vars, names(data))
  if (length(miss_cl)) stop("ivppmlhdfe: cluster variable(s) not found: ", paste(miss_cl, collapse = ", "))
  w_all <- .get_col(weights, data, "weight")
  has_w <- !is.null(w_all)
  if (has_w) w_all <- as.numeric(w_all)
  off_all <- .get_col(offset, data, "offset")
  exp_all <- .get_col(exposure, data, "exposure")
  if (!is.null(off_all) && !is.null(exp_all)) stop("ivppmlhdfe: cannot specify both offset and exposure")

  # ---- 2. Estimation sample (missing values, nonnegative y and weights, zero weights) ----
  touse <- !is.na(y_all) & stats::complete.cases(Xe, Xn, Zi)
  for (v in c(fe_vars, cl_vars)) touse <- touse & !is.na(data[[v]])
  if (has_w) touse <- touse & !is.na(w_all)
  if (!is.null(exp_all)) {
    exp_all <- as.numeric(exp_all)
    if (any(touse & !is.na(exp_all) & exp_all <= 0)) stop("ivppmlhdfe: exposure must be greater than zero")
    off_all <- log(exp_all)
  }
  has_off <- !is.null(off_all)
  if (has_off) { off_all <- as.numeric(off_all); touse <- touse & is.finite(off_all) }
  touse <- touse & is.finite(y_all)
  if (any(y_all[touse] < 0)) stop("ivppmlhdfe: ", yname, " must be nonnegative")
  if (has_w) {
    if (any(w_all[touse] < 0)) stop("ivppmlhdfe: negative weights encountered")
    touse <- touse & w_all > 0
  }
  if (any(touse)) {
    fin <- rowSums(!is.finite(cbind(Xe, Xn, Zi)[touse, , drop = FALSE])) == 0
    touse[which(touse)[!fin]] <- FALSE
  }
  N_full <- sum(touse)
  if (N_full == 0L) stop("ivppmlhdfe: no observations")
  if (sum(y_all[touse]) == 0) stop("ivppmlhdfe: ", yname, " is identically zero")
  rows <- which(touse)
  y <- y_all[rows]
  X <- cbind(Xe, Xn)[rows, , drop = FALSE]
  Zx <- Zi[rows, , drop = FALSE]
  w_user <- if (has_w) w_all[rows] else rep(1, length(rows))
  off <- if (has_off) off_all[rows] else rep(0, length(rows))
  if (length(inst)) {
    zc <- vapply(seq_len(ncol(Zx)), function(j) { r <- range(Zx[, j]); r[1L] == r[2L] }, NA)
    if (any(zc)) stop("ivppmlhdfe: instrument ", inst[which(zc)[1L]], " has zero variance (constant column)")
  }
  has_fe <- length(pf$fe) > 0L
  fes <- if (has_fe) lapply(pf$fe, function(vs) .interact(lapply(vs, function(v) data[[v]][rows])))
         else list(rep(1L, length(rows)))
  fe_labels <- if (has_fe) vapply(pf$fe, function(vs) paste(vs, collapse = "^"), "") else character(0)
  clusters <- NULL
  vtype <- vc$type
  if (vtype == "cluster") {
    clusters <- lapply(vc$vars, function(vs) .interact(lapply(vs, function(v) data[[v]][rows])))
    if (length(clusters) == 1L && max(clusters[[1L]]) <= 1L) {
      warning("ivppmlhdfe: only 1 cluster; switching to robust standard errors")
      vtype <- "robust"; clusters <- NULL
    }
  }
  subset_all <- function(keep) {
    rows <<- rows[keep]; y <<- y[keep]; X <<- X[keep, , drop = FALSE]; Zx <<- Zx[keep, , drop = FALSE]
    w_user <<- w_user[keep]; off <<- off[keep]
    fes <<- lapply(fes, function(g) .dense(g[keep]))
    if (!is.null(clusters)) clusters <<- lapply(clusters, function(g) g[keep])
  }

  # ---- 3. Singletons and fixed-effect separation (reghdfe5 loop, y as iweight) ----
  num_singletons <- 0L
  if (!keep_singletons) {
    ds <- .drop_singletons(fes, y, sep_fe = "fe" %in% sep)
    num_singletons <- ds$ndrop
    if (ds$ndrop > 0L) {
      subset_all(ds$keep)
      note(sprintf("(dropped %d observations that are either singletons or separated by a fixed effect)", ds$ndrop))
    }
  } else if (verbose >= 0L) {
    warning("ivppmlhdfe: keeping singleton groups keeps fixed effects that cause separation")
  }
  # sort by the fixed effect with the most levels (its groups become contiguous)
  big <- which.max(vapply(fes, max, 0L))
  o <- order(fes[[big]], method = "radix")
  if (is.unsorted(fes[[big]])) subset_all(o)

  # ---- 4. Collinearity (per block, after partialling out the FEs under the user weights) ----
  K_exog <- length(exog)
  t_coll <- 0
  S_last <- NULL
  coll_pass <- function() {
    tc0 <- proc.time()[["elapsed"]]
    on.exit(t_coll <<- t_coll + (proc.time()[["elapsed"]] - tc0))
    S0 <- .fe_setup(fes)
    S_last <<- S0
    .fe_weights(S0, w_user)
    raw <- cbind(X, Zx)
    part <- .fe_partial(S0, raw, target_tol)$E
    K <- ncol(X); Ke <- length(exog); Kn <- length(endog); L <- ncol(Zx)
    ke <- .select_not_collinear(part[, seq_len(Ke), drop = FALSE], raw[, seq_len(Ke), drop = FALSE], w_user)
    kn <- .select_not_collinear(part[, Ke + seq_len(Kn), drop = FALSE], raw[, Ke + seq_len(Kn), drop = FALSE], w_user)
    kz <- .select_not_collinear(part[, K + seq_len(L), drop = FALSE], raw[, K + seq_len(L), drop = FALSE], w_user)
    dropped <- c(exog[setdiff(seq_len(Ke), ke)], endog[setdiff(seq_len(Kn), kn)], inst[setdiff(seq_len(L), kz)])
    if (length(dropped)) {
      note("note: omitted because of collinearity: ", paste(dropped, collapse = " "))
      cols <- c(ke, Ke + kn)
      X <<- X[, cols, drop = FALSE]
      Zx <<- Zx[, kz, drop = FALSE]
      exog <<- exog[ke]; endog <<- endog[kn]; inst <<- inst[kz]
      if (length(inst) < length(endog)) stop("ivppmlhdfe: equation not identified after collinearity removal")
      if (length(exog) + length(endog) == 0L) stop("ivppmlhdfe: every regressor was dropped as collinear")
      part <- part[, c(cols, K + kz), drop = FALSE]
    }
    part
  }
  part <- coll_pass()

  # ---- 5. ReLU separation (regressors + FEs; instruments play no role) ----
  num_sep_relu <- 0L
  if ("relu" %in% sep && any(y == 0)) {
    sidx <- .relu_separation(y, X, fes, verbose = verbose)
    if (length(sidx)) {
      n_before <- length(y)
      keep <- rep(TRUE, n_before); keep[sidx] <- FALSE
      subset_all(keep)
      if (!keep_singletons) {
        ds <- .drop_singletons(fes, y, sep_fe = FALSE)
        if (ds$ndrop > 0L) subset_all(ds$keep)
      }
      num_sep_relu <- n_before - length(y)
      part <- coll_pass()
    }
  }
  K_exog <- length(exog)
  N <- length(y)
  K <- ncol(X); n_inst <- ncol(Zx)
  colnames(X) <- c(exog, endog)
  colnames(Zx) <- inst

  # ---- 6. IRLS-IV ----
  off_mean <- if (has_off) sum(w_user * off) / sum(w_user) else 0
  off_c <- off - off_mean
  S <- S_last      # solver structure of the final sample (built by the last coll_pass); counters cumulate
  fit <- .irls(y, X, Zx, K_exog, w_user, off_c, S, ws_init = if (acc$w) part else NULL,
               tol = tol, target_tol = target_tol, maxit = as.integer(maxit), acc = acc,
               sep_mu = "mu" %in% sep, verbose = verbose)
  b <- fit$b
  names(b) <- c(exog, endog)
  mu <- fit$mu
  eta <- fit$eta
  is_fw <- has_w && weight_type == "fweight"

  # ---- 7. Columns and weights for the VCE (refresh at the converged weights if accel) ----
  t_r0 <- proc.time()[["elapsed"]]
  if (acc$r) {
    w_b <- pmax(w_user * mu, 1e-20)
    z_f <- eta - off_c - 1 + y / mu
    if ("mu" %in% sep) { y0 <- which(y == 0); z_f[y0] <- eta[y0] - off_c[y0] - 1 }
    .fe_weights(S, w_b)
    D <- if (acc$w) { Dw <- fit$D; Dw[, 1L] <- Dw[, 1L] + (z_f - fit$z_last); Dw } else cbind(z_f, X, Zx)
    D <- .fe_partial(S, D, fit$inner_tol)$E
    Xdm <- D[, 1L + seq_len(K), drop = FALSE]
    Zdm <- cbind(Xdm[, seq_len(K_exog), drop = FALSE], D[, 1L + K + seq_len(n_inst), drop = FALSE])
    resid <- as.vector(D[, 1L] - Xdm %*% b)
  } else {
    w_b <- fit$irls_w
    Xdm <- fit$Xdm; Zdm <- fit$Zdm; resid <- fit$resid
  }
  t_refresh <- proc.time()[["elapsed"]] - t_r0
  w_m <- if (is_fw) sqrt(w_user) * mu else w_b
  Zw <- Zdm * w_b
  Xhat <- Zdm %*% solve(crossprod(Zw, Zdm), crossprod(Zw, Xdm))
  N_eff <- if (is_fw) sum(w_user) else N
  vv <- .vce(Xhat, Xdm, w_b, w_m, resid, N_eff, if (vtype == "cluster") clusters else NULL)
  V <- vv$V
  dimnames(V) <- list(names(b), names(b))
  se <- sqrt(abs(diag(V)))
  names(se) <- names(b)
  rank <- qr(crossprod(Xhat, Xdm * w_b))$rank

  # ---- 8. Constant, likelihoods, FE sum ----
  off_u <- off
  sw_b <- sum(w_b)
  b_cons <- sum(w_b * (eta - off_u)) / sw_b - sum(colSums(X * w_b) / sw_b * b)
  ll <- sum(w_user * (y * eta - mu - lgamma(y + 1)))
  if (has_off && sum(abs(off_u)) > 0) {
    c0 <- log(sum(w_user * y) / sum(w_user * exp(off_u)))
    ll_0 <- sum(w_user * (y * (c0 + off_u) - exp(c0 + off_u) - lgamma(y + 1)))
  } else {
    m0 <- sum(w_user * y) / sum(w_user)
    ll_0 <- sum(w_user * (y * log(m0) - m0 - lgamma(y + 1)))
  }
  d <- as.vector(eta - off_u - X %*% b - b_cons)
  d <- d - sum(w_b * d) / sw_b

  # ---- 9. First-stage statistics ----
  ff <- NULL
  if (isTRUE(ffirst) && length(endog) > 0L) {
    ff <- .ffirst(Xdm, Zdm, w_b, K_exog, endog, inst, N, N_eff, vtype,
                  if (vtype == "cluster") clusters else NULL, is_fw, w_user, mu, S, fes,
                  has_fe = has_fe)
  }

  esample <- logical(N0); esample[rows] <- TRUE
  out <- list(
    coefficients = b, vcov = V, se = se, constant = b_cons,
    N = N_eff, N_rows = N, N_full = N_full,
    num_singletons = num_singletons, num_separated = num_sep_relu + fit$num_sep_mu,
    num_sep_relu = num_sep_relu, num_sep_mu = fit$num_sep_mu,
    iterations = fit$iterations, converged = TRUE, deviance = fit$deviance,
    ll = ll, ll_0 = ll_0, r2_p = 1 - ll / ll_0, rank = rank,
    vcov_type = vtype, clustervars = if (vtype == "cluster") vc$labels else character(0),
    N_clust = vv$G,
    ffirst = ff,
    depvar = yname, exog = exog, endog = endog, instruments = inst, fe = fe_labels,
    weights = if (has_w) list(var = if (is.character(weights)) weights else "<vector>", type = weight_type) else NULL,
    offset = if (has_off) off_all else NULL, has_offset = has_off,
    accel = acc, tol = tol, itol = target_tol, separation = sep,
    esample = esample, obs = rows, eta = eta, mu = mu, d = d,
    inner = list(tol_trace = fit$tol_trace, cg_trace = fit$cg_trace, cg_total = S$cg_iters,
                 passes = S$passes, restarts = S$restarts),
    time = list(total = proc.time()[["elapsed"]] - t_start, partial = fit$t_partial,
                coll = t_coll, refresh = t_refresh),
    fml = fml, terms = list(exog = pf$exog, endog = pf$endog, inst = pf$inst),
    call = cl
  )
  class(out) <- "ivppmlhdfe"
  out
}
