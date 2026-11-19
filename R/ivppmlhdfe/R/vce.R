# Variance-covariance estimators (conventions of the Stata command ivppmlhdfe v0.9.4).
#   robust : (N/(N-1)) B M B,  B = (Xhat' W X)^{-1},  M = sum_i (xhat_i w_meat_i r_i)(...)'
#            N = number of rows, or sum of the frequency weights with fweights;
#            w_meat = w_user * mu (= the IRLS weight), or sqrt(w_user) * mu with fweights.
#   cluster: (G/(G-1)) B M_G B with M_G = sum_g s_g s_g', s_g = sum_{i in g} xhat_i w_i r_i.
#   multiway (Cameron, Gelbach and Miller 2011): inclusion-exclusion over the non-empty
#            subsets of the cluster variables, each term with its own G/(G-1); then, if the
#            smallest eigenvalue is negative, V - min_ev * I (as the Stata command).

.sym <- function(A) (A + t(A)) / 2

.bread <- function(Xhat, Xdm, w) {
  H <- crossprod(Xhat, Xdm * w)
  .sym(solve(.sym(H)))
}

.cluster_meat <- function(scores, codes) {
  sg <- rowsum(scores, codes, reorder = FALSE)
  list(meat = crossprod(sg), G = nrow(sg))
}

# clusters: list of integer code vectors (one per cluster variable), or NULL for robust.
.vce <- function(Xhat, Xdm, w_bread, w_meat, resid, N_eff, clusters = NULL) {
  K <- ncol(Xdm)
  if (K == 0L) return(list(V = matrix(0, 0, 0), G = integer(0)))
  B <- .bread(Xhat, Xdm, w_bread)
  if (is.null(clusters) || !length(clusters)) {
    sc <- Xhat * (w_meat * resid)
    V <- (N_eff / (N_eff - 1)) * (B %*% crossprod(sc) %*% B)
    return(list(V = .sym(V), G = integer(0)))
  }
  scores <- Xhat * (w_bread * resid)
  nC <- length(clusters)
  if (nC == 1L) {
    cm <- .cluster_meat(scores, clusters[[1L]])
    if (cm$G <= 1L) {
      sc <- Xhat * (w_meat * resid)
      V <- (N_eff / (N_eff - 1)) * (B %*% crossprod(sc) %*% B)
    } else {
      V <- (cm$G / (cm$G - 1)) * (B %*% cm$meat %*% B)
    }
    return(list(V = .sym(V), G = cm$G))
  }
  V <- matrix(0, K, K)
  Gs <- integer(nC)
  for (mask in seq_len(2L^nC - 1L)) {
    members <- which(bitwAnd(mask, 2L^(seq_len(nC) - 1L)) > 0L)
    cid <- if (length(members) == 1L) clusters[[members]] else .interact(clusters[members])
    cm <- .cluster_meat(scores, cid)
    sgn <- if (length(members) %% 2L == 1L) 1 else -1
    if (cm$G > 1L) V <- V + sgn * (cm$G / (cm$G - 1)) * (B %*% cm$meat %*% B)
    if (length(members) == 1L) Gs[members] <- cm$G
  }
  V <- .sym(V)
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 0) V <- V - min(ev) * diag(K)
  list(V = V, G = Gs)
}
