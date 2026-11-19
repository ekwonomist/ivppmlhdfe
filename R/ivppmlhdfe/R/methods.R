# S3 methods for "ivppmlhdfe" objects.

coef.ivppmlhdfe <- function(object, ...) object$coefficients

vcov.ivppmlhdfe <- function(object, ...) object$vcov

nobs.ivppmlhdfe <- function(object, ...) object$N

deviance.ivppmlhdfe <- function(object, ...) object$deviance

logLik.ivppmlhdfe <- function(object, ...) {
  structure(object$ll, df = object$rank, nobs = object$N, class = "logLik")
}

fitted.ivppmlhdfe <- function(object, ...) predict.ivppmlhdfe(object, type = "mu")

.coef_table <- function(object, level = 0.95) {
  b <- object$coefficients
  se <- object$se
  zst <- b / se
  p <- 2 * stats::pnorm(-abs(zst))
  q <- stats::qnorm(1 - (1 - level) / 2)
  tab <- cbind(Estimate = b, `Std. Error` = se, `z value` = zst, `Pr(>|z|)` = p,
               lower = b - q * se, upper = b + q * se)
  colnames(tab)[5:6] <- sprintf("%g%%", 100 * c((1 - level) / 2, 1 - (1 - level) / 2))
  tab
}

summary.ivppmlhdfe <- function(object, level = 0.95, ...) {
  out <- object
  out$coeftable <- .coef_table(object, level)
  class(out) <- c("summary.ivppmlhdfe", class(object))
  out
}

.header <- function(x) {
  cat("IV-PPML with high-dimensional fixed effects (ivppmlhdfe ", .IVPPML_VERSION, ")\n", sep = "")
  cat("Dependent variable: ", x$depvar, "\n", sep = "")
  if (length(x$endog))
    cat("Endogenous: ", paste(x$endog, collapse = " "), "   Instruments: ",
        paste(x$instruments, collapse = " "), "\n", sep = "")
  cat("Fixed effects: ", if (length(x$fe)) paste(x$fe, collapse = ", ") else "(none)", "\n", sep = "")
  cat(sprintf("Observations: %s", format(x$N, big.mark = ",")))
  if (x$num_singletons > 0) cat(sprintf("   (singletons/FE-separated dropped: %d)", x$num_singletons))
  if (x$num_separated > 0) cat(sprintf("   (separated: %d)", x$num_separated))
  cat("\n")
  se_lab <- if (x$vcov_type == "cluster")
    sprintf("clustered by %s (%s clusters)", paste(x$clustervars, collapse = ", "),
            paste(x$N_clust, collapse = ", ")) else "robust"
  cat("Standard errors: ", se_lab, "\n", sep = "")
  cat(sprintf("Iterations: %d   Deviance: %.6g   Log pseudo-likelihood: %.6g   Pseudo R2: %.4f\n",
              x$iterations, x$deviance, x$ll, x$r2_p))
}

.IVPPML_VERSION <- "1.0.0"

print.ivppmlhdfe <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  .header(x)
  tab <- .coef_table(x)[, 1:4, drop = FALSE]
  stats::printCoefmat(tab, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  if (!is.null(x$ffirst)) print_ffirst(x$ffirst, digits = digits)
  invisible(x)
}

print.summary.ivppmlhdfe <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  .header(x)
  stats::printCoefmat(x$coeftable[, 1:4, drop = FALSE], digits = digits, P.values = TRUE,
                      has.Pvalue = TRUE)
  cat(sprintf("_cons (FEs partialled out, no standard error): %.6g\n", x$constant))
  if (!is.null(x$ffirst)) print_ffirst(x$ffirst, digits = digits)
  invisible(x)
}

predict.ivppmlhdfe <- function(object, newdata = NULL, type = c("mu", "xb", "eta", "d"), ...) {
  type <- match.arg(type)
  b <- object$coefficients
  if (is.null(newdata)) {
    n0 <- length(object$esample)
    out <- rep(NA_real_, n0)
    eta <- object$eta
    val <- switch(type,
      mu = exp(eta),
      eta = eta,
      d = object$d,
      xb = eta - object$d - (if (object$has_offset) object$offset[object$obs] else 0))
    out[object$obs] <- val
    return(out)
  }
  newdata <- as.data.frame(newdata)
  X <- cbind(.part_matrix(object$terms$exog, newdata), .part_matrix(object$terms$endog, newdata))
  X <- X[, names(b), drop = FALSE]
  xb <- as.vector(X %*% b) + object$constant
  if (type == "xb") return(xb)
  if (length(object$fe))
    stop("predict: type = \"", type, "\" with newdata needs the fixed-effect estimates; ",
         "only type = \"xb\" is available for new data in models with fixed effects")
  if (type == "d") return(rep(0, length(xb)))
  if (object$has_offset)
    stop("predict: with newdata and an offset, compute exp(xb + offset) yourself")
  if (type == "eta") xb else exp(xb)
}
