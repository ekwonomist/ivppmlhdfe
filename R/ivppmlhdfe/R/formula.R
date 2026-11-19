# Formula grammar (fixest style):
#   y ~ x2 + x4 | id + year | x1 + x3 ~ z + z2 + z3
#   y ~ x2 | exp^year + imp^year                   (no IV part: plain PPML)
#   y ~ x2 | x1 ~ z                                 (IV without fixed effects)
# Parts are separated by top-level `|`; FE interactions use `^` (or `:`).

.split_op <- function(e, op) {
  if (is.call(e) && identical(e[[1L]], as.name(op)) && length(e) == 3L)
    c(.split_op(e[[2L]], op), list(e[[3L]]))
  else if (is.call(e) && identical(e[[1L]], as.name("(")))
    .split_op(e[[2L]], op)
  else list(e)
}

.is_tilde <- function(e) is.call(e) && identical(e[[1L]], as.name("~")) && length(e) == 3L

.parse_fml <- function(fml) {
  if (!inherits(fml, "formula") || length(fml) != 3L)
    stop("ivppmlhdfe: 'fml' must be a two-sided formula, e.g. y ~ x2 | id + year | x1 ~ z")
  lhs <- fml[[2L]]
  rhs <- fml[[3L]]
  iv <- FALSE
  inst <- NULL
  if (.is_tilde(lhs)) {             # `~` is left-associative: (y ~ a | b | c) ~ d
    y <- lhs[[2L]]; middle <- lhs[[3L]]; inst <- rhs; iv <- TRUE
  } else if (.is_tilde(rhs)) {
    y <- lhs; middle <- rhs[[2L]]; inst <- rhs[[3L]]; iv <- TRUE
  } else {
    y <- lhs; middle <- rhs
  }
  parts <- .split_op(middle, "|")
  exog <- parts[[1L]]
  fe <- NULL
  endog <- NULL
  if (iv) {
    if (length(parts) == 2L) {
      endog <- parts[[2L]]
    } else if (length(parts) == 3L) {
      fe <- parts[[2L]]; endog <- parts[[3L]]
    } else {
      stop("ivppmlhdfe: with an IV part the formula must be y ~ exog | fe | endog ~ instruments ",
           "or y ~ exog | endog ~ instruments")
    }
  } else {
    if (length(parts) == 2L) fe <- parts[[2L]]
    else if (length(parts) > 2L)
      stop("ivppmlhdfe: too many '|' parts; use y ~ exog | fe or y ~ exog | fe | endog ~ instruments")
  }
  fe_terms <- list()
  if (!is.null(fe) && !(is.numeric(fe) && fe %in% c(0, 1))) {
    for (t in .split_op(fe, "+")) {
      v <- all.vars(t)
      if (!length(v)) stop("ivppmlhdfe: cannot read the fixed-effect term ", deparse(t))
      if (is.call(t) && !(deparse(t[[1L]]) %in% c("^", ":")))
        stop("ivppmlhdfe: fixed-effect terms must be variables or interactions with '^' (got ",
             deparse(t), ")")
      fe_terms[[length(fe_terms) + 1L]] <- v
    }
  }
  list(y = y, exog = exog, fe = fe_terms, endog = endog, inst = inst, iv = iv)
}

# Numeric design columns of one formula part (no intercept). Rows are kept as in `data`
# (NA rows stay NA so the caller can drop them jointly).
.part_matrix <- function(part, data) {
  if (is.null(part)) return(matrix(0, nrow(data), 0L))
  if (is.numeric(part) && length(part) == 1L && part %in% c(0, 1))
    return(matrix(0, nrow(data), 0L))
  tt <- stats::terms(stats::as.formula(call("~", part)), data = data)
  if (!length(attr(tt, "term.labels"))) return(matrix(0, nrow(data), 0L))
  mf <- stats::model.frame(tt, data = data, na.action = stats::na.pass)
  mm <- stats::model.matrix(tt, mf)
  keep <- colnames(mm) != "(Intercept)"
  mm <- mm[, keep, drop = FALSE]
  attr(mm, "assign") <- NULL
  attr(mm, "contrasts") <- NULL
  if (nrow(mm) != nrow(data)) stop("internal: model.matrix dropped rows")
  storage.mode(mm) <- "double"
  mm
}

# A column given as a name in `data` or as a vector of length nrow(data).
.get_col <- function(x, data, role) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1L) {
    if (!x %in% names(data)) stop("ivppmlhdfe: ", role, " variable '", x, "' not found in data")
    return(data[[x]])
  }
  if (inherits(x, "formula")) {
    v <- all.vars(x)
    if (length(v) != 1L) stop("ivppmlhdfe: ", role, " must name one variable")
    return(data[[v]])
  }
  if (length(x) != nrow(data)) stop("ivppmlhdfe: ", role, " must have one value per row of data")
  x
}

# vcov argument -> list(type = "robust") or list(type = "cluster", vars = list of var sets)
.parse_vcov <- function(vcov) {
  if (is.null(vcov)) return(list(type = "robust", vars = list(), labels = character(0)))
  if (is.character(vcov)) {
    if (length(vcov) == 1L && tolower(vcov) %in% c("robust", "hetero", "hc1", "r", "iid"))
      return(list(type = "robust", vars = list(), labels = character(0)))
    vars <- lapply(vcov, function(s) strsplit(s, "[\\^:#]")[[1L]])
    return(list(type = "cluster", vars = vars, labels = vcov))
  }
  if (inherits(vcov, "formula")) {
    rhs <- vcov[[length(vcov)]]
    terms <- .split_op(rhs, "+")
    vars <- lapply(terms, all.vars)
    labels <- vapply(terms, function(t) paste(deparse(t), collapse = ""), "")
    return(list(type = "cluster", vars = vars, labels = labels))
  }
  stop("ivppmlhdfe: vcov must be \"robust\", a one-sided formula such as ~pair or ~exp + imp, ",
       "or a character vector of cluster variables")
}
