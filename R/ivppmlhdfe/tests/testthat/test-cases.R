# Every case of the shared cross-language test set (v010_build/testdata/cases.json) against the
# Stata v0.9.4 reference (results/stata_v094.json): coef and SE within 1e-6 relative, counts exact,
# for accel on and off. Runs only when IVPPML_TESTDATA (folder with cases.json and the CSVs) and
# IVPPML_REFERENCE (the reference JSON) are set; IVPPML_FFIRST_REFERENCE (the Stata ffirst
# reference, results/stata_ffirst_reference.json) adds the first-stage statistics at 1e-6.
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-3)

test_that("all cross-language cases match the Stata v0.9.4 reference", {
  tdir <- Sys.getenv("IVPPML_TESTDATA")
  ref_file <- Sys.getenv("IVPPML_REFERENCE")
  skip_if(!nzchar(tdir) || !file.exists(file.path(tdir, "cases.json")), "IVPPML_TESTDATA not set")
  skip_if(!nzchar(ref_file) || !file.exists(ref_file), "IVPPML_REFERENCE not set")
  skip_if_not_installed("jsonlite")
  cases <- jsonlite::fromJSON(file.path(tdir, "cases.json"), simplifyVector = FALSE)$cases
  ref <- jsonlite::fromJSON(ref_file, simplifyVector = FALSE)$cases
  ff_file <- Sys.getenv("IVPPML_FFIRST_REFERENCE")
  ffref <- if (nzchar(ff_file) && file.exists(ff_file))
    jsonlite::fromJSON(ff_file, simplifyVector = FALSE)$cases else NULL
  dat <- list()
  for (cs in cases) {
    if (is.null(dat[[cs$data]])) dat[[cs$data]] <- read.csv(file.path(tdir, cs$data))
    ex <- unlist(cs$exog); en <- unlist(cs$endog); iv <- unlist(cs$instruments)
    fe <- vapply(cs$absorb, function(a) paste(unlist(a), collapse = "^"), "")
    s <- paste(cs$depvar, "~", paste(ex, collapse = " + "), "|", paste(fe, collapse = " + "))
    if (length(en)) s <- paste(s, "|", paste(en, collapse = " + "), "~", paste(iv, collapse = " + "))
    vc <- if (identical(cs$vce$type, "cluster"))
      stats::as.formula(paste("~", paste(unlist(cs$vce$vars), collapse = " + "))) else "robust"
    rf <- ref[[cs$id]]
    for (acc in c(TRUE, FALSE)) {
      r <- ivppmlhdfe(stats::as.formula(s), dat[[cs$data]], vcov = vc, weights = cs$weights,
                      weight_type = if (is.null(cs$weight_type)) "pweight" else cs$weight_type,
                      offset = cs$offset, accel = acc, ffirst = length(en) > 0L, verbose = -1)
      lab <- paste(cs$id, if (acc) "accel" else "noaccel")
      if (!is.null(ffref) && length(en) && !is.null(ffref[[cs$id]]$ffirst)) {
        fr <- ffref[[cs$id]]$ffirst
        for (k in en) for (st in names(fr[[k]]))
          expect_lt(rel(r$ffirst$first[[k]][[st]], fr[[k]][[st]]), 1e-6, label = paste(lab, k, st))
        for (st in c("idstat", "iddf", "idp", "cdf", "widstat"))
          expect_lt(rel(r$ffirst[[st]], fr[[st]]), 1e-6, label = paste(lab, st))
      }
      bc <- unlist(rf$coef); bs <- unlist(rf$se)
      expect_lt(max(rel(coef(r)[names(bc)], bc)), 1e-6, label = paste(lab, "coef"))
      expect_lt(max(rel(r$se[names(bs)], bs)), 1e-6, label = paste(lab, "se"))
      expect_equal(r$N, rf$N, label = paste(lab, "N"))
      expect_equal(r$num_singletons, rf$num_singletons, label = paste(lab, "singletons"))
      expect_equal(r$num_separated, rf$num_separated, label = paste(lab, "separated"))
    }
  }
})
