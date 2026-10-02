#!/usr/bin/env Rscript
# Stage B colocalization backend. Called by stage_b/coloc_backend.py as
#   Rscript coloc_run.R request.json response.json
# Each task is coloc.abf or coloc.susie on two datasets that share the same snp vector.
# Nothing here chooses a threshold or a window; the Python side does that from PREREG.md.
#
# Every field is read with `[[ ]]`, which matches a name exactly. `$` matches partially: on a
# quantitative dataset sent without `s` and `sdY` (the eQTL / sQTL side), `d$s` returns `snp`. No
# `$` is used in this file, and the required fields of the request are checked right after parsing.
#
# Partial matching is an error in this process, wherever it happens: R reports a partial match by
# `$`, of an argument name or of an attribute name as a warning, and the handler below turns that
# warning into an error. A global calling handler runs below every tryCatch of this file, so the
# per-task handler does not catch the error: Rscript exits non-zero and writes no response.
options(warnPartialMatchDollar = TRUE, warnPartialMatchArgs = TRUE, warnPartialMatchAttr = TRUE)
globalCallingHandlers(warning = function(w) {
  if (grepl("partial (argument )?match of", conditionMessage(w))) {
    stop(sprintf("partial matching is an error in coloc_run.R: %s", conditionMessage(w)), call. = FALSE)
  }
})

suppressPackageStartupMessages({
  library(jsonlite)
  library(coloc)
})

TASK_FIELDS <- c("id", "method", "p1", "p2", "p12", "d1", "d2")
DATASET_FIELDS <- c("snp", "beta", "varbeta", "N", "type")
METHODS <- c("abf", "susie")

require_fields <- function(x, fields, what) {
  if (!is.list(x) || is.null(names(x))) stop(sprintf("%s is not a JSON object", what), call. = FALSE)
  missing <- setdiff(fields, names(x))
  if (length(missing) > 0) {
    stop(sprintf("%s lacks the required field(s): %s", what, paste(missing, collapse = ", ")), call. = FALSE)
  }
}

validate_request <- function(req) {
  require_fields(req, "tasks", "the request")
  tasks <- req[["tasks"]]
  if (!is.list(tasks) || !is.null(names(tasks))) stop("the request field `tasks` is not a JSON array", call. = FALSE)
  for (i in seq_along(tasks)) {
    what <- sprintf("task %d", i)
    require_fields(tasks[[i]], TASK_FIELDS, what)
    method <- tasks[[i]][["method"]]
    if (!is.character(method) || length(method) != 1 || !(method %in% METHODS)) {
      stop(sprintf("%s has a method other than abf or susie", what), call. = FALSE)
    }
    if (method == "susie") require_fields(tasks[[i]], "LD", what)
    for (side in c("d1", "d2")) require_fields(tasks[[i]][[side]], DATASET_FIELDS, sprintf("%s, dataset %s", what, side))
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: coloc_run.R request.json response.json")
req <- fromJSON(args[1], simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = TRUE)
validate_request(req)

as_dataset <- function(d) {
  out <- list(snp = as.character(d[["snp"]]), beta = as.numeric(d[["beta"]]), varbeta = as.numeric(d[["varbeta"]]),
              N = as.numeric(d[["N"]]), type = as.character(d[["type"]]))
  if (!is.null(d[["MAF"]])) out[["MAF"]] <- as.numeric(d[["MAF"]])
  if (!is.null(d[["sdY"]])) out[["sdY"]] <- as.numeric(d[["sdY"]])
  if (!is.null(d[["s"]])) out[["s"]] <- as.numeric(d[["s"]])
  out
}

run_abf <- function(t) {
  r <- coloc.abf(as_dataset(t[["d1"]]), as_dataset(t[["d2"]]), p1 = t[["p1"]], p2 = t[["p2"]], p12 = t[["p12"]])
  s <- r[["summary"]]
  results <- r[["results"]]
  list(id = t[["id"]], ok = TRUE, method = "abf",
       pp = unname(as.numeric(s[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])),
       nsnps = unname(as.numeric(s["nsnps"])),
       snp = as.character(results[["snp"]]), snp_pp_h4 = as.numeric(results[["SNP.PP.H4"]]))
}

run_susie <- function(t) {
  ld <- as.matrix(t[["LD"]])
  snp <- as.character(t[["d1"]][["snp"]])
  dimnames(ld) <- list(snp, snp)
  d1 <- as_dataset(t[["d1"]]); d1[["LD"]] <- ld
  d2 <- as_dataset(t[["d2"]]); d2[["LD"]] <- ld
  s1 <- runsusie(d1)
  s2 <- runsusie(d2)
  n1 <- length(s1[["sets"]][["cs"]])
  n2 <- length(s2[["sets"]][["cs"]])
  max_h4 <- NA_real_
  n_pairs <- 0
  if (n1 > 0 && n2 > 0) {
    r <- coloc.susie(s1, s2, p1 = t[["p1"]], p2 = t[["p2"]], p12 = t[["p12"]])
    pair_summary <- r[["summary"]]
    if (!is.null(pair_summary) && nrow(pair_summary) > 0) {
      max_h4 <- max(pair_summary[["PP.H4.abf"]])
      n_pairs <- nrow(pair_summary)
    }
  }
  list(id = t[["id"]], ok = TRUE, method = "susie", n_cs1 = n1, n_cs2 = n2, n_pairs = n_pairs,
       max_pp_h4 = max_h4)
}

out <- vector("list", 0)
for (t in req[["tasks"]]) {
  res <- tryCatch(
    if (t[["method"]] == "abf") run_abf(t) else run_susie(t),
    error = function(e) list(id = t[["id"]], ok = FALSE, method = t[["method"]], error = conditionMessage(e)))
  out[[length(out) + 1]] <- res
}
session <- list(R = R.version.string, coloc = as.character(packageVersion("coloc")),
                susieR = as.character(packageVersion("susieR")),
                jsonlite = as.character(packageVersion("jsonlite")))
writeLines(toJSON(list(tasks = out, session = session), auto_unbox = TRUE, digits = NA, na = "null"), args[2])
