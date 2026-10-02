#!/usr/bin/env Rscript
# Stage B colocalization backend. Called by stage_b/coloc_backend.py as
#   Rscript coloc_run.R request.json response.json
# Each task is coloc.abf or coloc.susie on two datasets that share the same snp vector.
# Nothing here chooses a threshold or a window; the Python side does that from PREREG.md.
suppressPackageStartupMessages({
  library(jsonlite)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: coloc_run.R request.json response.json")
req <- fromJSON(args[1], simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = TRUE)

as_dataset <- function(d) {
  out <- list(snp = as.character(d$snp), beta = as.numeric(d$beta), varbeta = as.numeric(d$varbeta),
              N = as.numeric(d$N), type = as.character(d$type))
  if (!is.null(d$MAF)) out$MAF <- as.numeric(d$MAF)
  if (!is.null(d$sdY)) out$sdY <- as.numeric(d$sdY)
  if (!is.null(d$s)) out$s <- as.numeric(d$s)
  out
}

run_abf <- function(t) {
  r <- coloc.abf(as_dataset(t$d1), as_dataset(t$d2), p1 = t$p1, p2 = t$p2, p12 = t$p12)
  s <- r$summary
  list(id = t$id, ok = TRUE, method = "abf",
       pp = unname(as.numeric(s[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])),
       nsnps = unname(as.numeric(s["nsnps"])),
       snp = as.character(r$results$snp), snp_pp_h4 = as.numeric(r$results$SNP.PP.H4))
}

run_susie <- function(t) {
  ld <- as.matrix(t$LD)
  snp <- as.character(t$d1$snp)
  dimnames(ld) <- list(snp, snp)
  d1 <- as_dataset(t$d1); d1$LD <- ld
  d2 <- as_dataset(t$d2); d2$LD <- ld
  s1 <- runsusie(d1)
  s2 <- runsusie(d2)
  n1 <- length(s1$sets$cs)
  n2 <- length(s2$sets$cs)
  max_h4 <- NA_real_
  n_pairs <- 0
  if (n1 > 0 && n2 > 0) {
    r <- coloc.susie(s1, s2, p1 = t$p1, p2 = t$p2, p12 = t$p12)
    if (!is.null(r$summary) && nrow(r$summary) > 0) {
      max_h4 <- max(r$summary$PP.H4.abf)
      n_pairs <- nrow(r$summary)
    }
  }
  list(id = t$id, ok = TRUE, method = "susie", n_cs1 = n1, n_cs2 = n2, n_pairs = n_pairs,
       max_pp_h4 = max_h4)
}

out <- vector("list", 0)
for (t in req$tasks) {
  res <- tryCatch(
    if (t$method == "abf") run_abf(t) else if (t$method == "susie") run_susie(t) else stop("unknown method"),
    error = function(e) list(id = t$id, ok = FALSE, method = t$method, error = conditionMessage(e)))
  out[[length(out) + 1]] <- res
}
session <- list(R = R.version.string, coloc = as.character(packageVersion("coloc")),
                susieR = as.character(packageVersion("susieR")),
                jsonlite = as.character(packageVersion("jsonlite")))
writeLines(toJSON(list(tasks = out, session = session), auto_unbox = TRUE, digits = NA, na = "null"), args[2])
