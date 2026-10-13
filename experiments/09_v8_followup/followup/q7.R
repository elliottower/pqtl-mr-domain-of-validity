# Q7 (PREREG §Statistical models, Q7): proportional-hazards test, Aalen-Johansen cumulative incidence
# and Fine-Gray, with R's survival package (Debian bookworm apt). Every output is a CSV.
#   Rscript q7.R <survival.csv> <terms_H1.txt> <terms_H2.txt> <out_dir>
library(survival)

args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(args[1], stringsAsFactors = FALSE)
terms <- list(H1 = readLines(args[2]), H2 = readLines(args[3]))
out <- args[4]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
d$fstat <- factor(d$status, levels = 0:2, labels = c("censor", "phase3", "competing"))

write.csv(data.frame(R = R.version.string, survival = as.character(packageVersion("survival"))),
          file.path(out, "environment.csv"), row.names = FALSE)

for (model in names(terms)) {
  sub <- if (model == "H1") d[d$cls %in% c("aligned", "blocking"), ] else d
  rhs <- paste(terms[[model]], collapse = " + ")

  # (a) Grambsch-Therneau test: scaled Schoenfeld residuals against rank of time, Efron ties
  fit <- coxph(as.formula(paste("Surv(time, status == 1) ~", rhs)), data = sub, ties = "efron")
  z <- cox.zph(fit, transform = "rank")
  write.csv(data.frame(term = rownames(z$table), z$table, row.names = NULL),
            file.path(out, paste0("ph_test_", model, ".csv")), row.names = FALSE)
  cf <- summary(fit)$coefficients
  write.csv(data.frame(term = rownames(cf), cf, row.names = NULL),
            file.path(out, paste0("cox_efron_", model, ".csv")), row.names = FALSE)

  # (c) Fine-Gray: finegray expansion, then a weighted Cox fit with gene-clustered robust variance
  fg <- finegray(Surv(time, fstat) ~ ., data = sub[, c("time", "fstat", terms[[model]], "gene_ensembl")],
                 etype = "phase3")
  if (!("fgwt" %in% names(fg))) stop("finegray output has no fgwt column")
  fgfit <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", rhs, "+ cluster(gene_ensembl)")),
                 data = fg, weights = fgwt, ties = "efron")
  s <- summary(fgfit)
  cf <- s$coefficients
  ci <- s$conf.int
  write.csv(data.frame(term = rownames(cf), log_shr = cf[, "coef"], robust_se = cf[, "robust se"],
                       shr = ci[, "exp(coef)"], shr_lower95 = ci[, "lower .95"], shr_upper95 = ci[, "upper .95"],
                       p = cf[, ncol(cf)], variance = "gene-clustered robust", row.names = NULL),
            file.path(out, paste0("fine_gray_", model, ".csv")), row.names = FALSE)
}

# (b) Aalen-Johansen cumulative incidence by class x S, all dated S1 hypotheses
aj <- survfit(Surv(time, fstat) ~ group, data = d)
years <- c(0, 2, 4, 6, 8, 10, 12, 15)
sm <- summary(aj, times = years * 365.25, extend = TRUE)
states <- colnames(sm$pstate)
write.csv(data.frame(group = sub("group=", "", as.character(sm$strata)), years = sm$time / 365.25,
                     n_risk = sm$n.risk[, 1],
                     cif_phase3 = sm$pstate[, which(states == "phase3")],
                     cif_competing = sm$pstate[, which(states == "competing")]),
          file.path(out, "aalen_johansen_risk_table.csv"), row.names = FALSE)
full <- summary(aj)
fstates <- colnames(full$pstate)
write.csv(data.frame(group = sub("group=", "", as.character(full$strata)), days = full$time,
                     cif_phase3 = full$pstate[, which(fstates == "phase3")],
                     cif_competing = full$pstate[, which(fstates == "competing")]),
          file.path(out, "aalen_johansen_curves.csv"), row.names = FALSE)
