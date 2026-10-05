"""Every R process a stage B test starts reads this profile (R_PROFILE_USER): partial matching by
`$` is reported (`warnPartialMatchDollar`) and that warning is turned into an error, raised outside
any handler of the script, so Rscript exits non-zero and the test fails. coloc_run.R reads every
field with `[[ ]]`; a `$` that matched a longer name would stop the R tests here. The exception is
the one coloc_run.R makes, by exact message: coloc 5.2.3's own reads of its suffixed p-value-form
columns (df$pvalues, df$MAF, df$N to pvalues.df<i>, MAF.df<i>, N.df<i>)."""
import pytest

STRICT_R_PROFILE = """options(warnPartialMatchDollar = TRUE)
COLOC_OWN <- as.vector(outer(c("pvalues", "MAF", "N"), c("df1", "df2"),
                             function(f, i) sprintf("partial match of '%s' to '%s.%s'", f, f, i)))
globalCallingHandlers(warning = function(w) {
  if (conditionMessage(w) %in% COLOC_OWN) invokeRestart("muffleWarning")
  if (grepl("partial match of", conditionMessage(w), fixed = TRUE)) {
    stop(paste("partial matching by `$`:", conditionMessage(w)), call. = FALSE)
  }
})
"""


@pytest.fixture(autouse=True)
def strict_r(tmp_path_factory, monkeypatch):
    profile = tmp_path_factory.mktemp("r_profile") / "strict.Rprofile"
    profile.write_text(STRICT_R_PROFILE)
    monkeypatch.setenv("R_PROFILE_USER", str(profile))
