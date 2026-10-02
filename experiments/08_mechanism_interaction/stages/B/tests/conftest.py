"""Every R process a stage B test starts reads this profile (R_PROFILE_USER): partial matching by
`$` is reported (`warnPartialMatchDollar`) and that warning is turned into an error, raised outside
any handler of the script, so Rscript exits non-zero and the test fails. coloc_run.R reads every
field with `[[ ]]`; a `$` that matched a longer name would stop the R tests here."""
import pytest

STRICT_R_PROFILE = """options(warnPartialMatchDollar = TRUE)
globalCallingHandlers(warning = function(w) {
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
