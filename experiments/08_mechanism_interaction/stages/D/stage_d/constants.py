"""Every number stage D takes from the frozen plan (PREREG.md, freeze b946087, OSF 9tzfk).

Section references are to PREREG.md. Nothing here is tuned; each value is quoted from the plan.
"""
from math import log

from pydantic import BaseModel, ConfigDict

PLAN_SHA256 = "6cbba0840b213f586ef455b33e569a84a94fc676103b74adef6c46e2e0e70040"
FREEZE_COMMIT = "b946087070e4"
OSF_REGISTRATION = "https://osf.io/9tzfk/"

SEED = 20260930

# §Inference criteria
PROB_THRESHOLD = 0.95
LOG_SESOI = log(1.5)
CREDIBLE_MASS = 0.90

# §Statistical models: priors
PRIOR_SD_FIXED = 1.5
PRIOR_SD_INTERCEPT = 2.5
PRIOR_SD_RANDOM = 1.0
PRIOR_SENSITIVITY_NORMAL_SD = 1.0
PRIOR_SENSITIVITY_T_NU = 3.0
PRIOR_SENSITIVITY_T_SCALE = 1.5

# §Statistical models, Identification
PRIOR_DOMINANCE_FRACTION = 0.8

# §Inference criteria, Estimability and reliability
RELIABILITY_MIN_SUPPORTIVE_GENES = 10
RELIABILITY_MIN_SUPPORTIVE_HYPOTHESES = 20

# §Statistical models, H3
H3_MIN_EACH = 20

# §Statistical models: sampler requirements
RHAT_MAX = 1.01
ESS_BULK_MIN = 400.0


class SamplerSettings(BaseModel):
    model_config = ConfigDict(frozen=True)
    chains: int
    draws: int
    tune: int
    target_accept: float


PRIMARY_SAMPLER = SamplerSettings(chains=4, draws=2000, tune=2000, target_accept=0.95)
RERUN_SAMPLER = SamplerSettings(chains=4, draws=4000, tune=2000, target_accept=0.99)

# §Statistical models, frequentist check; power/power_v9 null grid
N_BOOTSTRAP = 10_000
N_NULL_REPS = 1000
NULL_BASE_ADVANCE = 0.30
NULL_RE_SDS = (0.3, 0.7, 1.0)
NULL_SLOPE_SDS = (0.3, 0.7)
NULL_MAIN_ORS = (0.8, 1.0, 1.25)
NULL_REFERENCE_CELL = (0.7, 0.7, 1.0)
BOOTSTRAP_CHUNK = 1000
NULL_CHUNK = 100
MAX_ABS_LOGIT_COEF = 15.0   # power_v9_fast.run_rep_fast failure rule

# §Other planned analysis, table 11
TABLE11_MIN_HYPOTHESES = 5

# §Explanation of foreknowledge, item 5 (Karim et al. 2026 launched pQTL-supported targets)
KARIM_LAUNCHED_TARGETS = ("APOB", "CSF3", "CSF3R", "EGLN1", "F2", "FLT3", "IL12B", "IL4R", "KIT",
                          "PCSK9", "SERPINA1", "VWF")

# §Other planned analysis: forest plot sets
FOREST_MAIN_SETS = ("S1", "S2", "S3", "S4", "S6", "S7", "S8", "S13")
