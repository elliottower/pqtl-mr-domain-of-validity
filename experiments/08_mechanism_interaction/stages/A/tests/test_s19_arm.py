"""S19 arm (PREREG §Other planned analysis: "neutralizing biologics vs small-molecule blockers of
blood-secreted proteins"), through the whole stage A pipeline on a synthetic world. GA is
blood-secreted (HPA) with a UKB-PPP instrument; GB is not blood-secreted."""
import pandas as pd
import pytest

from stage_a.build import run_stage_a
from stage_a.mechanism import ROW1, ROW2, ROW3, ROW4, hypothesis_s19_arm, s19_arm
from tests.synth import CV, World, default_clinical, default_diseases, default_gwas_catalog, \
    default_mechanisms, default_molecules, disease_row, gcat_study, moa, world_inputs

# indication: (drugs as (id, molecule type, action type, mechanism text), expected class, expected arm)
CASES = {
    "MONDO_9101": ([("OL12", "Oligonucleotide", "RNAI INHIBITOR", "GA siRNA")], "aligned", ""),
    "MONDO_9102": ([("DG13", "Small molecule", "DEGRADER", "GA degrader")], "aligned", ""),
    "MONDO_9103": ([("AB14", "Antibody", "BINDING AGENT", "neutralizing antibody against GA")], "aligned",
                   "neutralizing_biologic"),
    "MONDO_9104": ([("SM15", "Small molecule", "INHIBITOR", "GA inhibitor")], "blocking", "small_molecule_blocker"),
    "MONDO_9105": ([("AB16", "Antibody", "INHIBITOR", "GA inhibitor"),
                    ("OL17", "Oligonucleotide", "ANTISENSE INHIBITOR", "GA antisense")], "aligned", ""),
    "MONDO_9106": ([("AB18", "Antibody", "ANTAGONIST", "GA antagonist"),
                    ("PR19", "Protein", "INHIBITOR", "GA decoy receptor")], "aligned", "neutralizing_biologic"),
}


@pytest.fixture(scope="module")
def hyp() -> pd.DataFrame:
    drugs = [d for ds, _, _ in CASES.values() for d in ds]
    world = World(
        diseases=default_diseases() + [disease_row(i, f"disease {i}", tas=[CV]) for i in CASES],
        mechanisms=default_mechanisms() + [moa(d, action, text, "GA") for d, _, action, text in drugs],
        molecules=default_molecules() + [{"id": d, "drugType": t, "parentId": None} for d, t, _, _ in drugs],
        clinical=default_clinical() + [(d[0], ind, "PHASE_2") for ind, (ds, _, _) in CASES.items() for d in ds],
        gwas_catalog=default_gwas_catalog() + [gcat_study(f"GCST{i[-4:]}", i, "EstBB", 4000, 40000) for i in CASES],
    )
    return run_stage_a(world_inputs(world)).hypotheses


def one(hyp: pd.DataFrame, gene: str, ind: str) -> pd.Series:
    rows = hyp[(hyp["gene_symbol"] == gene) & (hyp["indication_id"] == ind) & (hyp["variant"] == "primary")]
    assert len(rows) == 1, rows
    return rows.iloc[0]


@pytest.mark.parametrize("ind", list(CASES))
def test_arm_of_each_planted_hypothesis(hyp, ind):
    drugs, cls, arm = CASES[ind]
    r = one(hyp, "GA", ind)
    assert (r["mechanism_class"], r["s19_arm"]) == (cls, arm)
    assert r["drug_program_ids"] == ";".join(sorted(d[0] for d in drugs))
    assert r["in_s1"] and r["heldout"] and r["blood_secreted_hpa"]


def test_existing_rows_antibody_on_blood_secreted_and_small_molecule_on_non_secreted(hyp):
    assert one(hyp, "GA", "MONDO_1001")["s19_arm"] == "neutralizing_biologic"   # AB1 antibody inhibitor
    gb = hyp[(hyp["gene_symbol"] == "GB") & (hyp["indication_id"] == "MONDO_1001") & (hyp["variant"] == "primary")]
    blocking = gb[gb["mechanism_class"] == "blocking"].iloc[0]
    other = gb[gb["mechanism_class"] == "other"].iloc[0]
    assert blocking["s19_arm"] == ""          # SM2 blocks GB, which is not blood-secreted
    assert other["s19_arm"] == ""             # AB11: antibody against a non-secreted target is other


def test_arm_is_empty_on_every_other_row(hyp):
    other = hyp[hyp["mechanism_class"] == "other"]
    assert (other["s19_arm"] == "").all()
    arms = hyp[hyp["s19_arm"] != ""]
    assert (arms["blood_secreted_hpa"] | (arms["variant"] == "s9")).all()
    assert set(zip(arms["mechanism_class"], arms["s19_arm"])) <= {("aligned", "neutralizing_biologic"),
                                                                  ("blocking", "small_molecule_blocker")}


@pytest.mark.parametrize("reason,blood,arm", [
    (ROW1, True, "neutralizing_biologic"),
    (ROW2, True, ""),
    (ROW3, True, ""),
    (ROW4, True, "small_molecule_blocker"),
    (ROW4, False, ""),
    (f"{ROW1};{ROW3}", True, ""),
    ("other:biologic_target_not_blood_secreted", False, ""),
])
def test_s19_arm_of_a_program_target(reason, blood, arm):
    assert s19_arm(reason, blood) == arm


@pytest.mark.parametrize("programs,arm", [
    ({"P1": {"neutralizing_biologic"}, "P2": {"neutralizing_biologic"}}, "neutralizing_biologic"),
    ({"P1": {"neutralizing_biologic"}, "P2": {""}}, ""),
    ({"P1": {"small_molecule_blocker"}}, "small_molecule_blocker"),
    ({"P1": {"neutralizing_biologic", ""}}, ""),
    ({"P1": {"neutralizing_biologic"}, "P2": {"small_molecule_blocker"}}, ""),
])
def test_hypothesis_arm_requires_every_program_in_the_arm(programs, arm):
    assert hypothesis_s19_arm(programs) == arm
