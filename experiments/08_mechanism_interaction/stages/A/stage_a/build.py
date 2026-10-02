"""Stage A pipeline over loaded tables: hypotheses, eligibility, instruments, outcome-GWAS
selection, related-indication rule, sets, flags, funnel and descriptive tables.

Everything here is a pure function of `StageAInputs`; file reading is in `loaders.py`.
Phase enters only as the boolean `phase2plus` on clinical_indication rows (the loader reduces
`maxClinicalStage` to it on read) and leaves only as the boolean `eligible`.
"""
import hashlib
from collections import Counter
from datetime import date

import pandas as pd
from pydantic import BaseModel, ConfigDict

from stage_a.flags import earliest_publication
from stage_a.indications import NOT_IN_INDEX, PASS, TherapeuticAreaFlags, descendants_map, \
    indication_rule, specificity_rank, therapeutic_area_map
from stage_a.instruments import PLATFORM, DecodeStrictDiagnostic, InstrumentLists
from stage_a.mechanism import hypothesis_s19_arm, label_row, label_row_s8, resolve_program_targets, s19_arm
from stage_a.outcome_gwas import indication_tiers, outcome_trait_coding, select_outcome
from stage_a.related import related_indication_rule, s12_one_per_program
from stage_a.schemas import FunnelRow, HypothesisRow, MechanismCrosstabRow, OutcomeSelectionRow, OutcomeTraitCodingRow

KEY = ["gene", "indication_id", "direction", "mechanism_class"]
SINGLE_PROTEIN = "single protein"
ABSENT_TA = TherapeuticAreaFlags(neuro_psych=False, neuro_only=False, psych_only=False, oncology=False)


class StageAInputs(BaseModel):
    """Loaded inputs. Column contracts:
    clinical_indication: drugId, diseaseId, phase2plus (bool; no phase value).
    mechanisms: one row per (MoA row, drug, target): moa_drug, ensembl, mechanismOfAction,
        actionType, targetType.
    molecules: id, drugType, parentId.
    disease: id, name, exactSynonyms, descendants, ancestors, obsoleteTerms, therapeuticAreas,
        ontology, dbXRefs.
    candidates: outcome_gwas.combine_candidates output.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    clinical_indication: pd.DataFrame
    mechanisms: pd.DataFrame
    molecules: pd.DataFrame
    disease: pd.DataFrame
    ens2sym: dict[str, str]
    hpa_blood_ensembl: set[str]
    uniprot_secreted_symbols: set[str]
    strict_lists: InstrumentLists
    inclusive_lists: InstrumentLists
    candidates: pd.DataFrame
    pilot_keys: set[tuple[str, str]]
    pilot_indications: set[str]
    karim_keys: set[tuple[str, str]]
    publication_dates: dict[str, date]
    decode_diagnostic: DecodeStrictDiagnostic


class StageAOutputs(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    hypotheses: pd.DataFrame
    funnel: pd.DataFrame
    outcome_gwas_selection: pd.DataFrame
    mechanism_crosstab: pd.DataFrame
    outcome_trait_coding: pd.DataFrame


class DuplicateHypothesisId(ValueError):
    """Two hypothesis keys hashed to the same 16-hex id."""


# Sets whose hypotheses are re-formed: S8 (antibody antagonists of non-blood-secreted targets),
# S9 (UniProt Secreted), S21 (conflicted program-target rows restored, each row its own class).
VARIANTS = ("s8", "s9", "s21")


def hypothesis_id(gene_ensembl: str, indication_id: str, direction: str, mechanism_class: str,
                  variant: str = "primary") -> str:
    """sha256 of the key, first 16 hex. A row that exists only under the S8 or S9 class table
    appends `|s8` or `|s9`, so it never takes the id of a primary row with the same key but other
    programs."""
    key = f"{gene_ensembl}|{indication_id}|{direction}|{mechanism_class}"
    return hashlib.sha256((key if variant == "primary" else f"{key}|{variant}").encode()).hexdigest()[:16]


# ---- drug-target rows and hypotheses -----------------------------------------------------

def drug_target_rows(inp: StageAInputs) -> tuple[pd.DataFrame, dict[str, int]]:
    """Drug x target x indication rows with program, molecule type, localization and the
    per-row mechanism label. A drug's targets are its own mechanism rows, or its parent
    molecule's where it has none; the program is the parent molecule where present."""
    ci = inp.clinical_indication[["drugId", "diseaseId", "phase2plus"]]
    moa = inp.mechanisms[["moa_drug", "ensembl", "mechanismOfAction", "actionType", "targetType"]]
    moa = moa.dropna(subset=["moa_drug", "ensembl"]).drop_duplicates()  # exact duplicate rows collapse
    parent = {i: p for i, p in zip(inp.molecules["id"], inp.molecules["parentId"]) if isinstance(p, str)}
    dtype = dict(zip(inp.molecules["id"], inp.molecules["drugType"]))
    own = set(moa["moa_drug"])
    ci = ci.assign(moa_drug=[d if d in own else (parent.get(d) if parent.get(d) in own else None)
                             for d in ci["drugId"]])
    n_ci = len(ci)
    dt = ci.dropna(subset=["moa_drug"]).merge(moa, on="moa_drug", how="inner")
    dt["gene"] = dt["ensembl"].map(inp.ens2sym)
    dt = dt.dropna(subset=["gene"]).reset_index(drop=True)
    dt["drug_type"] = dt["drugId"].map(dtype).fillna(dt["moa_drug"].map(dtype)).fillna("Unknown")
    dt["program"] = [parent.get(d, d) for d in dt["drugId"]]
    dt["blood"] = dt["ensembl"].isin(inp.hpa_blood_ensembl)
    dt["uniprot_secreted"] = dt["gene"].isin(inp.uniprot_secreted_symbols)
    rows = list(zip(dt["drug_type"], dt["actionType"], dt["mechanismOfAction"], dt["blood"], dt["uniprot_secreted"]))
    # Primary table; S8 (antibody antagonists of non-blood-secreted targets -> blocking); S9 (the
    # same table with UniProt Secreted in place of HPA blood-secreted).
    primary = [label_row(t, a, m, b) for t, a, m, b, _ in rows]
    for prefix, labels in (("row", primary),
                           ("s8", [label_row_s8(t, a, m, b) for t, a, m, b, _ in rows]),
                           ("s9", [label_row(t, a, m, u) for t, a, m, _, u in rows])):
        dt[f"{prefix}_class"] = [x.mechanism_class for x in labels]
        dt[f"{prefix}_direction"] = [x.direction for x in labels]
        dt[f"{prefix}_reason"] = [x.reason for x in labels]
    dt["neutralizing_ba"] = [x.neutralizing_binding_agent for x in primary]
    dt = dt.rename(columns={"diseaseId": "indication_id"})
    counts = {"clinical_indication_rows": n_ci,
              "clinical_indication_rows_with_mechanism": int(ci["moa_drug"].notna().sum()),
              "drug_target_indication_rows": len(dt)}
    return dt, counts


def mechanism_row_table(dt: pd.DataFrame) -> pd.DataFrame:
    """Distinct mechanism rows per (program, gene)."""
    return dt[["program", "gene", "drug_type", "actionType", "mechanismOfAction", "targetType",
               "row_class", "row_direction", "row_reason"]].drop_duplicates()


def resolve_for_hypotheses(dt: pd.DataFrame) -> pd.DataFrame:
    """One class and direction per (program, gene), from the mechanism rows of the program-target's
    Phase II+ drug x indication rows: count_v4.py keeps only Phase II+ clinical_indication rows
    (line 595) before it derives the rows it resolves (lines 630-632), so a child molecule that
    reached Phase II+ for no indication adds no mechanism row to its parent's program. A
    program-target with no Phase II+ row at all forms no eligible hypothesis; it is resolved from
    its own rows, for the funnel's count of ineligible hypotheses only."""
    phase2 = resolve_program_targets(mechanism_row_table(dt[dt["phase2plus"]]))
    every = resolve_program_targets(mechanism_row_table(dt))
    have = set(zip(phase2["program"], phase2["gene"]))
    rest = every[[k not in have for k in zip(every["program"], every["gene"])]]
    return pd.concat([phase2, rest], ignore_index=True)


def relabelled(dt: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """The drug-target rows with the `<prefix>_*` labels (s8, s9) in place of the primary
    row_class / row_direction / row_reason, for re-forming hypotheses under that table."""
    return dt.assign(row_class=dt[f"{prefix}_class"], row_direction=dt[f"{prefix}_direction"],
                     row_reason=dt[f"{prefix}_reason"])


def assign_program_targets(dt: pd.DataFrame, resolved: pd.DataFrame, restore_conflicts: bool,
                           localization: str = "blood") -> pd.DataFrame:
    """Per-row class and direction for hypothesis formation. Default: the program-target's
    resolved class. `restore_conflicts` (S21): a conflicted program-target's rows each keep their
    own class and direction. `s19_arm` comes from the same reasons, with `localization` the
    column of the table's blood-secreted definition (S9: uniprot_secreted)."""
    d = dt.merge(resolved, on=["program", "gene"], how="left", validate="many_to_one")
    d["restored"] = d["conflict"] & restore_conflicts
    d["mechanism_class"] = d["row_class"].where(d["restored"], d["pt_class"])
    d["direction"] = d["row_direction"].where(d["restored"], d["pt_direction"])
    d["reason"] = d["row_reason"].where(d["restored"], d["pt_reason"])
    d["s19_arm"] = [s19_arm(r, bool(b)) for r, b in zip(d["reason"], d[localization])]
    return d


def build_hypotheses(assigned: pd.DataFrame) -> pd.DataFrame:
    """Hypotheses: gene + indication + direction + class. `eligible` if any drug row is Phase II+;
    programs, single-protein flag and class reason are taken over the Phase II+ rows."""
    rows = []
    for key, g in assigned.groupby(KEY, sort=True):
        p2 = g[g["phase2plus"]]
        basis = p2 if len(p2) else g
        rows.append({
            **dict(zip(KEY, key)),
            "gene_ensembl": basis["ensembl"].iloc[0],
            "eligible": bool(len(p2)),
            "programs": frozenset(p2["program"]),
            "single_protein_row": bool((p2["targetType"] == SINGLE_PROTEIN).any()),
            "class_reason": ";".join(sorted(set(basis["reason"]))),
            "restored": bool(basis["restored"].any()),
            "s19_arm": hypothesis_s19_arm({p: set(a) for p, a in p2.groupby("program")["s19_arm"]}) if len(p2) else "",
        })
    return pd.DataFrame(rows, columns=[*KEY, "gene_ensembl", "eligible", "programs", "single_protein_row",
                                       "class_reason", "restored", "s19_arm"])


def annotate(hyp: pd.DataFrame, inp: StageAInputs, ind_status: dict[str, str], rank: dict[str, int],
             ta: dict[str, TherapeuticAreaFlags], names: dict[str, str]) -> pd.DataFrame:
    h = hyp.copy()
    h["ind_rule"] = h["indication_id"].map(ind_status).fillna(NOT_IN_INDEX)
    h["indication_name"] = h["indication_id"].map(names).fillna("")
    h["specificity_rank"] = h["indication_id"].map(rank)
    keys = list(zip(h["gene"], h["indication_id"]))
    h["heldout"] = [k not in inp.pilot_keys for k in keys]
    h["pilot_indication"] = h["indication_id"].isin(inp.pilot_indications)
    h["karim_launched"] = [k in inp.karim_keys for k in keys]
    h["blood_secreted_hpa"] = h["gene_ensembl"].isin(inp.hpa_blood_ensembl)
    h["secreted_uniprot"] = h["gene"].isin(inp.uniprot_secreted_symbols)
    flags = [ta.get(i, ABSENT_TA) for i in h["indication_id"]]
    for f in ("neuro_psych", "neuro_only", "psych_only", "oncology"):
        h[f] = [getattr(x, f) for x in flags]
    return h


# ---- variant: instruments, selection, overlap, related rule -------------------------------

class VariantResult(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    lists: InstrumentLists
    retained: pd.DataFrame      # after the related-indication rule on the tier-1 / no pool (held out or not)
    s4_added: pd.DataFrame
    s10_added: pd.DataFrame
    funnel: list[FunnelRow]


def _reason(counter: Counter) -> str:
    return ";".join(f"{k}={v}" for k, v in sorted(counter.items(), key=lambda kv: str(kv[0])))


def run_variant(hyp: pd.DataFrame, lists: InstrumentLists, tiers: dict[str, dict[int, list[int]]],
                cands: pd.DataFrame, descendants: dict[str, set[str]], rank: dict[str, int],
                prefix: str) -> VariantResult:
    """`hyp`: eligible, annotated hypotheses. Applies the indication rule, instrument lists,
    joint outcome selection, tier and overlap steps and the related-indication rule, with the
    S4 and S10 pools re-run through the related rule."""
    funnel: list[FunnelRow] = []
    cur = hyp

    def step(name: str, mask: pd.Series, reason: str) -> None:
        nonlocal cur
        funnel.append(FunnelRow(step=f"{prefix}{name}", remaining=int(mask.sum()),
                                excluded=int((~mask).sum()), reason=reason))
        cur = cur[mask]

    step("indication_rule", cur["ind_rule"] == PASS, _reason(Counter(cur.loc[cur["ind_rule"] != PASS, "ind_rule"])))
    cur = cur.assign(sources=[lists.sources_for(g, e) for g, e in zip(cur["gene"], cur["gene_ensembl"])])
    step("cis_pqtl_instrument", cur["sources"].map(bool),
         "gene on no strict UKB-PPP, deCODE or INTERVAL list" if prefix == "" else "gene on no list of this variant")
    selections = {}
    for g, d, srcs in cur[["gene", "indication_id", "sources"]].drop_duplicates(["gene", "indication_id"]).itertuples(index=False):
        selections[(g, d)] = select_outcome(srcs, tiers[d], cands)
    sel = [selections[k] for k in zip(cur["gene"], cur["indication_id"])]
    cur = cur.assign(sel_source=[s.source for s in sel], sel_candidate=[s.candidate for s in sel],
                     sel_tier=[s.tier for s in sel], sel_overlap=[s.overlap for s in sel])
    after_instrument = cur
    tier_label = cur["sel_tier"].map(lambda t: f"tier{int(t)}" if pd.notna(t) else "no_candidate")
    step("outcome_gwas_tier1", cur["sel_tier"] == 1, _reason(Counter(tier_label[cur["sel_tier"] != 1])))
    after_tier1 = cur
    step("overlap_no", cur["sel_overlap"] == "no", _reason(Counter(cur.loc[cur["sel_overlap"] != "no", "sel_overlap"])))
    keep = related_indication_rule(cur, descendants, rank)
    step("related_indication_rule", keep,
         "ontology ancestor whose Phase II+ programs all sit on retained descendant hypotheses")
    retained = cur
    step("heldout", cur["heldout"], "gene+indication key in classification_v5 / classification_v5_1")

    pool4 = after_tier1[after_tier1["sel_overlap"].isin(["no", "unknown"])]
    k4 = related_indication_rule(pool4, descendants, rank)
    s4_added = pool4[k4 & (pool4["sel_overlap"] == "unknown") & pool4["heldout"]]
    funnel.append(FunnelRow(step=f"{prefix}S4_added_overlap_unknown", remaining=len(s4_added), excluded=0,
                            reason="tier-1 overlap-unknown hypotheses retained by the related rule on the no+unknown pool, held out"))

    pool10 = after_instrument[after_instrument["sel_tier"].isin([1, 2, 3]) & (after_instrument["sel_overlap"] == "no")]
    k10 = related_indication_rule(pool10, descendants, rank)
    s10_added = pool10[k10 & pool10["sel_tier"].isin([2, 3]) & pool10["heldout"]]
    funnel.append(FunnelRow(step=f"{prefix}S10_added_tier2_tier3", remaining=len(s10_added), excluded=0,
                            reason="tier-2/3 overlap-no hypotheses retained by the related rule on the tier-1/2/3 overlap-no pool, held out"))
    return VariantResult(lists=lists, retained=retained, s4_added=s4_added, s10_added=s10_added, funnel=funnel)


# ---- assembly ----------------------------------------------------------------------------

def _keys(df: pd.DataFrame) -> set[tuple]:
    return set(zip(*(df[c] for c in KEY)))


def _row(r: dict, lists: InstrumentLists, cands: pd.DataFrame, dates: dict[str, date], *,
         in_s1: bool, in_s4: bool, in_s10: bool, in_s12: bool, strict: bool, restored: bool,
         variant: str = "primary", in_var: frozenset[str] = frozenset()) -> HypothesisRow:
    c = cands.loc[int(r["sel_candidate"])]
    tier, overlap = int(r["sel_tier"]), r["sel_overlap"]
    assays = lists.source(r["sel_source"]).assays(r["gene"], r["gene_ensembl"])
    ukb_assays = lists.ukbppp.assays(r["gene"], r["gene_ensembl"])
    decode_assays = frozenset() if c["includes_iceland"] == "yes" else lists.decode.assays(r["gene"], r["gene_ensembl"])
    return HypothesisRow(
        hypothesis_id=hypothesis_id(r["gene_ensembl"], r["indication_id"], r["direction"], r["mechanism_class"], variant),
        gene_symbol=r["gene"], gene_ensembl=r["gene_ensembl"],
        indication_id=r["indication_id"], indication_name=r["indication_name"],
        direction=r["direction"], mechanism_class=r["mechanism_class"], class_reason=r["class_reason"],
        drug_program_ids=";".join(sorted(r["programs"])), eligible=bool(r["eligible"]),
        heldout=bool(r["heldout"]), specificity_rank=int(r["specificity_rank"]),
        instrument_source=r["sel_source"], instrument_assay_id=";".join(sorted(assays)),
        platform=PLATFORM[r["sel_source"]],
        outcome_accession=str(c["accession"]), outcome_source=c["source"], outcome_tier=tier, overlap=overlap,
        outcome_n_case=int(round(float(c["ncase"]))), outcome_n_control=int(round(float(c["ncontrol"]))),
        outcome_neff=float(c["neff"]),
        blood_secreted_hpa=bool(r["blood_secreted_hpa"]), secreted_uniprot=bool(r["secreted_uniprot"]),
        neuro_psych=bool(r["neuro_psych"]), neuro_only=bool(r["neuro_only"]), psych_only=bool(r["psych_only"]),
        oncology=bool(r["oncology"]), pilot_indication=bool(r["pilot_indication"]),
        karim_launched=bool(r["karim_launched"]), single_protein_row=bool(r["single_protein_row"]),
        conflicted_row_restored=restored, strict_instrument=strict,
        in_s1=in_s1, in_s4=in_s4, in_s10=in_s10, in_s12=in_s12,
        subtype_restricted=tier == 2, phenotype_broader=tier == 3,
        sample_overlap=overlap == "yes", overlap_unknown=overlap == "unknown",
        pre_pqtl_publication_date=earliest_publication(r["sources"], dates),
        variant=variant, in_s8="s8" in in_var, in_s9="s9" in in_var, in_s21="s21" in in_var,
        ukbppp_assay_ids=";".join(sorted(ukb_assays)), decode_assay_ids=";".join(sorted(decode_assays)),
        s19_arm=r["s19_arm"],
    )


def outcome_selection_table(hyp_rows: pd.DataFrame, tiers: dict[str, dict[int, list[int]]],
                            cands: pd.DataFrame, names: dict[str, str]) -> pd.DataFrame:
    """Per indication of hypotheses.csv, every qualifying candidate in tiers 1-3; `selected`
    when the candidate is the selected outcome GWAS of at least one hypothesis row."""
    chosen = set(zip(hyp_rows["indication_id"], hyp_rows["outcome_source"], hyp_rows["outcome_accession"]))
    rows = []
    for d in sorted(set(hyp_rows["indication_id"])):
        for tier in (1, 2, 3):
            for i in tiers[d][tier]:
                c = cands.loc[i]
                rows.append(OutcomeSelectionRow(
                    indication_id=d, indication_name=names.get(d, ""), outcome_accession=str(c["accession"]),
                    outcome_source=c["source"], trait=str(c["trait"] if c["trait"] is not None else ""),
                    tier=tier, n_case=int(round(float(c["ncase"]))), n_control=int(round(float(c["ncontrol"]))),
                    neff=float(c["neff"]), includes_ukb=c["includes_ukb"], includes_iceland=c["includes_iceland"],
                    includes_interval=c["includes_interval"],
                    selected=(d, c["source"], str(c["accession"])) in chosen).model_dump())
    return pd.DataFrame(rows, columns=list(OutcomeSelectionRow.model_fields))


def mechanism_crosstab(assigned: pd.DataFrame, s1_keys: set[tuple]) -> pd.DataFrame:
    """Descriptive table 4 over the Phase II+ mechanism rows of S1 hypotheses."""
    a = assigned[assigned["phase2plus"]]
    a = a[[k in s1_keys for k in zip(*(a[c] for c in KEY))]]
    group = ["drug_type", "actionType", "neutralizing_ba", "blood", "row_class", "row_direction", "row_reason",
             "pt_class", "pt_direction", "conflict"]
    rows = []
    for key, g in a.groupby(group, sort=True, dropna=False):
        k = dict(zip(group, key))
        rows.append(MechanismCrosstabRow(
            molecule_type=str(k["drug_type"]), action_type=str(k["actionType"]),
            neutralizing_binding_agent=bool(k["neutralizing_ba"]), blood_secreted_hpa=bool(k["blood"]),
            row_class=k["row_class"], row_direction=k["row_direction"], row_reason=k["row_reason"],
            program_target_class=k["pt_class"], program_target_direction=k["pt_direction"],
            conflicted=bool(k["conflict"]),
            n_mechanism_rows=int(len(g[["program", "gene", "drug_type", "actionType", "mechanismOfAction",
                                        "targetType"]].drop_duplicates())),
            n_program_targets=int(len(g[["program", "gene"]].drop_duplicates())),
            n_hypotheses_s1=int(len(g[KEY].drop_duplicates()))).model_dump())
    return pd.DataFrame(rows, columns=list(MechanismCrosstabRow.model_fields))


def decode_diagnostic_rows(d: DecodeStrictDiagnostic) -> list[FunnelRow]:
    """Outcome-blind diagnostic rows on the strict deCODE list. They count, and exclude nothing:
    the list is the registered one (instruments.decode_lists) whatever these rows say."""
    return [
        FunnelRow(step="diagnostic:decode_strict_seqids_st02_symbols_differ_from_st01_gene",
                  remaining=d.strict_seqids_st02_differs_from_st01, excluded=0,
                  reason=f"diagnostic, changes no list: of {d.strict_seqids} strict deCODE SeqIds (ST01: one Gene, one "
                         "UniProt, at most one Ensembl ID), those whose ST02 cis rows name gene symbols other than "
                         "exactly the ST01 Gene"),
        FunnelRow(step="diagnostic:decode_strict_symbols_only_through_st02_disagreement",
                  remaining=d.strict_symbols_only_through_disagreement, excluded=0,
                  reason=f"diagnostic, changes no list: of {d.strict_symbols} gene symbols on the strict deCODE list, "
                         "those that are the ST01 Gene of no strict SeqId and are listed because an ST02 cis row "
                         "names them"),
    ]


def run_stage_a(inp: StageAInputs) -> StageAOutputs:
    dates = inp.publication_dates
    dt, universe = drug_target_rows(inp)
    resolved = resolve_for_hypotheses(dt)
    assigned = assign_program_targets(dt, resolved, restore_conflicts=False)
    assigned_s21 = assign_program_targets(dt, resolved, restore_conflicts=True)
    all_hyp = build_hypotheses(assigned)
    all_hyp_s21 = build_hypotheses(assigned_s21)

    ind_status = indication_rule(inp.disease)
    rank = specificity_rank(inp.disease)
    descendants = descendants_map(inp.disease)
    ta = therapeutic_area_map(inp.disease)
    names = dict(zip(inp.disease["id"], inp.disease["name"]))

    funnel = [
        FunnelRow(step="universe_clinical_indication_rows", remaining=universe["clinical_indication_rows"],
                  excluded=0, reason="Open Targets 26.09 clinical_indication rows (drug x disease)"),
        FunnelRow(step="universe_rows_with_mechanism", remaining=universe["clinical_indication_rows_with_mechanism"],
                  excluded=universe["clinical_indication_rows"] - universe["clinical_indication_rows_with_mechanism"],
                  reason="drug and its parent molecule have no drug_mechanism_of_action row"),
        FunnelRow(step="drug_target_indication_rows", remaining=universe["drug_target_indication_rows"], excluded=0,
                  reason="unit: drug x HGNC-approved target gene x indication; targets without an HGNC symbol dropped"),
        FunnelRow(step="gene_indication_pairs", remaining=int(len(dt[["gene", "indication_id"]].drop_duplicates())),
                  excluded=0, reason="unit: gene x indication"),
        FunnelRow(step="hypotheses", remaining=len(all_hyp), excluded=0,
                  reason="unit from here: gene x indication x direction x mechanism class"),
        FunnelRow(step="eligible_phase2_or_later", remaining=int(all_hyp["eligible"].sum()),
                  excluded=int((~all_hyp["eligible"]).sum()), reason="no drug at Phase II or later for the indication"),
    ]
    hyp = annotate(all_hyp[all_hyp["eligible"]], inp, ind_status, rank, ta, names)
    hyp_s21 = annotate(all_hyp_s21[all_hyp_s21["eligible"]], inp, ind_status, rank, ta, names)
    # S8 and S9 re-form hypotheses from their own class tables, with the conflicting-rows rule
    # applied within each table exactly as in the primary.
    hyp_var: dict[str, pd.DataFrame] = {}
    hyp_var["s21"] = hyp_s21
    for v in ("s8", "s9"):
        dtv = relabelled(dt, v)
        av = build_hypotheses(assign_program_targets(dtv, resolve_for_hypotheses(dtv),
                                                     restore_conflicts=False,
                                                     localization="uniprot_secreted" if v == "s9" else "blood"))
        hyp_var[v] = annotate(av[av["eligible"]], inp, ind_status, rank, ta, names)

    candidates_needed = set()
    for h in (hyp, *hyp_var.values()):
        candidates_needed |= set(h.loc[h["ind_rule"] == PASS, "indication_id"])
    tiers = indication_tiers(candidates_needed, inp.candidates, inp.disease)

    main = run_variant(hyp, inp.strict_lists, tiers, inp.candidates, descendants, rank, prefix="")
    s20 = run_variant(hyp, inp.inclusive_lists, tiers, inp.candidates, descendants, rank, prefix="S20:")
    var = {v: run_variant(hyp_var[v], inp.strict_lists, tiers, inp.candidates, descendants, rank,
                          prefix=f"{v.upper()}:") for v in VARIANTS}

    # `in_s1` marks passing every S1 rule before the held-out filter, so the confirmatory S1 is
    # in_s1 & heldout and the pooled set S2 is in_s1; in_s4 and in_s10 are read the same way.
    # S12 is formed on the held-out S1 only.
    s1 = main.retained[main.retained["heldout"]]
    s1_keys = _keys(s1)
    s1_rule_keys = _keys(main.retained)
    s12_keys = _keys(s1[s12_one_per_program(s1, rank)])
    s4_keys = s1_rule_keys | _keys(main.s4_added)
    s10_keys = s1_rule_keys | _keys(main.s10_added)

    pending: dict[tuple, tuple[dict, InstrumentLists, bool, bool]] = {}
    selected: dict[tuple[str, str], int] = {}

    def add(frame: pd.DataFrame, lists: InstrumentLists, *, strict: bool, restored: bool) -> int:
        n = 0
        for r in frame.to_dict(orient="records"):
            k = tuple(r[c] for c in KEY)
            if k in pending:
                continue
            pending[k] = (r, lists, strict, restored)
            n += 1
        return n

    add(main.retained, main.lists, strict=True, restored=False)
    add(main.s4_added, main.lists, strict=True, restored=False)
    add(main.s10_added, main.lists, strict=True, restored=False)
    n20 = add(s20.retained[s20.retained["heldout"]], s20.lists, strict=False, restored=False)

    # A variant hypothesis identical to a written row (same key, programs, instrument source and
    # outcome selection, strict lists) is that row, flagged in_s8 / in_s9. Any other variant
    # hypothesis (a merged, split or reclassified one) is its own row, `variant` s8 or s9.
    shared: dict[str, set[tuple]] = {v: set() for v in VARIANTS}
    extra: dict[str, dict[tuple, dict]] = {v: {} for v in VARIANTS}
    for v in VARIANTS:
        for r in var[v].retained.to_dict(orient="records"):
            k = tuple(r[c] for c in KEY)
            prev = pending.get(k)
            if prev is not None and prev[2] and prev[0]["programs"] == r["programs"] \
                    and (prev[0]["sel_source"], prev[0]["sel_candidate"]) == (r["sel_source"], r["sel_candidate"]):
                shared[v].add(k)
            else:
                extra[v][k] = r

    rows: list[HypothesisRow] = []
    for k, (r, lists, strict, restored) in pending.items():
        rows.append(_row(r, lists, inp.candidates, dates, in_s1=k in s1_rule_keys, in_s4=k in s4_keys,
                         in_s10=k in s10_keys, in_s12=k in s12_keys, strict=strict, restored=restored,
                         in_var=frozenset(v for v in VARIANTS if k in shared[v])))
    for v in VARIANTS:
        for r in extra[v].values():
            rows.append(_row(r, var[v].lists, inp.candidates, dates, in_s1=False, in_s4=False, in_s10=False,
                             in_s12=False, strict=True, restored=v == "s21", variant=v, in_var=frozenset({v})))
    for r in [*(x[0] for x in pending.values()), *(r for v in VARIANTS for r in extra[v].values())]:
        c = inp.candidates.loc[int(r["sel_candidate"])]
        selected[(c["source"], str(c["accession"]))] = int(r["sel_candidate"])

    funnel += main.funnel
    funnel.append(FunnelRow(step="S12_one_hypothesis_per_program", remaining=len(s12_keys),
                            excluded=len(s1_keys) - len(s12_keys),
                            reason="within gene, a program already kept at a more specific indication"))
    funnel += s20.funnel
    funnel.append(FunnelRow(step="S20:rows_present_only_under_inclusive_lists", remaining=n20, excluded=0,
                            reason="written with strict_instrument = false"))
    for v in VARIANTS:
        funnel += var[v].funnel
        funnel.append(FunnelRow(step=f"{v.upper()}:rows_shared_with_primary", remaining=len(shared[v]), excluded=0,
                                reason="identical to a written row; flagged in_" + v))
        funnel.append(FunnelRow(step=f"{v.upper()}:rows_present_only_under_{v}", remaining=len(extra[v]), excluded=0,
                                reason=f"written with variant = {v}"))
    funnel.append(FunnelRow(step="hypotheses_csv_rows", remaining=len(rows), excluded=0,
                            reason="S1 and pilot keys surviving the V8 rules, plus S4, S10, S20, S21, S8, S9 additions"))
    funnel += decode_diagnostic_rows(inp.decode_diagnostic)

    hyp_rows = pd.DataFrame([r.model_dump() for r in rows], columns=list(HypothesisRow.model_fields))
    if hyp_rows["hypothesis_id"].duplicated().any():
        raise DuplicateHypothesisId(sorted(hyp_rows.loc[hyp_rows["hypothesis_id"].duplicated(), "hypothesis_id"]))
    hyp_rows = hyp_rows.sort_values("hypothesis_id").reset_index(drop=True)

    trait = pd.DataFrame([OutcomeTraitCodingRow(**r).model_dump()
                          for r in outcome_trait_coding(selected, inp.candidates).to_dict(orient="records")],
                         columns=list(OutcomeTraitCodingRow.model_fields))

    return StageAOutputs(
        hypotheses=hyp_rows,
        funnel=pd.DataFrame([f.model_dump() for f in funnel], columns=list(FunnelRow.model_fields)),
        outcome_gwas_selection=outcome_selection_table(hyp_rows, tiers, inp.candidates, names),
        mechanism_crosstab=mechanism_crosstab(assigned, s1_keys),
        outcome_trait_coding=trait,
    )
