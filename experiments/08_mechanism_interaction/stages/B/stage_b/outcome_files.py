"""GWAS Catalog outcome files without a harmonised copy: GWAS-SSF files and reviewed author formats.

The frozen plan asks the GWAS Catalog for "full summary statistics". Where an accession's study
directory has no `harmonised/` file, its summary statistics are the file the study directory holds:
a GWAS-SSF data file `<accession>.tsv[.gz]` or `<accession>_build<GRCh37|GRCh38>.tsv[.gz]` with its
`-meta.yaml` and the directory's `md5sum.txt` (GWAS-SSF v1.0.2 and v1.1.0,
github.com/gwas-catalog/summary-statistics-standard), or a file in the author's own format for
which stage_b/author_formats.py holds a reviewed column map. This module turns either into the
canonical outcome table the harmonised path produces (schemas.VARIANT_COLUMNS), with the registered
rules unchanged: variants are later matched by rsID (harmonize.py), no liftover is applied,
palindromes are handled at harmonization, and beta is log(OR) where only an odds ratio is given.

Collect (`resolve_ssf`, from the `-meta.yaml`, `md5sum.txt` and the data file's header line, before
any byte of data is downloaded). The `-meta.yaml` is decoded as UTF-8 (a leading byte-order mark is
dropped; any other encoding is refused) and parsed by PyYAML's SafeLoader, the loader of
`yaml.safe_load`, with a mapping that repeats a key refused and a plain scalar read as the string it
is written as, or null (`parse_meta_yaml`); its top level must be a mapping, and the four keys read
must be strings where present (`SsfMeta`). The file is
`unreadable`, with the reason, when:
- the directory holds more than one GWAS-SSF data file, or the data file has no `-meta.yaml`, or
  the directory no `md5sum.txt`;
- the `-meta.yaml` is not UTF-8, does not parse, repeats a key, or fails `SsfMeta`;
- `genome_assembly` is not GRCh37 or GRCh38; `coordinate_system` is neither `1-based` nor `0-based`
  (absent is read as 1-based, a convention for legacy files that the metadata does not state, as
  GWAS-SSF makes the key optional and permits either: positions only place a variant in or out of the ±1 Mb and ±500 kb
  windows, and variants are matched by rsID, so the reading can move only a variant at a window
  edge); `data_file_name` is not the data file; md5sum.txt gives no MD5 for the data file; or
  `data_file_md5sum` and md5sum.txt disagree;
- the header lacks `chromosome`, `base_pair_location`, `effect_allele`, `other_allele` or
  `standard_error`, has neither `beta` nor `odds_ratio` (`hazard_ratio` and the `z-score` fallback
  are not read), has both and column 4 (the standard's effect column) is neither, or repeats a
  column name.
A downloaded data file whose MD5 is not md5sum.txt's raises CollectError (refused, nothing
recorded; collect.collect_one). The MD5 md5sum.txt lists for the `-meta.yaml` itself is not a
condition: the GWAS Catalog rewrites a study's metadata in place without updating md5sum.txt (8 of
the 66 deposits resolved on 2026-10-03, every one a pre-GWAS-SSF or non-GWAS-SSF file whose
`date_metadata_last_modified` is January 2025, serve a `-meta.yaml` with another MD5, while in all
66 the data file's MD5 in `-meta.yaml` equals md5sum.txt's). A difference is written to the
record's `detail`, and the sha256 of the `-meta.yaml` read is in `meta_sha256`.

Columns (GWAS-SSF field names). chromosome (23 = X, 24 = Y), base_pair_location (+1 for a 0-based
file), effect_allele and other_allele (surrounding whitespace removed, upper case), beta, or
ln(odds_ratio) where the file gives an odds ratio and no beta (both: column 4 decides),
standard_error (read as the standard error of that beta), effect_allele_frequency (missing if
absent), p_value, else 10^-neg_log_10_p_value (missing if neither), n (missing if absent). Values
`#NA`, `NA` and empty are missing; every field is read with its surrounding whitespace removed.

Uncertainty mode. A file gives the uncertainty of its effect in one of three ways
(schemas.UncertaintyMode), chosen once for the whole file by the pre-analysis validation
(`choose_uncertainty_mode`, stage_b/validate.py), never row by row, from the header and a census of
the `standard_error` field over the rows with a position and a valid pair of alleles
(`RowReader.se_state`): `native_se` where it is present on every such row; where it is missing on
every such row, `or_ci_derived_se` for an odds-ratio file with `ci_lower` and `ci_upper`, and
`pvalue_coloc` for a beta file with a p-value and an effect-allele frequency and no confidence
limits; no mode, and the file unreadable, where it is present on some rows and missing on others
(the file fails closed) or where no route fits. An author format is read in `native_se` only.
`or_ci_derived_se`: SE of ln(OR) = (ln U - ln L) / (2 z), z the normal quantile of the interval's
level (`ci_z_for_level`; 1.96 at 95%, the level read for these files: the GWAS-SSF metadata schema
has no key stating another). `pvalue_coloc`: no standard error is read or formed; coloc.abf takes the
p-value, the MAF of the file's own effect-allele frequency, N and the case fraction
(coloc_backend.outcome_dataset), and beta gives the effect direction only.

Rows (`SsfReader`, `AuthorReader`). A data row is read, or not read for the first of these reasons
(REJECT_REASONS, in this order), and the pre-analysis validation counts each reason over the whole
file (stage_b/validate.py):
    row_width                    the row has another number of fields than the header
    no_position                  chromosome or position missing
    unparseable_number           a position or a numeric field that is not a number
    invalid_allele               effect or other allele (author formats: either allele of the pair)
                                 not [ACGT]+, the standard's accepted values
    equal_alleles                effect allele equal to other allele (author formats: the pair equal)
    effect_allele_not_in_pair    author formats: the effect allele is neither allele of the pair
    beta_nonfinite               a beta that is missing or infinite
    or_nonpositive_or_nonfinite  an odds ratio that is missing, infinite, zero or negative
    se_nonpositive_or_nonfinite  native_se: a standard error that is missing, infinite, zero or negative
    ci_nonpositive_or_nonfinite  or_ci_derived_se: a confidence limit missing, infinite, zero or negative
    ci_not_increasing            or_ci_derived_se: the lower limit not below the upper limit
    or_outside_ci                or_ci_derived_se: the odds ratio outside [lower, upper] by more than
                                 half a unit in the last printed decimal of the odds ratio plus half a
                                 unit in the last printed decimal of the limit it crosses
    p_outside_0_1                an ordinary p-value outside [0, 1]
    neg_log10_p_negative         a negative -log10 p-value
    eaf_outside_0_1              an effect-allele frequency outside [0, 1]
    p_missing_pvalue_coloc       pvalue_coloc: no p-value
    p_zero_pvalue_coloc          pvalue_coloc: p = 0 (coloc's p-value form cannot take it without a floor)
    eaf_not_inside_0_1_pvalue_coloc  pvalue_coloc: an effect-allele frequency missing, 0 or 1
A missing p-value is kept as missing outside `pvalue_coloc`, and p = 0 is read in `native_se` and
`or_ci_derived_se`. -log10 p is read as p = 10^-x; a value too large for a double gives p = 0. A
missing effect-allele frequency is kept as missing outside `pvalue_coloc`. In the analysis, rows
outside the window are skipped before their values are read.

CI-versus-p check of an `or_ci_derived_se` file (`CiPCheck`), fixed before it is run: on every row
read with 0 < p < 1, the two-sided p of the Wald statistic ln(OR) / SE, SE from the confidence
limits, is compared with the published p; the row agrees when |log10 p_implied - log10 p_published|
<= CI_P_LOG10_TOLERANCE (0.1), or when the gap is within the published p's printed rounding (half a
unit in its last printed digit, on the log10 scale) where that is looser. The file passes when at
least CI_P_MIN_AGREEMENT_PERCENT (95%) of the checkable rows agree; no checkable row fails it. The
check returns counts and pass/fail only.

rsID (`rsid_rule`). `column`: the `rsid` column, else `rs_id`; a value that is not rs<digits> gives
no rsID (the variant is dropped at harmonization, as on every other path). `ukbppp_map`: a file
with neither column takes rsIDs from the UKB-PPP rsID map of the variant's chromosome (Synapse
syn51396727, the map the plan registers for rsID matching, collected and hashed like any other
whole file): the map row whose position on the file's build (`POS19` for GRCh37, `POS38` for
GRCh38: the map's own columns, no liftover) equals the variant's position and whose {REF, ALT}
equals the variant's {effect_allele, other_allele}. A variant with no such row, or with rows
naming more than one rsID, has no rsID. `variant_id` is never read as an rsID: the standard defines
it as chromosome_position_reference_alternate.

Integrity of an odds-ratio author file (`OrCiCheck`, AuthorFormat.integrity "or_ci_se"; GCST008226,
whose read-me does not state the scale of Standard_error). Fixed before the file is read: on every
row whose odds ratio, confidence limits and standard error are finite, with OR > 0 and SE > 0, the
printed limits must equal exp(ln(OR) - 1.96 SE) and exp(ln(OR) + 1.96 SE) within a relative 1%, or
within the rounding the printed decimals imply where that is looser (half a unit in the last printed
digit of the limit, plus the error that half a unit in the last printed digit of OR and of SE
carries into the expected limit). The file passes when at least 99% of the rows checked agree; a
header without the two confidence-limit columns, no row to check, or fewer than 99% agreeing makes it
unreadable. The check returns counts and pass/fail only; it runs in the pre-analysis validation
(stage_b/validate.py), on Modal.
"""
import hashlib
import math
import re
from collections.abc import Hashable, Iterable, Mapping
from decimal import Decimal
from statistics import NormalDist
from typing import ClassVar

import pandas as pd
import yaml
from pydantic import BaseModel, ConfigDict, ValidationError

from stage_b.author_formats import AuthorFormat, check_header
from stage_b.parsers import MISSING, RowReader
from stage_b.schemas import InputContractError, SourceUnreadable, UncertaintyMode

SSF_DATA = re.compile(r"(?P<acc>GCST\d+)(?:_build(?:GRCh37|GRCh38))?\.tsv(?:\.gz)?")
MD5SUM_NAME = "md5sum.txt"
META_SUFFIX = "-meta.yaml"
BUILDS = ("GRCh37", "GRCh38")
COORDINATE_OFFSET = {"1-based": 0, "0-based": 1}
REQUIRED = ("chromosome", "base_pair_location", "effect_allele", "other_allele", "standard_error")
RSID = re.compile(r"rs[0-9]+")
ALLELE = re.compile(r"[ACGT]+")
MAP_POSITION = {"GRCh37": "POS19", "GRCh38": "POS38"}
MD5_LINE = re.compile(r"(?P<md5>[0-9a-fA-F]{32})[ \t]+\*?(?P<name>.+?)[ \t]*")
# Why a data row is not read, in the order the rules apply (module docstring); a row counts under the first.
REJECT_REASONS = ("row_width", "no_position", "unparseable_number", "invalid_allele", "equal_alleles",
                  "effect_allele_not_in_pair", "beta_nonfinite", "or_nonpositive_or_nonfinite", "se_nonpositive_or_nonfinite",
                  "ci_nonpositive_or_nonfinite", "ci_not_increasing", "or_outside_ci", "p_outside_0_1",
                  "neg_log10_p_negative", "eaf_outside_0_1", "p_missing_pvalue_coloc", "p_zero_pvalue_coloc",
                  "eaf_not_inside_0_1_pvalue_coloc")
# The reasons that apply in one uncertainty mode only; every other reason applies in all three.
MODE_REJECT_REASONS: dict[str, tuple[str, ...]] = {
    "native_se": ("se_nonpositive_or_nonfinite",),
    "or_ci_derived_se": ("ci_nonpositive_or_nonfinite", "ci_not_increasing", "or_outside_ci"),
    "pvalue_coloc": ("p_missing_pvalue_coloc", "p_zero_pvalue_coloc", "eaf_not_inside_0_1_pvalue_coloc")}
CI_Z = 1.96                       # the normal quantile of a 95% confidence interval (OrCiCheck, or_ci_derived_se)
CI_RELATIVE_TOLERANCE = 0.01      # OrCiCheck: agreement within 1% of the expected limit
CI_MIN_AGREEMENT_PERCENT = 99     # OrCiCheck: the file passes when at least 99% of the checked rows agree
# or_ci_derived_se, fixed before the route is run on any file:
CI_LEVEL = 0.95                   # the level ci_lower / ci_upper are read at; the GWAS-SSF metadata schema has no key
                                  # that states another, and the 95% reading stands only if CiPCheck passes
# OR inside its interval (`or_ci_standard_error`): OR may lie outside [L, U] only by the rounding the printed
# decimals imply, half a unit in the last printed decimal of OR plus half a unit in the last printed decimal of the
# limit it crosses.
CI_P_LOG10_TOLERANCE = 0.1        # CiPCheck: a row agrees when |log10 p_implied - log10 p_published| <= 0.1, or within
                                  # the published p's printed rounding (half a unit in its last digit) where looser
CI_P_MIN_AGREEMENT_PERCENT = 95   # CiPCheck: the file passes when at least 95% of the checkable rows (0 < p < 1) agree


class SsfColumns(BaseModel):
    """Where each canonical field is read from in a GWAS-SSF header."""
    model_config = ConfigDict(frozen=True)

    effect: str            # "beta" or "odds_ratio"
    p: str | None
    p_is_neg_log10: bool
    eaf: str | None
    n: str | None
    rsid: str | None       # "rsid" or "rs_id"; None: rsIDs come from the UKB-PPP map
    ci_lower: str | None = None    # the confidence limits, read in or_ci_derived_se only
    ci_upper: str | None = None


class SsfFile(BaseModel):
    """How one GWAS-SSF data file is collected and read; the last five fields go to its record."""
    model_config = ConfigDict(frozen=True)

    name: str
    md5: str
    build: str
    position_offset: int
    rsid_rule: str
    meta_sha256: str
    note: str = ""

    def record_fields(self) -> dict:
        return {"layout": "gwas_ssf", "build": self.build, "position_offset": self.position_offset,
                "rsid_rule": self.rsid_rule, "meta_sha256": self.meta_sha256, "detail": self.note}


def ssf_data_files(accession: str, names: Iterable[str]) -> list[str]:
    """The GWAS-SSF data files of `accession` among a study directory's file names."""
    return sorted(n for n in names if (m := SSF_DATA.fullmatch(n)) is not None and m.group("acc") == accession)


class SsfMeta(BaseModel):
    """The keys stage B reads from a GWAS-SSF -meta.yaml, top-level scalars in the standard's schema.
    A plain scalar is read as its string (_MetaLoader); a list or a mapping fails validation."""
    model_config = ConfigDict(frozen=True, extra="ignore", strict=True)

    genome_assembly: str
    data_file_name: str
    coordinate_system: str | None = None
    data_file_md5sum: str | None = None


class _MetaLoader(yaml.SafeLoader):
    """PyYAML's SafeLoader (what yaml.safe_load uses), refusing a mapping that repeats a key, and
    resolving a plain scalar only to null or to the string it is written as: SafeLoader would read
    an all-digit MD5 as an integer (and one with a leading zero as octal), losing the string."""

    yaml_implicit_resolvers: ClassVar[dict] = {first: [(tag, rx) for tag, rx in rules if tag == "tag:yaml.org,2002:null"]
                               for first, rules in yaml.SafeLoader.yaml_implicit_resolvers.items()}

    def construct_mapping(self, node, deep=False):
        seen = set()
        for key_node, _value in node.value:
            key = self.construct_object(key_node, deep=deep)
            if isinstance(key, Hashable):
                if key in seen:
                    raise yaml.constructor.ConstructorError(None, None, f"the key {key!r} is given twice",
                                                            key_node.start_mark)
                seen.add(key)
        return super().construct_mapping(node, deep=deep)


def parse_meta_yaml(meta_bytes: bytes, what: str) -> SsfMeta:
    """A GWAS-SSF -meta.yaml (module docstring); SourceUnreadable with the reason when it is not
    UTF-8, does not parse as YAML, repeats a key, is not a mapping, or fails SsfMeta."""
    try:
        text = meta_bytes.decode("utf-8")
    except UnicodeDecodeError:
        raise SourceUnreadable(f"{what}: the -meta.yaml is not UTF-8") from None
    try:
        doc = yaml.load(text.removeprefix("﻿"), Loader=_MetaLoader)  # a SafeLoader subclass
    except yaml.YAMLError as err:
        raise SourceUnreadable(f"{what}: the -meta.yaml does not parse as YAML: {getattr(err, 'problem', '') or type(err).__name__}") from None
    if not isinstance(doc, dict):
        raise SourceUnreadable(f"{what}: the -meta.yaml is not a mapping at its top level")
    try:
        return SsfMeta.model_validate(doc)
    except ValidationError as err:
        fields = sorted({".".join(str(x) for x in e["loc"]) for e in err.errors()})
        raise SourceUnreadable(f"{what}: the -meta.yaml has no string value for {', '.join(fields)}") from None


def parse_md5sum(text: str) -> dict[str, str]:
    """md5sum.txt: file name -> lower-case MD5. A name listed twice with two MD5s is left out."""
    seen: dict[str, set[str]] = {}
    for line in text.splitlines():
        m = MD5_LINE.fullmatch(line.rstrip("\r"))
        if m is not None:
            seen.setdefault(m.group("name"), set()).add(m.group("md5").lower())
    return {name: next(iter(md5s)) for name, md5s in seen.items() if len(md5s) == 1}


def md5_hex(data: bytes) -> str:
    return hashlib.md5(data, usedforsecurity=False).hexdigest()


def ssf_columns(header: list[str], what: str) -> SsfColumns:
    """The column map of a GWAS-SSF header, or SourceUnreadable naming what is missing."""
    if len(set(header)) != len(header):
        raise SourceUnreadable(f"{what}: the header repeats a column name")
    missing = [c for c in REQUIRED if c not in header]
    if missing:
        raise SourceUnreadable(f"{what}: the header has no {', '.join(missing)} column")
    has_beta, has_or = "beta" in header, "odds_ratio" in header
    if has_beta and has_or:
        effect = header[4] if len(header) > 4 and header[4] in ("beta", "odds_ratio") else ""
        if not effect:
            raise SourceUnreadable(f"{what}: the header has both beta and odds_ratio and column 4 is neither")
    elif has_beta or has_or:
        effect = "beta" if has_beta else "odds_ratio"
    else:
        raise SourceUnreadable(f"{what}: the header has neither a beta nor an odds_ratio column")
    p = "p_value" if "p_value" in header else "neg_log_10_p_value" if "neg_log_10_p_value" in header else None
    rsid = "rsid" if "rsid" in header else "rs_id" if "rs_id" in header else None
    return SsfColumns(effect=effect, p=p, p_is_neg_log10=p == "neg_log_10_p_value",
                      eaf="effect_allele_frequency" if "effect_allele_frequency" in header else None,
                      n="n" if "n" in header else None, rsid=rsid,
                      ci_lower="ci_lower" if "ci_lower" in header else None,
                      ci_upper="ci_upper" if "ci_upper" in header else None)


def choose_uncertainty_mode(cols: SsfColumns | None, missing: int, present: int) -> tuple[UncertaintyMode | None, str]:
    """The uncertainty mode of a whole file (module docstring), or None and the reason none fits.
    `cols` is the GWAS-SSF column map, None for an author format (read in native_se only);
    `missing` and `present` count the rows whose standard_error is missing or present
    (RowReader.se_state)."""
    if missing + present == 0:
        return None, "no data row with a position and a valid pair of alleles to choose the uncertainty mode from"
    if missing == 0:
        return "native_se", ""
    if present > 0:
        return None, (f"standard_error is missing in {missing} and present in {present} rows; the uncertainty mode is "
                      "chosen for the whole file, never row by row")
    if cols is not None and cols.effect == "odds_ratio" and cols.ci_lower and cols.ci_upper:
        return "or_ci_derived_se", ""
    if cols is not None and cols.effect == "beta" and cols.p and cols.eaf and not (cols.ci_lower or cols.ci_upper):
        return "pvalue_coloc", ""
    return None, ("standard_error is missing in every row and the header gives no other route (an odds ratio with "
                  "ci_lower and ci_upper, or a beta with a p-value and an effect-allele frequency and no confidence limits)")


def resolve_ssf(accession: str, name: str, meta_bytes: bytes, md5sum_text: str, header: list[str]) -> SsfFile:
    """The SsfFile of a study directory's GWAS-SSF data file `name`, from its -meta.yaml bytes,
    the directory's md5sum.txt and the data file's header line (module docstring)."""
    what = f"GWAS Catalog {accession}"
    md5s = parse_md5sum(md5sum_text)
    meta_name = name + META_SUFFIX
    note = (f"md5sum.txt lists another MD5 for {meta_name} than that of the file served"
            if meta_name in md5s and md5_hex(meta_bytes) != md5s[meta_name] else "")
    meta = parse_meta_yaml(meta_bytes, what)
    build = meta.genome_assembly
    if build not in BUILDS:
        raise SourceUnreadable(f"{what}: -meta.yaml gives genome_assembly {build!r}, not GRCh37 or GRCh38")
    coordinate = meta.coordinate_system or None
    if coordinate is not None and coordinate not in COORDINATE_OFFSET:
        raise SourceUnreadable(f"{what}: -meta.yaml gives coordinate_system {coordinate!r}, not 1-based or 0-based")
    if meta.data_file_name != name:
        raise SourceUnreadable(f"{what}: -meta.yaml names another data file than {name}")
    md5 = md5s.get(name)
    if md5 is None:
        raise SourceUnreadable(f"{what}: md5sum.txt gives no single MD5 for {name}")
    declared = meta.data_file_md5sum or None
    if declared is not None and declared.lower() != md5:
        raise SourceUnreadable(f"{what}: the MD5 of {name} in -meta.yaml and in md5sum.txt disagree")
    cols = ssf_columns(header, what)
    return SsfFile(name=name, md5=md5, build=build, position_offset=COORDINATE_OFFSET.get(coordinate or "1-based", 0),
                   rsid_rule="column" if cols.rsid else "ukbppp_map", meta_sha256=hashlib.sha256(meta_bytes).hexdigest(),
                   note=note)


# ---- reading: rows of the collected file -> the canonical outcome table ------------------------------

def _num(x: str) -> float:
    """A numeric field: NaN where missing; ValueError where it is not a number."""
    x = x.strip()
    return float("nan") if x in MISSING else float(x)


def _rsid(value: str) -> str:
    value = value.strip()
    return value if RSID.fullmatch(value) else ""


def p_value(value: float, neg_log10: bool) -> float | str:
    """The p-value of a row from the number its p column holds, or the reason the row is not read.
    Missing stays missing (NaN)."""
    if math.isnan(value):
        return float("nan")
    if neg_log10:
        if value < 0:
            return "neg_log10_p_negative"
        p = 10.0 ** (-value)          # 0.0 where -log10 p is too large for a double, or infinite
    else:
        if not 0.0 <= value <= 1.0:
            return "p_outside_0_1"
        p = value
    return p


def effect_beta(value: float, is_odds_ratio: bool) -> float | str:
    """beta as given, or ln(OR); a beta that is missing or infinite, or an odds ratio that is missing,
    infinite, zero or negative, is a reason."""
    if not is_odds_ratio:
        return value if math.isfinite(value) else "beta_nonfinite"
    if not (math.isfinite(value) and value > 0):
        return "or_nonpositive_or_nonfinite"
    return math.log(value)


def ci_z_for_level(level: float) -> float:
    """The normal quantile of a two-sided interval of `level`: CI_Z (1.96) at 95%, as the registered
    formula writes it; the exact quantile at any other level."""
    if not 0.0 < level < 1.0:
        raise ValueError(f"a confidence level must lie in (0, 1), not {level}")
    return CI_Z if level == 0.95 else NormalDist().inv_cdf(0.5 + level / 2)


def or_ci_standard_error(odds: float, lower: float, upper: float, texts: tuple[str, str, str], z: float) -> float | str:
    """SE of ln(OR) = (ln U - ln L) / (2 z) from a finite positive odds ratio and its limits, or the
    reason the row is not read (module docstring). `texts` are OR, L and U as printed, for the
    rounding tolerance of OR inside [L, U]."""
    if not (math.isfinite(lower) and math.isfinite(upper) and lower > 0 and upper > 0):
        return "ci_nonpositive_or_nonfinite"
    if not lower < upper:
        return "ci_not_increasing"
    if odds < lower and lower - odds > half_unit(texts[0]) + half_unit(texts[1]):
        return "or_outside_ci"
    if odds > upper and odds - upper > half_unit(texts[0]) + half_unit(texts[2]):
        return "or_outside_ci"
    return (math.log(upper) - math.log(lower)) / (2 * z)


def frequency_problem(eaf: float, pvalue_coloc: bool) -> str:
    """"" for an effect-allele frequency the mode reads; else the reason. Missing is read outside
    pvalue_coloc, where the frequency must lie strictly inside (0, 1)."""
    if math.isnan(eaf):
        return "eaf_not_inside_0_1_pvalue_coloc" if pvalue_coloc else ""
    if not 0.0 <= eaf <= 1.0:
        return "eaf_outside_0_1"
    if pvalue_coloc and not 0.0 < eaf < 1.0:
        return "eaf_not_inside_0_1_pvalue_coloc"
    return ""


def pvalue_coloc_problem(p: float) -> str:
    """"" for a p-value coloc's p-value form reads, p in (0, 1]; else the reason."""
    if math.isnan(p):
        return "p_missing_pvalue_coloc"
    return "p_zero_pvalue_coloc" if p == 0.0 else ""


def allele_problem(effect_allele: str, other_allele: str) -> str:
    """"" for a valid pair of alleles (already stripped and upper case); else the reason."""
    if ALLELE.fullmatch(effect_allele) is None or ALLELE.fullmatch(other_allele) is None:
        return "invalid_allele"
    if effect_allele == other_allele:
        return "equal_alleles"
    return ""


def _canonical(rsid: str, chrom: str, pos: int, ea: str, oa: str, eaf: float, beta: float, se: float, p: float,
               n: float) -> dict:
    return {"rsid": rsid, "chrom": chrom, "pos": pos, "ea": ea, "oa": oa, "eaf": eaf, "beta": beta, "se": se, "p": p, "n": n}


class SsfReader(RowReader):
    """Rows of a GWAS-SSF file under its column map (`ssf_columns`) and its uncertainty `mode`
    (module docstring); SourceUnreadable where the header has no column map, or none for the mode.
    Without an rsid column every row has an empty rsID here; `attach_map_rsids` fills it from the
    UKB-PPP map. In pvalue_coloc a row read has se NaN."""

    def __init__(self, header: list[str], position_offset: int, what: str, mode: UncertaintyMode = "native_se",
                 ci_level: float = CI_LEVEL):
        self.cols = cols = ssf_columns(header, what)
        if mode == "or_ci_derived_se" and not (cols.effect == "odds_ratio" and cols.ci_lower and cols.ci_upper):
            raise SourceUnreadable(f"{what}: or_ci_derived_se needs an odds_ratio effect with ci_lower and ci_upper")
        if mode == "pvalue_coloc" and not (cols.effect == "beta" and cols.p):
            raise SourceUnreadable(f"{what}: pvalue_coloc needs a beta effect and a p-value column")
        idx = {h: i for i, h in enumerate(header)}
        self.mode, self.ci_z = mode, ci_z_for_level(ci_level)
        self.width, self.position_offset = len(header), position_offset
        self.i_chrom, self.i_pos = idx["chromosome"], idx["base_pair_location"]
        self.i_ea, self.i_oa, self.i_se = idx["effect_allele"], idx["other_allele"], idx["standard_error"]
        self.i_effect = idx[cols.effect]
        self.i_p, self.i_eaf, self.i_n, self.i_rsid, self.i_lower, self.i_upper = (
            idx[c] if c else None for c in (cols.p, cols.eaf, cols.n, cols.rsid, cols.ci_lower, cols.ci_upper))

    def _alleles(self, r: list[str]) -> tuple[str, str] | str:
        ea, oa = r[self.i_ea].strip().upper(), r[self.i_oa].strip().upper()
        return allele_problem(ea, oa) or (ea, oa)

    def se_state(self, r: list[str]) -> str:
        """"missing" or "present": the standard_error field of a row with a position and a valid pair
        of alleles (the census the uncertainty mode is chosen from); "" for any other row."""
        if isinstance(self.position(r), str) or isinstance(self._alleles(r), str):
            return ""
        return "missing" if r[self.i_se].strip() in MISSING else "present"

    def values(self, r: list[str], chrom: str, pos: int) -> dict | str:
        alleles = self._alleles(r)
        if isinstance(alleles, str):
            return alleles
        ea, oa = alleles
        mode, nan = self.mode, float("nan")
        try:
            effect = _num(r[self.i_effect])
            se = _num(r[self.i_se]) if mode == "native_se" else nan
            lower, upper = (_num(r[self.i_lower]), _num(r[self.i_upper])) if mode == "or_ci_derived_se" else (nan, nan)
            pv = _num(r[self.i_p]) if self.i_p is not None else nan
            eaf = _num(r[self.i_eaf]) if self.i_eaf is not None else nan
            n = _num(r[self.i_n]) if self.i_n is not None else nan
        except ValueError:
            return "unparseable_number"
        beta = effect_beta(effect, self.cols.effect == "odds_ratio")
        if isinstance(beta, str):
            return beta
        if mode == "native_se" and not (math.isfinite(se) and se > 0):
            return "se_nonpositive_or_nonfinite"
        if mode == "or_ci_derived_se":
            se = or_ci_standard_error(effect, lower, upper, (r[self.i_effect], r[self.i_lower], r[self.i_upper]), self.ci_z)
            if isinstance(se, str):
                return se
        p = p_value(pv, self.cols.p_is_neg_log10)
        if isinstance(p, str):
            return p
        bad = frequency_problem(eaf, False)
        if not bad and mode == "pvalue_coloc":
            bad = pvalue_coloc_problem(p) or frequency_problem(eaf, True)
        if bad:
            return bad
        rsid = _rsid(r[self.i_rsid]) if self.i_rsid is not None else ""
        return _canonical(rsid, chrom, pos, ea, oa, eaf, beta, se, p, n)


class AuthorReader(RowReader):
    """Rows of a file in a reviewed author format (stage_b/author_formats.py); SourceUnreadable for
    an entry without a map, a header other than the reviewed one, or a mode other than native_se.
    The other allele is the allele of the pair that is not the effect allele; p is an ordinary p-value."""

    def __init__(self, fmt: AuthorFormat, header: list[str], mode: UncertaintyMode = "native_se"):
        if not fmt.readable:
            raise SourceUnreadable(f"GWAS Catalog {fmt.accession}: {fmt.reason}")
        check_header(fmt, header)
        if mode != "native_se":
            raise SourceUnreadable(f"GWAS Catalog {fmt.accession}: an author format is read in native_se only, not {mode}")
        self.fmt = fmt
        idx = {h: i for i, h in enumerate(header)}
        self.idx = idx
        self.width = len(header)
        self.i_chrom, self.i_pos = idx[fmt.chrom], idx[fmt.pos]

    def _alleles(self, r: list[str]) -> tuple[str, str] | str:
        fmt, idx = self.fmt, self.idx
        a1, a2 = fmt.allele_pair
        ea, x1, x2 = (r[idx[c]].strip().upper() for c in (fmt.effect_allele, a1, a2))
        if any(ALLELE.fullmatch(a) is None for a in (ea, x1, x2)):
            return "invalid_allele"
        if x1 == x2:
            return "equal_alleles"
        if ea not in (x1, x2):
            return "effect_allele_not_in_pair"
        return ea, (x2 if ea == x1 else x1)

    def se_state(self, r: list[str]) -> str:
        """As SsfReader.se_state, on the map's standard-error column."""
        if isinstance(self.position(r), str) or isinstance(self._alleles(r), str):
            return ""
        return "missing" if r[self.idx[self.fmt.se]].strip() in MISSING else "present"

    def values(self, r: list[str], chrom: str, pos: int) -> dict | str:
        fmt, idx = self.fmt, self.idx
        alleles = self._alleles(r)
        if isinstance(alleles, str):
            return alleles
        ea, oa = alleles
        try:
            effect, se = _num(r[idx[fmt.beta or fmt.odds_ratio]]), _num(r[idx[fmt.se]])
            pv = _num(r[idx[fmt.p]]) if fmt.p else float("nan")
            eaf = _num(r[idx[fmt.eaf]]) if fmt.eaf else float("nan")
            n = sum(_num(r[idx[col]]) for col in fmt.n_sum) if fmt.n_sum else float("nan")
        except ValueError:
            return "unparseable_number"
        beta = effect_beta(effect, not fmt.beta)
        if isinstance(beta, str):
            return beta
        if not (math.isfinite(se) and se > 0):
            return "se_nonpositive_or_nonfinite"
        p = p_value(pv, False)
        if isinstance(p, str):
            return p
        bad = frequency_problem(eaf, False)
        if bad:
            return bad
        return _canonical(_rsid(r[idx[fmt.rsid]]), chrom, pos, ea, oa, eaf, beta, se, p, n)


def filter_gwas_ssf(rows: Iterable[list[str]], header: list[str], chrom: str, center: int, half_width: int,
                    position_offset: int, what: str, mode: UncertaintyMode) -> pd.DataFrame:
    """Window rows of a GWAS-SSF file as the canonical table (SsfReader in the file's mode)."""
    return SsfReader(header, position_offset, what, mode).window(rows, chrom, center, half_width)


def filter_author(fmt: AuthorFormat, rows: Iterable[list[str]], header: list[str], chrom: str, center: int,
                  half_width: int, what: str, mode: UncertaintyMode) -> pd.DataFrame:
    """Window rows of a file in a reviewed author format (AuthorReader)."""
    return AuthorReader(fmt, header, mode).window(rows, chrom, center, half_width)


def half_unit(text: str) -> float:
    """Half a unit in the last digit of a number as printed: 0.005 for "1.23", 5e-07 for "1.2e-05"."""
    number = Decimal(text.strip())
    if not number.is_finite():
        raise ValueError(f"{text!r} is not a finite number")
    return 0.5 * 10.0 ** int(number.as_tuple().exponent)


class OrCiCheck:
    """The integrity rule of an author file whose effect is an odds ratio printed with its 95%
    confidence limits and standard error (module docstring). `add` takes every data row; `result`
    gives counts and pass/fail only."""

    def __init__(self, fmt: AuthorFormat, header: list[str]):
        missing = [c or "(none mapped)" for c in (fmt.ci_lower, fmt.ci_upper) if not c or c not in header]
        if missing:
            raise SourceUnreadable(f"GWAS Catalog {fmt.accession}: the header has no confidence-limit column "
                                   f"{', '.join(missing)}, so the OR/CI/SE integrity rule cannot be applied")
        idx = {h: i for i, h in enumerate(header)}
        self.width = len(header)
        self.cols = [idx[c] for c in (fmt.odds_ratio, fmt.ci_lower, fmt.ci_upper, fmt.se)]
        self.checked = self.agree = self.not_checkable = 0

    @staticmethod
    def _within(printed: float, printed_text: str, expected: float, rounding_of_expected: float) -> bool:
        tolerance = max(CI_RELATIVE_TOLERANCE * expected, half_unit(printed_text) + rounding_of_expected)
        return abs(printed - expected) <= tolerance

    def add(self, r: list[str]) -> None:
        if len(r) != self.width:
            self.not_checkable += 1
            return
        texts = [r[i].strip() for i in self.cols]
        try:
            odds, lower, upper, se = (float(t) for t in texts)
        except ValueError:
            self.not_checkable += 1
            return
        if not all(math.isfinite(v) for v in (odds, lower, upper, se)) or odds <= 0 or se <= 0:
            self.not_checkable += 1
            return
        self.checked += 1
        relative_rounding = half_unit(texts[0]) / odds + CI_Z * half_unit(texts[3])
        lo, hi = math.exp(math.log(odds) - CI_Z * se), math.exp(math.log(odds) + CI_Z * se)
        if self._within(lower, texts[1], lo, lo * relative_rounding) and self._within(upper, texts[2], hi, hi * relative_rounding):
            self.agree += 1

    def result(self) -> dict:
        passed = self.checked > 0 and 100 * self.agree >= CI_MIN_AGREEMENT_PERCENT * self.checked
        return {"rule": "or_ci_se", "z": CI_Z, "relative_tolerance": CI_RELATIVE_TOLERANCE,
                "min_agreement_percent": CI_MIN_AGREEMENT_PERCENT, "rows_checked": self.checked,
                "rows_agree": self.agree, "rows_disagree": self.checked - self.agree,
                "rows_not_checkable": self.not_checkable, "passed": passed}


SQRT2, LN10 = math.sqrt(2.0), math.log(10.0)


def log10_two_sided_p(z: float) -> float:
    """log10 of the two-sided normal p-value of `z`, erfc(|z|/sqrt 2); past the point where erfc
    leaves the doubles (|z| > 37), its asymptotic series ln erfc(x) = -x^2 - ln(x sqrt(pi))
    + ln(1 - 1/(2x^2) + 3/(4x^4) - 15/(8x^6)), whose error there is below 1e-9."""
    x = abs(z) / SQRT2
    tail = math.erfc(x)
    if tail > 1e-300:
        return math.log10(tail)
    x2 = x * x
    series = 1 - 1 / (2 * x2) + 3 / (4 * x2 ** 2) - 15 / (8 * x2 ** 3)
    return (-x2 - math.log(x * math.sqrt(math.pi)) + math.log(series)) / LN10


def p_rounding_log10(text: str, p: float, neg_log10: bool) -> float:
    """How far, on the log10 scale, a published p can lie from the value its printed digits give:
    half a unit in the last printed digit of -log10 p, or of p mapped to log10 (the larger side)."""
    h = half_unit(text)
    if neg_log10:
        return h
    return max(math.log10(p + h) - math.log10(p), math.log10(p) - math.log10(p - h)) if p > h else math.inf


class CiPCheck:
    """The frozen file-level CI-versus-p check of an or_ci_derived_se file (module docstring). `add`
    takes each row the reader read (with the derived SE) and the row as printed; `result` gives
    counts and pass/fail only."""

    def __init__(self, reader: SsfReader):
        if reader.mode != "or_ci_derived_se":
            raise ValueError(f"the CI-versus-p check applies to or_ci_derived_se, not {reader.mode}")
        self.i_p, self.neg_log10, self.z = reader.i_p, reader.cols.p_is_neg_log10, reader.ci_z
        self.checked = self.agree = self.not_checkable = 0

    def add(self, got: Mapping, r: list[str]) -> None:
        p = got["p"]
        if not 0.0 < p < 1.0:                       # p missing (NaN), 0 or 1: nothing to compare
            self.not_checkable += 1
            return
        self.checked += 1
        gap = abs(log10_two_sided_p(got["beta"] / got["se"]) - math.log10(p))
        if gap <= CI_P_LOG10_TOLERANCE or gap <= p_rounding_log10(r[self.i_p], p, self.neg_log10):
            self.agree += 1

    def result(self) -> dict:
        passed = self.checked > 0 and 100 * self.agree >= CI_P_MIN_AGREEMENT_PERCENT * self.checked
        return {"rule": "or_ci_vs_p", "ci_level": CI_LEVEL, "z": self.z, "log10_tolerance": CI_P_LOG10_TOLERANCE,
                "min_agreement_percent": CI_P_MIN_AGREEMENT_PERCENT, "rows_checked": self.checked,
                "rows_agree": self.agree, "rows_disagree": self.checked - self.agree,
                "rows_not_checkable": self.not_checkable, "passed": passed}


def map_rsids(rows: Iterable[list[str]], header: list[str], build: str, positions: set[int]) -> dict[tuple[int, frozenset], set[str]]:
    """UKB-PPP rsID map rows at `positions` on `build`: (position, {REF, ALT}) -> the rsIDs given."""
    need = ["rsid", "REF", "ALT", MAP_POSITION[build]]
    missing = [c for c in need if c not in header]
    if missing:
        raise InputContractError(f"UKB-PPP rsID map: header lacks {missing}; header was {header}")
    i_rs, i_ref, i_alt, i_pos = (header.index(c) for c in need)
    out: dict[tuple[int, frozenset], set[str]] = {}
    for r in rows:
        p = r[i_pos]
        if p in MISSING or not RSID.fullmatch(r[i_rs]):
            continue
        pos = int(float(p))
        if pos in positions:
            out.setdefault((pos, frozenset((r[i_ref].upper(), r[i_alt].upper()))), set()).add(r[i_rs])
    return out


def attach_map_rsids(df: pd.DataFrame, mapping: Mapping[tuple[int, frozenset], set[str]]) -> pd.DataFrame:
    """Each row's rsID from the map: the one rsID at its position with its two alleles; else ""."""
    def one(pos: int, ea: str, oa: str) -> str:
        found = mapping.get((int(pos), frozenset((ea, oa))), set())
        return next(iter(found)) if len(found) == 1 else ""
    out = df.copy()
    out["rsid"] = [one(p, a, b) for p, a, b in zip(df["pos"], df["ea"], df["oa"])]
    return out
