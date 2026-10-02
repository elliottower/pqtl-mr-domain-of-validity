"""Instrument lists (PREREG §Study design, "Instrument lists").

A gene has a cis-pQTL in a source when it is on the source's published cis list and the assay
maps to exactly one gene and one protein accession (strict, S1). The source-native inclusive
lists (every gene a source names, multi-gene assays split, the V3.4 EpiGraphDB list included)
form S20. Ported from `count_v4.py` (`ukbppp_lists`, `decode_lists`, `interval_lists`), with
the assay IDs kept so stage B knows which assay each gene's instrument uses.

The plan defines the lists "as implemented in `feasibility/v4_round3/count_v4.py`", so the port
keeps that script's rules as they are. For deCODE that means: strictness is judged on the SeqId's
ST01 row (one Gene, one UniProt, at most one Ensembl ID, so a blank Ensembl ID passes), and a
strict SeqId then puts on the list its ST01 gene and every gene symbol its ST02 cis rows name.
`decode_strict_diagnostic` counts where the two tables disagree; it changes no list.

Where count_v4.py does something that looks like a quirk, the port does the same and the line of
count_v4.py is cited beside it: `UniProt2` is not split on '_' (line 205); an empty ST02 gene cell
becomes the symbol "None" (line 248); an empty EpiGraphDB token is kept as the symbol "" (line
293); a falsy ST4 UniProt cell reads as no accession (line 280). No gene has the symbol "None" or
"", so these change no membership; they are kept so the lists are count_v4's, entry for entry.
"""
import re
from collections.abc import Iterable

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

SOURCE_ORDER = ("ukbppp", "decode", "interval")  # order of discovery N
PLATFORM = {"ukbppp": "Olink", "decode": "SomaScan", "interval": "SomaScan"}
UNIPROT_ACC = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")
EPIGRAPHDB_NO_ASSAY = "epigraphdb:no_assay"


class SourceList(BaseModel):
    """Genes (by HGNC symbol or Ensembl ID) with a cis-pQTL in one source, and their assays."""

    model_config = ConfigDict(frozen=True)

    by_symbol: dict[str, frozenset[str]] = Field(default_factory=dict)
    by_ensembl: dict[str, frozenset[str]] = Field(default_factory=dict)

    def has(self, symbol: str, ensembl: str) -> bool:
        return symbol in self.by_symbol or ensembl in self.by_ensembl

    def assays(self, symbol: str, ensembl: str) -> frozenset[str]:
        return self.by_symbol.get(symbol, frozenset()) | self.by_ensembl.get(ensembl, frozenset())


class InstrumentLists(BaseModel):
    model_config = ConfigDict(frozen=True)

    ukbppp: SourceList
    decode: SourceList
    interval: SourceList

    def source(self, name: str) -> SourceList:
        return getattr(self, name)

    def sources_for(self, symbol: str, ensembl: str) -> list[str]:
        """Sources with a cis-pQTL for the gene, in the registered order."""
        return [s for s in SOURCE_ORDER if self.source(s).has(symbol, ensembl)]


def _collect(pairs: Iterable[tuple[str, str]]) -> dict[str, frozenset[str]]:
    out: dict[str, set[str]] = {}
    for key, assay in pairs:
        if isinstance(key, str):   # "" is kept: count_v4.py line 293 keeps an empty EpiGraphDB token
            out.setdefault(key, set()).add(assay)
    return {k: frozenset(v) for k, v in out.items()}


def split_ids(x, sep: str = r"[ ,;|]+") -> list[str]:
    return [t for t in re.split(sep, str(x if x is not None else "")) if t and t not in ("None", "nan", "NA")]


def _accessions(values: Iterable) -> set[str]:
    """UniProt accessions of the `UniProt` column: '_'-joined values split, isoform suffix stripped,
    non-accession labels dropped (count_v4.py lines 202-203)."""
    return {a.split("-")[0] for v in values if isinstance(v, str) for a in v.split("_")
            if UNIPROT_ACC.match(a.split("-")[0])}


def _accessions_uniprot2(values: Iterable) -> set[str]:
    """UniProt accessions of the `UniProt2` column, read where `UniProt` holds a label. count_v4.py
    line 205 does not split this column on '_': a '_'-joined value is no accession and gives none."""
    return {a.split("-")[0] for a in values if isinstance(a, str) and UNIPROT_ACC.match(a.split("-")[0])}


def ukbppp_lists(cis_entries: list[str], olink: pd.DataFrame) -> tuple[SourceList, SourceList]:
    """(strict, inclusive). `cis_entries`: UKB-PPP cis-pQTL list (Olink assay names; multi-gene
    assays joined by '_'). `olink`: Olink map columns Assay, OlinkID, UniProt, UniProt2,
    HGNC.symbol, ensembl_id. Strict keeps a list entry that is an Olink assay with one HGNC
    symbol and one UniProt accession (UniProt2 used where UniProt holds a label)."""
    entries = set(cis_entries)
    ol = olink[olink["Assay"].isin(entries)]
    strict_sym, strict_ens = [], []
    for assay, g in ol.groupby("Assay", sort=True):
        syms = set(g["HGNC.symbol"].dropna())
        accs = _accessions(g["UniProt"].dropna())
        if not accs:
            accs = _accessions_uniprot2(g["UniProt2"].dropna())
        if "_" in assay or len(syms) != 1 or len(accs) != 1:
            continue
        oids = set(g["OlinkID"].dropna()) or {assay}
        for oid in oids:
            strict_sym += [(assay, oid)] + [(s, oid) for s in syms]
            strict_ens += [(e, oid) for e in g["ensembl_id"].dropna()]
    strict = SourceList(by_symbol=_collect(strict_sym), by_ensembl=_collect(strict_ens))

    incl_sym, incl_ens = [], []
    for entry in entries:
        rows = ol[ol["Assay"] == entry]
        oids = set(rows["OlinkID"].dropna()) or {entry}
        for oid in oids:
            incl_sym += [(s, oid) for s in entry.split("_")]
    for sym, ens, oid, assay in zip(ol["HGNC.symbol"], ol["ensembl_id"], ol["OlinkID"], ol["Assay"]):
        oid = oid if isinstance(oid, str) else assay
        if isinstance(sym, str):
            incl_sym.append((sym, oid))
        if isinstance(ens, str):
            incl_ens.append((ens, oid))
    inclusive = SourceList(by_symbol=_collect(incl_sym), by_ensembl=_collect(incl_ens))
    return strict, inclusive


NO_ST01_ROW: dict[str, set[str]] = {"genes": set(), "unis": set(), "ens": set()}


class DecodeStrictDiagnostic(BaseModel):
    """Outcome-blind counts on the strict deCODE list (reported in funnel.csv; no list reads them).
    A strict SeqId *disagrees* when the gene symbols its ST02 cis rows name are not exactly its
    ST01 Gene. A strict symbol is on the list *only through a disagreement* when no strict SeqId
    carries it as its ST01 Gene, so it is there because an ST02 row named it."""

    model_config = ConfigDict(frozen=True)

    strict_seqids: int
    strict_seqids_st02_differs_from_st01: int
    strict_symbols: int
    strict_symbols_only_through_disagreement: int


def decode_st01_ids(st01: pd.DataFrame) -> dict[str, dict[str, set[str]]]:
    """SeqId -> its ST01 gene symbols, UniProt accessions and Ensembl gene IDs, for aptamers
    included in the analysis."""
    seq: dict[str, dict[str, set[str]]] = {}
    for s, gene, uni, ens, inc in zip(st01["SeqId"], st01["Gene"], st01["UniProt"], st01["Ensembl"],
                                      st01["Included"]):
        if s is None or str(inc).strip() != "Yes":
            continue
        seq[str(s)] = {"genes": set(split_ids(gene)), "unis": set(split_ids(uni)),
                       "ens": {e for e in split_ids(ens) if e.startswith("ENSG")}}
    return seq


def st02_symbols(gene) -> set[str]:
    """The gene symbols of one ST02 `gene (prot.)` cell (multi-gene cells are joined by '.').
    count_v4.py line 248 drops only "NA", so an empty cell (None) gives the symbol "None"."""
    return {x for x in re.split(r"[ ,.;|]+", str(gene)) if x and x != "NA"}


def decode_is_strict(info: dict[str, set[str]]) -> bool:
    """count_v4.py: exactly one Gene, one UniProt and at most one Ensembl ID in ST01."""
    return len(info["genes"]) == 1 and len(info["unis"]) == 1 and len(info["ens"]) <= 1


def decode_lists(st01: pd.DataFrame, st02: pd.DataFrame) -> tuple[SourceList, SourceList]:
    """(strict, inclusive) from Ferkingstad 2021. `st01`: SeqId, Gene, UniProt, Ensembl,
    Included (aptamers; only 'Yes' rows count). `st02`: gene, SeqId, cis_trans. Strict keeps
    cis rows whose SeqId has, in ST01, exactly one Gene, one UniProt and at most one Ensembl ID."""
    seq = decode_st01_ids(st01)
    strict_sym, strict_ens, incl_sym, incl_ens = [], [], [], []
    for gene, s, ct in zip(st02["gene"], st02["SeqId"], st02["cis_trans"]):
        if ct != "cis":
            continue
        s = str(s)
        info = seq.get(s, NO_ST01_ROW)
        sym2 = st02_symbols(gene)
        incl_sym += [(g, s) for g in sym2 | info["genes"]]
        incl_ens += [(e, s) for e in info["ens"]]
        if decode_is_strict(info):
            strict_sym += [(g, s) for g in sym2 | info["genes"]]
            strict_ens += [(e, s) for e in info["ens"]]
    return (SourceList(by_symbol=_collect(strict_sym), by_ensembl=_collect(strict_ens)),
            SourceList(by_symbol=_collect(incl_sym), by_ensembl=_collect(incl_ens)))


def decode_strict_diagnostic(st01: pd.DataFrame, st02: pd.DataFrame) -> DecodeStrictDiagnostic:
    """Where ST01 and ST02 name different genes for a strict SeqId (class docstring). Written out
    beside `decode_lists`, not derived from its result, so the two can be compared."""
    seq = decode_st01_ids(st01)
    named: dict[str, set[str]] = {}
    for gene, s, ct in zip(st02["gene"], st02["SeqId"], st02["cis_trans"]):
        if ct == "cis" and decode_is_strict(seq.get(str(s), NO_ST01_ROW)):
            named.setdefault(str(s), set()).update(st02_symbols(gene))
    st01_genes = set().union(*(seq[s]["genes"] for s in named)) if named else set()
    st02_genes = set().union(*named.values()) if named else set()
    return DecodeStrictDiagnostic(
        strict_seqids=len(named),
        strict_seqids_st02_differs_from_st01=sum(1 for s, syms in named.items() if syms != seq[s]["genes"]),
        strict_symbols=len(st01_genes | st02_genes),
        strict_symbols_only_through_disagreement=len(st02_genes - st01_genes))


def interval_lists(st4: pd.DataFrame, uni2sym: dict[str, set[str]],
                   epigraphdb_genes: set[str]) -> tuple[SourceList, SourceList]:
    """(strict, inclusive) from Sun 2018 ST4. `st4`: somamer_id, UniProt, cis_trans. Strict keeps
    cis rows on a single-UniProt SOMAmer whose accession maps to exactly one HGNC symbol.
    Inclusive takes every symbol of every accession on cis rows plus the EpiGraphDB list, which
    carries no assay."""
    strict_sym, incl_sym = [], []
    for somamer, uni, ct in zip(st4["somamer_id"], st4["UniProt"], st4["cis_trans"]):
        if ct is None or str(ct).strip() != "cis":
            continue
        accs = [a for a in re.split(r"[ ,;|]+", str(uni or "")) if a]   # count_v4.py line 280
        syms_all = set().union(*[uni2sym.get(a, set()) for a in accs]) if accs else set()
        incl_sym += [(g, str(somamer)) for g in syms_all]
        if len(accs) == 1 and len(uni2sym.get(accs[0], set())) == 1:
            strict_sym += [(g, str(somamer)) for g in uni2sym[accs[0]]]
    incl_sym += [(g, EPIGRAPHDB_NO_ASSAY) for g in epigraphdb_genes]
    return SourceList(by_symbol=_collect(strict_sym)), SourceList(by_symbol=_collect(incl_sym))
