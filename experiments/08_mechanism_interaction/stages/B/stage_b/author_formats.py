"""Reviewed column maps for GWAS Catalog outcome files deposited in an author's own format, with
neither a harmonised copy nor a GWAS-SSF file (stage_b/outcome_files.py reads GWAS-SSF files).

Each entry is fixed by hand from the files the GWAS Catalog deposits beside the data file (the
read-me, README or sumstats.info) and from the data file's header line, read through an HTTP Range
request of its first 64 KB; no row of association data was read. Each entry holds the sha256 of
every document it was derived from and of the header line (its text without the line
terminator, UTF-8), so a deposit that changes is refused at collect (`check_header`), not read
under a map written for other columns.

A map is written only where the documents state every item the canonical outcome table needs: the
rsID column, the effect allele and the other allele, the effect (beta, or an odds ratio read as
log(OR)) and its standard error, and the genome build. An entry whose documents leave one of these
unstated is `readable=False` with the reason, and the study is recorded `unreadable`.

GCST008226 (Laskar et al. 2019, PMID 31231134, renal cell carcinoma, men). Source documents:
read-me_Laskar_31231134.txt (sha256 below). It lists 15 columns, in the order of the header line
and one to one with it; three of its names are spelled differently from the header (`Odda_ratio`
for `Odds_ratio`, `Number_Controls` for `Number_controls`, `effect_allele_frequency_contols` for
`effect_allele_frequency_controls`). What it states:
    Variant_ID          "Variant identification (rsID) for each SNP assigned by dbSNP"
    Chromosome          "Chromosome Number"
    base_pair_location  "Base Position of seach SNP mapped to hg19"                -> GRCh37
    Allele_1, Allele_2  "First allele or Reference allele given by dbSNP", "2nd allele or alternative allele"
    Effect_allele       "Effect alele (allele for which the effect (OR) is calculated)"
    Odds_ratio          "Odds ratio (OR, Effect size)"                              -> beta = ln(OR)
    Standard_error      "Standard Error"
    P_Value             "Significance"
    Number_controls, Number_cases                                                   -> n = their sum
    effect_allele_frequency_controls, effect_allele_frequency_cases
The other allele is the one of Allele_1 and Allele_2 that is not Effect_allele; a row where
Effect_allele is neither, or both alleles are equal, has no other allele and is dropped. The
read-me gives the effect-allele frequency in cases and in controls only, and no frequency over
the whole sample, so `eaf` is missing (a palindromic variant is then dropped at harmonization, and
coloc takes the MAF from the pQTL side, as for any outcome row without a frequency). The read-me
does not state the scale of Standard_error. It is read as the standard error of ln(OR), the reading
stage B gives every `standard_error` that accompanies an `odds_ratio` (harmonised and GWAS-SSF
files alike), and the reading is checked before the file is analyzed: the header's ci_lower and
ci_upper columns must equal exp(ln(OR) -/+ 1.96 SE) on at least 99% of the rows where all four
are finite (outcome_files.OrCiCheck, run by the pre-analysis validation, stage_b/validate.py),
else the study is recorded `unreadable`. The read-me does not state a coordinate system; positions are read as 1-based, the
reading stage B gives a GWAS-SSF file whose metadata states none. GWAS Catalog deposits no
md5sum.txt for this study, so no MD5 is checked; the collect record holds the file's size and
sha256.

GCST010514 (SPARK European ancestry with iPSYCH-PGC, ASD). The deposited `README copy.md` and
`sumstats.info copy` are saved Bitbucket web pages (the application shell of
bitbucket.org/steinlabunc/spark_asd_sumstats, src/master/README.md and sumstats.info): they hold
page code and no column definitions. The header line is a METAL layout (Chromosome, Position,
MarkerName, Allele1, Allele2, Effect, StdErr, P-value, Direction, TotalSampleSize). The deposited
documents state neither which allele the effect refers to, nor whether `Effect` is a log odds
ratio, nor whether MarkerName is an rsID, nor the genome build, so no map is written.
"""
import hashlib
from typing import Literal

from pydantic import BaseModel, ConfigDict

from stage_b.schemas import SourceUnreadable


class AuthorFormat(BaseModel):
    """One reviewed author format. Column fields name header columns; empty where not mapped."""
    model_config = ConfigDict(frozen=True)

    accession: str
    data_file: str
    documents: dict[str, str]          # deposited file name -> sha256 of the bytes reviewed
    header: tuple[str, ...]
    header_sha256: str
    readable: bool
    reason: str = ""
    build: str = ""
    rsid: str = ""
    chrom: str = ""
    pos: str = ""
    effect_allele: str = ""
    allele_pair: tuple[str, str] = ("", "")   # the two alleles; the other allele is the one that is not the effect allele
    odds_ratio: str = ""
    beta: str = ""
    se: str = ""
    p: str = ""
    eaf: str = ""
    n_sum: tuple[str, ...] = ()
    ci_lower: str = ""                 # the confidence limits of an odds ratio, read by the integrity rule only
    ci_upper: str = ""
    integrity: Literal["", "or_ci_se"] = ""   # outcome_files.OrCiCheck, applied in the pre-analysis validation


def header_sha256(header: list[str] | tuple[str, ...]) -> str:
    """sha256 of a header line as text: its columns joined by tabs, no line terminator, UTF-8."""
    return hashlib.sha256("\t".join(header).encode()).hexdigest()


LASKAR_HEADER = ("Variant_ID", "Chromosome", "base_pair_location", "Allele_1", "Allele_2", "Effect_allele", "Odds_ratio",
                 "ci_lower", "ci_upper", "Standard_error", "P_Value", "Number_controls", "Number_cases",
                 "effect_allele_frequency_controls", "effect_allele_frequency_cases")
SPARK_HEADER = ("Chromosome", "Position", "MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value", "Direction",
                "TotalSampleSize")

AUTHOR_FORMATS: dict[str, AuthorFormat] = {
    "GCST008226": AuthorFormat(
        accession="GCST008226", data_file="Laskar_31231134_Males.gz",
        documents={"read-me_Laskar_31231134.txt": "ac1f2a8b61c6447e3d052d40336de345771cd1217a189c953bef0be657ad4291"},
        header=LASKAR_HEADER, header_sha256="c666964a0bd48110c5473d401e087cca784e3d9f08d88574d46174ba8c4767fe",
        readable=True, build="GRCh37", rsid="Variant_ID", chrom="Chromosome", pos="base_pair_location",
        effect_allele="Effect_allele", allele_pair=("Allele_1", "Allele_2"), odds_ratio="Odds_ratio", se="Standard_error",
        p="P_Value", n_sum=("Number_controls", "Number_cases"), ci_lower="ci_lower", ci_upper="ci_upper",
        integrity="or_ci_se"),
    "GCST010514": AuthorFormat(
        accession="GCST010514", data_file="ASD_SPARK_EUR_iPSYCH_PGC.tsv",
        documents={"README copy.md": "24f8eb1c9d19dc87074eefc53c125aac3b4b6c25e2e3b03e867975f2832950da",
                   "sumstats.info copy": "8e6e8ac7986ebcb43eca8a9b3484f91671771e897a799f8bd95bbb8661729c9c",
                   "md5sum.txt": "d2ab54f162e394b817bc2bf26a8fd4a2937e76b28ba2bc267786ed940423811c"},
        header=SPARK_HEADER, header_sha256="cce9e53186302a4091e5a4aeae52c9875ea5f8aedc45079119c35b1c595db442",
        readable=False,
        reason=("author format; the deposited README copy.md and sumstats.info copy are saved Bitbucket web pages with no "
                "column definitions, so the effect allele, the scale of Effect, whether MarkerName is an rsID and the genome "
                "build are unstated")),
}


def check_header(fmt: AuthorFormat, header: list[str]) -> None:
    """Raises SourceUnreadable unless `header` is the header line the map was reviewed against."""
    if tuple(header) != fmt.header or header_sha256(header) != fmt.header_sha256:
        raise SourceUnreadable(f"GWAS Catalog {fmt.accession}: the header of {fmt.data_file} is not the one its reviewed "
                               f"author-format map was written for")

