# The Domain of Validity of *Cis*-pQTL Mendelian Randomization for Drug-Target Prioritization

A pre-registered, outcome-blind evaluation of a frozen, zero-parameter *cis*-pQTL MR classifier against 195 Phase III drug-target pairs across 24 diseases.

**Finding:** *Cis*-pQTL MR has a definable domain of validity set by drug mechanism. It is uninformative for activity-blocking targets (BA = 0.511) and suggestively informative for abundance-modulating targets (BA = 0.590, AUC = 0.595, permutation *p* = 0.038). The pooled BA of 0.559 is a mixture of these two strata.

**Translational rule:** Trust a significant *cis*-pQTL MR result for an abundance-modulating target; treat MR silence on an activity-blocking target as uninformative, not disqualifying.

## Repository structure

```
paper/              Manuscript (LaTeX source + PDF)
protocol/           Pre-registration documents (V1-V4), SHA-256 hashes, deviation log
  v1/               Initial 3-class taxonomy
  v2/               6-class taxonomy, binary endpoint
  v3/               cis-pQTL classifier, 24 diseases, n=195 evaluable
  v4/               LoF burden second instrument (Genebass)
code/               Classification and analysis scripts
data/               Per-pair classification table (138 analyzed pairs) and frozen candidate set
results/            V4 LoF burden results
```

## Pre-registration

The protocol evolved across six publicly versioned iterations (V1-V3.4), each SHA-256 hashed before execution. All amendments were made before outcome adjudication. The V4 extension (LoF burden instrument) was separately pre-registered. See `protocol/` for all specifications and hashes.

## Data sources

This study uses publicly available summary-level data from:

- **EpiGraphDB** *cis*-pQTL MR catalog (Zheng et al., 2020)
- **UK Biobank Pharma Proteomics Project** (UKB-PPP) *cis*-pQTLs
- **Open Targets Platform** (Ochoa et al., 2021) for Phase III drug-target pairs
- **Genebass** (Karczewski et al., 2022) for pLoF burden tests
- Outcome GWAS from FinnGen R10, CARDIoGRAMplusC4D, IMSGC, PGC, IIBDGC, ILCCO/TRICL, and disease-specific consortia

No individual-level data were used. All data sources are publicly accessible.

## License

MIT License. See [LICENSE](LICENSE).
