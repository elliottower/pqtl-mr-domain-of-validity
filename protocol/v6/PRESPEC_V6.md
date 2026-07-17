# Pre-Specification V6: Expanded Mechanism Taxonomy, Open Targets Blind Spot, and ACE Case Study

**Author:** Elliot Tower
**Date:** 2026-07-17
**Status:** FROZEN (pending SHA-256 hash)
**Depends on:** V3.4 (Zenodo Version 2), V4 (LoF Extension), V5 (pQTL Instrument Expansion)

---

## Motivation

The V3.4 analysis established a mechanism-dependent domain of validity for cis-pQTL MR: BA=0.511 for activity-blocking targets (uninformative) and BA=0.590 for abundance-modulating targets (trends positive). The binary classification (abundance_modulating vs activity_blocking) was sufficient to demonstrate the boundary but is too coarse to characterize the *gradient* of MR informativeness predicted by the underlying theory.

The paper's causal framework predicts that MR informativeness depends on how closely the drug's mechanism aligns with the pQTL instrument's causal axis (circulating protein abundance). This alignment is not binary---it varies continuously with how directly the drug acts on the same circulating molecule the pQTL instruments. A finer mechanism taxonomy should reveal this gradient: drugs that directly neutralize the same soluble protein the pQTL measures should produce the strongest MR signal, while drugs that block enzymatic activity or receptor signaling should produce none.

Three new analyses test predictions derived from this gradient theory. All predictions are stated before any results are computed.

---

## Analysis 5: Expanded Mechanism Taxonomy

### 5.1 Rationale

The binary abundance/activity split collapses heterogeneous drug mechanisms into two bins. Within the abundance-modulating stratum, a biologic that neutralizes the exact circulating protein the pQTL instruments (such as anti-TNF) has near-perfect mechanistic alignment, while a small molecule MMP inhibitor that blocks enzymatic activity of a secreted protease has much weaker alignment despite being classified as abundance-modulating. Similarly, within activity-blocking, an enzyme inhibitor and a receptor tyrosine kinase antagonist share the same structural mismatch with the pQTL instrument, but the pharmacological mechanisms differ.

The paper's theory---that MR informativeness tracks the degree to which drug and instrument share a causal axis---predicts an ordered gradient across a finer taxonomy.

### 5.2 Five-Way Taxonomy

Each gene in classification_v34.csv is assigned to exactly one of five subcategories based on the drug's pharmacological mechanism and its relationship to the circulating protein the pQTL instruments.

#### Within abundance_modulating:

**A1. Soluble ligand neutralization.** The drug directly binds, neutralizes, or IS the same circulating soluble protein the pQTL instruments. This includes monoclonal antibodies against soluble cytokines, ligand traps, and recombinant protein replacement therapies. The drug-pQTL alignment is maximal: both act on the identical molecular species in circulation.

**A2. Surface target depletion.** The drug targets a cell-surface protein (checkpoint receptor, adhesion molecule, integrin). The pQTL instruments the shed or soluble form of that protein circulating in blood. Alignment is moderate: the soluble form measured by the pQTL correlates with surface expression, but shedding kinetics, surface-to-soluble ratios, and tissue-specific expression introduce noise.

**A3. Indirect abundance modulation.** The drug modulates protein levels or activity indirectly. This includes antisense oligonucleotides (ASOs) that reduce mRNA rather than directly binding protein, small molecule inhibitors of secreted enzymes (where level and activity are correlated but not identical), and drugs that target a different molecule in the same pathway. The pQTL instruments protein level, but the drug acts on a related-but-distinct quantity.

#### Within activity_blocking:

**B1. Enzyme inhibitor.** The drug blocks an enzymatic active site. The protein is present at normal or elevated circulating levels, but its catalytic function is abolished. The pQTL instruments protein level, which is causally unrelated to enzymatic activity under drug treatment. Examples: ACE inhibitors, cholinesterase inhibitors, PDE5 inhibitors, Factor Xa inhibitors.

**B2. Receptor/kinase antagonist.** The drug blocks receptor binding or kinase phosphorylation without changing the target protein's circulating level. The pQTL instruments protein level or a shed soluble form, while the drug blocks signaling. Examples: EGFR tyrosine kinase inhibitors, FGFR inhibitors, multi-kinase inhibitors.

### 5.3 Collapsed Three-Way Taxonomy

For analyses requiring larger strata, the five subcategories collapse to three:

- **Direct neutralization** = A1 (soluble ligand) + A2 (surface target). Drug directly engages the same protein the pQTL instruments.
- **Indirect abundance** = A3 (indirect). Drug-pQTL alignment is real but attenuated.
- **Function blocking** = B1 (enzyme) + B2 (receptor/kinase). Drug acts on protein function; pQTL acts on protein level. Structural mismatch.

### 5.4 Frozen Gene-to-Subcategory Classification

The complete classification is frozen below. Every gene from classification_v34.csv is assigned exactly once. The mixed-class genes (APP, IGF1R) are excluded from subcategory analysis.

#### A1. Soluble Ligand Neutralization (14 unique genes, 22 pairs)

| Gene | Representative Drug | Rationale |
|------|-------------------|-----------|
| ANGPT1 | Trebananib (peptibody trap) | Neutralizes circulating angiopoietin-1 |
| ANGPT2 | Trebananib (peptibody trap) | Neutralizes circulating angiopoietin-2 |
| C5 | Eculizumab/ravulizumab | Neutralizes circulating complement C5 |
| CSF2 | Mavrilimumab/otilimab | Neutralizes circulating GM-CSF |
| IL12B | Ustekinumab/guselkumab | Neutralizes soluble IL-12/IL-23 p40 subunit |
| IL6R | Tocilizumab/sarilumab | Targets soluble IL-6R (paper Section 2 rationale) |
| KLKB1 | Lanadelumab | mAb neutralizes circulating plasma kallikrein |
| PLAT | Alteplase/tenecteplase | Drug IS the circulating protein (tPA replacement) |
| PLAU | Urokinase | Drug IS the circulating protein (urokinase thrombolytic) |
| PLG | Plasminogen concentrate | Drug IS the circulating protein (plasminogen replacement) |
| SOST | Romosozumab | Neutralizes circulating sclerostin |
| TNF | Infliximab/adalimumab | Neutralizes circulating soluble TNF |
| TNFSF11 | Denosumab | Neutralizes circulating RANKL |
| VEGFA | Bevacizumab/aflibercept | Neutralizes circulating VEGF-A |

Gene-disease pairs in this subcategory:
ANGPT1 x Ovarian cancer; ANGPT2 x Ovarian cancer; C5 x ALS, MI; CSF2 x RA; IL12B x Crohn's, RA, SLE, UC; IL6R x JIA, RA; KLKB1 x Alzheimer's; PLAT x MI; PLAU x MI; PLG x MI; SOST x Anorexia; TNF x IBD; TNFSF11 x Anorexia; VEGFA x Glioma, Lung cancer, Ovarian cancer, Pancreatic cancer.

#### A2. Surface Target Depletion (13 unique genes, 20 pairs)

| Gene | Representative Drug | Rationale |
|------|-------------------|-----------|
| CD22 | Epratuzumab | Targets CD22 on B cell surface |
| CD274 | Atezolizumab/durvalumab | Targets PD-L1 on cell surface |
| CD80 | Abatacept (CTLA4-Ig) | Blocks CD80 on antigen-presenting cell surface |
| CSF2RB | Targeting antibody | Surface GM-CSF receptor beta chain |
| CSF3R | G-CSF biologics | Surface G-CSF receptor |
| FOLR1 | Mirvetuximab soravtansine | ADC targeting surface folate receptor alpha |
| ICAM1 | Lifitegrast | Targets surface ICAM-1 adhesion molecule |
| IFNAR1 | Anifrolumab | Targets surface type I interferon receptor |
| IL1R1 | Anakinra/rilonacept | Targets surface IL-1 receptor |
| IL6ST | Olamkicept | Targets surface gp130 co-receptor |
| ITGAV | Abituzumab/intetumumab | Targets surface integrin alpha-V |
| ITGB7 | Etrolizumab | Targets surface integrin beta-7 |
| TIGIT | Tiragolumab | Targets surface TIGIT checkpoint receptor |

Gene-disease pairs in this subcategory:
CD22 x SLE; CD274 x Lung cancer, Ovarian cancer; CD80 x CKD, Crohn's, RA, UC; CSF2RB x Neuroblastoma, Pancreatic cancer; CSF3R x MI, Ovarian cancer; FOLR1 x Ovarian cancer; ICAM1 x Crohn's; IFNAR1 x MS, SLE; IL1R1 x RA; IL6ST x RA; ITGAV x MI; ITGB7 x IBD; TIGIT x Lung cancer.

#### A3. Indirect Abundance Modulation (12 unique genes, 17 pairs)

| Gene | Representative Drug | Rationale |
|------|-------------------|-----------|
| CLU | Custirsen (OGX-011) | ASO reduces clusterin mRNA, indirectly lowering protein |
| ELN | Elastase inhibitors | Structural protein; drug relationship to circulating level is indirect |
| GAS6 | AXL kinase inhibitors | Drug targets receptor (AXL), not the GAS6 ligand itself |
| MMP1 | Marimastat/batimastat | Small molecule blocks MMP active site; level and activity correlated for secreted enzymes but not identical |
| MMP3 | Marimastat | Same rationale as MMP1 |
| MMP7 | Marimastat | Same rationale as MMP1 |
| MMP8 | Marimastat | Same rationale as MMP1 |
| MMP9 | Andecaliximab/marimastat | MMP inhibitor; mAb or small molecule |
| MMP12 | MMP-12 inhibitor | Small molecule MMP inhibitor |
| PLA2G2A | Varespladib | Small molecule sPLA2 inhibitor blocks active site |
| SERPINC1 | LMWH/heparin | Drug potentiates antithrombin activity; does not change circulating level |
| SOD1 | Tofersen | ASO reduces SOD1 mRNA, indirectly lowering protein |

Gene-disease pairs in this subcategory:
CLU x Lung cancer; ELN x CKD; GAS6 x Ovarian cancer; MMP1 x Alzheimer's, Lung cancer; MMP3 x Lung cancer; MMP7 x Alzheimer's, Lung cancer; MMP8 x Alzheimer's; MMP9 x Lung cancer, UC; MMP12 x Lung cancer; PLA2G2A x MI, RA; SERPINC1 x Lung cancer, Pancreatic cancer; SOD1 x ALS.

#### B1. Enzyme Inhibitor (20 unique genes, 38 pairs)

| Gene | Representative Drug | Rationale |
|------|-------------------|-----------|
| ACE | Enalapril/ramipril/lisinopril | ACE inhibitor blocks angiotensin-converting enzyme active site |
| ACHE | Donepezil/rivastigmine/galantamine | Cholinesterase inhibitor |
| ALDH5A1 | Vigabatrin/valproate | GABA aminotransferase pathway inhibitor |
| BCHE | Rivastigmine | Butyrylcholinesterase inhibitor |
| CA2 | Acetazolamide | Carbonic anhydrase inhibitor |
| COMT | Entacapone/tolcapone | Catechol-O-methyltransferase inhibitor |
| DDC | Carbidopa | DOPA decarboxylase inhibitor |
| DPP4 | Sitagliptin/saxagliptin | Dipeptidyl peptidase-4 inhibitor |
| EGLN1 | Roxadustat/daprodustat | Prolyl hydroxylase (HIF-PHD) inhibitor |
| F10 | Rivaroxaban/apixaban | Factor Xa inhibitor |
| F2 | Dabigatran | Thrombin (Factor IIa) inhibitor |
| GSR | Glutathione reductase inhibitor | Oxidoreductase inhibitor |
| IMPA1 | Lithium | Inositol monophosphatase inhibitor |
| MME | Sacubitril | Neprilysin inhibitor |
| NOS1 | NOS inhibitor | Neuronal nitric oxide synthase inhibitor |
| NOS2 | iNOS inhibitor | Inducible nitric oxide synthase inhibitor |
| PARP1 | Olaparib/niraparib/rucaparib | Poly(ADP-ribose) polymerase inhibitor |
| PDE5A | Sildenafil/tadalafil | Phosphodiesterase 5 inhibitor |
| SNAP25 | Botulinum toxin type A | Protease (BoNT) cleaves SNAP25 protein |
| TOP2B | Etoposide/doxorubicin | Topoisomerase II inhibitor/poison |

Gene-disease pairs in this subcategory:
ACE x CKD, MI; ACHE x Alzheimer's, MDD, MS, Parkinson's, Schizophrenia; ALDH5A1 x ALS, Bipolar, MDD, Schizophrenia; BCHE x Alzheimer's, Parkinson's; CA2 x Bipolar; COMT x Parkinson's; DDC x Parkinson's; DPP4 x CKD, MI; EGLN1 x CKD; F10 x MI; F2 x MI; GSR x Glioma; IMPA1 x ALS, ASD, Bipolar, MDD, Schizophrenia; MME x CKD, MI; NOS1 x MI; NOS2 x MI; PARP1 x Ovarian cancer; PDE5A x Alzheimer's, CKD, MI; SNAP25 x MS; TOP2B x Lung cancer, Ovarian cancer.

#### B2. Receptor/Kinase Antagonist (19 unique genes, 38 pairs)

| Gene | Representative Drug | Rationale |
|------|-------------------|-----------|
| AGER | Azeliragon (TTP488) | RAGE receptor antagonist |
| CACNA2D3 | Pregabalin/gabapentin | Alpha-2-delta calcium channel subunit modulator |
| CACNB3 | Calcium channel blockers | Calcium channel beta-3 subunit modulator |
| EGFR | Erlotinib/gefitinib/osimertinib | EGFR tyrosine kinase inhibitor |
| EPHB2 | Kinase inhibitor | EphB2 receptor tyrosine kinase inhibitor |
| F2R | Vorapaxar | PAR1 (protease-activated receptor 1) antagonist |
| FGFR2 | Erdafitinib/pemigatinib | FGFR2 kinase inhibitor |
| FGFR3 | Erdafitinib/infigratinib | FGFR3 kinase inhibitor |
| FGFR4 | Futibatinib | FGFR4-selective kinase inhibitor |
| FLT3 | Gilteritinib/midostaurin | FLT3 kinase inhibitor |
| FLT4 | Lenvatinib/axitinib | VEGFR3/multi-kinase inhibitor |
| GRIK2 | Glutamate receptor modulator | Kainate receptor (GluK2) antagonist |
| KDR | Sunitinib/sorafenib | VEGFR2 kinase inhibitor |
| KIT | Imatinib/sunitinib | KIT kinase inhibitor |
| MAP2K1 | Trametinib/selumetinib | MEK1 kinase inhibitor |
| MET | Capmatinib/tepotinib | MET kinase inhibitor |
| PDGFRA | Imatinib/sunitinib | PDGFRA kinase inhibitor |
| PDGFRB | Imatinib | PDGFRB kinase inhibitor |
| RET | Selpercatinib/pralsetinib | RET kinase inhibitor |

Gene-disease pairs in this subcategory:
AGER x Alzheimer's; CACNA2D3 x MS; CACNB3 x MS; EGFR x Lung cancer, Pancreatic cancer; EPHB2 x Thyroid cancer; F2R x MI; FGFR2 x Ovarian cancer; FGFR3 x Alzheimer's, Ovarian cancer, Pancreatic cancer; FGFR4 x Ovarian cancer; FLT3 x Alzheimer's, Pancreatic cancer; FLT4 x Ovarian cancer, Pancreatic cancer; GRIK2 x Bipolar; KDR x Ovarian cancer, Pancreatic cancer; KIT x Alzheimer's, ALS, MS, Ovarian cancer, Pancreatic cancer; MAP2K1 x Glioma, Lung cancer; MET x Lung cancer; PDGFRA x Alzheimer's, ALS, MS, Ovarian cancer, Pancreatic cancer; PDGFRB x Alzheimer's, ALS, MS, Ovarian cancer, Pancreatic cancer; RET x Pancreatic cancer.

#### Excluded from subcategorization

| Gene | Mechanism Class | Reason |
|------|----------------|--------|
| APP | mixed | Both abundance and activity mechanisms apply |
| IGF1R | mixed | Both abundance and activity mechanisms apply |

3 gene-disease pairs excluded (APP x Alzheimer's, IGF1R x ALS, IGF1R x Anorexia).

#### Pair Count Summary

| Subcategory | Unique Genes | Gene-Disease Pairs |
|------------|-------------|-------------------|
| A1. Soluble ligand neutralization | 14 | 22 |
| A2. Surface target depletion | 13 | 20 |
| A3. Indirect abundance modulation | 12 | 17 |
| B1. Enzyme inhibitor | 20 | 38 |
| B2. Receptor/kinase antagonist | 19 | 38 |
| Mixed (excluded) | 2 | 3 |
| **Total** | **80** | **138** |

Collapsed three-way:

| Collapsed Category | Subcategories | Pairs |
|-------------------|--------------|-------|
| Direct neutralization | A1 + A2 | 42 |
| Indirect abundance | A3 | 17 |
| Function blocking | B1 + B2 | 76 |

### 5.5 Classification Edge Cases and Rationale

Several assignments require justification:

1. **IL6R in A1 (soluble) rather than A2 (surface).** IL-6 receptor exists in both membrane-bound and soluble forms. The paper explicitly classifies tocilizumab as abundance-modulating because "the cis-pQTL instruments soluble IL-6 receptor level, and tocilizumab blocks a soluble ligand-receptor axis." Tocilizumab binds both forms, but the pQTL alignment is through the soluble receptor, placing IL6R in A1.

2. **KLKB1 in A1 (soluble) rather than B1 (enzyme).** Plasma kallikrein is a serine protease (enzyme), and lanadelumab inhibits its enzymatic activity. However, lanadelumab is a monoclonal antibody that directly binds and sequesters circulating kallikrein protein. The drug acts on the same circulating molecule the pQTL instruments. The V3.4 classification as abundance_modulating is consistent with this interpretation.

3. **PLAT/PLAU/PLG in A1 (soluble) rather than A3 (indirect).** These are thrombolytics where the drug IS the circulating protein (tPA, urokinase, plasminogen). The drug-pQTL alignment is maximal because the genetic instrument directly models the drug's mechanism: a pQTL variant that increases endogenous tPA level mimics the effect of exogenous tPA administration. The direction is augmentation rather than neutralization, but the causal axis alignment is identical.

4. **MMPs (MMP1/3/7/8/9/12) in A3 (indirect) rather than A1 (soluble).** MMPs are secreted enzymes circulating in blood. The Phase III MMP inhibitors (marimastat, batimastat, tanomastat) are small molecule active-site inhibitors that block enzymatic activity without changing circulating protein levels. The pQTL instruments MMP level; the drug blocks MMP activity. For secreted enzymes, level and activity are correlated (more enzyme = more activity), which justifies the V3.4 abundance_modulating classification. But the alignment is indirect compared to anti-TNF (where the antibody directly removes the protein), placing MMPs in A3.

5. **GAS6 in A3 (indirect) rather than A1 (soluble).** GAS6 is a soluble ligand, but the drugs targeting this pathway (AXL kinase inhibitors) act on the receptor, not the ligand. The pQTL instruments GAS6 level, while the drug blocks AXL signaling. This is an indirect relationship.

6. **SERPINC1 in A3 (indirect) rather than A1 (soluble).** Antithrombin (SERPINC1) is a circulating protein, and the drugs include heparin/LMWH which potentiate antithrombin's anticoagulant activity through a conformational mechanism. The drug does not change circulating antithrombin level; it changes the protein's activity state. This is the inverse of the activity-blocking mismatch: the pQTL instruments level, and the drug modulates activity.

### 5.6 Hypotheses

**H5.1 (Gradient hypothesis).** Balanced accuracy decreases monotonically across the five subcategories in the predicted order:

BA(A1, soluble) > BA(A2, surface) > BA(A3, indirect) >= BA(B1, enzyme) ~ BA(B2, receptor/kinase)

Predicted BA ranges (from theory, before seeing results):
- A1 (soluble ligand, n=22): BA in [0.65, 0.75]. The drug and pQTL instrument act on the identical circulating molecule. This is where MR should perform best.
- A2 (surface target, n=20): BA in [0.55, 0.65]. Moderate alignment; shed soluble form correlates with surface expression but imperfectly.
- A3 (indirect, n=17): BA in [0.50, 0.55]. Weak alignment; most pairs in this category involve small molecule enzyme inhibitors of secreted proteases, where level and activity are correlated but the drug does not change level.
- B1 (enzyme, n=38): BA in [0.48, 0.53]. Null. The pQTL instruments level; the drug blocks activity. Complete axis mismatch.
- B2 (receptor/kinase, n=38): BA in [0.48, 0.53]. Null. Same structural mismatch as B1.

**H5.2 (True positive concentration).** At least 60% of all abundance-modulating true positives (pairs where MR p<0.05 and outcome is SUCCESS) will be in the A1 (soluble ligand neutralization) subcategory. This follows from the theory: the strongest drug-pQTL alignment should produce the most true positive MR predictions.

**H5.3 (Three-way gradient).** Under the collapsed taxonomy:

BA(direct neutralization) > BA(indirect abundance) > BA(function blocking)

Predicted:
- Direct neutralization (n=42): BA in [0.60, 0.70]
- Indirect abundance (n=17): BA in [0.50, 0.55]
- Function blocking (n=76): BA in [0.50, 0.52]

**H5.4 (Sensitivity concentration).** Sensitivity (TP rate) in the A1 subcategory will be at least twice the sensitivity in the combined A2+A3 subcategories. The A1 subcategory should concentrate the true signal, while A2 and A3 contribute mainly true negatives and false negatives.

### 5.7 Statistical Plan

All analyses use the V3.4 analyzed set (n=138 with MR results). Mixed-class pairs (n=3) are excluded.

1. **Per-subcategory BA.** For each of the five subcategories, compute balanced accuracy (average of sensitivity and specificity) with 10,000-iteration bootstrap 90% CI.

2. **Jonckheere-Terpstra trend test.** Test for an ordered trend in MR informativeness across the five subcategories (ordered A1 > A2 > A3 > B1 > B2). The test statistic is the Jonckheere-Terpstra J, computed as the sum of Mann-Whitney U statistics for all ordered pairs of groups. Significance is assessed by 10,000 label permutations (permuting subcategory labels while preserving within-group sizes). The outcome variable for the trend test is the binary MR-correct indicator (1 if the MR prediction matches the trial outcome, 0 otherwise).

3. **Pairwise permutation tests.** For each adjacent pair of subcategories (A1 vs A2, A2 vs A3, A3 vs B1, B1 vs B2), test the difference in BA by 10,000 label permutations. Report raw and Holm-corrected p-values for 4 comparisons.

4. **Three-way analysis.** Repeat the trend test and pairwise permutation tests for the collapsed three-way taxonomy (direct neutralization, indirect abundance, function blocking). Two pairwise comparisons, no multiplicity adjustment needed for two tests.

5. **True positive distribution.** Tabulate the number of true positives and false positives in each subcategory. Compute the fraction of abundance-modulating TPs in A1. Test H5.2 against a null of proportional distribution (expected A1 fraction = 22/59 = 0.373 under the null that TPs are distributed proportionally to subcategory size) using a one-sided exact binomial test.

6. **Sensitivity by subcategory.** Compute sensitivity (TP / [TP + FN]) for each subcategory. Test H5.4 by comparing A1 sensitivity to (A2+A3) pooled sensitivity via Fisher exact test on the 2x2 table of (TP, FN) x (A1, A2+A3).

### 5.8 Decision Rules

- **Gradient confirmed:** If the Jonckheere-Terpstra trend test yields p < 0.10 AND the A1 subcategory BA 90% CI lower bound exceeds 0.55, report the gradient as the primary finding of Analysis 5.
- **Gradient consistent but underpowered:** If the point estimates follow the predicted monotonic order (BA_A1 > BA_A2 > BA_A3 >= BA_B1 ~ BA_B2) but the trend test p >= 0.10, report as "consistent with the gradient hypothesis but underpowered for formal confirmation" and state the sample size needed for 80% power.
- **Gradient violated:** If any non-adjacent inversion occurs (e.g., BA_A3 > BA_A1, or BA_B1 > BA_A2), report the deviation, examine whether it is driven by a specific disease or a few pairs, and discuss possible explanations.
- **True positive concentration confirmed:** If >= 60% of abundance TPs fall in A1 AND the binomial test p < 0.10, report as confirmed.
- **True positive concentration failed:** If < 50% of abundance TPs fall in A1, the concentration hypothesis is rejected. Examine which subcategory captures the unexpected TPs and whether the mechanism classification is correct for those genes.

---

## Analysis 6: Open Targets Pipeline Blind Spot

### 6.1 Rationale

Drug target validation pipelines such as Open Targets aggregate genetic evidence to prioritize targets. The genetic associations datasource (GWAS credible sets from GWAS Catalog, UK Biobank, FinnGen) provides a key evidence channel for target-disease associations. If pQTL-based genetic evidence is structurally uninformative for activity-blocking targets (as established in the V3.4 primary result), then genetic evidence pipelines that rely heavily on this channel should systematically underweight activity-blocking targets.

The prediction extends beyond pQTL-MR. For diseases where the therapeutic mechanism involves blocking protein *function* (enzyme activity, receptor signaling), common genetic variants near the gene that affect protein *level* (through expression or splicing) may also fail to associate with disease, because the causal path runs through function, not level. GWAS credible sets map common-variant associations to genes through locus-to-gene (L2G) scoring, but if the functional variants that matter are rare coding variants affecting protein function (not common variants affecting protein level), the GWAS evidence channel will miss them.

This predicts a systematic pipeline blind spot: activity-blocking targets with proven clinical efficacy should have lower genetic evidence scores than abundance-modulating targets on Open Targets.

### 6.2 Data Collection

For each of the 138 gene-disease pairs in the V3.4 analyzed set (plus the 57 missing pairs if possible), query the Open Targets Platform GraphQL API for:

1. **Overall association score** for the gene-disease pair
2. **Genetic associations datasource score** (specifically the GWAS credible set evidence type)
3. **Number of GWAS credible set evidence entries** linking the gene to the disease

The query uses the Open Targets GraphQL endpoint (`https://api.platform.opentargets.org/api/v4/graphql`). Scores are normalized to [0, 1] by Open Targets.

### 6.3 Hypotheses

**H6.1 (Activity-blocking GWAS gap).** Among activity-blocking gene-disease pairs with clinical SUCCESS (approved drugs), at least 40% will have zero GWAS credible set evidence score on Open Targets. These are targets where the drug works (RCT evidence) but the genetic evidence pipeline cannot detect the target-disease link because the relevant variants affect function, not level.

**H6.2 (Mechanism-differential genetic evidence).** The mean GWAS credible set score will be higher for abundance-modulating pairs than for activity-blocking pairs. Predicted difference: >= 0.05 on the Open Targets normalized scale (0-1). The direction follows from the theory: variants affecting protein level (which GWAS preferentially detects) are informative for diseases where level is causal (abundance targets) and uninformative where function is causal (activity targets).

**H6.3 (Null MR + null GWAS co-occurrence).** Among activity-blocking pairs where MR p >= 0.05 (MR null), at least 50% will also have GWAS credible set score = 0. This tests whether MR nullity and GWAS absence are correlated for activity-blocking targets, as the theory predicts (both fail for the same structural reason: common variants near the gene affect level, and level is not the causal axis).

**H6.4 (MR-GWAS correlation by stratum).** The Spearman correlation between -log10(MR p-value) and GWAS credible set score will be positive and significant within the abundance-modulating stratum (both capture the same level-mediated biology) and null within the activity-blocking stratum (both are uninformative, so neither is correlated with the other).

**H6.5 (Clinical success without genetic evidence).** Among all 138 analyzed pairs with clinical SUCCESS, the fraction with zero GWAS evidence will be higher for activity-blocking successes than for abundance-modulating successes. Predicted: >= 50% of activity-blocking successes have zero GWAS evidence vs <= 30% of abundance-modulating successes.

### 6.4 Statistical Plan

1. **Proportions with zero GWAS evidence.** For each mechanism stratum (abundance, activity), compute the fraction of pairs with GWAS credible set score = 0. Compare proportions by Fisher exact test.

2. **Distributional comparison.** Compare the distribution of GWAS scores between abundance and activity strata using Mann-Whitney U test (two-sided).

3. **MR-GWAS correlation.** Compute Spearman rho between -log10(MR p) and GWAS score, overall and separately for abundance and activity strata. Report p-values from permutation (10,000 iterations).

4. **Success-specific analysis.** Restrict to SUCCESS pairs only. Compare GWAS evidence rates between mechanism strata via Fisher exact test.

5. **Subcategory analysis.** Repeat the above analyses using the five-way subcategory taxonomy from Analysis 5. Test whether the gradient (A1 > A2 > A3 > B1 ~ B2) holds for genetic evidence scores as well as MR performance.

### 6.5 Decision Rules

- **Blind spot confirmed:** If H6.1 (>= 40% zero-evidence activity successes) AND H6.2 (significant Mann-Whitney, abundance > activity) both hold, report the Open Targets blind spot as confirmed.
- **Partial support:** If H6.2 holds but H6.1 does not reach 40%, report the differential as confirmed but the absolute gap as smaller than predicted.
- **Blind spot rejected:** If GWAS evidence scores are similar across mechanism strata (Mann-Whitney p > 0.10 and proportion difference < 5 percentage points), the blind spot hypothesis is rejected. Discuss why the theory's prediction failed (e.g., GWAS detects coding variants that affect function, not just level-altering eQTL/pQTL variants).

---

## Analysis 7: ACE Inhibitor Case Study

### 7.1 Selection Rationale

ACE (angiotensin-converting enzyme) is pre-selected as the case study to illustrate the domain-of-validity framework. The selection is based on five criteria, all evaluable without reference to MR results:

1. **Activity-blocking mechanism (B1, enzyme inhibitor).** ACE inhibitors (enalapril, ramipril, lisinopril, captopril, perindopril) block the active site of angiotensin-converting enzyme. The drug inhibits ACE enzymatic activity---specifically, the conversion of angiotensin I to angiotensin II---without removing ACE protein from circulation. ACE protein levels are typically *increased* by ACE inhibitor treatment (compensatory upregulation via renin-angiotensin feedback), while enzymatic activity is decreased. This is the clearest possible illustration of the level-vs-function mismatch.

2. **Clinical SUCCESS.** ACE inhibitors have regulatory approval for hypertension, heart failure, diabetic nephropathy, post-myocardial infarction, and chronic kidney disease. They have been validated in multiple landmark Phase III RCTs:
   - CONSENSUS (1987): enalapril reduced mortality 40% in severe heart failure
   - SOLVD (1991): enalapril reduced mortality/hospitalization in heart failure
   - HOPE (2000): ramipril reduced cardiovascular events 22% in high-risk patients
   - ALLHAT (2002): lisinopril confirmed as effective first-line antihypertensive
   - RENAAL/IDNT (2001): ACE inhibitors/ARBs slowed diabetic nephropathy progression

3. **Null MR prediction.** The cis-pQTL for ACE instruments circulating ACE protein level. ACE inhibitors block enzymatic activity at the active site while leaving the protein intact and measurable. A genetic variant that lowers ACE protein level is NOT equivalent to a drug that blocks ACE activity: less ACE protein means proportionally less of both angiotensin-converting activity and other ACE substrates (such as bradykinin), while an ACE inhibitor selectively blocks the active site while leaving structural/non-catalytic functions intact. The MR estimate answers "what happens with lifelong lower ACE protein?" when the clinically relevant question is "what happens when ACE enzymatic activity is pharmacologically blocked?"

4. **High prescription volume.** ACE inhibitors are among the most widely prescribed drug classes globally: approximately 160 million prescriptions per year in the United States alone. They are first-line therapy for hypertension in multiple national guidelines (JNC, ESC/ESH, NICE).

5. **Mechanistic clarity.** The renin-angiotensin-aldosterone system is one of the best-characterized drug target pathways in medicine. The disconnect between "less ACE protein" (pQTL model) and "blocked ACE activity" (drug mechanism) is pharmacologically unambiguous.

### 7.2 Case Study Structure

The case study will present:

1. **The mechanistic disconnect.** A cis-pQTL for ACE instruments circulating ACE protein concentration (measured by aptamer/antibody binding to the ACE protein, regardless of enzymatic state). An ACE inhibitor blocks the catalytic site of ACE protein that is present at normal or elevated levels. The pQTL models "what if you had less ACE protein throughout life?" The drug does something categorically different: "what if the ACE protein you have cannot convert angiotensin I to angiotensin II?"

2. **The MR result.** ACE appears in the dataset as ACE x CKD and ACE x MI, both classified as activity_blocking. The zero-parameter MR classifier is applied to both pairs. The predicted result is MR p >= 0.05 for both (null), because the pQTL perturbation (less ACE protein) does not model the drug perturbation (blocked ACE activity).

3. **The clinical evidence.** Present the RCT evidence demonstrating unambiguous clinical efficacy of ACE inhibitors for both MI prevention and CKD progression.

4. **The pipeline implication.** Under an undifferentiated pQTL-MR screen (as used in Open Targets and other genetic evidence pipelines), ACE would receive no genetic support from the pQTL channel for either indication. A target validation pipeline that interprets MR silence as evidence against the target would deprioritize ACE---one of the most validated drug targets in cardiovascular medicine.

### 7.3 Hypotheses

**H7.1 (MR null for ACE).** Both ACE x MI and ACE x CKD will have MR p >= 0.05. This is predicted by the activity-blocking mechanism: the pQTL instruments ACE level, and the drug blocks ACE activity.

**H7.2 (Open Targets genetic evidence gap for ACE).** ACE will have low or zero GWAS credible set evidence for MI and/or CKD on Open Targets through the genetic associations channel, despite being a validated target for both diseases.

**H7.3 (Direction paradox).** If ACE has a nominally non-null MR result (p < 0.20), the effect direction may be *paradoxical*: lower ACE protein (the pQTL-lowering allele) might not associate with lower cardiovascular risk, because reduced ACE protein also reduces bradykinin degradation, and bradykinin has cardioprotective effects that complicate the directionality.

### 7.4 Statistical Plan

The case study is illustrative rather than hypothesis-testing. It demonstrates the domain-of-validity concept with a concrete, well-known drug class. No novel statistical tests are required beyond:

1. Report the MR p-value and effect estimate (beta, SE) for ACE x MI and ACE x CKD.
2. Report the Open Targets GWAS credible set evidence score for ACE-MI and ACE-CKD.
3. Tabulate the landmark RCT evidence (trial name, year, intervention, primary endpoint, result, N) for ACE inhibitors in cardiovascular disease and CKD.
4. Present the contrast: null MR + null genetic evidence + proven clinical efficacy.

### 7.5 Decision Rules

- **Case study supports framework:** If H7.1 is confirmed (MR null for both ACE pairs), the case study is presented as illustrative of the activity-blocking null.
- **Case study complications:** If ACE x MI or ACE x CKD unexpectedly reaches MR p < 0.05, investigate whether the pQTL variant is in LD with a coding variant that affects ACE enzymatic activity (not just level), which would mean the instrument captures function as well as level. Report this as a boundary case where the level-function distinction is blurred by genetics.
- **Open Targets gap not present:** If ACE has strong GWAS evidence for MI or CKD on Open Targets, investigate the evidence source. It may reflect a disease GWAS signal near the ACE locus that is driven by variants affecting ACE function (coding variants, splice variants) rather than ACE level, which would represent a different genetic evidence channel than pQTL-MR.

---

## Classification Audit Trail

The five-way subcategory classification was assigned based on:

1. The drug's pharmacological mechanism of action (how the drug physically interacts with the target protein)
2. Whether the target protein circulates as a soluble molecule, is shed from cell surfaces, or is membrane-bound
3. Whether the drug directly engages the protein the pQTL instruments (neutralization, depletion) or acts indirectly (ASO targeting mRNA, inhibitor of a secreted enzyme's active site, drug targeting a different molecule in the pathway)

The V3.4 binary mechanism_class was not altered. Every gene classified as abundance_modulating in V3.4 is assigned to one of A1/A2/A3. Every gene classified as activity_blocking in V3.4 is assigned to one of B1/B2. The subcategories are a strict refinement of the binary classification.

Representative drug assignments were determined from the author's pharmacological knowledge of Phase III drug-target pairs in the Open Targets database. Where multiple drugs target the same gene with different mechanisms (e.g., anti-MMP9 antibody vs small molecule MMP inhibitor), the predominant Phase III mechanism was used.

---

## Integrity Statement

This document specifies three analyses (Analysis 5: Expanded Taxonomy, Analysis 6: Open Targets Blind Spot, Analysis 7: ACE Case Study) with frozen subcategory classifications, testable hypotheses, and quantitative decision rules. All predictions are derived from the paper's causal framework (pQTL instruments circulating protein abundance; MR is informative when drug and instrument share a causal axis) applied to each gene's pharmacological mechanism. No results from the expanded taxonomy analysis, Open Targets queries, or ACE case study have been inspected before writing this document.

The gene-to-subcategory classification table is frozen as written above. Any post-hoc changes to classifications must be documented in the deviation log with justification.

This document will be SHA-256 hashed before any Analysis 5, 6, or 7 results are computed. The hash and timestamp will be recorded in `protocol/v6/PRESPEC_V6_sha256.txt`.
