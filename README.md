# When It Jumped?

**Candidate transposon transfer in baculoviruses: a combined exploratory and targeted phylogenetic analysis**

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![R 4.3+](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/)
[![Status: complete](https://img.shields.io/badge/status-complete-brightgreen.svg)]()

**Authors:** Polina Malysheva, Viktoria Orlova

**Course:** Microbial Omics / Bioinformatics, ITMO University, October 2026

---

## Abstract

Horizontal transfer (HT) of transposable elements (TEs) between hosts and viruses is increasingly recognized as a driver of genome evolution. A 2026 *Science* study (Wint3r et al., DOI: 10.1126/science.aeh3441) reported two insertions of the human BC200 non-coding RNA gene in molluscum contagiosum virus (MCV), likely acquired through LINE-1-mediated retrotransposition. We asked whether similar signals are detectable in **baculoviruses** — dsDNA viruses of Lepidoptera — using a dual approach:

1. **Exploratory** — de novo repeat annotation (RepeatModeler + RepeatMasker vs Dfam) to survey TE families across 4 baculovirus genomes.
2. **Targeted** — BLAST and phylogenetic analysis of the **piggyBac** transposon across 5 viral and 4 host genomes.

The exploratory annotation identified candidate TE-related sequences, but most long candidates were tandem repeats, and Dfam is poorly curated for Lepidoptera. Targeted BLAST identified **90 piggyBac homologs**. Phylogenetic analysis of **85 unique taxa** (GTR+F+R3, 1000 UFBoot) recovered a well-supported viral clade (**100% UFBoot**) that is sister to a sub-clade of *Spodoptera frugiperda* sequences (**93.6% UFBoot**). These results identify candidate viral piggyBac homologs that warrant further validation; the direction and timing of transfer cannot be established from this analysis alone.

---

## 1. Introduction

### 1.1. Background

Transposable elements are mobile genetic elements that can move within and occasionally between genomes. A 2026 *Science* paper (Wint3r et al.) reported two insertions of the human **BC200** gene in molluscum contagiosum virus (MCV), a poxvirus. BC200 is a functional non-coding RNA gene that retains mobilization capacity through LINE-1-encoded machinery; it is not a canonical retrotransposon, and the MCV insertions are best described as **LINE-1-mediated retrotransposition events**.

### 1.2. Why baculoviruses?

Baculoviruses are dsDNA viruses with narrow host ranges (mostly Lepidoptera) that replicate in the host nucleus. Importantly, the canonical **piggyBac** element (IFP2) was itself first characterized as an insertion in the AcMNPV genome (Fraser et al., 1983), so the presence of piggyBac-related sequences in baculoviruses is not a new observation. What remains unclear is whether viral and host piggyBac copies share a **recent evolutionary history**, or whether they represent **independent, ancient acquisitions**.

### 1.3. Hypotheses

- **H₀:** Viral piggyBac homologs are evolutionarily distinct from host copies (consistent with ancient acquisition or independent evolution).
- **H₁:** Viral piggyBac homologs share a recent common ancestor with host copies (consistent with horizontal transfer).

> **Note on scope:** These hypotheses can be tested with phylogenetics, but the direction and timing of any transfer cannot be established without additional evidence (molecular clock, target-site duplications, population data).

---

## 2. Data

### 2.1. Viral genomes (5)

| Virus | Accession | Length | Host |
|---|---|---|---|
| AcMNPV | NC_001623.1 | 133,894 bp | *Autographa californica* |
| BmNPV | NC_001962.1 | 128,413 bp | *Bombyx mori* |
| CpGV | NC_002816.1 | 123,500 bp | *Cydia pomonella* |
| NeseNPV | NC_002169.1 | 135,611 bp | *Neodiprion sertifer* |
| LdMNPV | NC_001973.1 | 161,046 bp | *Lymantria dispar* |

Additional genomes used in Part I: SfMNPV (NC_009011.2), HearNPV (NC_002654.1).

### 2.2. Host genomes (4)

| Host | Accession | Assembly level |
|---|---|---|
| *Bombyx mori* | GCF_014905235.1 | chromosome |
| *Spodoptera frugiperda* | GCF_023101765.2 | chromosome |
| *Helicoverpa armigera* | GCF_030705265.1 | chromosome |
| *Cydia pomonella* | GCF_033807575.1 | chromosome |

### 2.3. Query

piggyBac reference from *Trichoplusia ni* (**J04364.2**, ~2,475 bp; Fraser et al., 1983) — `data/query/piggyBac_IFP2_J04364.2.fasta`.

Host genomes (~350–620 MB each) are not committed to the repository. Re-download via `COMMANDS.md`, section 1.

---

## 3. Repository structure

```
When-it-jumped/
├── README.md                       ← you are here
├── REPORT.md                       ← full report with discussion
├── COMMANDS.md                     ← all commands used in the pipeline
├── LICENSE                         ← MIT
│
├── data/
│   ├── viruses/                    ← 7 baculovirus genomes (RefSeq)
│   └── query/                      ← piggyBac IFP2 reference (J04364.2)
│
├── scripts/
│   ├── fig1_ml_tree.R              ← Figure 1 (ggtree)
│   ├── fig2_heatmap.R              ← Figure 2 (pheatmap)
│   ├── fig3_summary.R              ← Figure 3 (patchwork + ggplot2)
│   ├── fig4_mds.R                  ← Figure 4 (ggplot2 + ggrepel)
│   ├── check_tandem.py             ← tandem-repeat detector (k=3..6)
│   ├── classifier.py               ← k-mer Random Forest classifier
│   ├── nj_tree.py                  ← NJ tree (Bio.Phylo)
│   ├── mds_heatmap.py              ← Python version of heatmap + MDS
│   ├── extract_hits.py             ← BLAST hit extraction
│   ├── unique_ids.py               ← FASTA header cleanup
│   ├── locus_filter.py             ← locus-specificity filter
│   └── plot_tree.R                 ← quick R tree visualization
│
├── results/
│   ├── te_annotation/              ← 23 TE loci + 6 long candidates
│   │   ├── all_TE_filtered.bed
│   │   ├── TE_long.bed
│   │   └── TE_long.fasta
│   ├── blast/                      ← 4 filtered BLAST tables (.tsv.gz)
│   └── phylogeny/                  ← Part I alignment + trees
│       ├── combined_aligned.fasta
│       ├── combined_iqtree.treefile
│       ├── combined_tree_NJ.nwk
│       └── combined_iqtree.iqtree
│
├── figures/                        ← publication-ready (PDF + PNG 300 dpi)
│   ├── fig1_ml_tree.pdf / .png     ← Figure 1
│   ├── fig2_heatmap.pdf / .png     ← Figure 2
│   ├── fig3_summary.pdf / .png     ← Figure 3
│   └── fig4_mds.pdf / .png         ← Figure 4
│
└── piggyBac_analysis/              ← Part II (collaborator)
    ├── README.md
    ├── data/
    ├── scripts/
    ├── results/
    └── figures/
```

---

## 4. Reproducibility

### 4.1. Software

| Tool | Version | Purpose |
|---|---|---|
| RepeatModeler | 2.0.8 | de novo TE library |
| RepeatMasker | 4.2.3 | TE annotation vs Dfam |
| BLAST+ | 2.14.1 | targeted search |
| MAFFT | v7.526 | multiple alignment |
| trimAl | 1.x | alignment trimming (Part II) |
| IQ-TREE | 3.1.4 (Part I) / 2.x (Part II) | ML phylogeny |
| samtools | 1.20+ | FASTA indexing |
| bedtools | 2.31+ | sequence extraction |
| Python | 3.11 | pipeline scripting |
| R | 4.3.3 | figures (ape, ggtree, pheatmap, patchwork, ggrepel) |

Install everything via conda/mamba:

```bash
mamba install -c bioconda -c conda-forge -y \
    blast repeatmasker repeatmodeler bedtools samtools mafft iqtree trimAl
mamba install -c conda-forge -y \
    pandas seaborn matplotlib biopython scikit-learn \
    r-ape r-ggplot2 r-ggrepel r-patchwork r-pheatmap
mamba install -c bioconda -c conda-forge -y \
    bioconductor-ggtree bioconductor-treeio
pip install pyfaidx
```

### 4.2. How to reproduce

**Step 1 — download data** (see `COMMANDS.md` section 1):

```bash
# Viruses
cd data/viruses
for acc in NC_001623.1 NC_001962.1 NC_002816.1 NC_002654.1 \
           NC_009011.2 NC_002169.1 NC_001973.1; do
    efetch -db nuccore -id "$acc" -format fasta > "${acc}.fasta"
done

# Hosts (large; not in repo) — follow COMMANDS.md, section 1.2

# piggyBac reference
efetch -db nuccore -id J04364.2 -format fasta > data/query/piggyBac_IFP2_J04364.2.fasta
```

**Step 2 — run pipeline** (see `COMMANDS.md` sections 2–5):

```bash
# Exploratory
RepeatModeler -database baculo_db -threads 4
RepeatMasker -lib Dfam-RepeatMasker.lib -engine ncbi -no_is -norna -gff ...
python3 scripts/check_tandem.py

# BLAST vs hosts (sensitive)
blastn -task blastn -word_size 7 -evalue 10 -dust no ...

# Phylogeny
mafft --auto combined.fasta > combined_aligned.fasta
iqtree -s combined_aligned.fasta -m MFP -bb 1000 -nt 4 -pre tree_piggyBac
```

**Step 3 — regenerate figures**:

```bash
for s in fig1_ml_tree fig2_heatmap fig3_summary fig4_mds; do
    Rscript "scripts/$s.R"
done
```

Full command log: **[COMMANDS.md](COMMANDS.md)**.

---

## 5. Results

### 5.1. Part I — Exploratory annotation

| Run | Total loci | Interspersed TE | TE fraction |
|---|---|---|---|
| `-nolow` | 28 | 28 | **100%** |
| without `-nolow` | 148 | **23** | **16%** |

**The `-nolow` flag produced ~80% false-positive TE annotations.** After filtering `Simple_repeat`, `Low_complexity`, `Satellite`, `tRNA`, `rRNA`, `snRNA`, 23 interspersed TE loci remained. Of 6 candidates ≥100 bp, **4 were tandem repeats** (GATGA/TCATC microsatellites). The only retained candidate — a **hAT-Tip100 fragment in BmNPV** — was placed **outside** all Dfam clusters at **87% UFBoot**, consistent with an ancient, divergent origin.

> ⚠️ A single locus cannot establish a "locked" TE compartment. Absence of a close host homolog in Dfam is not absence in nature.

**Figures:** [Figure 1](#figure-1) (ML tree), [Figure 2](#figure-2) (distance heatmap), [Figure 3](#figure-3) (exploratory summary), [Figure 4](#figure-4) (MDS).

### 5.2. Part II — Targeted piggyBac search

**BLAST results: 90 hits total.**

| Genome | Hits |
|---|---|
| *S. frugiperda* | 54 |
| *B. mori* | 18 |
| *H. armigera* | 14 |
| AcMNPV | 2 |
| BmNPV | 2 |
| *C. pomonella*, CpGV, NeseNPV, LdMNPV | **0** |

**Phylogenetic dataset:** 91 sequences (1 reference + 90 hits) → **85 unique taxa**, **2,472 sites**, model **GTR+F+R3**.

| Clade | UFBoot / SH-aLRT | Observation |
|---|---|---|
| Viral clade (AcMNPV + BmNPV) | **100 / 100** | Monophyletic |
| Viral clade + *S. frugiperda* sub-clade | **93.6 / 99** | Sister to subset of *Spodoptera* copies |
| Host clades (*Bombyx*, *Spodoptera*, *Helicoverpa*) | 90–100 | Sort by species |

**Long branch:** viral clade branch = **3.62 subs/site** — may reflect rapid divergence, incomplete sampling, or rate heterogeneity.

### 5.3. What the tree does NOT show

- **Direction** — sister relationship is consistent with host→virus, virus→host, or shared unsampled source.
- **Timing** — no molecular clock applied; "recent" cannot be inferred from topology alone.
- **Source species** — *S. frugiperda* is a candidate, not a proven donor.

---

## 6. Figures

### Figure 1
**Maximum-likelihood phylogeny of viral TE candidates and Dfam consensus sequences.**
IQ-TREE 3.1.4, HKY+F+G4, UFBoot 1000, 13 taxa, 2,222 sites. The viral hAT candidate (red) is sister to all other sequences (87% bootstrap). `Academ` and `LINE` (CpGV) form a clade with branch lengths ≈ 0 — they are the same GATGA tandem repeat mis-annotated twice. Dfam references (blue) form a separate clade.

[![Figure 1](figures/fig1_ml_tree.png)](figures/fig1_ml_tree.pdf)

### Figure 2
**Pairwise identity-distance heatmap.**
Near-zero block between `Academ` and `LINE` confirms their tandem-repeat identity. Viral candidates cluster separately from Dfam consensus.

[![Figure 2](figures/fig2_heatmap.png)](figures/fig2_heatmap.pdf)

### Figure 3
**Exploratory annotation summary.**
A — `-nolow` vs no `-nolow` (28 vs 23 interspersed loci; ~80% false positives with `-nolow`). B — of 6 long candidates, 4 are tandem repeats and 2 retained.

[![Figure 3](figures/fig3_summary.png)](figures/fig3_summary.pdf)

### Figure 4
**MDS of pairwise identity distances.**
Viral candidates (red) separate from Dfam consensus (blue). `hAT` is a borderline case at MDS1 ≈ −0.3.

[![Figure 4](figures/fig4_mds.png)](figures/fig4_mds.pdf)

---

## 7. Comparison with MCV/BC200

| Feature | MCV/BC200 (Science 2026) | This study (piggyBac) |
|---|---|---|
| Element | BC200 non-coding RNA (LINE-1-mobilized) | piggyBac DNA transposon |
| Evidence | 2 insertions, ~74% identity over 143 bp, TSDs, phylogenetic nesting | Viral clade nested within host diversity (93.6% UFBoot) |
| Direction | Human → MCV (inferred) | **Not established** |
| Timing | Modern human history (inferred) | **Not established** |
| **Conclusion** | **Recent HT** | **Candidate HGT warranting validation** |

---

## 8. Discussion

### 8.1. Strengths of the targeted approach

Blind RepeatMasker/Dfam annotation in Lepidoptera is unreliable (~80% false positives with `-nolow`; most long candidates are tandem repeats). Targeted BLAST with a curated query (J04364.2) is more robust for identifying specific TE families and allowed recovery of divergent, full-length viral homologs that consensus libraries would miss.

### 8.2. What the piggyBac tree adds

The monophyletic viral clade nested within the host tree is consistent with a shared evolutionary history between viral and host piggyBac copies. Follow-up work should include:

- **Target-site duplication (TSD)** analysis at viral insertion sites.
- **Molecular clock dating** (BEAST) to estimate divergence time.
- **Expanded taxon sampling** (more baculoviruses and Lepidoptera) to identify the donor lineage.
- **Terminal inverted repeat (TIR)** analysis to test whether viral copies retain a functional transposon structure.

### 8.3. Reconciling exploratory and targeted results

The two analyses answer different questions. The exploratory search identifies candidate TE-related sequences across families; absence of recent host homologs for a single hAT locus cannot establish a "locked compartment." The targeted search shows that at least one family (**piggyBac**) has viral copies phylogenetically nested within host diversity. Together, they indicate that baculovirus genomes contain a **mix of ancient, divergent TE fragments** and at least one family with a **plausible shared evolutionary history with hosts** — but do not establish direction, timing, or a general "recent HGT" phenomenon.

---

## 9. Limitations

1. **No molecular clock** — timing not estimated.
2. **No TSD analysis** — recent piggyBac activity not independently verified.
3. **Limited taxon sampling** — 5 viruses and 4 hosts; the true donor lineage may be unsampled.
4. **Single-locus exploratory analysis** — cannot generalize about "locked" vs. "active" TE compartments.
5. **Model selection** — results should be compared across alternative models (GTR+F+G4, HKY+F+R3) as sensitivity check.
6. **Alignment** — `--adjustdirection` and `--auto` used; alternative trimming strategies (`--gappyout`, `--strict`) may be informative.
7. **BC200 framing** — described as a "retrotransposon" in earlier drafts; corrected to "LINE-1-mobilized non-coding RNA gene."

---

## 10. Data availability

All sequences are publicly available from NCBI:

- Viral genomes (RefSeq): NC_001623.1, NC_001962.1, NC_002816.1, NC_002654.1, NC_009011.2, NC_002169.1, NC_001973.1
- Host genomes (RefSeq): GCF_014905235.1, GCF_023101765.2, GCF_030705265.1, GCF_033807575.1
- piggyBac reference: J04364.2 (IFP2, Fraser et al., 1983)

Dfam consensus library: https://www.dfam.org/releases/current/families/Dfam-RepeatMasker.lib.gz

---

## 11. References

1. Pu Gao et al., Escape of the BC200 gene to a human poxvirus reveals its persistent transposition in primates.Science393,1342-1346(2026).DOI:10.1126/science.aeh3441
2. Fraser, M. J., Smith, G. E., & Summers, M. D. (1983). Acquisition of host cell DNA sequences by baculoviruses. *Journal of Virology*, 47(2), 287–300.
3. Smit, A. F. A., Hubley, R., & Green, P. (2013–2015). RepeatMasker Open-4.0.
4. Flynn, J. M. et al. (2020). RepeatModeler2 for automated genomic discovery of transposable element families. *PNAS*, 117(17), 9451–9457.
5. Storer, J. et al. (2021). The Dfam community resource of transposable element families. *Mobile DNA*, 12, 2.
6. Minh, B. Q. et al. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference. *Molecular Biology and Evolution*, 37(5), 1530–1534.
7. Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7. *Molecular Biology and Evolution*, 30(4), 772–780.
8. Yu, G. et al. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees. *Methods in Ecology and Evolution*, 8(1), 28–36.

---

## 12. License

MIT License — see [LICENSE](LICENSE).
