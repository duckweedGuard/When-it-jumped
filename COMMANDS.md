# Commands Log — When It Jumped Project

All commands used in the project, grouped by stage, with explanatory comments.

## 0. Environment setup

```bash
conda create -n metagenome_project python=3.11 r-base=4.3 -y
conda activate metagenome_project
conda install -c bioconda -c conda-forge -y blast repeatmasker repeatmodeler bedtools hmmer mafft iqtree samtools
conda install -c conda-forge -y pandas seaborn matplotlib biopython scikit-learn r-ape
pip install pyfaidx
```

NCBI Datasets CLI:
```bash
cd ~/ITMO/omics/tools
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod +x datasets
```

Entrez Direct (efetch/esearch/xtract):
```bash
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
```

## 1. Data acquisition

Viral genomes (NCBI RefSeq, 7 baculoviruses; 4 used in analysis):
```bash
cd ~/ITMO/omics/data/viruses
for acc in NC_001623.1 NC_001962.1 NC_002816.1 NC_002654.1 NC_009011.2 NC_002169.1 NC_001973.1; do
    efetch -db nuccore -id "$acc" -format fasta > "${acc}.fasta"
done
```

Host genomes (RefSeq chromosome-level assemblies):
```bash
cd ~/ITMO/omics/data/hosts
for acc in GCF_014905235.1 GCF_023101765.2 GCF_030705265.1 GCF_033807575.1; do
    ~/ITMO/omics/tools/datasets download genome accession "$acc" \
        --filename "${acc}.zip" --include genome --no-progressbar
    unzip -o -q "${acc}.zip" -d "${acc}"
    fna=$(find "${acc}" -name "*.fna" | head -n 1)
    mv "$fna" "${acc}.fasta"
done
rm -rf GCF_*.zip GCF_*/
```

Dfam consensus library (15.8 MB compressed; must be decompressed for RepeatMasker):
```bash
mkdir -p ~/ITMO/omics/tools/dfam && cd ~/ITMO/omics/tools/dfam
wget https://www.dfam.org/releases/current/families/Dfam-RepeatMasker.lib.gz
gzip -t Dfam-RepeatMasker.lib.gz && echo "VALID"
gzip -dk Dfam-RepeatMasker.lib.gz
```

Reference pipeline (from Science paper):
```bash
cd ~/ITMO/omics
git clone https://github.com/FelWint3r/BC200-MCV-analysis.git
```

## 2. TE discovery in baculoviruses

Concatenate target viruses and index:
```bash
cd ~/ITMO/omics/data/viruses
cat NC_001962.1.fasta NC_009011.2.fasta NC_002654.1.fasta NC_002816.1.fasta > baculo4.fasta
samtools faidx baculo4.fasta
```

De novo RepeatModeler (v2.0.8: no -engine flag):
```bash
mkdir -p rm_work && cd rm_work
BuildDatabase -name baculo4_db ../baculo4.fasta
RepeatModeler -database baculo4_db -threads 4 > repeatmodeler.log 2>&1
grep -c '^>' baculo4_db-families.fa
```

First RepeatMasker run WITH -nolow (wrong; kept for comparison):
```bash
cd ~/ITMO/omics/data/viruses && mkdir -p rm_known
# -nolow skips simple-repeat masking -> Dfam over-annotates tandem repeats as TEs
# -bed NOT available in RepeatMasker 4.2.3; use -gff + RM2Bed.py
for f in NC_001962.1.fasta NC_009011.2.fasta NC_002654.1.fasta NC_002816.1.fasta; do
    RepeatMasker -lib ~/ITMO/omics/tools/dfam/Dfam-RepeatMasker.lib \
        -engine ncbi -no_is -nolow -norna -gff -dir rm_known "$f"
done
cd rm_known
for f in *.fasta.out; do
    python3 ~/miniconda3/envs/metagenome_project/bin/RM2Bed.py "$f"
done
```

Second RepeatMasker run WITHOUT -nolow (correct):
```bash
cd ~/ITMO/omics/data/viruses && mkdir -p rm_correct
for f in NC_001962.1.fasta NC_009011.2.fasta NC_002654.1.fasta NC_002816.1.fasta; do
    RepeatMasker -lib ~/ITMO/omics/tools/dfam/Dfam-RepeatMasker.lib \
        -engine ncbi -no_is -norna -gff -dir rm_correct "$f"
done
cd rm_correct
for f in *.fasta.out; do
    python3 ~/miniconda3/envs/metagenome_project/bin/RM2Bed.py "$f"
done
```

Filter interspersed TE loci and extract >=100 bp:
```bash
cd ~/ITMO/omics/data/viruses/rm_correct
> all_TE_filtered.bed
for f in *_rm.bed; do
    base=$(basename "$f" _rm.bed)
    awk -v OFS='\t' -v name="$base" \
        '$7 != "Simple_repeat" && $7 != "Low_complexity" && $7 != "Satellite" \
         && $7 != "tRNA" && $7 != "rRNA" && $7 != "snRNA" {
            print $1, $2, $3, name"|"$4, $5, $6, $7, $8
        }' "$f" >> all_TE_filtered.bed
done
awk -F'\t' '$3-$2 >= 100' all_TE_filtered.bed > TE_long.bed
bedtools getfasta -fi ../baculo4.fasta -bed TE_long.bed -s -name -fo TE_long.fasta
```

Tandem-repeat detector (script `check_tandem.py`): counts max k-mer fraction for k=3..6 (skips k=1,2 to avoid homopolymer false positives); verdict TANDEM if fraction >= 0.30.

## 3. BLAST against host genomes

Build DBs:
```bash
cd ~/ITMO/omics/data/hosts && mkdir -p blastdb
for acc in GCF_014905235.1 GCF_023101765.2 GCF_030705265.1 GCF_033807575.1; do
    makeblastdb -in "${acc}.fasta" -dbtype nucl -out "blastdb/${acc}" -title "${acc}" -parse_seqids
done
```

Strict BLAST (megablast, 1e-5):
```bash
cd ~/ITMO/omics/data/viruses && mkdir -p blast_hosts
for nsq in ~/ITMO/omics/data/hosts/blastdb/GCF_*.nsq; do
    dbname=$(basename "$nsq" .nsq)
    blastn -query viral_TE_copies.fasta -db "$HOME/ITMO/omics/data/hosts/blastdb/${dbname}" \
        -out "blast_hosts/${dbname}.tsv" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
        -evalue 1e-5 -max_target_seqs 50 -num_threads 8
done
```

Sensitive BLAST (blastn, word_size 7, e-value 10):
```bash
cd ~/ITMO/omics/data/viruses && mkdir -p blast_hosts_sensitive
for nsq in ~/ITMO/omics/data/hosts/blastdb/GCF_*.nsq; do
    dbname=$(basename "$nsq" .nsq)
    blastn -query viral_TE_copies.fasta -db "$HOME/ITMO/omics/data/hosts/blastdb/${dbname}" \
        -out "blast_hosts_sensitive/${dbname}.tsv" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
        -task blastn -word_size 7 -evalue 10 -dust no \
        -max_target_seqs 100 -num_threads 8
done
```

Self-BLAST control (positive control for host DB): 500 bp from B. mori chr1 vs B. mori DB -> expect ~100% identity, e-value 0.

Locus-specificity filter (`locus_filter.py`): counts unique 10-kb windows per (query, host); <=10 loci = LOCUS-SPECIFIC (candidate HT); 11-50 = moderate; >50 = REPETITIVE.

## 4. Phylogenetics (ML + NJ)

Build dataset (6 viral + 7 Dfam consensus):
```bash
cd ~/ITMO/omics/data/viruses && mkdir -p phylo
zcat ~/ITMO/omics/tools/dfam/Dfam-RepeatMasker.lib.gz \
  | awk 'BEGIN{RS=">"} /hAT-Tip100|Chlamys|Academ-1|L1-Tx1|Gypsy|TcMar/ {print ">"$0}' \
  | head -100 > phylo/Dfam_refs.fasta
```

Clean viral FASTA with unique IDs (`unique_ids.py`): IQ-TREE requires unique names; append coordinates to disambiguate the two Chlamys fragments.
```bash
python3 unique_ids.py
cat phylo/viral_TE_clean.fasta phylo/Dfam_refs.fasta > phylo/combined.fasta
```

MAFFT alignment:
```bash
mafft --auto phylo/combined.fasta > phylo/combined_aligned.fasta 2> phylo/mafft.log
```

IQ-TREE ML (ModelFinder + UFBoot 1000). Binary is `iqtree` (v3.x); `iqtree2` not present in this env:
```bash
rm -f phylo/combined_iqtree.*
iqtree -s phylo/combined_aligned.fasta -m MFP -bb 1000 -nt 4 -pre phylo/combined_iqtree -quiet
grep -E "Best-fit model|Log-likelihood|BIC" phylo/combined_iqtree.iqtree
```

NJ tree (Bio.Phylo, identity distance) -- `nj_tree.py`.

## 5. Visualization and classifier

Distance heatmap + MDS -- `mds_heatmap.py` (seaborn heatmap, sklearn MDS):
```bash
python3 mds_heatmap.py
```

Random Forest classifier on 3-mer composition -- `classifier.py`. Training: 7 Dfam TE (positive) + 8 synthetic tandems (negative). Vocabulary built ONLY on train (105 3-mers); test uses same vocabulary (fixes "X has 62 features" error). Output: CV accuracy + per-candidate P(TE).

R tree visualization -- `plot_tree.R` (requires `ape`):
```bash
conda install -c conda-forge r-ape -y
Rscript plot_tree.R
```

## 6. Outputs

```
├── REPORT.md
├── COMMANDS.md
└── data/viruses/
    ├── rm_correct/
    │   ├── all_TE_filtered.bed
    │   ├── TE_long.bed
    │   └── TE_long.fasta
    ├── blast_hosts_sensitive/
    └── phylo/
        ├── combined_aligned.fasta
        ├── combined_iqtree.treefile
        ├── combined_tree_NJ.nwk
        ├── heatmap_distances.png
        ├── mds_plot.png
        ├── tree_ML.pdf
        └── tree_NJ.pdf
```
