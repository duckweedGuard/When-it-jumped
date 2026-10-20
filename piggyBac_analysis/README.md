# piggyBac analysis

# 1. Data downloading
```bash
wget -O viruses/NC_001623.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001623.1&rettype=fasta&retmode=text"
wget -O viruses/NC_001962.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001962.1&rettype=fasta&retmode=text"
wget -O viruses/NC_002816.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_002816.1&rettype=fasta&retmode=text"
wget -O viruses/NC_002169.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_002169.1&rettype=fasta&retmode=text"
wget -O viruses/NC_001973.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001973.1&rettype=fasta&retmode=text"
#hosts
wget -O hosts/GCF_014905235.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=GCF_014905235.1&rettype=fasta&retmode=text"
wget -O hosts/GCF_023101765.2.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=GCF_023101765.2&rettype=fasta&retmode=text"
wget -O hosts/GCF_030705265.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=GCF_030705265.1&rettype=fasta&retmode=text"
wget -O hosts/GCF_033807575.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=GCF_033807575.1&rettype=fasta&retmode=text"
#piggyBac(IFP2)
wget -O transposon_queries/piggyBac_clean.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=J04364.2&rettype=fasta&retmode=text"
```
# 2. BLAST db preparation&search piggyBac
```bash
for fasta in hosts/*.fasta viruses/*.fasta; do
    makeblastdb -in "$fasta" -dbtype nucl -out "${fasta%.fasta}_db"
done
```

```bash
for db in hosts/*_db.nsq viruses/*_db.nsq; do
    db_name=${db%.nsq}
    blastn -query transposon_queries/piggyBac_clean.fasta -db "$db_name" \
           -out "blast_results/$(basename $db_name)_hits.tsv" \
           -outfmt "6 qseqid sseqid pident length evalue bitscore sstart send sstrand" \
           -evalue 1e-5 -num_threads 8
done
```

# 3. Run samtools

```bash
for f in hosts/*.fasta viruses/*.fasta; do samtools faidx "$f"; done
```
## extract_hits.sh
```bash
cat << 'EOF' > extract_hits.sh
#!/bin/bash
for blast_file in blast_results/*.tsv; do
    [ -s "$blast_file" ] || continue
    genome=$(basename "$blast_file" _hits.tsv | sed 's/.*_vs_//')
    [[ "$genome" == GCF* ]] && genome_path="hosts/${genome}.fasta" || genome_path="viruses/${genome}.fasta"
    
    while IFS=$'\t' read -r qseqid sseqid pident length evalue bitscore sstart send sstrand; do
        sseqid_clean=$(echo "$sseqid" | sed 's/^ref|//; s/|$//')
        hit_name="${qseqid}_${genome}_${sstart}-${send}_${sstrand}"
        
        if [ "$sstart" -gt "$send" ]; then
            start=$(( send > 100 ? send - 100 : 1 )); end=$(( sstart + 100 ))
        else
            start=$(( sstart > 100 ? sstart - 100 : 1 )); end=$(( send + 100 ))
        fi
        
        samtools faidx "$genome_path" "${sseqid_clean}:${start}-${end}" 2>/dev/null | \
        sed "1s/.*/>${hit_name}/" >> "extracted_hits/piggyBac_hits.fasta"
    done < "$blast_file"
done
EOF
```

```bash
chmod +x extract_hits.sh
./extract_hits.sh
```
# 4. MAFFT + trimAl

```bash
cat transposon_queries/piggyBac_clean.fasta extracted_hits/piggyBac_hits.fasta > combined_piggyBac.fasta

mafft --adjustdirection --auto combined_piggyBac.fasta > alignment_piggyBac.fasta

trimal -in alignment_piggyBac.fasta -out trimmed_piggyBac.fasta -automated1
```
# 5. IQ-TREE. ModelFinder + UFBoot 1000

```bash
iqtree2 -s trimmed_piggyBac.fasta -m MFP -B 1000 -alrt 1000 -T 16 --prefix tree_piggyBac
```

# 6. Visualisation at iTOL.embl.de
