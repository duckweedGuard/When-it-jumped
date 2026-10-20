from pyfaidx import Fasta
from Bio import SeqIO

hosts = {
    "GCF_014905235.1": "/home/user16/ITMO/omics/data/hosts/GCF_014905235.1.fasta",
    "GCF_023101765.2": "/home/user16/ITMO/omics/data/hosts/GCF_023101765.2.fasta",
    "GCF_030705265.1": "/home/user16/ITMO/omics/data/hosts/GCF_030705265.1.fasta",
    "GCF_033807575.1": "/home/user16/ITMO/omics/data/hosts/GCF_033807575.1.fasta",
}
candidate = "Acol_gypsy_Ele46::NC_002816.1:20385-20583(+)"
out_fa = "phylo/Ele46_dataset.fasta"

def clean_id(sseq):
    # 'ref|NC_051358.1|' -> 'NC_051358.1'
    s = sseq
    if s.startswith("ref|"):
        s = s[len("ref|"):]
    s = s.rstrip("|")
    return s

seqs = []
for acc, fa in hosts.items():
    genome = Fasta(fa, as_raw=False)
    fn = f"blast_hosts_sensitive/{acc}.tsv"
    hits = []
    for line in open(fn):
        c = line.rstrip("\n").split("\t")
        if c[0] != candidate:
            continue
        pident = float(c[2]); length = int(c[3])
        if pident < 70 or length < 100:
            continue
        sseq = clean_id(c[1])
        sstart = int(c[8]); send = int(c[9]); bs = float(c[11])
        start, end = sorted([sstart, send])
        hits.append((bs, sseq, start, end, pident, length))
    hits.sort(reverse=True)
    seen, taken = set(), 0
    for bs, sseq, start, end, pid, ln in hits:
        locus = (sseq, start // 10000)
        if locus in seen:
            continue
        seen.add(locus)
        try:
            seq = str(genome[sseq][start-1:end])
        except KeyError:
            continue
        name = f"{acc}|{sseq}:{start}-{end}|{pid:.0f}pct|{ln}bp"
        seqs.append((name, seq))
        taken += 1
        if taken >= 5:
            break

for rec in SeqIO.parse("viral_TE_copies.fasta", "fasta"):
    if rec.id == candidate:
        seqs.append((f"VIRAL_CpGV|{rec.id}", str(rec.seq)))
        break

with open(out_fa, "w") as f:
    for name, seq in seqs:
        f.write(f">{name}\n{seq}\n")
print(f"Wrote {len(seqs)} sequences to {out_fa}")
