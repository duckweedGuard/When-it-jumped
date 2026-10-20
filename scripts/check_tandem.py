from Bio import SeqIO
from collections import Counter

print(f"{'LEN':>5}  {'FRAC':>5}  {'KMER':>8}  {'VERDICT':>10}  ID")
for rec in SeqIO.parse("rm_correct/TE_long.fasta", "fasta"):
    s = str(rec.seq).upper()
    n = len(s)
    best_frac, best_kmer = 0.0, None
    for k in range(3, 7):
        if n < 3 * k:
            continue
        cnt = Counter(s[i:i+k] for i in range(n - k + 1))
        kmer, c = cnt.most_common(1)[0]
        frac = c * k / n
        if frac > best_frac:
            best_frac, best_kmer = frac, kmer
    verdict = "TANDEM" if best_frac >= 0.30 else "OK"
    print(f"{n:>5}  {best_frac:>5.2f}  {str(best_kmer):>8}  {verdict:>10}  {rec.id}")
