from Bio import SeqIO
from collections import Counter
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

def kmer_counts(seq, k=3):
    s = str(seq).upper()
    return Counter(s[i:i+k] for i in range(len(s)-k+1))

def build_vocab(seqs, k=3):
    vocab = set()
    for s in seqs:
        s = str(s).upper()
        for i in range(len(s)-k+1):
            vocab.add(s[i:i+k])
    return sorted(vocab)

def vectorize(seqs, vocab, k=3):
    X = np.zeros((len(seqs), len(vocab)))
    idx = {km: j for j, km in enumerate(vocab)}
    for i, s in enumerate(seqs):
        cnt = kmer_counts(s, k)
        total = sum(cnt.values()) or 1
        for km, c in cnt.items():
            if km in idx:
                X[i, idx[km]] = c / total
    return X

# Training set: Dfam consensus (TE=1) + synthetic tandems (TANDEM=0)
dfam = [r.seq for r in SeqIO.parse("phylo/Dfam_refs.fasta", "fasta")]
tandems = [(m * 100)[:300] for m in ["GAT","TCAT","ATT","AT","A","T","GATGAT","TAT"]]
train_seqs = dfam + tandems
y_train = np.array([1]*len(dfam) + [0]*len(tandems))

vocab = build_vocab(train_seqs, k=3)
X_train = vectorize(train_seqs, vocab, k=3)

clf = RandomForestClassifier(n_estimators=200, random_state=42)
cv = cross_val_score(clf, X_train, y_train, cv=3)
clf.fit(X_train, y_train)

# Test set: viral candidates
viral = list(SeqIO.parse("phylo/viral_TE_clean.fasta", "fasta"))
X_test = vectorize([r.seq for r in viral], vocab, k=3)

pred = clf.predict(X_test)
prob = clf.predict_proba(X_test)

print(f"Training: {len(dfam)} Dfam TE + {len(tandems)} synthetic tandems")
print(f"Feature vocabulary: {len(vocab)} 3-mers")
print(f"CV accuracy: {cv.mean():.3f} +/- {cv.std():.3f}\n")
print(f"{'ID':<58} {'pred':>8} {'P(TE)':>8}")
for rec, p, pr in zip(viral, pred, prob):
    label = "TE" if p == 1 else "TANDEM"
    print(f"{rec.id:<58} {label:>8} {pr[1]:>8.3f}")
