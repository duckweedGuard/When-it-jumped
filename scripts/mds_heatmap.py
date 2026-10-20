import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator

aln = AlignIO.read("phylo/combined_aligned.fasta", "fasta")
calc = DistanceCalculator("identity")
dm = calc.get_distance(aln)

names = [r.id for r in aln]
n = len(names)
mat = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        mat[i, j] = dm[i, j]

# --- Heatmap ---
plt.figure(figsize=(12, 10))
sns.heatmap(mat, xticklabels=names, yticklabels=names, cmap="viridis", square=True)
plt.title("Distance matrix: viral TE candidates + Dfam consensus sequences")
plt.xticks(rotation=90, fontsize=7)
plt.yticks(fontsize=7)
plt.tight_layout()
plt.savefig("phylo/heatmap_distances.png", dpi=150)
plt.close()

# --- MDS (PCoA) ---
from sklearn.manifold import MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42, normalized_stress="auto")
coords = mds.fit_transform(mat)

plt.figure(figsize=(12, 8))
colors = ["red" if x.startswith("viral_") else "steelblue" for x in names]
plt.scatter(coords[:, 0], coords[:, 1], c=colors, s=80)
for i, name in enumerate(names):
    plt.annotate(name, (coords[i, 0], coords[i, 1]), fontsize=7, alpha=0.85)
plt.title("MDS of TE candidates vs Dfam consensus (red = viral, blue = Dfam)")
plt.xlabel("MDS1")
plt.ylabel("MDS2")
plt.tight_layout()
plt.savefig("phylo/mds_plot.png", dpi=150)
plt.close()

print("Saved: phylo/heatmap_distances.png, phylo/mds_plot.png")
