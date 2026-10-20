from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

aln = AlignIO.read("phylo/combined_aligned.fasta", "fasta")
calc = DistanceCalculator("identity")
dm = calc.get_distance(aln)
tree = DistanceTreeConstructor().nj(dm)
Phylo.write(tree, "phylo/combined_tree_NJ.nwk", "newick")
print(f"NJ tree: {len(aln)} taxa, written to phylo/combined_tree_NJ.nwk")
