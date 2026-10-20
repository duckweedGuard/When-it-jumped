library(ape)
ml <- read.tree("phylo/combined_iqtree.treefile")
nj <- read.tree("phylo/combined_tree_NJ.nwk")
cols <- ifelse(grepl("^viral_", ml$tip.label), "red", "steelblue")

pdf("phylo/tree_ML.pdf", width=14, height=10)
plot(ml, cex=0.7, main="ML tree: viral TE candidates (red) + Dfam consensus (blue)")
tiplabels(cex=0.5, bg=cols, frame="none")
dev.off()

pdf("phylo/tree_NJ.pdf", width=14, height=10)
plot(nj, cex=0.7, main="NJ tree: viral TE candidates + Dfam consensus")
dev.off()
cat("Trees plotted\n")
