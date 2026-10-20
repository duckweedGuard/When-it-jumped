# Figure 2 — Pairwise identity-distance heatmap
# Output: figures/fig2_heatmap.pdf, figures/fig2_heatmap.png

suppressPackageStartupMessages({ library(ape); library(pheatmap) })

aln_file <- "results/phylogeny/combined_aligned.fasta"
out_pdf  <- "figures/fig2_heatmap.pdf"
out_png  <- "figures/fig2_heatmap.png"

aln <- read.dna(aln_file, format = "fasta")
dm  <- dist.dna(aln, model = "N", pairwise.deletion = TRUE)
m   <- as.matrix(dm)

relabel <- function(x) {
  if (grepl("^viral_", x)) {
    if (grepl("hAT", x))                 return("hAT (BmNPV)")
    if (grepl("Chlamys-1_bBra_3712", x)) return("Chlamys a (HearNPV)")
    if (grepl("Chlamys-1_bBra_4094", x)) return("Chlamys b (HearNPV)")
    if (grepl("Academ", x))              return("Academ (CpGV)")
    if (grepl("LINE_M", x))              return("LINE (CpGV)")
    if (grepl("Acol_gypsy_Ele46", x))    return("Gypsy (CpGV)")
    return(x)
  }
  name <- sub("#.*", "", x); fam <- sub(".*#", "", x); fam <- sub(" .*", "", fam)
  fam  <- sub("^(DNA|LTR|LINE|SINE|RC)/", "", fam)
  paste0(name, " (", fam, ")")
}

nm <- sapply(rownames(m), relabel)
rownames(m) <- colnames(m) <- nm

ann <- data.frame(
  Origin = ifelse(grepl("\\(BmNPV\\)|\\(HearNPV\\)|\\(CpGV\\)", nm), "Viral", "Dfam"),
  row.names = nm
)
ann_colors <- list(Origin = c(Viral = "#D7263D", Dfam = "#1B98E0"))

common_args <- list(
  mat = m,
  clustering_distance_rows = as.dist(m),
  clustering_distance_cols = as.dist(m),
  clustering_method = "average",
  annotation_row = ann, annotation_col = ann,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("white", "orange", "firebrick"))(50),
  border_color = "grey80",
  fontsize = 9,
  main = "Pairwise identity distance: viral TE candidates vs Dfam consensus",
  width = 8, height = 7
)

do.call(pheatmap, c(common_args, list(filename = out_pdf)))
do.call(pheatmap, c(common_args, list(filename = out_png)))
message("Wrote: ", out_pdf, " and ", out_png)
