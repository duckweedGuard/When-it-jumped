# Figure 4 — MDS of pairwise identity distances
# Output: figures/fig4_mds.pdf, figures/fig4_mds.png

suppressPackageStartupMessages({
  library(ape); library(ggplot2); library(ggrepel)
})

aln_file <- "results/phylogeny/combined_aligned.fasta"
out_pdf  <- "figures/fig4_mds.pdf"
out_png  <- "figures/fig4_mds.png"

aln <- read.dna(aln_file, format = "fasta")
dm  <- dist.dna(aln, model = "N", pairwise.deletion = TRUE)

mds <- cmdscale(dm, k = 2, eig = TRUE)
var_explained <- round(100 * mds$eig[1:2] / sum(mds$eig[mds$eig > 0]), 1)

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

df <- data.frame(
  label = sapply(rownames(as.matrix(dm)), relabel),
  MDS1 = mds$points[, 1],
  MDS2 = mds$points[, 2],
  stringsAsFactors = FALSE
)
df$origin <- ifelse(grepl("\\(BmNPV\\)|\\(HearNPV\\)|\\(CpGV\\)", df$label),
                    "Viral candidate", "Dfam consensus")

p <- ggplot(df, aes(MDS1, MDS2, color = origin)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 30, show.legend = FALSE) +
  scale_color_manual(values = c("Viral candidate" = "#D7263D",
                                "Dfam consensus"  = "#1B98E0"), name = NULL) +
  labs(title    = "MDS of pairwise identity distances",
       subtitle = "Classical MDS \u00b7 2,222 aligned sites",
       x = paste0("MDS1 (", var_explained[1], "%)"),
       y = paste0("MDS2 (", var_explained[2], "%)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(out_pdf, p, width = 8, height = 6)
ggsave(out_png, p, width = 8, height = 6, dpi = 300)
message("Wrote: ", out_pdf, " and ", out_png)
