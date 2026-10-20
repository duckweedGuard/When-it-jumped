# Figure 1 — ML phylogeny (Part I): viral TE candidates + Dfam consensus
# Output: figures/fig1_ml_tree.pdf, figures/fig1_ml_tree.png

suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2)
})

tree_file <- "results/phylogeny/combined_iqtree.treefile"
out_pdf   <- "figures/fig1_ml_tree.pdf"
out_png   <- "figures/fig1_ml_tree.png"

tree <- read.tree(tree_file)
tree <- ladderize(tree, right = FALSE)

relabel <- function(x) {
  if (grepl("^viral_", x)) {
    if (grepl("hAT", x))                 return("hAT \u00b7 BmNPV")
    if (grepl("Chlamys-1_bBra_3712", x)) return("Chlamys a \u00b7 HearNPV")
    if (grepl("Chlamys-1_bBra_4094", x)) return("Chlamys b \u00b7 HearNPV")
    if (grepl("Academ", x))              return("Academ \u00b7 CpGV")
    if (grepl("LINE_M", x))              return("LINE \u00b7 CpGV")
    if (grepl("Acol_gypsy_Ele46", x))    return("Gypsy \u00b7 CpGV")
    return(x)
  }
  name <- sub("#.*", "", x)
  fam  <- sub(".*#", "", x); fam <- sub(" .*", "", fam)
  fam  <- sub("^(DNA|LTR|LINE|SINE|RC)/", "", fam)
  paste0(name, " \u00b7 ", fam)
}

origin <- ifelse(grepl("^viral_", tree$tip.label), "Viral candidate", "Dfam consensus")
tree$tip.label <- sapply(tree$tip.label, relabel)
tipdata <- data.frame(label = tree$tip.label, origin = origin)

p <- ggtree(tree, size = 0.5) %<+% tipdata +
  geom_tippoint(aes(color = origin), size = 3) +
  geom_tiplab(aes(label = label, color = origin), size = 3.4, offset = 0.02) +
  geom_nodelab(aes(subset = !isTip & label != "" &
                     suppressWarnings(as.numeric(label)) >= 50, label = label),
               size = 2.8, color = "grey30", nudge_x = -0.01) +
  scale_color_manual(values = c("Viral candidate" = "#D7263D",
                                "Dfam consensus"  = "#1B98E0"), name = NULL) +
  theme_tree2(base_size = 11) +
  theme(legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(10, 60, 10, 10)) +
  labs(title    = "ML phylogeny of viral TE candidates and Dfam consensus sequences",
       subtitle = "IQ-TREE 3.1.4 \u00b7 HKY+F+G4 \u00b7 UFBoot 1000 \u00b7 13 taxa \u00b7 2,222 sites",
       x = "substitutions per site")

ggsave(out_pdf, p, width = 10, height = 7)
ggsave(out_png, p, width = 10, height = 7, dpi = 300)
message("Wrote: ", out_pdf, " and ", out_png)
