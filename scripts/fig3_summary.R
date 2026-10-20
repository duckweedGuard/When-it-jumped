# Figure 3 — Exploratory annotation summary (-nolow vs no -nolow)
# Output: figures/fig3_summary.pdf, figures/fig3_summary.png

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})

out_pdf <- "figures/fig3_summary.pdf"
out_png <- "figures/fig3_summary.png"

df_a <- data.frame(
  run   = factor(rep(c("With -nolow", "Without -nolow"), each = 2),
                 levels = c("With -nolow", "Without -nolow")),
  class = rep(c("All annotated loci", "Interspersed TE only"), 2),
  n     = c(28, 28, 148, 23)
)

p_a <- ggplot(df_a, aes(run, n, fill = class)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(0.7), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("All annotated loci" = "grey70",
                               "Interspersed TE only" = "#D7263D")) +
  labs(x = NULL, y = "Number of loci", fill = NULL,
       title = "A. RepeatMasker vs Dfam") +
  theme_classic(base_size = 11) + theme(legend.position = "top")

df_b <- data.frame(
  status = factor(c("Tandem repeat", "Retained"),
                  levels = c("Tandem repeat", "Retained")),
  n = c(4, 2)
)

p_b <- ggplot(df_b, aes(status, n, fill = status)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Tandem repeat" = "grey50", "Retained" = "#1B98E0")) +
  labs(x = NULL, y = "Candidates \u2265100 bp", title = "B. Long candidates (n = 6)") +
  theme_classic(base_size = 11)

p <- (p_a | p_b) +
  plot_annotation(
    title    = "Exploratory annotation: -nolow inflates TE calls by ~80%",
    subtitle = "4 baculovirus genomes (BmNPV, SfMNPV, HearNPV, CpGV)",
    theme    = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(out_pdf, p, width = 11, height = 4.5)
ggsave(out_png, p, width = 11, height = 4.5, dpi = 300)
message("Wrote: ", out_pdf, " and ", out_png)
