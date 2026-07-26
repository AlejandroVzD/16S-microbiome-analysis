# =====================================================================
#  16S rRNA MICROBIOME ANALYSIS PIPELINE
#  ---------------------------------------------------------------
#  Reproducible diversity analysis on a public 16S dataset
#  (GlobalPatterns, bundled with phyloseq). Runs in seconds.
#
#  Outputs:
#    - pcoa_public.png    Bray-Curtis PCoA + PERMANOVA
#    - alpha_public.png   Alpha diversity (Observed, Shannon) + Kruskal-Wallis
#
#  Dependencies: phyloseq, vegan, ggplot2
# =====================================================================

library(phyloseq)
library(vegan)
library(ggplot2)

# Journal-style palette (NPG), defined inline so no extra package is needed
npg <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F",
         "#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")

# ---------------------------------------------------------------
# Load public data and convert counts to relative abundance
# (deterministic; avoids random rarefaction)
# ---------------------------------------------------------------
data(GlobalPatterns)
ps     <- GlobalPatterns
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# ---------------------------------------------------------------
# FIGURE 1 — Beta diversity: Bray-Curtis PCoA + PERMANOVA
# ---------------------------------------------------------------
set.seed(123)
ord <- ordinate(ps_rel, method = "PCoA", distance = "bray")

dm  <- phyloseq::distance(ps_rel, method = "bray")
grp <- data.frame(sample_data(ps_rel))
per <- adonis2(dm ~ SampleType, data = grp, permutations = 999)

subt <- paste0("Bray-Curtis  |  PERMANOVA  R\u00b2 = ", round(per$R2[1], 3),
               ",  p = ", format.pval(per$`Pr(>F)`[1], digits = 2))

p1 <- plot_ordination(ps_rel, ord, color = "SampleType") +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = npg) +
  labs(title = "Microbiome Ordination by Sample Type",
       subtitle = subt, color = "Sample type") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave("pcoa_public.png", p1, width = 12.8, height = 7.69, dpi = 200, bg = "white")

# ---------------------------------------------------------------
# FIGURE 2 — Alpha diversity: Observed + Shannon, Kruskal-Wallis
# ---------------------------------------------------------------
rich <- estimate_richness(ps, measures = c("Observed", "Shannon"))
rich$SampleType <- sample_data(ps)$SampleType

long <- rbind(
  data.frame(SampleType = rich$SampleType, Metric = "Observed (ASVs)",     Value = rich$Observed),
  data.frame(SampleType = rich$SampleType, Metric = "Shannon (diversity)", Value = rich$Shannon)
)

kw_obs <- kruskal.test(Observed ~ SampleType, data = rich)
kw_sha <- kruskal.test(Shannon  ~ SampleType, data = rich)
cap <- paste0("Kruskal-Wallis  |  Observed p = ", format.pval(kw_obs$p.value, digits = 2),
              "   \u00b7   Shannon p = ", format.pval(kw_sha$p.value, digits = 2))

p2 <- ggplot(long, aes(SampleType, Value, fill = SampleType)) +
  geom_violin(alpha = 0.55, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = npg) +
  labs(title = "Alpha Diversity by Sample Type", caption = cap, x = NULL, y = "Value") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1),
        strip.background = element_rect(fill = "#3C5488"),
        strip.text = element_text(colour = "white", face = "bold"))

ggsave("alpha_public.png", p2, width = 12.8, height = 7.69, dpi = 200, bg = "white")

cat("\nDone. Created pcoa_public.png and alpha_public.png\n")
