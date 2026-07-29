# =============================================================================
# 02_downstream_analysis.R
# From an abundance table to diversity, statistics and biomarkers
# -----------------------------------------------------------------------------
# Part 2 of a reproducible 16S rRNA pipeline. Part 1 (01_upstream_dada2.R)
# demonstrates the raw-reads -> ASV-table step on a small tutorial dataset.
# This part runs the full community analysis on a real two-group study.
#
# Data: kostic_crc, bundled with the microbiomeMarker package -- 16S gut
# microbiome, Healthy vs Tumor (colorectal carcinoma).
# Reference: Kostic et al. (2012) Genome Research 22:292-298.
#
# Workflow: QC/filtering -> taxonomic composition -> alpha diversity
#           -> beta diversity (PERMANOVA + betadisper) -> DESeq2 -> LEfSe.
# Output:   01_composition.png ... 05_lefse.png (written to the working directory)
# =============================================================================

library(phyloseq); library(vegan); library(ggplot2)
library(DESeq2);   library(microbiomeMarker)

theme_set(theme_classic(base_size = 14))
npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
         "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")


# ---- 1. LOAD + CHECK GROUPS (run first and confirm the output) -------------
data(kostic_crc)
ps <- kostic_crc
cat("Sample variables:\n"); print(sample_variables(ps))
cat("\nDIAGNOSIS values:\n"); print(table(sample_data(ps)$DIAGNOSIS, useNA = "ifany"))


# ---- 2. QC / FILTERING -----------------------------------------------------
sample_data(ps)$Group <- factor(sample_data(ps)$DIAGNOSIS,
                                levels = c("Healthy", "Tumor"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
# Prevalence filter: keep taxa with >= 3 counts in at least 5% of samples.
keep <- genefilter_sample(ps, filterfun_sample(function(x) x >= 3),
                          A = 0.05 * nsamples(ps))
ps_f   <- prune_taxa(keep, ps)
ps_rel <- transform_sample_counts(ps_f, function(x) x / sum(x))
cat("Taxa after filtering:", ntaxa(ps_f), "| Samples:", nsamples(ps_f), "\n")


# ---- 3. TAXONOMIC COMPOSITION (top phyla, relative abundance) --------------
ps_phy  <- tax_glom(ps_rel, taxrank = "Phylum")
top_phy <- names(sort(taxa_sums(ps_phy), decreasing = TRUE))[1:8]
ps_top  <- prune_taxa(top_phy, ps_phy)
dfc <- psmelt(ps_top)
dfc$Phylum <- as.character(dfc$Phylum)

p_comp <- ggplot(dfc, aes(Sample, Abundance, fill = Phylum)) +
  geom_col(position = "fill", width = 1) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = npg) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Taxonomic Composition (top phyla)",
       x = "Samples", y = "Relative abundance") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "#3C5488"),
        strip.text = element_text(colour = "white", face = "bold"))
ggsave("01_composition.png", p_comp, width = 11, height = 6,
       dpi = 150, bg = "white")


# ---- 4. ALPHA DIVERSITY (points + boxplot + Kruskal-Wallis) ----------------
rich <- estimate_richness(ps_f, measures = c("Observed", "Shannon"))
rich$Group <- sample_data(ps_f)$Group
long <- rbind(
  data.frame(Group = rich$Group, Metric = "Observed", Value = rich$Observed),
  data.frame(Group = rich$Group, Metric = "Shannon",  Value = rich$Shannon))

p_alpha <- ggplot(long, aes(Group, Value, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.5) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = npg) +
  labs(title = "Alpha Diversity", x = NULL, y = "Value") +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#3C5488"),
        strip.text = element_text(colour = "white", face = "bold"))
ggsave("02_alpha_diversity.png", p_alpha, width = 9, height = 6,
       dpi = 150, bg = "white")
cat("Alpha KW  Observed p =", signif(kruskal.test(Observed ~ Group, rich)$p.value, 3),
    "| Shannon p =", signif(kruskal.test(Shannon ~ Group, rich)$p.value, 3), "\n")


# ---- 5. BETA DIVERSITY (Bray-Curtis PCoA + PERMANOVA + betadisper) ---------
set.seed(123)   # reproducible permutation tests (PERMANOVA / betadisper)
ord  <- ordinate(ps_rel, method = "PCoA", distance = "bray")
dm   <- phyloseq::distance(ps_rel, method = "bray")
meta <- data.frame(sample_data(ps_rel))
per  <- adonis2(dm ~ Group, data = meta, permutations = 999)
bd   <- betadisper(dm, meta$Group); bd_p <- permutest(bd)$tab$`Pr(>F)`[1]
sub  <- paste0("PERMANOVA R2 = ", round(per$R2[1], 3), ", p = ", per$`Pr(>F)`[1],
               "  |  betadisper p = ", round(bd_p, 3))

p_beta <- plot_ordination(ps_rel, ord, color = "Group") +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Group), type = "norm", linetype = 2) +
  scale_color_manual(values = npg) +
  labs(title = "Beta Diversity (Bray-Curtis PCoA)", subtitle = sub)
ggsave("03_beta_pcoa.png", p_beta, width = 9, height = 6.5,
       dpi = 150, bg = "white")


# ---- 6. DIFFERENTIAL ABUNDANCE (DESeq2, genus-labeled only) ----------------
dds <- phyloseq_to_deseq2(ps_f, ~ Group)
dds <- DESeq(dds, test = "Wald", fitType = "parametric", sfType = "poscounts")
res <- as.data.frame(results(dds))
res$Genus <- as.character(tax_table(ps_f)[rownames(res), "Genus"])
sig <- subset(res, padj < 0.05 & !is.na(padj) & !is.na(Genus) & Genus != "")
sig <- head(sig[order(-abs(sig$log2FoldChange)), ], 25)
if (nrow(sig) > 0) {
  sig$Genus <- make.unique(sig$Genus)
  sig$Genus <- factor(sig$Genus, levels = sig$Genus[order(sig$log2FoldChange)])
  sig$Dir <- ifelse(sig$log2FoldChange > 0, "Enriched in Tumor", "Enriched in Healthy")
  p_deseq <- ggplot(sig, aes(log2FoldChange, Genus, color = Dir)) +
    geom_segment(aes(x = 0, xend = log2FoldChange, yend = Genus), color = "grey70") +
    geom_point(size = 3.5) +
    scale_color_manual(values = c("Enriched in Healthy" = "#0072B5",
                                  "Enriched in Tumor"   = "#BC3C29")) +
    labs(title = "Differential Abundance (DESeq2)",
         x = "log2 Fold Change (Tumor vs Healthy)", y = NULL, color = NULL) +
    theme(legend.position = "top")
  ggsave("04_deseq2.png", p_deseq, width = 9, height = 8,
         dpi = 150, bg = "white")
  cat("DESeq2 genus-labeled significant taxa:", nrow(sig), "\n")
}


# ---- 7. BIOMARKERS (LEfSe, genus level, top by LDA) ------------------------
ps_g <- tax_glom(ps_f, taxrank = "Genus")
mm <- run_lefse(ps_g, group = "Group", taxa_rank = "Genus",
                norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2)
mt <- marker_table(mm)
cat("LEfSe markers (genus):", nrow(mt), "\n")
if (nrow(mt) > 0) {
  p_lefse <- plot_ef_bar(mm, label_level = 1) +
    scale_fill_manual(values = c("#BC3C29", "#0072B5")) +
    labs(title = "LEfSe Biomarkers (genus, LDA > 2)")
  ggsave("05_lefse.png", p_lefse, width = 9, height = 8,
         dpi = 150, bg = "white")
}

cat("\nDownstream analysis complete. Figures 01-05 saved.\n")
# sessionInfo()   # uncomment to record exact package versions for reproducibility
