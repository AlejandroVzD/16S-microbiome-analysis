# =============================================================================
# 01_upstream_dada2.R
# From raw reads to an ASV table with DADA2
# -----------------------------------------------------------------------------
# Part 1 of a reproducible 16S rRNA pipeline. Takes paired-end Illumina reads
# and produces an amplicon sequence variant (ASV) table plus taxonomy, ready
# for downstream diversity and differential-abundance analysis.
#
# Demo data: DADA2 tutorial set (mouse gut, 16S V4, 2x250 MiSeq) -- a small,
# canonical, fully public dataset chosen so the whole workflow runs in minutes
# and is 100% reproducible. Point `path` at your own demultiplexed, primer-free
# fastq files to run the exact same pipeline on real data.
#
# Reference: Callahan et al. (2016) Nature Methods 13:581-583 (DADA2).
# =============================================================================

## ---- 0. Setup --------------------------------------------------------------
# install.packages("BiocManager")
# BiocManager::install(c("dada2", "phyloseq", "Biostrings"))
library(dada2)        # v1.26+
library(phyloseq)
library(Biostrings)

set.seed(100)         # reproducible error learning and sample inference

# Folder with demultiplexed, primer-free fastq files.
# Demo data (run once): download and unzip
#   https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip
# then set `path` to the resulting MiSeq_SOP folder.
path <- "data/MiSeq_SOP"
list.files(path)

## ---- 1. Read files & parse sample names ------------------------------------
# Forward = R1, reverse = R2. Adjust the pattern to your own file naming.
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## ---- 2. Inspect read quality -----------------------------------------------
# Use these plots to choose the truncation lengths in step 3.
plotQualityProfile(fnFs[1:2])   # forward
plotQualityProfile(fnRs[1:2])   # reverse (quality usually drops earlier)

## ---- 3. Filter & trim ------------------------------------------------------
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# truncLen: read it off the quality profiles (values below fit the demo V4 data).
# maxEE: max expected errors -- the main quality filter. Relax if you lose too
# many reads. Keep enough overlap between F and R for merging (>= ~20 bp).
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(240, 160),
                     maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)

## ---- 4. Learn error rates --------------------------------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)   # observed points should track the fitted line

## ---- 5. Sample inference (denoising) ---------------------------------------
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

## ---- 6. Merge pairs & build the ASV table ----------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
seqtab  <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))  # amplicon length distribution

## ---- 7. Remove chimeras ----------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
sum(seqtab.nochim) / sum(seqtab)    # fraction of reads kept (usually high)

## ---- 8. Track reads through the pipeline (QC sanity check) -----------------
getN  <- function(x) sum(getUniques(x))
track <- cbind(out,
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample.names
head(track)   # no single step should drop the majority of reads

## ---- 9. Assign taxonomy (SILVA v138.1) -------------------------------------
# Download the training set once (place it in data/):
#   https://zenodo.org/record/4587955  ->  silva_nr99_v138.1_train_set.fa.gz
taxa <- assignTaxonomy(seqtab.nochim,
                       "data/silva_nr99_v138.1_train_set.fa.gz",
                       multithread = TRUE)
# taxa <- addSpecies(taxa, "data/silva_species_assignment_v138.1.fa.gz")

## ---- 10. Hand off to phyloseq ----------------------------------------------
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(taxa))

# Store full ASV sequences in refseq() and rename to short IDs (ASV1, ASV2, ...)
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Save so 02_downstream_analysis.R can pick up from here.
saveRDS(ps, "data/ps_dada2_demo.rds")

# =============================================================================
# Output: data/ps_dada2_demo.rds  (phyloseq object: ASV table + taxonomy)
# Next:   02_downstream_analysis.R (diversity, PERMANOVA, LEfSe, DESeq2)
# =============================================================================
