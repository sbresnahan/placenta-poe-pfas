# =============================================================================
# Stage 02 | infer_ancestry.R
# Assign every trio to Chinese / Indian / Malay ancestry. PCA on LD-pruned
# autosomal SNPs in infants; infants with reported maternal ethnicity (the
# RNA-seq subset) train a 3-NN classifier on PC1/PC2 that labels the rest.
# Each infant's parents inherit the infant's label.
# In : TRIO_VCF_HG38, INFANT_META_CSV (SubjectID, maternal ethnicity),
#      infant .fam (first PLINK_PREFIXES entry), TRIO_PED
# Out: ANCESTRY_DIR/samples_<group>.txt (VCF sample IDs: infants, mothers, fathers)
#      ANCESTRY_DIR/ancestry_assignments.tsv; figures in FIG_DIR
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages({
  library(SNPRelate); library(data.table); library(ggplot2); library(class)
})

# ---- 1. Inputs ----
vcf_file   <- P("TRIO_VCF_HG38")
infant_fam <- file.path(P("TRIOS_DIR"), paste0(words("PLINK_PREFIXES")[1], ".fam"))
out_dir    <- P("ANCESTRY_DIR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Reported maternal ethnicity for the RNA-seq infants
meta <- fread(P("INFANT_META_CSV"), colClasses = "character")
setnames(meta, 1:2, c("sample.id", "maternal_ethnicity"))
meta <- meta[, .(sample.id, maternal_ethnicity)]

# PLINK writes VCF sample names as FID_IID
infant_ids <- fread(infant_fam, select = 1:2, header = FALSE, colClasses = "character",
                    col.names = c("sample.id", "iid"))
infant_ids[, sample.id.vcf := paste0(sample.id, "_", iid)]
meta <- merge(meta, infant_ids[, .(sample.id, sample.id.vcf)], by = "sample.id")

# Pedigree -> VCF IDs of each infant's parents (FID is the IID minus its role letter)
peds <- fread(P("TRIO_PED"), header = FALSE, colClasses = "character", select = 2:4,
              col.names = c("infant_id", "father_id", "mother_id"))
peds <- peds[!(father_id == "0" & mother_id == "0")]
to_vcf_id <- function(x) ifelse(x == "0", NA_character_, paste0(sub("^.", "", x), "_", x))
peds[, (names(peds)) := lapply(.SD, to_vcf_id)]

# ---- 2. VCF -> GDS ----
gds_file <- file.path(P("WORK_DIR"), "trios_hg38.gds")
snpgdsVCF2GDS(vcf.fn = vcf_file, out.fn = gds_file, method = "biallelic.only",
              ignore.chr.prefix = "chr", verbose = FALSE)
genofile <- snpgdsOpen(gds_file)
on.exit(snpgdsClose(genofile), add = TRUE)

# ---- 3. PCA on infants ----
set.seed(1)
snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2,
                          sample.id = infant_ids$sample.id.vcf,
                          autosome.only = TRUE, num.thread = NCPU,
                          maf = 0.01, missing.rate = 0.05,
                          remove.monosnp = FALSE)
pruned <- unlist(snpset, use.names = FALSE)
message("LD-pruned SNPs: ", length(pruned))

# MAF distribution of pruned SNPs on chr1
snp_chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp_ids <- read.gdsn(index.gdsn(genofile, "snp.id"))
chr1 <- pruned[snp_chr[match(pruned, snp_ids)] == "1"]
maf_df <- data.frame(MAF = snpgdsSNPRateFreq(genofile, snp.id = chr1,
                                             sample.id = infant_ids$sample.id.vcf)$MinorFreq)
g <- ggplot(maf_df, aes(MAF)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
  labs(title = "MAF of LD-pruned chr1 SNPs in infants", x = "Minor allele frequency", y = "SNPs") +
  theme_minimal()
ggsave(fig("infant_chr1_MAF.png"), g, width = 5.5, height = 5.5, dpi = 300, bg = "white")

pca <- snpgdsPCA(genofile, sample.id = infant_ids$sample.id.vcf, snp.id = pruned,
                 autosome.only = TRUE, num.thread = NCPU, verbose = FALSE)
pv <- round(pca$varprop[1:2] * 100, 2)
pc.df <- data.table(sample.id.vcf = pca$sample.id,
                    PC1 = pca$eigenvect[, 1], PC2 = pca$eigenvect[, 2])

# ---- 4. Reported maternal ethnicity ----
pc.df <- merge(pc.df, meta[, .(sample.id.vcf, maternal_ethnicity)],
               by = "sample.id.vcf", all.x = TRUE)

pca_plot <- function(dt, col, legend, subtitle) {
  ggplot(dt, aes(PC1, PC2, color = .data[[col]])) +
    geom_point(size = 2) + theme_minimal() +
    labs(title = sprintf("PCA of infant genotypes (%s LD-pruned autosomal SNPs)",
                         format(length(pruned), big.mark = ",")),
         subtitle = subtitle,
         x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"), color = legend)
}
ggsave(fig("PCA_infants_reported_ethnicity.png"),
       pca_plot(pc.df, "maternal_ethnicity", "Maternal ethnicity", "Reported maternal ethnicity (RNA-seq infants)"),
       width = 6.5, height = 4.5, dpi = 300, bg = "white")

# ---- 5. 3-NN assignment on PC1/PC2 ----
train <- !is.na(pc.df$maternal_ethnicity)
pc.df[, ancestry_group := maternal_ethnicity]
pc.df[!train, ancestry_group := as.character(
  knn(train = pc.df[train, .(PC1, PC2)], test = pc.df[!train, .(PC1, PC2)],
      cl = pc.df$maternal_ethnicity[train], k = 3))]
pc.df[, ancestry_group := tolower(gsub("[ /]", "_", ancestry_group))]

ggsave(fig("PCA_infants_inferred_ancestry.png"),
       pca_plot(pc.df, "ancestry_group", "Inferred ancestry", "3-NN on PC1/PC2 trained on reported maternal ethnicity"),
       width = 6.5, height = 4.5, dpi = 300, bg = "white")

# ---- 6. Ancestry-stratified sample lists (infant + parents) ----
groups <- words("ANCESTRY_GROUPS")
unexpected <- setdiff(unique(pc.df$ancestry_group), groups)
if (length(unexpected)) stop("ancestry labels not in ANCESTRY_GROUPS: ", paste(unexpected, collapse = ", "))

assign_tab <- merge(pc.df[, .(infant_id = sample.id.vcf, ancestry_group, PC1, PC2,
                              reported = !is.na(maternal_ethnicity))],
                    peds, by = "infant_id", all.x = TRUE)
fwrite(assign_tab, file.path(out_dir, "ancestry_assignments.tsv"), sep = "\t")

for (grp in groups) {
  rows <- assign_tab[ancestry_group == grp]
  ids <- unique(na.omit(c(rows$infant_id, rows$mother_id, rows$father_id)))
  writeLines(ids, file.path(out_dir, paste0("samples_", grp, ".txt")))
  message(sprintf("%-8s infants=%d  total samples=%d", grp, nrow(rows), length(ids)))
}
