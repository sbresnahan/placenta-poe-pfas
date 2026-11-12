# =============================================================================
# Stage 06 | ornaments_summarizeToGene.R
# Build maternal / paternal allele-count matrices from Ornaments output.
#
#   Rscript ornaments_summarizeToGene.R --level gene|transcript
#
# Per sample (in parallel):
#   1. TPMs -> within-gene transcript usage weights w = TPM_t / sum_{t in g} TPM_t
#      (equal weights when every transcript of the gene has TPM 0)
#   2. allele counts (transcript, tx_pos, L, R) -> maternal / paternal using the
#      sample's haplotype origin on that chromosome (HAP_ORIGIN, stage 05);
#      haplotypes of undefined origin contribute 0
#   3a. transcript level: feature = transcript:tx_pos
#   3b. gene level: transcript positions mapped to genomic SNPs (TXVCF_MERGED);
#       feature = gene:chr:pos with count = sum over the gene's transcripts
#       carrying the SNP of w_t * count_t (union of SNPs across isoforms)
# Rows = union of features over samples (absent = 0);
# columns = <sample>_maternal ..., <sample>_paternal ...
#
# In : QUANT_DIR/<sample>/{allele_counts.txt,tpms.txt}, HAP_ORIGIN, TX_CLASSIF, TXVCF_MERGED
# Out: POE_GENE_COUNTS or POE_TX_COUNTS (first column feature_id)
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages({ library(data.table); library(optparse); library(parallel) })
setDTthreads(1)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--level", type = "character", default = "gene", help = "gene | transcript"))))
stopifnot(opt$level %in% c("gene", "transcript"))
gene_level <- opt$level == "gene"

# ---- 1. transcript usage weights ----
# tpms.txt is JSON-like: "transcript_id": value, -- IDs may contain colons, so
# split each line on its last colon.
read_ornaments_tpms <- function(path) {
  x <- trimws(readLines(path, warn = FALSE))
  x <- sub(",\\s*$", "", x[nzchar(x) & !grepl("^[{}]$", x)])
  at <- regexpr(":[^:]*$", x)
  keep <- at > 0
  key <- gsub('^"|"$', "", trimws(substr(x[keep], 1, at[keep] - 1)))
  val <- suppressWarnings(as.numeric(trimws(substring(x[keep], at[keep] + 1))))
  val[!is.finite(val) | val < 0] <- 0
  data.table(transcript_id = key, tpm = val)
}

tx_usage_weights <- function(tpm_dt, tx2gene) {
  strip <- function(x) sub("\\..*", "", x)          # ignore transcript version
  W <- merge(tpm_dt[, .(tx_norm = strip(transcript_id), tpm)],
             tx2gene[, .(tx_norm = strip(transcript_id), transcript_id, gene_id)], by = "tx_norm")
  W[, w := if (sum(tpm) > 0) tpm / sum(tpm) else rep(1 / .N, .N), by = gene_id]
  unique(W[, .(transcript_id, gene_id, w)])
}

# ---- 2. allele counts -> parent of origin ----
# allele_counts.txt: V1 transcript, V2 tx position, V6 L-allele count, V7 R-allele count
read_allele_counts <- function(path) {
  dt <- fread(path, header = FALSE, sep = "\t", select = c(1, 2, 6, 7), quote = "",
              na.strings = c("NA", ".", ""), col.names = c("transcript_id", "tx_pos", "L", "R"))
  dt[, tx_pos := as.integer(tx_pos)]
  dt[!is.finite(L) | L < 0, L := 0][!is.finite(R) | R < 0, R := 0]
  dt[, .(L = sum(L), R = sum(R)), by = .(transcript_id, tx_pos)]
}

assign_parent_of_origin <- function(ac, sid, hap, tx_chrom) {
  ac <- merge(ac, tx_chrom, by = "transcript_id")
  h <- hap[sample_id == sid, .(chrom, oriL, oriR)]
  ac <- merge(ac, h, by = "chrom", all.x = TRUE)
  ac[, `:=`(maternal = (oriL %in% "maternal") * L + (oriR %in% "maternal") * R,
            paternal = (oriL %in% "paternal") * L + (oriR %in% "paternal") * R)]
  ac[, .(transcript_id, tx_pos, maternal, paternal)]
}

# ---- 3. transcript position -> genomic SNP ----
read_txomic_map <- function(vcf_gz) {
  m <- fread(cmd = sprintf("zcat %s | grep -v '^#' | cut -f1-3", shQuote(vcf_gz)),
             header = FALSE, sep = "\t", col.names = c("transcript_id", "tx_pos", "genomic_id"))
  p <- tstrsplit(m$genomic_id, ":", fixed = TRUE, keep = 1:2)
  unique(m[, .(transcript_id, tx_pos = as.integer(tx_pos), snp_id = paste0(p[[1]], ":", p[[2]]))])
}

summarize_to_gene <- function(poe_tx, tx_map, w_df) {
  DT <- merge(poe_tx, tx_map, by = c("transcript_id", "tx_pos"))
  DT <- merge(DT, w_df, by = "transcript_id")
  DT[, .(maternal = sum(w * maternal), paternal = sum(w * paternal)),
     by = .(feature_id = paste0(gene_id, ":", snp_id))]
}

# ---- inputs ----
cl <- fread(P("TX_CLASSIF"), sep = "\t", select = c(1, 2, 7),
            col.names = c("transcript_id", "chrom", "gene_id"))
tx2gene  <- unique(cl[, .(transcript_id, gene_id)])
tx_chrom <- unique(cl[, .(transcript_id, chrom = sub("^chr", "", chrom))])

hap_wide <- fread(P("HAP_ORIGIN"), sep = "\t", check.names = FALSE)
hap <- melt(hap_wide, id.vars = "chr", variable.name = "col", value.name = "origin", variable.factor = FALSE)
hap[, `:=`(sample_id = sub("_[LR]$", "", col), hap = sub("^.*_", "", col), origin = tolower(trimws(origin)))]
hap <- dcast(hap, chr + sample_id ~ hap, value.var = "origin")
setnames(hap, c("chr", "L", "R"), c("chrom", "oriL", "oriR"))
hap[, chrom := sub("^chr", "", chrom)]

tx_map <- if (gene_level) read_txomic_map(P("TXVCF_MERGED")) else NULL

sample_dirs <- list.dirs(P("QUANT_DIR"), full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[file.exists(file.path(sample_dirs, "allele_counts.txt")) &
                           file.exists(file.path(sample_dirs, "tpms.txt"))]
in_hap <- basename(sample_dirs) %in% hap$sample_id
if (any(!in_hap)) message("no haplotype origin for: ", paste(basename(sample_dirs)[!in_hap], collapse = ", "))
sample_dirs <- sample_dirs[in_hap]
if (!length(sample_dirs)) stop("no Ornaments quantifications with haplotype origin in ", P("QUANT_DIR"))
message(sprintf("[%s] %d samples, %d workers", opt$level, length(sample_dirs), NCPU))

# ---- per sample ----
process_one <- function(d) {
  sid <- basename(d)
  poe <- assign_parent_of_origin(read_allele_counts(file.path(d, "allele_counts.txt")), sid, hap, tx_chrom)
  out <- if (gene_level) {
    summarize_to_gene(poe, tx_map, tx_usage_weights(read_ornaments_tpms(file.path(d, "tpms.txt")), tx2gene))
  } else {
    poe[, .(maternal = sum(maternal), paternal = sum(paternal)),
        by = .(feature_id = paste0(transcript_id, ":", tx_pos))]
  }
  out[, sample_id := sid][]
}
res <- mclapply(sample_dirs, process_one, mc.cores = NCPU, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop("failed samples: ", paste(basename(sample_dirs)[bad], collapse = ", "), "\n", res[[which(bad)[1]]])
res <- rbindlist(res)

# ---- matrices ----
features <- sort(unique(res$feature_id))
samples  <- sort(unique(res$sample_id))
i <- match(res$feature_id, features); j <- match(res$sample_id, samples)
M <- P_ <- matrix(0, length(features), length(samples))
M[cbind(i, j)]  <- res$maternal
P_[cbind(i, j)] <- res$paternal
colnames(M)  <- paste0(samples, "_maternal")
colnames(P_) <- paste0(samples, "_paternal")

out_file <- P(if (gene_level) "POE_GENE_COUNTS" else "POE_TX_COUNTS")
fwrite(data.table(feature_id = features, M, P_), out_file, sep = "\t")
message(sprintf("[%s] wrote %d features x %d samples to %s", opt$level, length(features), length(samples), out_file))
