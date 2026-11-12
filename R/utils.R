# =============================================================================
# R/utils.R -- input readers and small helpers shared by stages 07-09
# =============================================================================

#' Keep features (gene or transcript) with at least N informative SNP rows.
#' Row names are <feature>:<chr>:<pos> or <feature>:<pos>; the ":chr..:"
#' segment is dropped so every row becomes <feature>:<pos>.
filter_features_by_snps <- function(x, N = 3) {
  stopifnot(!is.null(rownames(x)))
  rn <- sub(":chr[^:]*:", ":", rownames(x), ignore.case = TRUE)
  rownames(x) <- rn
  ids <- sub(":[^:]*$", "", rn)
  n <- table(ids)
  x[ids %in% names(n)[n >= N], , drop = FALSE]
}

#' Read a stage-06 count matrix (first column feature_id; then
#' <sample>_maternal ... <sample>_paternal) and apply the SNP filter
read_poe_counts <- function(file, min_snps = 3) {
  dt <- data.table::fread(file, sep = "\t", check.names = FALSE)
  m <- as.matrix(dt[, -1, with = FALSE])
  rownames(m) <- dt[[1]]
  colnames(m) <- gsub("[-.]", "_", colnames(m))
  filter_features_by_snps(m, N = min_snps)
}

#' Transcript -> gene map from the SQANTI3 classification file.
#' Transcript IDs carry no chromosome segment (ESPRESSO:chrN:a:b -> ESPRESSO:a:b),
#' matching the feature IDs produced by filter_features_by_snps().
read_tx2g <- function(file = P("TX_CLASSIF")) {
  cl <- data.table::fread(file, sep = "\t", select = c("isoform", "associated_gene"))
  data.table::setnames(cl, c("transcript_id", "gene_id"))
  cl[grepl("ESPRESSO", transcript_id), transcript_id := sub("chr[0-9XYM]+:", "", transcript_id)]
  unique(cl, by = "transcript_id")
}

#' Ensembl gene ID -> symbol table (stage 00 output)
read_ensembl <- function(file = P("PLACENTA_DB")) {
  e <- data.table::fread(file, sep = "\t")
  e <- unique(data.frame(gene_id = e[[1]], symbol = e[[2]], stringsAsFactors = FALSE))
  if (!"ENSG00000124019" %in% e$gene_id)
    e <- rbind(e, data.frame(gene_id = "ENSG00000124019", symbol = "FAM124B"))
  e
}

#' |beta| corresponding to a majority-allele ratio r
ratio_to_abs_beta <- function(r) abs(qlogis(r))

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Feature IDs -> gene symbols (falls back to the input ID)
convert_to_symbols <- function(feature_ids, ensembl, tx2g = NULL, level = "gene") {
  feature_ids <- as.character(feature_ids)
  gene_id <- if (level == "transcript") {
    as.character(tx2g$gene_id[match(feature_ids, tx2g$transcript_id)])
  } else feature_ids
  g <- sub("\\.\\d+$", "", gene_id)
  out <- feature_ids
  sym <- ensembl$symbol[match(g, ensembl$gene_id)]
  ok <- !is.na(g) & !is.na(sym) & nzchar(sym)
  out[ok] <- sym[ok]
  # novelGene_<ENSG>... : substitute the symbol for the embedded Ensembl ID
  nov <- !is.na(g) & grepl("^novelGene_.*ENSG[0-9]+", g)
  if (any(nov)) {
    ensg <- sub("^.*?(ENSG[0-9]+).*$", "\\1", g[nov])
    s2 <- ensembl$symbol[match(ensg, ensembl$gene_id)]
    hit <- !is.na(s2) & nzchar(s2)
    out[which(nov)[hit]] <- mapply(function(x, s) sub("ENSG[0-9]+", s, x), g[nov][hit], s2[hit],
                                   USE.NAMES = FALSE)
  }
  out
}

#' Sample IDs listed in EXCLUDE_SAMPLES_FILE (one per line; empty if unset)
excluded_samples <- function() {
  f <- Sys.getenv("EXCLUDE_SAMPLES_FILE")
  if (!nzchar(f)) return(character(0))
  if (!file.exists(f)) stop("EXCLUDE_SAMPLES_FILE not found: ", f)
  x <- trimws(readLines(f))
  x[nzchar(x) & !startsWith(x, "#")]
}
