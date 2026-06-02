library(dada2); packageVersion("dada2")

args <- commandArgs(trailingOnly = TRUE)
amp <- args[1]

if (is.na(amp) || amp == "") {
  stop("Usage: Rscript scripts/merge_runs_remove_chimera.R <amplicon>")
}

seqtab_files <- sort(Sys.glob(file.path("dada2", "*", amp, "seqtab.rds")))

if (length(seqtab_files) == 0) {
  stop(paste("No seqtab.rds files found for amplicon:", amp))
}

cat("Found", length(seqtab_files), "sequence tables for", amp, "\n")
print(seqtab_files)

outdir <- file.path("dada2_merged", amp)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Read and merge
seqtabs <- lapply(seqtab_files, readRDS)
st.all <- do.call(mergeSequenceTables, seqtabs)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method = "consensus", multithread = TRUE)

# Save core R objects
saveRDS(st.all, file.path(outdir, "seqtab_merged.rds"))
saveRDS(seqtab.nochim, file.path(outdir, "seqtab_merged_nochim.rds"))

# Convert to ASV x sample orientation
asv_mat <- t(seqtab.nochim)
seqs <- rownames(asv_mat)
asv_ids <- paste0("ASV_", seq_len(nrow(asv_mat)))

# Replace rownames with ASV IDs
rownames(asv_mat) <- asv_ids

# Lookup table
asv_lookup <- data.frame(
  ASV_ID = asv_ids,
  Sequence = seqs,
  stringsAsFactors = FALSE
)

# Counts table with IDs only
write.csv(asv_mat, file.path(outdir, "asv_counts.csv"), quote = FALSE)

# Lookup table
write.csv(asv_lookup, file.path(outdir, "asv_lookup.csv"), row.names = FALSE, quote = FALSE)

# Summary
track <- data.frame(
  step = c("merged_before_chimera", "nonchimera"),
  reads = c(sum(st.all), sum(seqtab.nochim))
)

write.csv(track, file.path(outdir, "merge_chimera_summary.csv"),
          row.names = FALSE, quote = FALSE)

cat("Merged table dimensions (samples x ASVs):", dim(st.all)[1], "x", dim(st.all)[2], "\n")
cat("Nonchimera table dimensions (samples x ASVs):", dim(seqtab.nochim)[1], "x", dim(seqtab.nochim)[2], "\n")
cat("Done.\n")
