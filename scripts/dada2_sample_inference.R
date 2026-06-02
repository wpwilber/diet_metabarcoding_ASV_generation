library(dada2); packageVersion("dada2")

# Define arguments
args <- commandArgs(trailingOnly = TRUE)
project <- args[1]
amp <- args[2]

# Warnings
if (is.na(project) || project == "") {
  stop("Project name must be supplied")
}
if (is.na(amp) || amp == "") {
  stop("Amplicon name must be supplied")
}

# Create output directory
top_outdir <- file.path("dada2", project, amp)
if (!dir.exists(top_outdir)) {
  dir.create(top_outdir, recursive = TRUE)
}

# Discover forward read paths
path <- file.path("trim_clean_qc", project, "trimmed")
pattern_f <- paste0("_", amp, "_R1\\.primertrim\\.fastq\\.gz$")

filtFs <- sort(list.files(
  path,
  pattern = pattern_f,
  full.names = TRUE
))

if (length(filtFs) == 0) {
  stop("No forward reads found.")
}

# Derive reverse read paths
filtRs <- gsub("_R1.primertrim.fastq.gz$", "_R2.primertrim.fastq.gz", filtFs)

# Extract file names
sample_names <- sub(paste0("_", amp, "_R1\\.primertrim\\.fastq\\.gz$"), "", basename(filtFs))
sample_namesR <- sub(paste0("_", amp, "_R2\\.primertrim\\.fastq\\.gz$"), "", basename(filtRs))

# Validate
if (!identical(sample_names, sample_namesR)) {
  stop("Forward and reverse files do not match.")
}

names(filtFs) <- sample_names
names(filtRs) <- sample_names

set.seed(100)

# Learn error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Visualize estimated error rates
pdf(file.path(top_outdir, paste0(amp, "_pooled_error_model_forward.pdf")))
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf(file.path(top_outdir, paste0(amp, "_pooled_error_model_reverse.pdf")))
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Helper
getN <- function(x) sum(getUniques(x))

# Read retention tracking
input_counts <- setNames(numeric(length(sample_names)), sample_names)
denoisedF_counts <- setNames(numeric(length(sample_names)), sample_names)
denoisedR_counts <- setNames(numeric(length(sample_names)), sample_names)
merged_counts <- setNames(numeric(length(sample_names)), sample_names)

# Keep only mergers for final sequence table
mergers <- vector("list", length(sample_names))
names(mergers) <- sample_names

# Sample inference and merge paired-end reads
for (sam in sample_names) {
  cat("Processing:", sam, "\n")
  
  derepF <- derepFastq(filtFs[[sam]])
  derepR <- derepFastq(filtRs[[sam]])
  
  input_counts[sam] <- sum(derepF$uniques)
  
  ddF <- dada(derepF, err = errF, multithread = TRUE)
  ddR <- dada(derepR, err = errR, multithread = TRUE)
  
  denoisedF_counts[sam] <- getN(ddF)
  denoisedR_counts[sam] <- getN(ddR)
  
  mergers[[sam]] <- mergePairs(ddF, derepF, ddR, derepR)
  merged_counts[sam] <- getN(mergers[[sam]])
  
  rm(derepF, derepR, ddF, ddR)
  gc(verbose = FALSE)
}

# Build read retention table
track <- data.frame(
  input     = input_counts,
  denoisedF = denoisedF_counts,
  denoisedR = denoisedR_counts,
  merged    = merged_counts,
  check.names = FALSE
)

write.csv(track, file.path(top_outdir, "read_retention.csv"), row.names = TRUE)

# Build sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(top_outdir, "seqtab.rds"))

# Optional cleanup
rm(mergers)
gc(verbose = FALSE)
