library(ArchR)
library(SummarizedExperiment)
library(Rsamtools)
library(Matrix)

# 1. Setup
addArchRThreads(threads = snakemake@threads)
addArchRGenome("mm10")
input_path  <- as.character(snakemake@input[["rse"]])
output_dir  <- as.character(snakemake@output[["arrow"]])
sample_name <- "benchmark_sample"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Since we can't easily run the archr implementation, we'll have to create an
# archr project with arrow files and all... to this end, we first generate
# fragments from the count matrix:

# 2. Fast Pseudo-Fragment Creation
rse <- readRDS(input_path)
peaks <- rowRanges(rse)
counts <- as(assay(rse), "dgCMatrix") # Ensure sparse for speed

message("Generating pseudo-fragments...")
# Get indices of non-zero counts
nz <- which(counts > 0, arr.ind = TRUE)

# Create a data frame of coordinates (one line per non-zero peak-cell pair)
# We use the start of the peak as a point-fragment
frag_df <- data.frame(
  chr    = as.character(seqnames(peaks)[nz[,1]]),
  start  = start(peaks)[nz[,1]]+1L,
  end    = end(peaks)[nz[,1]]-1L,
  barcode = colnames(rse)[nz[,2]]
)
frag_df$barcode <- gsub("^.*\\.", "", frag_df$barcode)

chr_levels <- seqlevels(peaks)
frag_df <- frag_df[order(match(frag_df$chr, chr_levels), frag_df$start), ]

# 3. Quick Zip and Index
# This is much faster than expanding counts and creates a tiny file
tmp_frag <- file.path(output_dir, "pseudo_frags.tsv")
write.table(frag_df, tmp_frag, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
zipped_frag <- bgzip(tmp_frag, overwrite = TRUE)
indexTabix(zipped_frag, format = "bed")
unlink(tmp_frag)

# 4. Create Arrow File
arrow_file <- createArrowFiles(
    inputFiles = zipped_frag, 
    sampleNames = sample_name,
    outputNames = sample_name,
    minTSS = 0, 
    minFrags = 0, 
    minFragSize = 1,
    maxFragSize = 2000,
    addTileMat = FALSE,
    addGeneScoreMat = FALSE,    
    force = TRUE
)

# 5. Initialize Project and Build the Matrix
# ArchR will now see one fragment for every peak that had a count > 0
proj <- ArchRProject(
    ArrowFiles = arrow_file, 
    outputDirectory = output_dir,
    copyArrows = TRUE
)

# Crucial: Add the peak set and let ArchR count the "fragments" 
# into a standard PeakMatrix structure.
proj <- addPeakSet(proj, peakSet = peaks, force = TRUE)
proj <- addPeakMatrix(proj, force = TRUE)

saveArchRProject(proj, load = FALSE)
