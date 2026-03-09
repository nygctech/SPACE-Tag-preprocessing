library(Signac)
library(GenomicRanges)
library(Matrix)
library(stringr)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
frag_file <- args[1]
peaks_file <- args[2]
output_dir <- args[3]

# Get fragments file
fragments <- CreateFragmentObject(
  path = frag_file,
  cells = NULL,
  validate.fragments = FALSE
)

# Read peaks from macs2 output and covert to granges
peaks <- read.table(peaks_file, sep='\t', header=F)
features <- makeGRangesFromDataFrame(peaks, seqnames.field='V1', start.field='V2', end.field='V3')
features <- keepStandardChromosomes(features, pruning.mode = "coarse")

mtx <- FeatureMatrix(
  fragments,
  features,
  verbose = TRUE
)

# Write output
# barcodes.tsv
write.table(colnames(mtx), paste(output_dir, '/barcodes.tsv', sep=''), quote=F, sep='\t', row.names=F, col.names=F)
# peaks.bed
write.table(str_split_fixed(rownames(mtx), '-', 3), paste(output_dir, '/peaks.bed', sep=''), quote=F, sep='\t', row.names=F, col.names=F)
# matrix.mtx
writeMM(mtx, paste(output_dir, '/matrix.mtx', sep=''))