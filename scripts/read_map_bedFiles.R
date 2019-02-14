# source("https://bioconductor.org/biocLite.R")
# biocLite("rtracklayer")
library(rtracklayer)

# H3K4me1 modification for H1 cells
hist_mod <- import("../data/GSM409307_UCSD.H1.H3K4me1.LL228.bed.gz", format="bed")
hist_mod <- reduce(unique(hist_mod))
# hist_mod

# All human genes
# human_genes <- import("../data/hg.bed", format="bed")
# chrom_names <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
#                  "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
#                  "chrX", "chrY")
# human_genes <- human_genes[seqnames(human_genes) %in% chrom_names, ]    # get only 'normal' chromosome genes
# human_genes <- reduce(unique(human_genes))      # merge overlapping transcripts
# human_genes <- flank(human_genes, width = 1000)
# human_genes

# brca1 gene (chr17:43044295-43125483) GRCh38 assembly
# brca1_GRCh38 <- import("Data/brca1_gene_GRCh38.txt", format = "bed")

# brca1 gene (chr17:41196312-41277500) GRCh37 assembly
brca1_GRCh37 <- import("../data/brca1_gene_GRCh37.txt", format = "bed")

sum(countOverlaps(brca1_GRCh37, hist_mod)) # counts for H3K4me1 modifications in brca1 gene in H1 cells
# sum(countOverlaps(brca1_GRCh37, hist_mod)) 

# # H3K4me1 modification for liver Cells
# hist_mod <- import("../data/GSE18927_RAW/GSM1027296_UW.CD19_Primary_Cells.H3K4me1.RO_01679.Histone.DS22584.bed.gz", format="bed")
# hist_mod <- reduce(unique(hist_mod))
# hist_mod
# sum(countOverlaps(brca1_GRCh37, hist_mod)) # counts for H3K4me1 modifications in brca1 gene in H1 cells

# H3K4me3 modification for H1 cells
hist_mod <- import("../data/GSM409308_UCSD.H1.H3K4me3.LL227.bed.gz", format="bed")
hist_mod <- reduce(unique(hist_mod))
sum(countOverlaps(brca1_GRCh37, hist_mod))
