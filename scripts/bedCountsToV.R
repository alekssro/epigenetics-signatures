#!/usr/bin/Rscript

###################################################################################################
# Script used to read the bed files with the counts for each bin and merge them into a V matrix ###
# which can be used to perform a V matrix analysis. Takes the processed bed files path as input ###
# and outputs a file with the V matrix in csv format                                            ###
###################################################################################################
# Done by Alejandro Roca (alekss.ro@gmail.com)                                                  ###
###################################################################################################

# TODO (Maybe): add an option to normalize data

# set working directory
setwd("/mnt/DataHDD/alekssro/Nextcloud/Universidad/Bioinformatics/MasterThesis/mscthesis/")

# Get arguments from terminal call
args <- commandArgs(trailingOnly = T)       # trailingOnly = T, gets only arguments not the call

# Load required packages
# source("https://bioconductor.org/biocLite.R")
# biocLite("rtracklayer")
suppressPackageStartupMessages(require(rtracklayer, quietly = T, warn.conflicts = F))
# require(IRanges, verbose = F, warn.conflicts = F)
# require(rtracklayer, verbose = F, warn.conflicts = F)

input.datafile <- "results/genomic_survey/input_data_files.txt" # args[1]           # path for raw counts files
datasheetsamples.file <- "data/DatasetInfoFile.tsv" # args[2]    # datasets used in the analysis (used to get names/structure)
binnedgenome.file <- "data/hg19binned.200bp.bed" # args[3]        # binned genome bed file
mappabilitytrack.file <- "data/wgEncodeDukeMapabilityUniqueness35bp.uniqueMapRegions.bedGraph" # args[4]

###################################
# Import datasheet samples file ###
###################################
cat("Import datasheet samples file", basename(datasheetsamples.file), "\n")
datasheetsamples <- read.table(datasheetsamples.file, header = F, sep = "\t", stringsAsFactors = F, fill = T)
colnames(datasheetsamples) <- c("cell_line" , "data_type", "assay", "sampleID", "alignfile", "map_type" ,
                                "replicate" ,"lab" ,"dwnld_link")

#########################################
# Import counts per bin files path ######
#########################################
input.data <- scan(input.datafile, what = "character", sep = "\t")

##############################
# Import mappability track ###
##############################
cat("Import uniqueness mappability track", "\n", sep=" ")
mappabilitytrack <- import.bed(mappabilitytrack.file)
grmap <- mappabilitytrack

##############################
# Import genomic bins file ###
##############################
cat("Import binned genome file", "\n", sep=" ")
# binned.genome <- import.bed(binnedgenome.file)
binned.genome = read.table(binnedgenome.file, header=F, sep="\t", stringsAsFactors=F, fill=TRUE)
colnames(binned.genome) <- c("chrom" , "chromStart" ,"chromEnd")
binned.genome <- data.frame(binned.genome, counts=rep(0, nrow(binned.genome)), stringsAsFactors=FALSE)
binned.genome <- data.frame(binned.genome, name=paste(binned.genome$chrom, binned.genome$chromStart, sep="_"), stringsAsFactors=F)

#################################################################
### Import raw-counts and save their counts into the V matrix ###
#################################################################

# Calculate mappable genome size:
mappabilityIntervalSizes <- width(ranges(mappabilitytrack));
mappabilityGenomeSize <- sum(as.numeric(mappabilityIntervalSizes));


# Import raw counts files and get average counts
for (i in c(1:length(input.data))) {
    
    filepath <- input.data[i]
    cat("import", basename(filepath), "for estimating normalized signal, observed score, expected score ", "\n", sep=" ")
    
    # Import single rawCount file (only bins with >= 1 tag)
    bincounts <- read.table(file=filepath, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    colnames(bincounts) <- c("chrom", "chromStart", "chromEnd", "counts" )
    bincounts <- data.frame(bincounts, name=paste(bincounts$chrom, bincounts$chromStart, sep="_"), stringsAsFactors=F)
    
    # Add empty bins to the bincounts table (genomic bins with no tags)
    binned.genome2 <- binned.genome
    binned.genome2[match(bincounts$name, binned.genome2$name), "counts"] <- bincounts$counts
    bincounts <- binned.genome2
    rm(binned.genome2)
    
    # 5.4 Convert bincounts in GRanges object!
    ### NB: this 'GRanges' command is good for all inputs that are already sorted by chromosome.
    grbins <- GRanges(
        seqnames=( Rle(names(table(bincounts$chrom)[unique(bincounts$chrom)]), as.numeric(table(bincounts$chrom)[unique(bincounts$chrom)]) ) ),
        ranges=IRanges(bincounts$chromStart, bincounts$chromEnd),
        strand=Rle('*', nrow(bincounts)),
        score=bincounts$counts
    )
    
    
    # 5.5 Find overlap between bins and uniqueness mappability frames
    cat("examine overlap between genomic bins and mappability regions", "\n", sep=" ")
    bin2mapHits <- findOverlaps( query=grbins, subject=grmap, ignore.strand = TRUE)
    
    # Determine sizes of overlap regions for each query with any match
    w <- width(overlapsRanges(ranges(grbins), ranges(grmap)))
    
    
    # Build a three-column table; queryIndexes, subjectIndexes, subjectItemSizes
    queryToSubject <- cbind( bin2mapHits@queryHits, bin2mapHits@subjectHits, w)
    colnames(queryToSubject) <- c("queryidx", "subjectidx", "size")
    # Estiamte number of unique-mappability positions within each genomic bin
    queryToSubject <- as.data.frame(queryToSubject, stringsAsFactors=F)
    uqmapp <- as.numeric(tapply(queryToSubject$size, queryToSubject$queryidx, sum) )
    # Estimate the saling factor on sample size
    ScalingFactorLibrarySize = round( (totMappedSamples[i]/mean(totMappedAllExp)), 3) # totMappedSamples[i] = tot.mapped tags in sample [i]
    # Estimate the normalized expected fragment count in each bin
    expCounts= round( (uqmapp * (totMappedSamples[i] /mappabilityGenomeSize ) )* ScalingFactorLibrarySize, 3)
    # Estimate the normalized observed count in each bin
    obsCounts= round ( (score(grbins)[unique(bin2mapHits@queryHits)] ) , 3)
    NormalizedSignals <- round((obsCounts/expCounts), 3)
    
    
    
    
    # 5.6 Export:
    filepathexport <- gsub( ".rawCounts.bed", ".normalizedSignal.bed", filepath)
    
    cat("exporting normalizedSignal bed file for", filepath, "\n", sep=" ")
    write.table( data.frame( chrom=as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
                             chromStart=start(ranges(grbins))[unique(bin2mapHits@queryHits)],
                             chromEnd=end(ranges(grbins))[unique(bin2mapHits@queryHits)],
                             #name=paste(as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
                             #			start(ranges(grbins))[unique(bin2mapHits@queryHits)], sep="_"),
                             score=NormalizedSignals,
                             #strand=rep("*", length(unique(bin2mapHits@queryHits))),
                             stringsAsFactors=FALSE),
                 file=filepathexport, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);
    
    
    
    
    
    # 5.7 Export only expCount values
    filepathexportExp <- gsub( ".rawCounts.bed", ".expectedScore.bg", filepath)
    cat("exporting expected signal bed file for", filepath, "\n", sep=" ")
    write.table( data.frame( chrom=as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
                             chromStart=start(ranges(grbins))[unique(bin2mapHits@queryHits)],
                             chromEnd=end(ranges(grbins))[unique(bin2mapHits@queryHits)],
                             score=expCounts,
                             stringsAsFactors=FALSE),
                 file=filepathexportExp, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);
    
    
    
    
    
    # 5.8 Export only obsCount values
    filepathexportObs <- gsub( ".rawCounts.bed", ".observedScore.bg", filepath)
    cat("exporting observed signal bed file for", filepath, "\n", sep=" ")
    write.table( data.frame( chrom=as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
                             chromStart=start(ranges(grbins))[unique(bin2mapHits@queryHits)],
                             chromEnd=end(ranges(grbins))[unique(bin2mapHits@queryHits)],
                             score=obsCounts,
                             stringsAsFactors=FALSE),
                 file=filepathexportObs, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);
    
}

# 
# # H3K4me1 modification for H1 cells
# hist_mod <- import("../data/GSM409307_UCSD.H1.H3K4me1.LL228.bed.gz", format="bed")
# 
# # hist_mod
# 
# # All human genes
# # human_genes <- import("../data/hg.bed", format="bed")
# # chrom_names <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
# #                  "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
# #                  "chrX", "chrY")
# # human_genes <- human_genes[seqnames(human_genes) %in% chrom_names, ]    # get only 'normal' chromosome genes
# # human_genes <- reduce(unique(human_genes))      # merge overlapping transcripts
# # human_genes <- flank(human_genes, width = 1000)
# # human_genes
# 
# # brca1 gene (chr17:43044295-43125483) GRCh38 assembly
# # brca1_GRCh38 <- import("Data/brca1_gene_GRCh38.txt", format = "bed")
# 
# # brca1 gene (chr17:41196312-41277500) GRCh37 assembly
# brca1_GRCh37 <- import("../data/brca1_gene_GRCh37.txt", format = "bed")
# 
# sum(countOverlaps(brca1_GRCh37, hist_mod)) # counts for H3K4me1 modifications in brca1 gene in H1 cells
# # sum(countOverlaps(brca1_GRCh37, hist_mod)) 
# 
# 
# # H3K4me3 modification for H1 cells
# hist_mod <- import("../data/GSM409308_UCSD.H1.H3K4me3.LL227.bed.gz", format="bed")
# hist_mod <- reduce(unique(hist_mod))
# sum(countOverlaps(brca1_GRCh37, hist_mod))
