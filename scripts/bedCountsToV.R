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
setwd("/home/alekssro/mscthesis")

# Get arguments from terminal call
args <- commandArgs(trailingOnly = T)       # trailingOnly = T, gets only arguments not the call

# Load required packages
# source("https://bioconductor.org/biocLite.R")
# biocLite("rtracklayer")
suppressPackageStartupMessages(require(rtracklayer, quietly = T, warn.conflicts = F))
# require(IRanges, verbose = F, warn.conflicts = F)
# require(rtracklayer, verbose = F, warn.conflicts = F)

input.datafile <-  args[1] #  #    "results/genomic_survey/input_data_files.txt" #        # path for raw counts files
datasheetsamples.file <- args[2] # "data/DatasetInfoFile.tsv" #   # datasets used in the analysis (used to get names/structure)
binnedgenome.file <- args[3]    #"data/hg19binned.200bp.bed" #     # binned genome bed file
mappabilitytrack.file <- args[4] # "data/wgEncodeDukeMapabilityUniqueness35bp.uniqueMapRegions.bedGraph" #  "data/wgEncodeDukeMapabilityUniqueness35bp.uniqueMapRegions.bedGraph" # args[4]
totMappedAllExpFile <- args[5] # "results/genomic_survey/A549/CTCF/A549_CTCF_allExp_totals.txt" #
outputVfile <- args[6]  # "results/genomic_survey/A549/V_matrix.csv"

print(c(outputVfile, totMappedAllExpFile))
# Retrieve the number of informative reads mapped for each sample Bi in the group P  (cell-linX/Signal-trackY)
totMappedAllExp <- scan(totMappedAllExpFile, what="character")
totMappedAllExp <- as.numeric(unlist(strsplit(totMappedAllExp, split=" ")))

###################################
# Import datasheet samples file ###
###################################
cat("Import datasheet samples file", basename(datasheetsamples.file), "\n")
datasheetsamples <- read.table(datasheetsamples.file, header = F, sep = "\t", stringsAsFactors = F, fill = T)
colnames(datasheetsamples) <- c("cell_line" , "data_type", "assay", "sampleID", "alignfile", "map_type" ,
                                "replicate" ,"lab" ,"dwnld_link")

#################################
# Get epigenetic mark name ######
#################################
input.data <- scan(input.datafile, what = "character", sep = "\t")
path.sections <- unlist(strsplit(input.data[1], split="[/]"))
epigen.mark <- path.sections[length(path.sections) - 1]

##############################
# Import mappability track ###
##############################
cat("Import uniqueness mappability track", "\n", sep=" ")
grmap <- import.bed(mappabilitytrack.file)
# Calculate mappable genome size:
mappabilityIntervalSizes <- width(ranges(grmap));
mappabilityGenomeSize <- sum(as.numeric(mappabilityIntervalSizes));

##############################
# Import genomic bins file ###
##############################
cat("Import binned genome file", "\n", sep=" ")
# binned.genome <- import.bed(binnedgenome.file)
binned.genome <- read.table(binnedgenome.file, header=F, sep="\t", stringsAsFactors=F, fill=TRUE)
colnames(binned.genome) <- c("chrom" , "chromStart" ,"chromEnd")
binned.genome <- data.frame(binned.genome, counts=rep(0, nrow(binned.genome)), stringsAsFactors=FALSE)
binned.genome.id <- paste(binned.genome$chrom, binned.genome$chromStart, sep="_")
# binned.genome <- data.frame(binned.genome, name=paste(binned.genome$chrom, binned.genome$chromStart, sep="_"), stringsAsFactors=F)

#################################################################
### Import raw-counts and save their counts into the V matrix ###
#################################################################

# Create list to store counts for the different replicates
bincountslist <- list()

# Import raw counts files and get average counts
for (i in c(1:length(input.data))) {

    filepath <- input.data[i]
    cat("import", basename(filepath), "for estimating normalized signal, observed score, expected score ", "\n", sep=" ")

    # Import single rawCount file (only bins with >= 1 tag)
    bincounts <- read.table(file=filepath, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    colnames(bincounts) <- c("chrom", "chromStart", "chromEnd", "counts")
    bincounts.id <- paste(bincounts$chrom, bincounts$chromStart, sep="_")
    # bincounts <- data.frame(bincounts, name=paste(bincounts$chrom, bincounts$chromStart, sep="_"), stringsAsFactors=F)

    # Add empty bins to the bincounts table (genomic bins with no tags)
    binned.genome2 <- binned.genome
    binned.genome2[binned.genome.id %in% bincounts.id, "counts"] <- bincounts$counts
    bincounts <- binned.genome2

    rm(binned.genome2)

    # 5.4 Convert bincounts in GRanges object!
    ## NB: this 'GRanges' command is good for all inputs that are already sorted by chromosome.
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
    overlaps <- pintersect(grbins[queryHits(bin2mapHits)], grmap[subjectHits(bin2mapHits)])
    w <- width(overlaps)
    # intersect.ranges <- import.bed(intersect.file)
    # w <- width(intersect.ranges)
    # overlaps <- overlapsRanges(ranges(grbins), ranges(grmap))
    # w <- width(overlaps)

    length(bin2mapHits)
    length(w)

    # Build a three-column table; queryIndexes, subjectIndexes, subjectItemSizes
    queryToSubject <- cbind( queryHits(bin2mapHits), subjectHits(bin2mapHits), w)
    colnames(queryToSubject) <- c("queryidx", "subjectidx", "size")
    # Estimate number of unique-mappability positions within each genomic bin
    queryToSubject <- as.data.frame(queryToSubject, stringsAsFactors=F)
    uqmapp <- as.numeric(tapply(queryToSubject$size, queryToSubject$queryidx, sum) )
    # Estimate the scaling factor on sample size
    ScalingFactorLibrarySize = round((length(bincounts.id)/mean(totMappedAllExp)), 3)
    # Estimate the normalized expected fragment count in each bin
    expCounts= round( (uqmapp * (length(bincounts.id) / mappabilityGenomeSize ) ) * ScalingFactorLibrarySize, 3)
    # Estimate the normalized observed count in each bin
    obsCounts= round ( (score(grbins)[unique(queryHits(bin2mapHits))] ) , 3)
    NormalizedSignals <- round((obsCounts/expCounts), 3)
    length(NormalizedSignals)
    head(NormalizedSignals)

    bincountslist[[i]] <- NormalizedSignals


    # # 5.6 Export:
    # filepathexport <- gsub( ".rawCounts.bed", ".normalizedSignal.bed", filepath)
    #
    # cat("exporting normalizedSignal bed file for", filepath, "\n", sep=" ")
    # write.table( data.frame( chrom=as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
    #                          chromStart=start(ranges(grbins))[unique(bin2mapHits@queryHits)],
    #                          chromEnd=end(ranges(grbins))[unique(bin2mapHits@queryHits)],
    #                          #name=paste(as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
    #                          #			start(ranges(grbins))[unique(bin2mapHits@queryHits)], sep="_"),
    #                          score=NormalizedSignals,
    #                          #strand=rep("*", length(unique(bin2mapHits@queryHits))),
    #                          stringsAsFactors=FALSE),
    #              file=filepathexport, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);
    #
    #
    #
    #
    #
    # # 5.7 Export only expCount values
    # filepathexportExp <- gsub( ".rawCounts.bed", ".expectedScore.bg", filepath)
    # cat("exporting expected signal bed file for", filepath, "\n", sep=" ")
    # write.table( data.frame( chrom=as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
    #                          chromStart=start(ranges(grbins))[unique(bin2mapHits@queryHits)],
    #                          chromEnd=end(ranges(grbins))[unique(bin2mapHits@queryHits)],
    #                          score=expCounts,
    #                          stringsAsFactors=FALSE),
    #              file=filepathexportExp, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);
    #
    #
    #
    #
    #
    # # 5.8 Export only obsCount values
    # filepathexportObs <- gsub( ".rawCounts.bed", ".observedScore.bg", filepath)
    # cat("exporting observed signal bed file for", filepath, "\n", sep=" ")
    # write.table( data.frame( chrom=as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
    #                          chromStart=start(ranges(grbins))[unique(bin2mapHits@queryHits)],
    #                          chromEnd=end(ranges(grbins))[unique(bin2mapHits@queryHits)],
    #                          score=obsCounts,
    #                          stringsAsFactors=FALSE),
    #              file=filepathexportObs, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);

}

# List to matrix
bincountsmatrix <- do.call(cbind, bincountslist)

cat("Appending combined counts for the epigenetic mark to the V matrix file \n")
count2V <- paste(c(epigen.mark, ceiling(rowMeans(bincountsmatrix))), collapse = ",")
write(x = count2V, file = outputVfile, append = T)
