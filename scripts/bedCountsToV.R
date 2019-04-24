#!/usr/bin/Rscript

###################################################################################################
# Script used to read the bed files with the counts for each bin and merge them into a V matrix ###
# which can be used to perform a V matrix analysis. Takes the processed bed files path as input ###
# and outputs a file with the V matrix in csv format, containing normalized counts in each bin. ###
###################################################################################################
# Done by Alejandro Roca (alekss.ro@gmail.com)                                                  ###
###################################################################################################

# Beginning of the script
cat("Initiating normalization of counts for each bin to add a column in the V matrix\n")

# set working directory
setwd("/home/alekssro/mscthesis")

# Get arguments from terminal call
args <- commandArgs(trailingOnly = T)       # trailingOnly = T, gets only arguments not the call

# Load required packages
# source("https://bioconductor.org/biocLite.R")
# biocLite("rtracklayer")
suppressPackageStartupMessages(require(rtracklayer, quietly = T, warn.conflicts = F))

input.datafile <- args[1]           # "results/genomic_survey/input_data_files.txt" # raw counts files #  
datasheetsamples.file <- args[2]    # "data/DatasetInfoFile.tsv" # datasets (used to get names/structure) #  
binnedgenome.file <- args[3]        # "data/hg19binned.200bp.bed" # binned genome bed file #  
mappabilitytrack.file <- args[4]    # "data/wgEncodeDukeMapabilityUniqueness35bp.uniqueMapRegions.bedGraph" # 
totMappedAllExpFile <- args[5]      # "results/genomic_survey/A549/CTCF/A549_CTCF_allExp_totals.txt" # 
outputVfile <- args[6]              # "results/genomic_survey/A549/V_matrix.csv" # 


#################################
# Get epigenetic mark name ######
#################################
input.data <- scan(input.datafile, what = "character", sep = "\t")
path.sections <- unlist(strsplit(input.data[1], split="[/]"))
epigen.mark <- path.sections[length(path.sections) - 1]
cat("Epignetic mark: ", epigen.mark, "\n")


###################################
# Import datasheet samples file ###
###################################
datasheetsamples <- read.table(datasheetsamples.file, header = F, sep = "\t",
                               stringsAsFactors = F, fill = T)
colnames(datasheetsamples) <- c("cell_line" , "data_type", "assay", "sampleID", "alignfile", "map_type" ,
                                "replicate" ,"lab" ,"dwnld_link")

# Retrieve the number of informative reads mapped for each replicate of the epigenetic mark
totMappedAllExp <- scan(totMappedAllExpFile, what="character")
totMappedAllExp <- as.numeric(unlist(strsplit(totMappedAllExp, split=" ")))


##############################
# Import mappability track ###
##############################
cat("  Loading uniqueness mappability track and calculating mappable genome size...", "\n", sep=" ")
grmap <- import.bed(mappabilitytrack.file)
# Calculate mappable genome size:
mappabilityIntervalSizes <- width(ranges(grmap));
mappabilityGenomeSize <- sum(as.numeric(mappabilityIntervalSizes));


##############################
# Import genomic bins file ###
##############################
cat("  Loading 200-bp binned genome file...", "\n", sep=" ")
# binned.genome <- import.bed(binnedgenome.file)
binned.genome <- read.table(binnedgenome.file, header=F, sep="\t", stringsAsFactors=F, fill=TRUE)
colnames(binned.genome) <- c("chrom" , "chromStart" ,"chromEnd")
binned.genome <- data.frame(binned.genome, counts=rep(0, nrow(binned.genome)), stringsAsFactors=FALSE)
binned.genome.id <- paste(binned.genome$chrom, binned.genome$chromStart, sep="_")


#################################################################
### Import raw-counts and save their counts into the V matrix ###
#################################################################

# Create list to store counts for the different replicates
bincountslist <- list()

# Import raw counts files and get average counts
for (i in c(1:length(input.data))) {

    filepath <- input.data[i]
    cat("  Loading", basename(filepath), "for estimating normalized signal...", "\n", sep=" ")

    # Import single rawCount file (only bins with >= 1 tag)
    bincounts <- read.table(file=filepath, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    colnames(bincounts) <- c("chrom", "chromStart", "chromEnd", "counts")
    bincounts.id <- paste(bincounts$chrom, bincounts$chromStart, sep="_")

    # Add empty bins to the bincounts table (genomic bins with no tags)
    binned.genome2 <- binned.genome
    binned.genome2[binned.genome.id %in% bincounts.id, "counts"] <- bincounts$counts
    bincounts <- binned.genome2
    rm(binned.genome2)

    # Convert bincounts in GRanges object!
    grbins <- GRanges(
        seqnames=( Rle(names(table(bincounts$chrom)[unique(bincounts$chrom)]),
                       as.numeric(table(bincounts$chrom)[unique(bincounts$chrom)]) ) ),
        ranges=IRanges(bincounts$chromStart, bincounts$chromEnd),
        strand=Rle('*', nrow(bincounts)),
        score=bincounts$counts
    )


    # 5.5 Find overlap between bins and uniqueness mappability frames
    cat("    Estimate expected counts by finding overlaps between bins and mappability regions", "\n", sep=" ")
    bin2mapHits <- findOverlaps(query=grbins, subject=grmap, ignore.strand = TRUE)

    # Determine sizes of overlap regions for each query with any match
    overlaps <- pintersect(grbins[queryHits(bin2mapHits)], grmap[subjectHits(bin2mapHits)])
    w <- width(overlaps)

    # Build a three-column table; queryIndexes, subjectIndexes, subjectItemSizes
    queryToSubject <- cbind(queryHits(bin2mapHits), subjectHits(bin2mapHits), w)
    colnames(queryToSubject) <- c("queryidx", "subjectidx", "size")

    # Estimate number of unique-mappability positions within each genomic bin
    queryToSubject <- as.data.frame(queryToSubject, stringsAsFactors=F)
    uniqPos <- tapply(queryToSubject$size, queryToSubject$queryidx, sum)
    allPos <- rep(0, length(bincounts.id))
    allPos[as.numeric(names(uniqPos))] <- uniqPos
    uqmapp <- allPos
    # uqmapp <- as.numeric(tapply(queryToSubject$size, queryToSubject$queryidx, sum))

    # Estimate the scaling factor on sample size
    ScalingFactorLibrarySize <- round((length(bincounts.id)/mean(totMappedAllExp)), 3)

    # Estimate the normalized expected fragment count in each bin
    expCounts <- (uqmapp * (length(bincounts.id) / mappabilityGenomeSize ) ) * ScalingFactorLibrarySize

    # Estimate the normalized observed count in each bin
    obsCounts= score(grbins)
    NormalizedSignals <- obsCounts/expCounts
    NormalizedSignals[!is.numeric(NormalizedSignals) | is.na(NormalizedSignals) | is.infinite(NormalizedSignals)] <- 0

    cat("    Normalized counts for ", basename(filepath), " calculated\n")
    bincountslist[[i]] <- NormalizedSignals
}

# List to matrix
bincountsmatrix <- do.call(cbind, bincountslist)

cat("  Appending combined counts for the epigenetic mark", epigen.mark, " to the V matrix file \n\n")
count2V <- paste(c(epigen.mark, ceiling(rowMeans(bincountsmatrix))), collapse = ",")
write(x = count2V, file = outputVfile, append = T)
# write(binned.genome.id, file = "results/genomic_survey/HepG2/bin_names.txt", sep = ",")
