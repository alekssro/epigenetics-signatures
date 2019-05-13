#!/usr/bin/Rscript

###################################################################################################
# Script used to analyse and summarise the data in the reads matrix V                           ###
# 
###################################################################################################
# Done by Alejandro Roca (alekss.ro@gmail.com)                                                  ###
###################################################################################################

# Set working directory (may need to be changed)
# setwd("/mnt/DataHDD/alekssro/Nextcloud/Universidad/Bioinformatics/MasterThesis/mscthesis")

# Load necessary packages (tidyverse, reshape2)
suppressPackageStartupMessages(require(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(reshape2, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(MASS, quietly = T, warn.conflicts = F))

# Get arguments from terminal call
args <- commandArgs(trailingOnly = T)       # trailingOnly = T, gets only arguments not the call
infile <- args[1]
outfile <- args[2]                          
makePlots <- args[3]                        # if "-p" then makePlots will be True 
dirPlots <- args[4]
makePlots <- ifelse(is.na(makePlots), F, ifelse(makePlots == "-p", T, F))

# Load V matrix and chromosome info
cat("Loading V matrix...\n")
V <- read.csv(infile, header = T)
chrmsL <- read.table("data/chromInfo.txt")

# Assign chrm names and rebuild bin names
cat("Rebuilding bin names...\n")
chrmsNames <- rep("", times=nrow(V))
binNames <- rep(0, times=nrow(V))
start <- 1
for (i in 1:nrow(chrmsL)) {
    
    nBins <- ceiling(chrmsL[i, 2]/200.0)
    
    chr <- rep(as.character(chrmsL[i, 1]), nBins)
    bin <- as.numeric(seq(from = 1, to = chrmsL[i, 2], by=200))
    
    end <- start + nBins - 1
    # print(c(end - start, nBins, length(chr), length(bin)))
    chrmsNames[start:end] <- chr
    binNames[start:end] <- bin
    start <- end + 1
}

V$chr <- chrmsNames
V$binStart <- binNames


## Filter bins in V matrix
# combination function
cat("Filtering bins with higher probabilities to 2.5% in at least one epigenetic mark...\n")

# Filter out those bins which doesn't have at least 1 mark in higher 0.025 quantile
numV <- V[,1:(ncol(V) - 2)]     # exclude chrs names and bin names
thrshld <- apply(numV, 2, function(x) quantile(x, probs = 0.975))   # get total counts for each epigen mark

bins2keep <- apply(numV, 1, function(x) sum(x >= thrshld) > 1)

filteredV <- V[bins2keep, ]

cat("Number of bins after filtering: ", nrow(filteredV), "\n")

cat("Saving filtered V matrix (with chromosome and bin name)...\n")
write.csv(filteredV, file = outfile, quote = F, row.names = F)

# Plot reads coverage by bins along each chromosome
chrs <- unique(V$chr)
if (makePlots) {
    cat("Creating coverage plots for every chromosome, using filtered data...\n")
    for (i in 1:length(chrs)) {
        
        d <- filteredV %>% filter(chr == chrs[i])
        d <- melt(d, id.vars = "binStart", variable.name = "mark", value.name = "counts")
        
        p <- d %>% filter(mark != "chr") %>%
            ggplot(aes(x=binStart, y=as.numeric(counts), col=mark)) + 
            geom_point(aes(color=factor(mark))) +
            xlab("Normalized Reads") + ylab("Bin Position")
        
        ggsave(filename = paste(dirPlots, chrs[i], "_readsDistr.png", sep=""), plot = p)
        
    }
}



