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
# outfile <- args[1]                          # uncomment to add path in command line
outfile <- "results/genomic_survey/HepG2/filteredV.csv"
makePlots <- args[2]                        # if "-p" then makePlots will be True 
makePlots <- ifelse(is.na(makePlots), F, ifelse(makePlots == "-p", T, F))

# Load V matrix and chromosome info
cat("Loading V matrix...\n")
V <- read.csv("results/genomic_survey/HepG2/V_matrix_transposed.csv", header = T)
chrmsL <- read.table("data/chromInfo.txt")

# Assign chrm names and rebuild bin names
cat("Rebuilding bin names...\n")
chrmsNames <- rep("", times=length(V[,1]))
binNames <- rep(0, times=length(V[,1]))
start <- 1
for (i in 1:length(chrmsL[,1])) {
    
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

# Filter out those bins which doesn't have at least 1% of the read frequency in any of the marks
numV <- V[,1:8]
thrshld <- apply(numV, 2, function(x) quantile(x, probs = 0.975))   # get total counts for each epigen mark

bins2keep <- apply(numV, 1, function(x) sum(x >= thrshld) > 1)

sum(bins2keep)

filteredV <- V[bins2keep, ]

cat("Number of bins after filtering: ", dim(filteredV), "\n")

cat("Saving filtered V matrix (with chromosome and bin name)...\n")
write.csv(filteredV, file = outfile, quote = F, row.names = F)

# Plot reads coverage by bins along each chromosome
chrs <- unique(V$chr)
if (makePlots) {
    cat("Creating coverage plots for every chromosome, using filtered data...\n")
    for (i in 1:length(chrs)) {
        
        d <- filteredV %>% filter(chr == chrs[i])
        d <- melt(d, id.vars = "binStart")
        
        p <- d %>% filter(variable != "chr") %>%
            ggplot(aes(x=binStart, y=as.numeric(value), col=variable)) + 
            geom_point(aes(color=factor(variable))) +
            xlab("Normalized Reads") + ylab("Bin Position")
        
        ggsave(filename = paste("plots/", chrs[i], "_readsDistr.png", sep=""), plot = p)
        
    }
}



