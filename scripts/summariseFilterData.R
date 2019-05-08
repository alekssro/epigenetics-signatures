#!/usr/bin/Rscript

###################################################################################################
# Script used to analyse and summarise the data in the reads matrix V                           ###
# 
###################################################################################################
# Done by Alejandro Roca (alekss.ro@gmail.com)                                                  ###
###################################################################################################

# #Resampling 1 set of counts using a multinomial distribution (here 3 categories)
# 
# ncounts <- c(100, 230, 45)
# 
# tot <- sum(ncounts)
# 
# vecp <- ncounts/tot
# 
# rmultinom(n = 1,size = tot,prob = vecp)

# Set working directory (may need to be changed)
# setwd("/mnt/DataHDD/alekssro/Nextcloud/Universidad/Bioinformatics/MasterThesis/mscthesis")

# Load necessary packages (tidyverse, reshape2)
suppressPackageStartupMessages(require(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(reshape2, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(MASS, quietly = T, warn.conflicts = F))

# Get arguments from terminal call
args <- commandArgs(trailingOnly = T)       # trailingOnly = T, gets only arguments not the call
makePlots <- args[1]           # "results/genomic_survey/input_data_files.txt" # raw counts files #  
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
cat("Filtering bins with higher probabilities to 1% in at least one epigenetic mark...\n")
comb = function(n, x) {
    factorial(n) / factorial(n-x) / factorial(x)
}

# Fit gamma distribution parameters for each epigenetic mark
cat("    Obtaining gamma distribution parameters for each epigenetic mark, based on reads...\n")
pars <- matrix(0, ncol = 2, nrow = 8)
for (j in 1:8) {
    subEpigMark <- sample(V[,j] + 1, size = length(V[,j])*0.01)     # pseudocount and 1% sample
    gammaDistr <- fitdistr(subEpigMark, densfun = "gamma", lower = 0.001)
    
    pars[j, ] <- gammaDistr$estimate
}

# Generate Q probability matrix for each element in V
cat("    Generating Q matrix of probabilities for each read in V matrix...\n")
Q <- matrix(0, ncol = ncol(V[,1:8]), nrow = nrow(V[,1:8]))      # init matrix
for (j in 1:8) {
    a <- pars[j, 1]
    b <- pars[j, 2]
    for (i in 1:nrow(V)){

        n <- V[i, j] + 1    # add pseudocount
        
        Q[i, j] <- comb(n + b, n) * (a/(a+1))^b * (1/(a+1))^n
        print(i)
        print(Q[i, j])
        if (i > 100) {
            break
        }
    }
    break
}

warnings()

print(V[50:100, ])
print(Q[50:100, ])

# bins2keep <- apply(Q, 1, function(x) sum(x > 0.01) > 1)
# cat("Number of bins before filtering: ", dim(V)[1], "\n")
# cat("Number of bins after filtering: ", sum(bins2keep), "\n")


# row_means <- rowMeans(V[,1:8])
# 
# noHits_i <- row_means < max(row_means) * 0.001
# sum(noHits_i)
# 
# filteredV <- V[!noHits_i,]
# dim(filteredV)
# 
# chrs <- unique(V$chr)

# Plot reads coverage by bins along each chromosome
if (makePlots) {
    for (i in 1:length(chrs)) {
        
        d <- filteredV %>% filter(chr == chrs[i])
        d <- melt(d, id.vars = "binStart")
        
        p <- d %>% filter(variable != "chr") %>%
            ggplot(aes(x=binStart, y=as.numeric(value), col=variable)) + 
            geom_point(aes(color=factor(variable)))
        
        ggsave(filename = paste("plots/", chrs[i], "_readsDistr.png", sep=""), plot = p)
        
    }
}



