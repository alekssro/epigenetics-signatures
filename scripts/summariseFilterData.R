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

# Load V matrix and chromosome info
V <- read.csv("results/genomic_survey/HepG2/V_matrix_transposed.csv", header = T)
chrmsL <- read.table("data/chromInfo.txt")

# Assign chrm names and rebuild bin names
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
comb = function(n, x) {
    factorial(n) / factorial(n-x) / factorial(x)
}

# Fit gamma distribution parameters for each epigenetic mark
pars <- matrix(0, ncol = 2, nrow = 8)
for (j in 1:8) {
    subEpigMark <- sample(V[,j] + 1, size = length(V[,j])*0.1)     # pseudocount and 10% sample
    gammaDistr <- fitdistr(subEpigMark, densfun = "gamma", lower = 0.001)
    
    pars[j, ] <- gammaDistr$estimate
}

# Generate Q probability matrix for each element in V

Q <- matrix(1, ncol = ncol(V), nrow = nrow(V))      # init matrix
for (j in 1:8) {
    for (i in 1:10){
        a <- pars[j, 1]
        b <- pars[j, 2]
        n <- V[i, j]
        
        Q[i, j] <- comb(n + b, n) * (a/(a+1))^b * (1/(a+1))^n
    }
}

bins2keep <- apply(head(Q), 1, function(x) sum(x > 0.01) > 1)
print(dim(V))
print(sum(bins2keep))


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
if (makePlots == "-p") {
    for (i in 1:length(chrs)) {
        
        d <- filteredV %>% filter(chr == chrs[i])
        d <- melt(d, id.vars = "binStart")
        
        p <- d %>% filter(variable != "chr") %>%
            ggplot(aes(x=binStart, y=as.numeric(value), col=variable)) + 
            geom_point(aes(color=factor(variable)))
        
        ggsave(filename = paste("plots/", chrs[i], "_readsDistr.png", sep=""), plot = p)
        
    }
}



