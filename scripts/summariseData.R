# #Resampling 1 set of counts using a multinomial distribution (here 3 categories)
# 
# ncounts <- c(100, 230, 45)
# 
# tot <- sum(ncounts)
# 
# vecp <- ncounts/tot
# 
# rmultinom(n = 1,size = tot,prob = vecp)

library(tidyverse)
library(reshape2)

setwd("/mnt/DataHDD/alekssro/Nextcloud/Universidad/Bioinformatics/MasterThesis/mscthesis")

V <- read.csv("results/genomic_survey/HepG2/V_matrix_transposed.csv", header = T)
chrmsL <- read.table("data/chromInfo.txt")

# assign chrm names and rebuild bin names
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

row_means <- rowMeans(V[,1:8])

noHits_i <- row_means < max(row_means) * 0.001
sum(noHits_i)

filteredV <- V[!noHits_i,]
dim(filteredV)

chrs <- unique(V$chr)

plots <- list()
for (i in 1:length(chrs)) {
    d <- V %>% filter(chr == chrs[i])
    d <- melt(d, id.vars = "binStart")
    
    p <- d %>% filter(variable != "chr") %>%
        ggplot(aes(x=binStart, y=value, col=variable)) + geom_point()
        
}


ggplot(V, aes(binStart)) +
    facet_wrap()
    geom_point() +
    facet_wrap(~chr)

tot_counts <- sum(V)

row_sums <- apply(V, 1, sum)

plot(V[,1])
