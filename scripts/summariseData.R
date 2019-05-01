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

for (i in 1:length(chrs)) {
    d <- filteredV %>% filter(chr == chrs[i])
    d <- melt(d, id.vars = "binStart")
    
    p <- d %>% filter(variable != "chr") %>%
        ggplot(aes(x=binStart, y=as.numeric(value), col=variable)) + 
        geom_point(aes(color=factor(variable)))

    ggsave(filename = paste("plots/", chrs[i], "_readsDistr.png", sep=""), plot = p)
        
}

source(file = "scripts/more/MajorizeMinimizeNQP.R")

epimarks <- colnames(filteredV[,1:8])
chr1V <- filteredV %>% filter(chr == "chr1")
V <- as.matrix(chr1V[,1:8])
n_bins <- nrow(V)
n_epimarks <- ncol(V)
nIter <- 30   # From a run of 100 iters, RSS doesn't change after 25-30 iterations
n_signatures <- 3   # stomach, skin and background
s <- 1e-3      # span needed for convergence



# dim(W) 10 x 3
# dim(H) 3 x 6

MM.RSS <- c()
MM.RSS.n <- matrix(data = 0, nrow = nIter, ncol = 8)
for (n in 1:8) {
    
    n_signatures <- n
    # Initialize W and load H
    W <- matrix(runif(n_bins*n_signatures, 0, 1), nrow = n_bins, ncol = n_signatures)
    H <- matrix(runif(n_signatures*n_epimarks, 0, 1), nrow = n_signatures, ncol = n_epimarks)
    
    for (iter in 1:nIter) {
        
        # Update rows in W; every patient
        for (i in 1:n_bins) {  # nrow(W)
            v <- V[i, ]     # we use rows of V
            
            # For A and b, opposite dimensions as in column update in H
            A = 2 * H %*% t(H)     
            b = as.vector(-2 * H %*% v)
            
            MMres <- MajorizeMinimizeNQP(A = A, b = b, init = W[i,], maxIter = 10, tol = 1e-10)
            
            # Replace old by updated row
            W[i, ] <- MMres$x
        }
        
        # Calculate RSS for each iteration of the algorithm
        MM.RSS[iter] <- sum((V - W %*% H) ^ 2)
        
        MM.RSS.n[iter, n] <- MM.RSS[iter]
        
        # # Keep track of the convergence condition evolution
        # if (iter %% 10 == 0) {cat("Iteration ", iter, ", log RSS decrease: ", 
        #                           log(MM.RSS[iter - 1]) - log(MM.RSS[iter]), "\n", 
        #             sep = '')}
        
        # Convergence conditioned on log(RSS) increase
        if ( (iter > 1) && ((log(MM.RSS[iter - 1]) - log(MM.RSS[iter])) < s) ) {break}
        
    }
    
    print(paste("N:", n, "sum(V - W * H)^2 =", MM.RSS[iter]))
}

iters <- apply(MM.RSS.n, 2, function(x) sum(x > 0))
meanRSS <- colSums(MM.RSS.n) / iters
plot(x = 1:8, y = meanRSS)

ggplot(NULL) + 
    geom_line(mapping = aes(x = 1:n_bins, W[1:length(W[,1]),1], color = "Signature 1")) + 
    geom_line(mapping = aes(x = 1:n_bins, W[1:length(W[,1]),2], color = "Signature 2")) + 
    geom_line(mapping = aes(x = 1:n_bins, W[1:length(W[,1]),3], color = "Signature 3")) +
    theme_minimal() +
    xlab("Patient") + ylab('Mutation Load') +
    ggtitle('Mutation Load for all patients towards each signature')

