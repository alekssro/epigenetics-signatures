library(tidyverse)

setwd("/mnt/DataHDD/alekssro/Nextcloud/Universidad/Bioinformatics/MasterThesis/mscthesis/")

source(file = "scripts/MajorizeMinimizeNQP.R")

V <- read.csv("results/genomic_survey/A549/V_matrix.csv")
epimarks <- colnames(V)
V <- as.matrix(V)
n_bins <- nrow(V)
n_epimarks <- ncol(V)
nIter <- 3   # From a run of 100 iters, RSS doesn't change after 25-30 iterations
n_signatures <- 2   # stomach, skin and background
s <- 1e-5      # span needed for convergence

# Initialize W and load H
W <- matrix(runif(n_bins*n_signatures, 0, 1), nrow = n_bins, ncol = n_signatures)
H <- matrix(runif(n_signatures*n_epimarks, 0, 1), nrow = n_signatures, ncol = n_epimarks)

# dim(W) 10 x 3
# dim(H) 3 x 6

MM.RSS <- c()
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
    
    # # Keep track of the convergence condition evolution
    # if (iter %% 10 == 0) {cat("Iteration ", iter, ", log RSS decrease: ", 
    #                           log(MM.RSS[iter - 1]) - log(MM.RSS[iter]), "\n", 
    #             sep = '')}
    
    # Convergence conditioned on log(RSS) increase
    if ( (iter > 1) && ((log(MM.RSS[iter - 1]) - log(MM.RSS[iter])) < s) ) {break}
    
}

ggplot(NULL) + 
    geom_line(mapping = aes(x = 1:10, W[1:length(W[,1]),1], color = "Signature 1")) + 
    geom_line(mapping = aes(x = 1:10, W[1:length(W[,1]),2], color = "Signature 2")) + 
    theme_minimal() +
    xlab("Patient") + ylab('Mutation Load') +
    ggtitle('Mutation Load for all patients towards each signature')