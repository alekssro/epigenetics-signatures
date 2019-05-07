#!/usr/bin/Rscript

###################################################################################################
# Script used to perform NMF analysis in the filtered V matrix data                             ###
# 
###################################################################################################
# Done by Alejandro Roca (alekss.ro@gmail.com)                                                  ###
###################################################################################################


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