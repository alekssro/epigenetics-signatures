#!/usr/bin/Rscript

###################################################################################################
# Script used to perform NMF analysis in the filtered V matrix data                             ###
#
###################################################################################################
# Done by Alejandro Roca (alekss.ro@gmail.com)                                                  ###
###################################################################################################

suppressPackageStartupMessages(require(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(NMF, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(reshape2, quietly = T, warn.conflicts = F))

# Get arguments from terminal call
args <- commandArgs(trailingOnly = T)       # trailingOnly = T, gets only arguments not the call
infile <- args[1]
outfile <- args[2]               

# Load V matrix and chromosome info
cat("Loading V matrix...\n")
filteredV <- read.csv(infile, header = T)
numeric_cols <- 1:(ncol(V) - 2)
V <- as.matrix(filteredV[, numeric_cols])

epimarks <- colnames(filteredV[, numeric_cols])

n_bins <- nrow(V)
n_epimarks <- length(epimarks)
n_reps <- 1

res_reps <- list()
for (i in 1:n_reps) {

    randomV <- randomize(V)

    resids <- matrix(0, ncol = 2, nrow = n_epimarks)
    colnames(resids) <- c("real", "random")
    for (n in 1:n_epimarks) {

        nmf_res_real <- nmf(V, n)
        nmf_res_rand <- nmf(randomV, n)


        print(paste("N:", n, "--> Residual error: [REAL]", nmf_res_real@residuals, "[RANDOM]", nmf_res_rand@residuals))
        resids[n, ] <- c(nmf_res_real@residuals, nmf_res_rand@residuals)

    }
    res_reps[[i]] <- as.data.frame(resids)
}

all_reps <- do.call(rbind, res_reps)
all_reps$n <- rep(1:n_epimarks, length.out=nrow(all_reps))
all_reps$rep <- rep(1:n_reps, each=n_epimarks)

mean_resids_real <- {all_reps %>% group_by(n) %>% summarise(mean = mean(real))}$mean
mean_resids_random <- {all_reps %>% group_by(n) %>% summarise(mean = mean(random))}$mean

p <- ggplot(NULL) +
    geom_line(mapping = aes(x = 1:n_epimarks, y = mean_resids_real, color = "Real Data")) +
    geom_line(mapping = aes(x = 1:n_epimarks, y = mean_resids_random, color = "Random Data")) +
    geom_boxplot(data = all_reps, mapping = aes(x = n, y = real, group=n)) +
    geom_boxplot(data = all_reps, mapping = aes(x = n, y = random, group=n)) +
    xlab("N") + ylab("Reconstruction Error") +
    theme_minimal()

ggsave(outfile, plot = p, height=5, width=7, units='in', dpi=600)

#
# ggplot(NULL) +
#     geom_line(mapping = aes(x = 1:n_bins, W[1:length(W[,1]),1], color = "Signature 1")) +
#     geom_line(mapping = aes(x = 1:n_bins, W[1:length(W[,1]),2], color = "Signature 2")) +
#     geom_line(mapping = aes(x = 1:n_bins, W[1:length(W[,1]),3], color = "Signature 3")) +
#     theme_minimal() +
#     xlab("Patient") + ylab('Mutation Load') +
#     ggtitle('Mutation Load for all patients towards each signature')
