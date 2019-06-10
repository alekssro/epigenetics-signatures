# Create performance plots

suppressPackageStartupMessages(require(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(NMF, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(reshape2, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(scales, quietly = T, warn.conflicts = F))

# Load creating V matrix times
times <- read.table("plots/performance/times.txt", header = T)

times <- rbind(times, c(10000000, 2845.473))

times["perRead"] <- times$bedToNorm_time / times$N


p <- ggplot(times, aes(log(N), perRead)) +
    geom_point() +
    ggtitle("Pre-processing performance by read") +
    xlab("Input size (log(N))") + ylab("Time per each read (seconds)") +
    theme_minimal()

ggsave(filename = "plots/performance/performance_per_read.png", plot = p)



# Study nmf performance:

infile <- "results/genomic_survey/Hela-S3/filteredV.csv"

V_mat <- read.csv(infile, header = T)
numeric_cols <- 1:(ncol(V_mat) - 2)
V <- V_mat[, numeric_cols]

nmf_times <- c()
nmf_resids <- c()
for (i in 2:9) {
    
    nmf_HelaS3 <- nmf(V, i)
    
    nmf_times <- c(nmf_times, nmf_HelaS3@runtime[1])
    nmf_resids <- c(nmf_resids, nmf_HelaS3@residuals)
}

nmf_perform <- data.frame(K = 2:9, residuals = nmf_resids, runtimes = nmf_times)

p <- ggplot(nmf_perform, aes(K, runtimes)) +
    geom_point() +
    ggtitle("NMF performance for each k value") +
    xlab("k value") + ylab("Time (seconds)") +
    theme_minimal()

ggsave(filename = "plots/performance/performance_nmf.png", plot = p)
