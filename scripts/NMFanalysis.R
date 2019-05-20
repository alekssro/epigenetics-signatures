#!/usr/bin/Rscript

###################################################################################################
# Script used to perform NMF analysis in the filtered V matrix data                             ###
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
n <- args[3]
n <- 7

order <- c("H3K27me3", "H3K9me3", "CTCF", "H3K27ac", "POLR2A", "H3K36me3", "H3K4me1", "H2A.Z", "H3K4me3", "H3K9ac", "EP300")

##########################################################
# NFM on A549 V matrix
##########################################################

# infile <- "results/genomic_survey/A549/filteredV.csv"
# 
# V_mat <- read.csv(infile, header = T)
# numeric_cols <- 1:(ncol(V_mat) - 2)
# V <- V_mat[, numeric_cols]
# 
# nmf_A549 <- nmf(V, 7)
# 
# H_A549 <- nmf_A549@fit@H
# W_A549 <- nmf_A549@fit@W

##########################################################
# NFM on Hela-S3 V matrix
##########################################################

infile <- "results/genomic_survey/Hela-S3//filteredV.csv"

V_mat <- read.csv(infile, header = T)
numeric_cols <- 1:(ncol(V_mat) - 2)
V <- V_mat[, numeric_cols]

nmf_HelaS3 <- nmf(V, 7)

H_HelaS3 <- nmf_HelaS3@fit@H
W_HelaS3 <- nmf_HelaS3@fit@W

##########################################################
# NFM on HepG2 V matrix
##########################################################

infile <- "results/genomic_survey/HepG2//filteredV.csv"

V_mat <- read.csv(infile, header = T)
numeric_cols <- 1:(ncol(V_mat) - 2)
V <- V_mat[, numeric_cols]

nmf_HepG2 <- nmf(V, 7)

H_HepG2 <- nmf_HepG2@fit@H
W_HepG2 <- nmf_HepG2@fit@W

##########################################################
# NFM on K562 V matrix
##########################################################

infile <- "results/genomic_survey/K562//filteredV.csv"

V_mat <- read.csv(infile, header = T)
numeric_cols <- 1:(ncol(V_mat) - 2)
V <- V_mat[, numeric_cols]

nmf_K562 <- nmf(V, 7)

H_K562 <- nmf_K562@fit@H
W_K562 <- nmf_K562@fit@W

#########################################################
###     PLOTS                                       #####
#########################################################


###############
# Heatmap   ###
###############

H <- H_HelaS3[, order]

heatmap(H, scale = "column", Rowv = NA, Colv = NA)

cormat <- round(cor(H),2)
head(cormat)

melted_cormat <- melt(cormat)
head(melted_cormat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue")

ggplot(melted_cormat, aes(variable, Name)) + 
    geom_tile(aes(fill = rescale), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue")



#######################################
#   PLOTS FROM PROJECT              ###
#######################################

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
# dim(V) # 21x96
# rownames(W) <- rownames(V) # 21x4
# colnames(H) <- colnames(V) # 4x96

H_plot <- data.frame(matrix(0, nrow = ncol(H), ncol = nrow(H)+2))
colnames(H_plot) <- c("Signature1", "Signature2", "Signature3", "Signature4",
                      "MutType", "MutGroup")
for (i in 1:nrow(H)) {     # Transform to percentages
    H_plot[,i] <- H[i,]/sum(H[i,]) * 100
}

# Associate mutation with color by group
colors <- c(rep("C>A", 16), rep("C>G", 16), rep("C>T", 16),
            rep("T>A", 16), rep("T>C", 16), rep("T>G", 16))
H_plot[,6] <- colors

# Order by mutation type
H_plot[,5] <- reorder(factor(colnames(V)), seq(1, 96))
# H_plot <- H_plot %>% arrange(MutType)

labels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
p1 <- ggplot(H_plot) +
    geom_col(mapping = aes (x = MutType, y = Signature1, fill = MutGroup)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(c(0, 30)) +
    xlab("Mutation Type") + ylab("Mutation type Count") +
    ggtitle("Signature 1") +
    annotate("text", x = seq(from = 8, to = 88, by = 16), y = 25, label = labels,
             colour = gg_color_hue(6)) +
    guides(fill=FALSE)

p2 <- ggplot(H_plot) +
    geom_col(mapping = aes (x = MutType, y = Signature2, fill = MutGroup)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(c(0, 30)) +
    xlab("Mutation Type") + ylab("Mutation type Count") +
    ggtitle("Signature 2") +
    annotate("text", x = seq(from = 8, to = 88, by = 16), y = 25, label = labels,
             colour = gg_color_hue(6)) +
    guides(fill=FALSE)

p3 <- ggplot(H_plot) +
    geom_col(mapping = aes (x = MutType, y = Signature3, fill = MutGroup)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(c(0, 30)) +
    xlab("Mutation Type") + ylab("Mutation type Count") +
    ggtitle("Signature 3") +
    annotate("text", x = seq(from = 8, to = 88, by = 16), y = 25, label = labels,
             colour = gg_color_hue(6)) +
    guides(fill=FALSE)

p4 <- ggplot(H_plot) +
    geom_col(mapping = aes (x = MutType, y = Signature4, fill = MutGroup)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(c(0, 30)) +
    xlab("Mutation Type") + ylab("Mutation type Count") +
    ggtitle("Signature 4") +
    annotate("text", x = seq(from = 8, to = 88, by = 16), y = 25, label = labels,
             colour = gg_color_hue(6)) +
    guides(fill=FALSE)
# 
# p1
# p2
# p3
# p4


# ggplot(NULL) +
#     geom_point(mapping = aes(x = rowSums(V), y = rowSums(W %*% H))) +
#     geom_abline(slope = 1) +
#     ggtitle('Comparison of mutation counts for each patient') + 
#     ylab('Expected Total Mutations (WxH)') + xlab('Observed Total Mutations (V)')
# 
# ggplot(NULL) +
#     geom_point(aes(x = rowSums(V)[1:nrow(V)-1], y = rowSums(W %*% H)[1:nrow(V)-1])) +
#     geom_abline(slope = 1) +
#     ggtitle('Comparison of mutation counts for each patient (excluding patient 21)') + 
#     ylab('Expected Total Mutations (WxH)') + xlab('Observed Total Mutations (V)')


# ggplot(NULL) + 
#     geom_line(mapping = aes(x = 1:20, W[1:length(W[,1])-1,1], color = "Signature 1")) + 
#     geom_line(mapping = aes(x = 1:20, W[1:length(W[,1])-1,2], color = "Signature 2")) + 
#     geom_line(mapping = aes(x = 1:20, W[1:length(W[,1])-1,3], color = "Signature 3")) + 
#     geom_line(mapping = aes(x = 1:20, W[1:length(W[,1])-1,4], color = "Signature 4")) + 
#     theme_minimal() +
#     xlab("Patient") + ylab('Mutation Load') +
#     ggtitle('Mutation Load for all patients towards each signature')