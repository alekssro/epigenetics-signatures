#!/usr/bin/Rscript

###################################################################################################
# Script used to perform NMF analysis in the filtered V matrix data                             ###
###################################################################################################
# Done by Alejandro Roca (alekss.ro@gmail.com)                                                  ###
###################################################################################################

suppressPackageStartupMessages(require(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(NMF, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(reshape2, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(scales, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(karyoploteR, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(Homo.sapiens, quietly = T, warn.conflicts = F))

# Get arguments from terminal call
args <- commandArgs(trailingOnly = T)       # trailingOnly = T, gets only arguments not the call
infile <- args[1]
outfile <- args[2]
n <- args[3]
n <- 7
set.seed(1234)

order <- c("H3K27me3", "H3K9me3", "CTCF", "H3K27ac", "POLR2A", "H3K36me3", "H3K4me1", "H2A.Z", "H3K4me3", "H3K9ac", "EP300")
signature_names <- c("Signature 1", "Signature 2", "Signature 3", "Signature 4",
                                              "Signature 5", "Signature 6", "Signature 7")

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

infile <- "results/genomic_survey/Hela-S3/filteredV.csv"

V_mat <- read.csv(infile, header = T)
numeric_cols <- 1:(ncol(V_mat) - 2)
V <- V_mat[, numeric_cols]

nmf_HelaS3 <- nmf(V, 7)

H_HelaS3 <- nmf_HelaS3@fit@H
W_HelaS3 <- nmf_HelaS3@fit@W

# ##########################################################
# # NFM on HepG2 V matrix
# ##########################################################
# 
# infile <- "results/genomic_survey/HepG2/filteredV.csv"
# 
# V_mat <- read.csv(infile, header = T)
# numeric_cols <- 1:(ncol(V_mat) - 2)
# V <- V_mat[, numeric_cols]
# 
# nmf_HepG2 <- nmf(V, 7)
# 
# H_HepG2 <- nmf_HepG2@fit@H
# W_HepG2 <- nmf_HepG2@fit@W
# 
# ##########################################################
# # NFM on K562 V matrix
# ##########################################################
# 
# infile <- "results/genomic_survey/K562/filteredV.csv"
# 
# V_mat <- read.csv(infile, header = T)
# numeric_cols <- 1:(ncol(V_mat) - 2)
# V <- V_mat[, numeric_cols]
# 
# nmf_K562 <- nmf(V, 7)
# 
# H_K562 <- nmf_K562@fit@H
# W_K562 <- nmf_K562@fit@W

#########################################################
###     PLOTS                                       #####
#########################################################


###############
# Heatmap   ###
###############

# generate_H_heatmap <- function(H, cell_type, epimark_order, sign_order) {
#     
#     signature_names <- c("Signature 1", "Signature 2", "Signature 3", "Signature 4",
#                          "Signature 5", "Signature 6", "Signature 7")
#     
#     H <- H[, epimark_order]
#     rownames(H) <- signature_names
#     
#     H_m <- melt(H[sign_order,], varnames = c("Signature", "Mark"))
#     
#     p <- ggplot(data = H_m, aes(x=Mark, y=Signature, fill=value)) + 
#         geom_tile(aes(fill = value), colour = "white") +
#         ggtitle(paste(cell_type, "H matrix", sep = " ")) +
#         theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1)) +
#         scale_fill_gradient(low = "white", high = "steelblue")
#     
#     p
#     
#     ggsave(paste("plots/NMF/", cell_type, "_H_heatmap.png", sep = ""), p, width = 10, height = 6)
#     
# }
# 
# # H heatmap
# generate_H_heatmap(H = H_HelaS3, cell_type = "Hela-S3", epimark_order = order, sign_order = rev(c(7,4,6,5,3,1,2)))
# generate_H_heatmap(H = H_K562, cell_type = "K562", epimark_order = order, sign_order = rev(c(6,4,1,5,2,7,3)))
# generate_H_heatmap(H = H_HepG2, cell_type = "HepG2", epimark_order = order, sign_order = rev(c(4,5,7,3,1,6,2)))



# V_mat to GenomicRanges
# 

V_mat["end"] <- V_mat$binStart + 200
V_mat["strand"] <- rep("*", nrow(V_mat))
V_mat["max_per_bin"] <- apply(W_HelaS3, 1, which.max) # Get maximum load signature for each bin




gr_HelaS3 <- makeGRangesFromDataFrame(V_mat[,12:16],
                                    seqnames.field = "chr",
                                    start.field = "binStart",
                                    end.field = "end", 
                                    keep.extra.columns = T)

# plot density of the epigenetic modifications
par(oma=c(4,1,1,1), xpd=TRUE)
kp <- plotKaryotype(plot.type = 1, chromosomes = "chr4", bty='L')
kpPlotDensity(kp, data=gr_HelaS3[mcols(gr_HelaS3)$max_per_bin == 1], r0=0, r1 = 0.2, col = "blue")
kpPlotDensity(kp, data=gr_HelaS3[mcols(gr_HelaS3)$max_per_bin == 2], r0=0.2, r1 = 0.4, col = "red")
kpPlotDensity(kp, data=gr_HelaS3[mcols(gr_HelaS3)$max_per_bin == 3], r0=0.4, r1 = 0.6, col = "green")
kpPlotDensity(kp, data=gr_HelaS3[mcols(gr_HelaS3)$max_per_bin == 4], r0=0.6, r1 = 0.8, col = "pink")
kpPlotDensity(kp, data=gr_HelaS3[mcols(gr_HelaS3)$max_per_bin == 5], r0=0.8, r1 = 1, col = "purple")
kpPlotDensity(kp, data=gr_HelaS3[mcols(gr_HelaS3)$max_per_bin == 6], r0=1, r1 = 1.2, col = "grey")
kpPlotDensity(kp, data=gr_HelaS3[mcols(gr_HelaS3)$max_per_bin == 7], r0=1.2, r1 = 1.4, col = "black")
legend("bottom", legend = signature_names, fill = c("blue", "red", "green", "pink", "purple", "grey", "black"),
       xpd = TRUE, horiz = TRUE, inset = c(0.3,-0.5), col = 1:4, cex = 0.7, bty = "n")




### ENRICHMENT STUDY

human_genes_gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
geneid2symbol <- as.data.frame(org.Hs.egSYMBOL)

get_genes_per_signature <- function(gr, human, tissue){
    
    for (i in 1:7) {
        sign_gr <- gr[mcols(gr)$max_per_bin == i]
        countOverlaps(human_genes_gr, sign_gr)
        
    }
    
    overlaps <- subsetByOverlaps(human, gr)
    
    
    
}

sign1_gr <- gr_HelaS3[mcols(gr_HelaS3)$max_per_bin == 1]


overlaps <- subsetByOverlaps(human_genes_gr, sign1)

geneids <- names(countOverlaps(human_genes_gr, sign1_gr)[countOverlaps(human_genes_gr, sign1_gr) > 100])

symbol_id <- as.data.frame(mapIds(org.Hs.eg.db, geneids, 'SYMBOL', 'ENTREZID'))

write.table(symbol_id, file="results/genomic_survey/Hela-S3/enrichment_sign1.txt", quote = F, row.names = F, col.names = F)

present_genes_symbol <- geneid2symbol[mcols(overlaps)$gene_id %in% geneid2symbol$gene_id,]




W <- W_HelaS3
colnames(W) <- c("Signature 1", "Signature 2", "Signature 3", "Signature 4",
                 "Signature 5", "Signature 6", "Signature 7")

W_m <- melt(W, varnames = c("Bin", "Signature"))

ggplot(data = W_m, aes(x=Bin, y=Signature, fill=value)) + 
    geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue")

############
# Lines  ###
############

# Get maximum load signature for each bin
max_per_bin <- apply(W_HelaS3[chr1, ], 1, which.max)
grouped_bins <- matrix(max_per_bin, nrow = ceiling(length(max_per_bin)/10))
max_per_group <- apply(grouped_bins, 1, function(x){
    as.numeric(names(which.max(table(x))))})

p <- ggplot(NULL) + 
    geom_step(mapping = aes(x = 1:length(max_per_group), y = max_per_group, color = "Signature 1")) + 
    # geom_line(mapping = aes(x = 1:nrow(W_HelaS3), W_HelaS3[1:length(W_HelaS3[,1]),2], color = "Signature 2")) + 
    theme_minimal() +
    xlab("Bin Position") + ylab('Predominant Epigenetic Load') +
    xlim(c(0, 10000)) +
    ggtitle('Mutation Load for all patients towards each signature')

p

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