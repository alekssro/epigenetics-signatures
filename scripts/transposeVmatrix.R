require(data.table)

# Get arguments from terminal call
args <- commandArgs(trailingOnly = T)       # trailingOnly = T, gets only arguments not the call

V <- fread(args[1], sep = ",", header = F)

cat("Transposing V matrix...\n")
if (nrow(V) < ncol(V)) {
    V_mat <- t(V)
    
    V_matrix <- matrix(as.numeric(V_mat[-1, ]), nrow = nrow(V_mat) - 1, ncol = ncol(V_mat))
    
    colnames(V_matrix) <- V_mat[1, ]
    
    write.csv(V_matrix, file = args[2], quote = F, row.names = F)
}



