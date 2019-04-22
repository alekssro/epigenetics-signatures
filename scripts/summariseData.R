# #Resampling 1 set of counts using a multinomial distribution (here 3 categories)
# 
# ncounts <- c(100, 230, 45)
# 
# tot <- sum(ncounts)
# 
# vecp <- ncounts/tot
# 
# rmultinom(n = 1,size = tot,prob = vecp)

V <- read.csv("results/genomic_survey/HepG2/V_matrix_transposed.csv", header = T)
chrmsL <- read.table("data/chromInfo.txt")

# assign chrm names and rebuild bin names
chrmsNames <- rep("", times=length(V[,1]))
start <- 1
for (i in 1:length(chrmsL[,1])) {
    nBins <- ceiling(chrmsL[i, 2]/200.0)
    # print(c(chrmsL[i, 2], nBins))
    chr <- rep(as.character(chrmsL[i, 1]), nBins)
    
    end <- start + nBins - 1
    print(c(end - start, nBins, length(chr)))
    chrmsNames[start:end] <- chr
    start <- end + 1
}

V$chr <- chrmsNames

tot_counts <- sum(V)

row_sums <- apply(V, 1, sum)

plot(V[,1])
