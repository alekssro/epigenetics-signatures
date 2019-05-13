#!/usr/bin/Rscript

###################################################################################################
# Script used to perform NMF analysis in the filtered V matrix data                             ###
###################################################################################################
# Done by Alejandro Roca (alekss.ro@gmail.com)                                                  ###
###################################################################################################

suppressPackageStartupMessages(require(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(require(NMF, quietly = T, warn.conflicts = F))

# Get arguments from terminal call
args <- commandArgs(trailingOnly = T)       # trailingOnly = T, gets only arguments not the call
infile <- args[1]
outfile <- args[2]
n <- args[3]