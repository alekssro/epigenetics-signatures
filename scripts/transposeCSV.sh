#!/bin/bash

# Script to transpose a csv file uses more memory, faster
# Input:
#   - First argument: file to be transposed
#   - Second argument: output file
# Warning: reads whole file into memory. Use lowMem if file is too big.

awk '
BEGIN { FS=OFS="," }
{
    for (rowNr=1;rowNr<=NF;rowNr++) {
        cell[rowNr,NR] = $rowNr
    }
    maxRows = (NF > maxRows ? NF : maxRows)
    maxCols = NR
}
END {
    for (rowNr=1;rowNr<=maxRows;rowNr++) {
        for (colNr=1;colNr<=maxCols;colNr++) {
            printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
        }
    }
}' $1 > $1.temp

mv $1.temp $2
