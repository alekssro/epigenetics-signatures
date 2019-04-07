#!/bin/bash

# Script to transpose a csv file. Uses less memory, slower.
# Input:
#   - First argument: file to be transposed
#   - Second argument: output file
# Gets the file to be transposed as argument

awk '
BEGIN { FS=OFS="," }
{ printf "%s%s", (FNR>1 ? OFS : ""), $ARGIND }
ENDFILE {
    print ""
    if (ARGIND < NF) {
        ARGV[ARGC] = FILENAME
        ARGC++
    }
}' $1 > $1.temp

mv $1.temp $2
