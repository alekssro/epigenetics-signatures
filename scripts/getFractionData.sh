#!/bin/bash

####################################################################################################
### Script to create bed files with a fraction of the data included in the original file.      	####
### This bed files are used to check the performance of the analysis.							####
####################################################################################################

# 0. define root directorties
PROJECT_DIR=`pwd`
DATA_DIR=${PROJECT_DIR}/data;
DATASHEET_SAMPLES_FILE=$1	# file containing info for each sample
# FRACTION_DATA=$2
# cd ${DATA_DIR} || exit;		# change to data directory or exit


# 1 read the datasheet-tab-delimited file line by line
while read LINE; do			# Read ${DATASHEET_SAMPLES_FILE}

	# get values in single fields;
	CELL_LINE=$(echo ${LINE} | awk '{split($0,a," ");print a[1]}');
	# ASSAY_TYPE=$(echo ${LINE} | awk '{split($0,a," ");print a[2]}');
	ASSAY=$(echo ${LINE} | awk '{split($0,a," ");print a[3]}');
	SAMPLE_ID=$(echo ${LINE} | awk '{split($0,a," ");print a[4]}');
	ALIGN_FILENAME=$(echo ${LINE} | awk '{split($0,a," ");print a[5]}');


	#progress
	echo "Reducing sample ${SAMPLE_ID} assay ${ASSAY} cell-line ${CELL_LINE}"
	CELL_LINE_DIR=${DATA_DIR}/${CELL_LINE};
	SAMPLE_DIR=${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID};
	ALIGN_FILENAME=${SAMPLE_ID}.processed.bed;
	cd ${SAMPLE_DIR} || exit;

	# get number of reads
	N_READS=`wc -l $ALIGN_FILENAME`
	echo "  original bed file ${SAMPLE_ID}.processed.bed, Nº Reads: ${N_READS}";

	# Save fraction of the file
	echo "Saving 100, 1000, 10000, 100000 reads"
	head -n 100 $ALIGN_FILENAME > ${SAMPLE_ID}.100.processed.bed
	head -n 1000 $ALIGN_FILENAME > ${SAMPLE_ID}.1000.processed.bed
	head -n 10000 $ALIGN_FILENAME > ${SAMPLE_ID}.10000.processed.bed
	head -n 100000 $ALIGN_FILENAME > ${SAMPLE_ID}.100000.processed.bed


done<${DATASHEET_SAMPLES_FILE};


###################
### End script ####
##################
