#!/bin/bash

####################################################################################################
### Script to prepare the BAM or BED files into processed BED files ready for further analysis #####
####################################################################################################
## Script done based on the processing pipeline used in “A computational approach for the functional
## classification of the epigenome” by Francesco Gandolfi and Anna Tramontano
## (https://doi.org/10.1186/s13072-017-0131-7)
####################################################################################################

# Process raw BAM files into filtered BED files
# 	- duplicates removal
# 	- quality alignment filtering
# 	- 200bp-tag extension


#################################################
## phase 1: prepare bed files from .bam files####
#################################################

# keep only uniquely mapped reads in .bam file;;
# remove duplicate reads from .bam file;;
# sort and index .bam file;;
# convert each .bam into .bed file;;


# 0. define root directorties
PROJECT_DIR=`pwd`
DATA_DIR=${PROJECT_DIR}/data;
DATASHEET_SAMPLES_FILE=$1	# file containing info for each sample
# cd ${DATA_DIR} || exit;		# change to data directory or exit


# 1 read the datasheet-tab-delimited file line by line
while read LINE; do			# Read ${DATASHEET_SAMPLES_FILE}

	# get values in single fields;
	CELL_LINE=$(echo ${LINE} | awk '{split($0,a," ");print a[1]}');
	ASSAY_TYPE=$(echo ${LINE} | awk '{split($0,a," ");print a[2]}');
	ASSAY=$(echo ${LINE} | awk '{split($0,a," ");print a[3]}');
	SAMPLE_ID=$(echo ${LINE} | awk '{split($0,a," ");print a[4]}');
	ALIGN_FILENAME=$(echo ${LINE} | awk '{split($0,a," ");print a[5]}');
	# ALIGN_FILETYPE=$(echo ${LINE} | awk '{split($0,a," ");print a[6]}');

	# 1.2 check if the alignment file is available or not
	if [ -z "$ALIGN_FILENAME" ]
	then
	    continue
	fi

	#progress
	echo "Working on sample ${SAMPLE_ID} assay ${ASSAY} cell-line ${CELL_LINE}"

	# @bam files or .bed not-processed@:
	echo "  raw bam file $ALIGN_FILENAME";
	CELL_LINE_DIR=${DATA_DIR}/${CELL_LINE};
	SAMPLE_DIR=${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID};
	# get the path to the alignment file
	ALIGN_FILEPATH=${SAMPLE_DIR}/${ALIGN_FILENAME};

	# take note of every bam file contained;
	cd ${SAMPLE_DIR} || exit;

	# If there is a .bed file, convert it into .bam
	bedcount=`ls -1 *.bed 2>/dev/null | wc -l`
	if [ "$bedcount" -gt 1 ]
	then
	    echo "  Converting bed to bam..";
	    # define the output file name just generated (same name, change extension .bed to .bam )
	    BAM_FILE=${ALIGN_FILENAME/".bed"/".bam"};
	    bedtools bedtobam -i ${ALIGN_FILEPATH} -g ${DATA_DIR}/chromInfo.txt > ${BAM_FILE};
	else
	    BAM_FILE=${ALIGN_FILENAME};
	fi

	# Retrieve only the sample name:
	BAM_NAME=${BAM_FILE/".bam"/""};

	# sort and index bam file;
	echo "  Sort and index ${BAM_NAME}";
    samtools sort ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_FILE} -o ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.sorted.bam
    samtools index ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.sorted.bam

	# remove duplicate reads
	echo "  Mark and remove duplicates in ${BAM_NAME}";
	samtools rmdup ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.sorted.bam ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.nodup.bam

    # take only reliable alignments ( only matches with 99% probability of being real matches are kept )
	echo "  Make a new bam-index and keep only high-quality alignments in ${BAM_NAME}";
	samtools index ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.nodup.bam
	samtools view -q 20 -b ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.nodup.bam > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.filtered.bam
	# Alternatively, if Bowtie was used, you select uniquematch reads:
	# /usr/bin/samtools view -f 0x2 -q 255 -b tumor1.bam -o unique.bam

	# convert .bam to .bed (bedtools bamtobed);
	echo "  Convert ${BAM_NAME} into .bed format";
	bamToBed -i ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.filtered.bam  > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed

	# For TFs and HistMod signals apply tag extension of 200bp
	# 	- split by strand +/- ;
	# 	- apply tag extension of 200bp;
	# 	- merge extended bed files

	cd ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID} || exit;

	if [ "$ASSAY_TYPE" != "openchromatin" ]
	then
	    echo "  Apply tag extension on ${BAM_NAME} and saving into processed bed file";
	    awk -F "\t" '$6=="+" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.pos.bed
	    awk -F "\t" '$6=="-" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.neg.bed
	    awk -F "\t" '{print $1,$2,$2+199,$4,$5,$6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.pos.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.extpos.bed
	    awk -F "\t" '{print $1,$3-199,$3,$4,$5,$6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.neg.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.extneg.bed
	    cat ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.extneg.bed ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.extpos.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed

		echo "  Sort by chromosome,pos $BAM_NAME";
	    sort -k1,1 -k2,2n ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed -o ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed
	    # remove intermediate files
	    rm ${BAM_NAME}.ext*; rm ${BAM_NAME}.pos.bed;rm ${BAM_NAME}.neg.bed;

	else
		# otherwise sort the .bed file only
	    echo "  Sort by strand,chromosome,pos $BAM_NAME";
	    sort -k6,6 -k1,1 -k2,2n ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed -o ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed
	fi

	# introduce tab-delimited separator
	awk '$1=$1' FS=" " OFS="\t" ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/out.bed

	# Remove anomalous position
	awk -F "\t" '$2 >= 0 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/out.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed

	# remove intermediate files
	rm ${BAM_NAME}.nodup*; rm ${BAM_NAME}.sorted*; rm ${BAM_NAME}.filtered.bam; rm out.bed;


done<${DATASHEET_SAMPLES_FILE};


###################
### End script ####
##################
