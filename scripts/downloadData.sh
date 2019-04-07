#!/bin/bash
# Script to download the datasets defined in the data/DatasetInfoFile.tsv file

# 1 define root directorties
PROJECT_DIR=`pwd`   # should be run from project folder
DATA_DIR=${PROJECT_DIR}/data;

# 2 Sample-datasheet path
DATASHEET_SAMPLES_FILE=${DATA_DIR}/DatasetInfoFile.tsv;

# Every line contains:
# CELL_LINE MOD_TYPE    SIGNAL_TRACK    SAMPLE_ID   FILENAME    raw REPLICATE_NUM   LAB DOWNLOAD_LINK

# 3 Download every file into a correspondent directory
while read LINE;do

    CELL_LINE=$(echo ${LINE} | awk '{split($0,a," ");print a[1]}');
    SIGNAL_TRACK=$(echo ${LINE} | awk '{split($0,a," ");print a[3]}');
    SAMPLE_ID=$(echo ${LINE} | awk '{split($0,a," ");print a[4]}');
    DOWNLOAD_LINK=$(echo ${LINE} | awk '{split($0,a," ");print a[9]}');
    echo $CELL_LINE $SIGNAL_TRACK $SAMPLE_ID $DOWNLOAD_LINK

    # Create directory for each file (this file tree is needed for later steps) and download
    mkdir -p ${DATA_DIR}/${CELL_LINE}/${SIGNAL_TRACK}/${SAMPLE_ID}

    wget -P ${DATA_DIR}/${CELL_LINE}/${SIGNAL_TRACK}/${SAMPLE_ID} $DOWNLOAD_LINK

done < ${DATASHEET_SAMPLES_FILE};
