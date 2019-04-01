#!/bin/bash

# 0.1 define root directorties
PROJECT_DIR=`pwd`   # should be run from project folder
DATA_DIR=${PROJECT_DIR}/data;

# 0.3 Sample-datasheet
DATASHEET_SAMPLES_FILE=${DATA_DIR}/DatasetInfoFile.tsv;

while read LINE;do

    CELL_LINE=$(echo ${LINE} | awk '{split($0,a," ");print a[1]}');
    SIGNAL_TRACK=$(echo ${LINE} | awk '{split($0,a," ");print a[3]}');
    SAMPLE_ID=$(echo ${LINE} | awk '{split($0,a," ");print a[4]}');
    DOWNLOAD_LINK=$(echo ${LINE} | awk '{split($0,a," ");print a[9]}');
    echo $CELL_LINE $SIGNAL_TRACK $SAMPLE_ID $DOWNLOAD_LINK

    mkdir -p ${DATA_DIR}/${CELL_LINE}/${SIGNAL_TRACK}/${SAMPLE_ID}

    wget -P ${DATA_DIR}/${CELL_LINE}/${SIGNAL_TRACK}/${SAMPLE_ID} $DOWNLOAD_LINK


done < ${DATASHEET_SAMPLES_FILE};
