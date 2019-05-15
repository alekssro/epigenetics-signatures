#!/bin/bash

####################################################################################################
### Script to convert BED files into Normalized Signal tracks 								########
####################################################################################################
## Script done based on the processing pipeline done in “A computational approach for the functional
## classification of the epigenome” by Francesco Gandolfi and Anna Tramontano
## (https://doi.org/10.1186/s13072-017-0131-7)
####################################################################################################

# Convert processed BED files into Normalized Signal tracks.
# - human genome segmentation into 200bp genomic intervals (commented out as it is available)
# - assignment of filtered reads to the segmented genome (obtain read count distribution).
# - estimate both observed and expected read count distribution for each sample in a given
# 		epigenetic mark.
# - Estimation of the normalized coverage signal for each chromatin mark.


# 0.1 define root directorties
PROJECT_DIR=`pwd`   # should be run from project folder
DATA_DIR=${PROJECT_DIR}/data;
SCRIPT_DIR=${PROJECT_DIR}/scripts;

# 0.2 Analysis directory (to store results)
ANALYSIS_DIR=${PROJECT_DIR}/results/genomic_survey;
mkdir -p ${ANALYSIS_DIR};       # make directory if does not exist already

# 0.3 Sample-datasheet
DATASHEET_SAMPLES_FILE=${DATA_DIR}/DatasetInfoFile.tsv;

# 0.4 Uniqueness mappability track
MAPPABILITY_TRACK=${DATA_DIR}/wgEncodeDukeMapabilityUniqueness35bp.uniqueMapRegions.bedGraph

# 0.5 Chrominfo file (Be sure that you already have the k-binned genome file (k=bin size))
# CHROMINFO_FILEPATH=${DATA_DIR}/chromInfo.txt 	# used if needed to generated k-binned genome file
# bedtools makewindows -g ${CHROMINFO_FILEPATH} -w 200 > ${DATA_DIR}/hg19binned.200bp.bed

# 0.6 Genomic bin size
# GENOME_BINSIZE=200;

# 0.7 Define Binned hg19 genome path
BINNED_GENOME=${DATA_DIR}/hg19binned.200bp.bed


##############################
####### LOOP-1 ###############
##############################

# 1. Calculate the total number of informative reads mapped across all the datasets
# Take note of the number of informative reads in that sample to estimate total number of
#   tags over all experiments in a cell line

# Define cell types of interest
# CELL_LINES=(A549 Hela-S3)     # (HepG2 K562 A549 Hela-S3);
CELL_LINES=("$@")     # get cell lines from command line
# Define signal-types
SIGNAL_TRACKS=(POLR2A CTCF H2A.Z EP300 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9ac H3K9me3)

echo "Counting reads:"
for CELL_LINE in "${CELL_LINES[@]}";do

    mkdir -p ${ANALYSIS_DIR}/${CELL_LINE}

    for SIGNAL_TRACK in "${SIGNAL_TRACKS[@]}";do

        # Progress
        echo "  cell line -> $CELL_LINE	signal track -> $SIGNAL_TRACK"

        # Define the output directory (and create if needed)
        OUT_DIR=${ANALYSIS_DIR}/${CELL_LINE}/${SIGNAL_TRACK};
        # echo $OUT_DIR
        mkdir -p ${OUT_DIR};

        # Save the lines of DatasetInfoFile corresponding to the specific epigenetic mark
        awk -F "\t" '$1=="'${CELL_LINE}'" && $3=="'${SIGNAL_TRACK}'" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' ${DATASHEET_SAMPLES_FILE} > ${ANALYSIS_DIR}/ENCODE_project_datasheet_sub.tsv;


        TAG_TOTALS=();  # should be the total of informative reads across all ChIP-seq samples in a given sample group (cell line/signal-track)
        sep=' ';

        while read LINE;do

            SAMPLE_ID=$(echo ${LINE} | awk '{split($0,a," ");print a[4]}');
            # Identify and retrieve processed.bed file and its path
            PROCESSED_FILEDIR=${DATA_DIR}/${CELL_LINE}/${SIGNAL_TRACK}/${SAMPLE_ID};
            if [ -d "${PROCESSED_FILEDIR}" ]; then
                cd ${PROCESSED_FILEDIR} || continue
            else
                continue
            fi
            PROCESSED_FILE="*.processed.bed"
            filelines=$(wc -l ${PROCESSED_FILE})
            N_INFORMATIVE_TAGS=$(echo ${filelines} | awk '{split($0,a," ");print a[1]}');
            echo "      Increase the total count of informative reads across all ChIP-seq assay samples in ${CELL_LINE}: $SIGNAL_TRACK";
            TAG_TOTALS+=($N_INFORMATIVE_TAGS);
            TAG_TOTALS+=($sep);

        done < ${DATASHEET_SAMPLES_FILE};

        # Save read counts into a txt file for all replicates corresponding to a epigen mark
        echo ${TAG_TOTALS[*]} > ${OUT_DIR}/${CELL_LINE}_${SIGNAL_TRACK}_allExp_totals.txt
    done;

done;

echo ""

##############################
####### LOOP-2 ###############
##############################


# 2. Once the cell-type and the signal-type are defined, extract all samples from the cell-type X and
#		signal-type Y from the ENCODE_project_datasheet_samples.tsv then loop row by row on the subset

for CELL_LINE in "${CELL_LINES[@]}";do

    for SIGNAL_TRACK in "${SIGNAL_TRACKS[@]}";do

        # 2.1 Regenerate OUTPUT directory:
        OUT_DIR=${ANALYSIS_DIR}/${CELL_LINE}/${SIGNAL_TRACK};

        # 2.2 Create an empty file with the list of files used as input
        :>${ANALYSIS_DIR}/input_data_files.txt

        # 2.3 Select from the full-datasheet file only row corresponding to CELL_LINE and SIGNAL_TRACK
        awk -F "\t" '$1=="'${CELL_LINE}'" && $3=="'${SIGNAL_TRACK}'" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' ${DATASHEET_SAMPLES_FILE} > ${ANALYSIS_DIR}/ENCODE_project_datasheet_sub.tsv

        # Progress
        echo ""
        echo "Calculating counts per genomic bin for Cell line ${CELL_LINE}, Mark ${SIGNAL_TRACK}."

        # Read through the datasheet subset (foreach replicate Bi belonging to the signal-type Y...)
        while read LINE;do

            # Get values in single fields
            SAMPLE_ID=$(echo ${LINE} | awk '{split($0,a," ");print a[4]}')
            ALIGN_FILENAME=$(echo ${LINE} | awk '{split($0,a," ");print a[5]}')

            # Check if the alignment file is available or not
            if [ -z "$ALIGN_FILENAME" ]
            then
                continue
            fi

            # Reconstruct path and retrieve processed.bed file
            PROCESSED_FILEDIR=${DATA_DIR}/${CELL_LINE}/${SIGNAL_TRACK}/${SAMPLE_ID};
            cd ${PROCESSED_FILEDIR} || continue;
            PROCESSED_FILE="*processed.bed"
            PROCESSED_FILEPATH=${PROCESSED_FILEDIR}/${PROCESSED_FILE};
            SORTED_PROCESS_FILEPATH=${PROCESSED_FILEDIR}/${SAMPLE_ID}.sorted.processed.bed

            # Get the raw counts in each genomic bin
            echo "  Intersect genomic bins with processed bed ${SAMPLE_ID}.processed.bed to get rawCounts"
            bedtools sort -k1,1 -k2,2n -k3,3n ${PROCESSED_FILEPATH} > ${SORTED_PROCESS_FILEPATH}
            if [ -s "$SORTED_PROCESS_FILEPATH" ]
            then
                rm ${PROCESSED_FILEPATH}
                mv ${SORTED_PROCESS_FILEPATH} ${PROCESSED_FILEPATH}
            fi
            #
            bedtools intersect -sorted -a ${BINNED_GENOME} -b ${PROCESSED_FILEDIR} -c > ${OUT_DIR}/${SAMPLE_ID}.rawCounts.bed

            # Extract only informative bins in 'temp.bed' , and replace 'rawCounts' by this one.
            awk -F "\t" '$4 > 0 {print $1 "\t" $2 "\t" $3 "\t" $4}' ${OUT_DIR}/${SAMPLE_ID}.rawCounts.bed > ${OUT_DIR}/temp.bed
            mv ${OUT_DIR}/temp.bed ${OUT_DIR}/${SAMPLE_ID}.rawCounts.bed

            # Put the rawCount file name in a list for normalized signal estimation step
            echo "${OUT_DIR}/${SAMPLE_ID}.rawCounts.bed" >> ${ANALYSIS_DIR}/input_data_files.txt

        done < ${ANALYSIS_DIR}/ENCODE_project_datasheet_sub.tsv;



        # 2.4 Process rawCounts; the R script takes the list of rawCount files just generated.
        #   INPUT: rawCounts files for a given group (cell-lineX/signaltrackY)
        #   OUTPUT: normalizedSignal files for the group (cell-lineX/signalrackY)

        # Arguments:
        # 1) input_data_files.txt file containing path for raw counts files
        # 2) datasheet_samples_file.tsv
        # 3) binned genome path
        # 4) Mappability track file
        # 5) total read counts per cell type.txt
        # 6) output file path for V matrix

        Rscript ${SCRIPT_DIR}/bedCountsToV.R ${ANALYSIS_DIR}/input_data_files.txt ${DATASHEET_SAMPLES_FILE} ${BINNED_GENOME} ${MAPPABILITY_TRACK} ${OUT_DIR}/${CELL_LINE}_${SIGNAL_TRACK}_allExp_totals.txt ${ANALYSIS_DIR}/${CELL_LINE}/V_matrix.csv

    done

    # # Transpose V matrix
    Rscript ${SCRIPT_DIR}/transposeVmatrix.R ${ANALYSIS_DIR}/${CELL_LINE}/V_matrix.csv ${ANALYSIS_DIR}/${CELL_LINE}/V_matrix_transpose.csv
    # ${SCRIPT_DIR}/transposeCSV.lowMem.sh ${ANALYSIS_DIR}/${CELL_LINE}/V_matrix.csv ${ANALYSIS_DIR}/${CELL_LINE}/V_matrix_transpose.csv
    # mv ${ANALYSIS_DIR}/${CELL_LINE}/V_matrix_transpose.csv ${ANALYSIS_DIR}/${CELL_LINE}/V_matrix.csv

done



### End script ####
