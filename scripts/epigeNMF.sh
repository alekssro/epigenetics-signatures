#!/bin/bash

# command line tool to perform the analysis described in the thesis

 # ./epigeNMF.sh --download [.tsv file path]    # Download data
 # ./epigeNMF.sh --process [.tsv file path]     # Process BAM data into BED alignment files
 # ./epigeNMF.sh --counts                       # Obtain normalized counts from BED files
 # ./epigeNMF.sh --filter [-p]                  # filter V matrix and obtain coverage plots
 # ./epigeNMF.sh --chooseN [n_reps]             # test multiple 'n' values for NMF, rep n_reps times
 # ./epigeNMF.sh --NMF [n]                      # perform NMF analysis with 'n' signatures
 # ./epigeNMF.sh --enrichment                   # enrichment/association study

DATASHEET_SAMPLES_FILE=data/DatasetInfoFile.tsv
N_REPEATS=1
N_SIGNATURES=7

while test $# -gt 0; do   # check arguments one by one until there is none
  case "$1" in    # Perform actions for different cases

    -d|--download)
    shift
    if test $# -gt 0; then    # get nº arguments after -d option (DATASHEET_SAMPLES_FILE)
      echo "Download files defined in $1 into data directory"
      echo ""
      export DATASHEET_SAMPLES_FILE=$1
      ./scripts/downloadData.sh $DATASHEET_SAMPLES_FILE
      shift
      # exit 0
    else
      echo "No .tsv file specified with data info, using $DATASHEET_SAMPLES_FILE"
      echo ""
      ./scripts/downloadData.sh $DATASHEET_SAMPLES_FILE
      exit 0
    fi
    ;;

    -p|--process)
    shift
    if test $# -gt 0; then    # get nº arguments after -p option (DATASHEET_SAMPLES_FILE)
        echo "Processing BAM files in $1 directory and saving processed BED files"
        echo ""
        export DATASHEET_SAMPLES_FILE=$1
        ./scripts/prepareBedAlignment.sh $DATASHEET_SAMPLES_FILE
        shift
    else
        echo "No .tsv file specified with data info, using $DATASHEET_SAMPLES_FILE"
        echo ""
        ./scripts/prepareBedAlignment.sh $DATASHEET_SAMPLES_FILE
    fi
    ;;

    -c|--counts)
    shift         # remove first argument ($2 now is $1 if any) (rm "-l")
    echo "Calculating V matrix of normalized counts from data in 'data/' directory "
    echo "Saving results in 'results/' directory"
    echo ""

    ./scripts/bedToNormCounts.sh
    ;;

    -f|--filter)
    shift         # remove first argument ($2 now is $1 if any) (rm "-l")
    if test $# -eq 2; then
      export IN_FILE=$1
      export OUT_FILE=$2
      shift 2
      echo "Filtering bins in V matrix (remove noise)"
      echo "  INPUT: V matrix path: $IN_FILE"
      echo "  OUTPUT: filtered V matrix path: $OUT_FILE"
      echo ""
      Rscript summariseFilterData.R $IN_FILE $OUT_FILE

    else
      if test $# -eq 4; then
          export IN_FILE=$1
          export OUT_FILE=$2
          export PLOT_DIR=$4
          shift 4
          echo "Filtering bins in V matrix (remove noise) and plotting coverage"
          echo "  INPUT: V matrix path: $IN_FILE"
          echo "  OUTPUT: filtered V matrix path: $OUT_FILE"
          echo "          plots directory: $PLOT_DIR"
          echo ""
          Rscript summariseFilterData.R $IN_FILE $OUT_FILE -p $PLOT_DIR

      fi
      echo "Incorrect number of arguments specified: 2 or 4 arguments expected"
      echo ""
      echo "Usage: $0 INPUT_FILE OUTPUT_FILE [-p PLOTS_DIR]"
      exit 1
    fi
    ;;

    -k|--chooseN)
    shift   # Remove flag (-k|--chooseN)
    if test $# -eq 3; then
      echo "Test multiple 'n' values for NMF in order to choose the optimal one from the plot:"
      echo "  Number of repetitions: $3"
      echo ""
      export IN_FILE=$1
      export OUT_FILE=$2
      export N_REPEATS=$3
      shift 3
      Rscript chooseN.R $IN_FILE $OUT_FILE $N_REPEATS
      # exit 0
    else
      if test $# -eq 2; then
        echo "Test multiple 'n' values for NMF in order to choose the optimal one from the plot:"
        echo "  Number of repetitions: $N_REPEATS"
        echo ""
        export IN_FILE=$1
        export OUT_FILE=$2
        shift 2
        Rscript chooseN.R $IN_FILE $OUT_FILE $N_REPEATS
      fi
      echo "Incorrect number of arguments specified: needs at least input and output arguments"
      echo ""
      echo "Usage: $0 INPUT_FILE OUTPUT_FILE [N_REPEATS]"
      exit 1
    fi
    ;;

    -n|--nmf)
    shift
    if test $# -eq 3; then
      echo "Performing an NMF analysis on V matrix $1:"
      echo "  Number of signatures: $3"
      echo ""
      export IN_FILE=$1
      export OUT_FILE=$2
      export N_SIGNATURES=$3
      shift 3
      Rscript NMFanalysis.R $IN_FILE $OUT_FILE $N_SIGNATURES
      # exit 0
    else
      if test $# -eq 2; then
        echo "Test multiple 'n' values for NMF in order to choose the optimal one from the plot:"
        echo "  Number of signatures: $N_SIGNATURES"
        echo ""
        export IN_FILE=$1
        export OUT_FILE=$2
        shift 2
        Rscript NMFanalysis.R $IN_FILE $OUT_FILE $N_SIGNATURES
      fi
      echo "Incorrect number of arguments specified: needs at least input and output arguments"
      echo ""
      echo "Usage: $0 INPUT_FILE OUTPUT_FILE [N_SIGNATURES]"
      exit 1
    fi
    ;;

    -h|--help)
    echo "$0:    NMF ANALYSIS OF MULTIPLE EPIGENETIC MARKS COUNTS PER BIN"
    echo ""
    echo "You should follow the next pipeline for the analysis:"
    echo "    1) Download the data (option '-d' or '--download')"
    echo "    2) Prepare BED alignment files (option '-p' or '--process')"
    echo "    3) Create V matrix of counts for marks (cols) by bins (rows) (option '-c' or '--counts')"
    echo "    4) Filter V matrix bins to remove noise (option '-f' or '--filter')"
    echo "    5) Choose the optimal 'n' number of signatures for NMF (option '-k' or '--chooseN')"
    echo "    6) NMF analysis (option '-n' or '--nmf')"
    echo ""
    echo " "
    echo "$0 [options] [arguments]"
    echo " "
    echo "Options:"
    echo ""
    echo "  -c, --counts            obtain normalized counts V matrix from BED files"
    echo "  -d, --download          download data defined in a tsv file with the corresponding format"
    echo "  -f, --filter            filter out bins without informative data (remove noise)"
    echo "  -h, --help              display this help"
    echo "  -k, --chooseN           test multiple 'n' values for NMF, choose 'n' based on generated plot"
    echo "  -n, --nmf               perform NMF on a V matrix"
    echo "                              save results in .RData file and generate plots."
    echo "  -p, --process           process BAM data files into BED alignment files"
    echo "                                  Several arguments can be added."
    echo "                                  ITERATIONS should be input as argument"
    echo ""
    echo "Arguments:"
    echo ""
    echo "  -d | -p options:"
    echo "      SAMPLES_INFO.tsv    tsv file containing the following information for each sample:"
    echo ""
    echo "    CELL_TYPE MARK_TYPE MARK SAMPLE_ID SAMPLE_FILENAME STATUS REPLICATE LAB DOWNLOAD_LINK"
    echo ""
    echo ""
    echo "  -c option:"
    echo "      CELL_TYPES          define cell types of interest, separated by space"
    echo ""
    echo "  -f | -k | -n options:"
    echo "      IN_FILE             input file, should be a V matrix"
    echo "      OUTPUT_FILE         output name"
    echo "                              -f: filtered V matrix csv file name"
    echo "                              -k: reconstruction error png file name"
    echo "                              -n: NMF plots directory name"
    echo " Optional arguments:"
    echo ""
    echo "  -f option:"
    echo "      [-p]                whether to produce coverage plots or not (argument missing)"
    echo "      [PLOTS_DIR]         in case -p argument is specified, directory to save plots in."
    echo ""
    echo "  -k option:"
    echo "      [N_REPS]            number of times each reconstruction error is calculated"
    echo "                              default: 1"
    echo ""
    echo "  -n option:"
    echo "      [N]                 number of signatures used in the NMF analysis"
    echo "                              default: 7"
    echo ""
    exit 0
    ;;

    *)    # other arguments
    echo ""
    echo "Yo no comprendo :("
    echo ""
    echo "Please use '$0 --help' for usage info"
    exit 1
    ;;
  esac

done
