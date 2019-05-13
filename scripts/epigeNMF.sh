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
    # -a|--all)
    # shift         # remove first argument ($2 now is $1 if any) (remove flag (-a|--all))
    # if test $# -gt 0; then    # get nº arguments after -l option (should be one -> CSVDIRECTORY)
    #   echo ""
    #   echo "Running all analysis with default values:"
    #   echo ""
    #   export CSVDIRECTORY=$1
    #   $0 -l $CSVDIRECTORY -t $ITERATIONS -m $METHOD $ITERATIONS $N_FEATURES $TUNE_LENGTH $CV_REPEATS
    #   $0 -p $METHOD $CSVFILE
    #   exit 0
    # else
    #   echo "No directory specified for loading CSV data"
    #   exit 1
    # fi
    # ;;

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
    if test $# -eq 2; then    # if nº arguments after -m option > 0
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
    if test $# -eq 3; then    # get nº arguments after -t option (should be one -> ITERATIONS)
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
    if test $# -eq 3; then    # get nº arguments after -t option (should be one -> ITERATIONS)
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

    # -h|--help)
    # echo "$0:    ANTIBIOTIC RESISTANT GENES MODEL"
    # echo ""
    # echo "You should follow the next pipeline for the analysis:"
    # echo "    1) Load the data (option '-l' or '--load')"
    # echo "    2) Find significantly different features (option '-t' or '--test')"
    # echo "    3) Create model and test predictions (option '-m' or '--model')"
    # echo ""
    # echo "OR use '$0 -a [CSVDIRECTORY]' command to run whole pipeline with default parameters"
    # echo ""
    # echo " "
    # echo "$0 [options] [arguments]"
    # echo " "
    # echo "Options:"
    # echo ""
    # echo "  -a, --all                 run all analysis with default values"
    # echo "  -c, --convert             convert fasta files to csv files"
    # echo "  -h, --help                display this help"
    # echo "  -l, --load                load csv files to .RData."
    # echo "                              CSVDIRECTORY should be input as argument"
    # echo "  -m, --model               create model and give test results back."
    # echo "                              Several arguments can be added."
    # echo "                              ITERATIONS should be input as argument"
    # # echo "  -o, --output-dir=DIR      specify a directory to store output in"
    # echo "  -p, --predict                predict antibiotic resistance genes from csv file"
    # echo "  -s, --standardize         standardize fasta files for analysis"
    # echo "  -t, --test                test for significantly different features"
    # echo "                              between resistant and non-resistant genes."
    # echo "                              ITERATIONS should be input as argument"
    # echo ""
    # echo "Arguments:"
    # echo ""
    # echo "  -c | -s option:"
    # echo "      FASTADIRECTORY"
    # echo ""
    # echo "  -l | -a options:"
    # echo "      CSVDIRECTORY            directory with the CSV files"
    # echo ""
    # echo "  -t option:"
    # echo "      ITERATIONS              number of iterations used in resampling"
    # echo ""
    # echo "  -m option:"
    # echo "      METHOD                 statistical method used to fit the model (KNN/GLM)"
    # echo "      ITERATIONS             number of iterations used to fit and test the model"
    # echo "      N_FEATURES             number of (most significant) features used to fit the model"
    # echo "      TUNE_LENGTH            KNN tune length"
    # echo "      CV_REPEATS             number of repeats used in cross validation when method = KNN"
    # echo "  *arguments for -m should be input in the next order:"
    # echo "    -m [METHOD] [ITERATIONS] [N_FEATURES] [TUNE_LENGTH] [CV_REPEATS]"
    # echo ""
    # echo "  -p option:"
    # echo "      METHOD                 statistical method used to fit the model (KNN/GLM)"
    # echo "      CSVFILE                path to CSV file"
    # exit 0
    # ;;

    *)    # other arguments
    echo ""
    echo "Yo no comprendo :("
    echo ""
    echo "Please use '$0 --help' for usage info"
    exit 1
    ;;
  esac

done
