%----------------------------------------------------------------------------------------
%	SECTION 3
%----------------------------------------------------------------------------------------

\section{Pipeline}

In order to integrate the different scripts used in the project, a bash command line tool was created. The command line requires of a conda environment which can be set up using the \texttt{createCondaEnv.sh} script. This command line includes options to follow the next pipeline for the analysis:

\begin{enumerate}
    \item Download the data (option \texttt{-d} or \texttt{--download}).
    \item Prepare BED alignment files (option \texttt{-p} or \texttt{--process}).
    \item Create V matrix of counts for marks (cols) by bins (rows) (option \texttt{-c} or \texttt{--counts}).
    \item Filter V matrix bins to remove noise (option \texttt{-f} or \texttt{--filter}).
    \item Choose the optimal `n' number of signatures for NMF (option \texttt{-k} or \texttt{--chooseN}).
    \item NMF analysis (option \texttt{-n} or \texttt{--nmf}).
\end{enumerate}

All the scripts used for the analysis as well as for creating the report are available as a \href{https://bitbucket.org/alekssro/mscthesis/src/master/}{Bitbucket repository}.

\subsection{Download Data}

Command: \texttt{scripts/epigeNMF.sh -d [Datasets.tsv]}

\smallskip

Makes use of \texttt{downloadData.sh} script. It requires the file with the datasets info as explained in the Data section \ref{data}. For each sample, the dataset will be downloaded and saved in the \path{data} directory using the next file tree: \path{data/[cell_line]/[epigenetic_mark_name]/[sample_id]/}.

\subsection{Processing BAM files}

Command: \texttt{scripts/epigeNMF.sh -p [Datasets.tsv]}

\smallskip

Makes use of \texttt{prepareBedAlignment.sh} script, based on the processing pipeline used in \cite{Gandolfi2017}, following the scheme and updating outdated methods. The script process raw BAM files into filtered BED files by (1) removing duplicate reads, (2) filter reads by quality, (3) add a 200-bp tag extension in both 3' ends to transcription factors and histone modification signals (resembling half of the average ChIP-seq fragments) and (4) convert BAM files into BED files. It outputs processed BED files to the same path as the BAM file.

\subsection{Generate V matrix}

Command: \texttt{scripts/epigeNMF.sh -c [cell-lines]}

\smallskip

Calls \texttt{bedToNormCounts.sh} script, based on the processing pipeline used in \cite{Gandolfi2017}, following the scheme and updating outdated methods. Takes the cell lines desired to generate the \(V\) matrix from their samples. It requires the datasets info file, human genome segmented into 200-bp fragments as well as a uniqueness mappability track file (using the reference \textit{Duke Uniqueness Regions}) to be in the data directory. For each sample in each cell line, the script assigns the reads to a bin in the segmented genome and calls \texttt{bedCountsToV.R} script which adds the columns of the corresponding cell line \(V\) matrix.

\medskip

The bin counts are normalized with an scaling factor based on the uniqueness mappability positions present in each bin and the total number of reads in the sample. The technical replicates available for a epigenetic sample were combined by taking the mean of the normalized counts for each bin. The result of the execution is a \(V\) matrix csv file for each of the cell lines of interest, saved in the file tree directory \path{results/genomic_survey/[cell_line]/V_matrix.csv}.

\subsection{Filter V matrix}

Command: \texttt{scripts/epigeNMF.sh -f [input-V.csv] [output-V.csv]}

\smallskip

Uses \texttt{summariseFilterData.R} script. Takes as arguments the input and output directories, being the input an unfiltered \(V\) matrix in CSV format and the output a \(V\) filtered matrix also in CSV format. The bins are filtered taking only those for which at least one of the epigenetic marks fall into the top 2.5\% percentile. The threshold was set to the 0.975 quantile in order to filter out the counts considered as noise for our analysis.

\subsection{Choose number of signatures}

Command: \texttt{scripts/epigeNMF.sh -k [input] [output] [n-repeats]}

\smallskip

With \texttt{chooseN.R} script used to choose the optimal \(k\) number of signatures for the NMF analysis in the filtered \(V\) matrix data. Generates a plot comparing real and random data reconstruction errors (choose \(k\) based on that). It also includes an option to define the number of times NMF and reconstruction error are calculated for each \(k\) value. Outputs a residual error plot as the one described in \ref{chooseK} section.

\subsection{NMF analysis}

Command: \texttt{scripts/epigeNMF.sh -n [input-file] [output-file] [n-signatures]}

\smallskip

Calls \texttt{NMFanalysis.R} with a filtered \(V\) matrix in CSV file as input, and produces several informative plots describing the \(H\) and \(W\) matrices generated, saving them to the output directory. The number of signatures can be passed as an argument, with a default of \(k = 7\). 