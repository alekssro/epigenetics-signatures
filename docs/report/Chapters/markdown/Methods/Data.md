%----------------------------------------------------------------------------------------
%	SECTION 3
%----------------------------------------------------------------------------------------

\section{Data} \label{data}

The data used in this study was obtained entirely from the ENCODE project database \cite{Feingold2004}, where multiple epigenetic marks and reference epigenomes are available. Data sets from 11 types of epigenetic marks were collected, including histone modifications (H3K4me3, H3K4me1, H3K27ac, H3K36me3, H3K9me3, H3K27me3 and H3K9ac), chromatin remodeling proteins (EP300 and H2A.Z) and transcription regulation factors (CTCF, POLR2A). Replicate samples for each of the 11 epigenetic marks were used in the three cell lines of study, which includes a human liver cancer cell line (HepG2) \cite{Aden1979}, myelogenous leukemia cell line (K562) \cite{Andersson1979} and cells derived from HeLa cancerous cervical tumor line (HeLa-S3) \cite{Douglas1973,Chen2008}.

\medskip

In order to standardize the input and facilitate the data sets processing, the information for all the samples was arranged in a tab separated file with the following fields as columns:

\begin{enumerate}
    \item \textbf{Cell line type}. Either HepG2, K562 or Hela-S3 for this analysis.
    \item \textbf{Epigenetic modification category}. Where does the modification apply. Either histone modification (`histonemod'), chromatin modulation (`openchromatin') or transcription factor (`TFBS').
    \item \textbf{Epigenetic modification name}. The label of the mark which the sample corresponds to, as mentioned above.
    \item \textbf{Accession ID}. The accession name for the sample in the database.
    \item \textbf{File name}.
    \item \textbf{Processing status}. Either `raw' or `process'. Raw samples need to be pre-processed.
    \item \textbf{Replicate number}. Integer enumerating the various replicates for a particular epigenetic mark.
    \item \textbf{Database name}. In the present case, `ENCODE'.
    \item \textbf{Download link}. Used to download the sample reads file.
\end{enumerate}

As a matter of convenience, `\path{downloadData.sh} \path{[DataInfo.tsv]}' or `\path{epigeNMF.sh} \path{-d}  \path{[DataInfo.tsv]}' can be used to automatically download the data sets into the appropriate directory tree: \path{CELL_LINE/SIGNAL_TRACK/SAMPLE_ID}.