%----------------------------------------------------------------------------------------
%	SECTION 3
%----------------------------------------------------------------------------------------

\section{Objectives}

The main aim of this master thesis was to develop a computational pipeline capable of finding combinatorial patterns in the epigenomic data.

The absence of mapping information for the different epigenetic marks implied that a large proportion of the work consisted of processing the data. First, tools were developed for downloading epigenetic raw data and preparing the data for the analysis. The final employed data contained per-bin-counts for every epigenetic modification included in the analysis. Data from three cell types was used, hence three per-bin-counts matrices were produced.

With the correct data, NMF algorithm was used on the assumption that we are able to find combinatorial patters or ``signatures''. The underlying idea of the present analysis is the ability of these signatures to summarize the epigenetic modifications states in two angles: (1) the interaction between the several epigenetic modification types, comprised by the \(H\) matrix, and (2) the interaction between the genomic regions, held by the \(W\) matrix. This been fulfilled, we can seek to relate the signatures with biological processes, as well as compare their differences among tissues.

Over the present analysis, aspects from previous related analysis have been adopted, reproducing the arrangement of the data information for the processing step \cite{Gandolfi2017}, or segmenting the genome into 200-bp bins. Results show that there are similar chromatin signatures found for these tissues, suggesting a possible generalization of the chromatin profiles in normal and cancer cells. Moreover, some differences were found and are discussed in the upcoming sections.