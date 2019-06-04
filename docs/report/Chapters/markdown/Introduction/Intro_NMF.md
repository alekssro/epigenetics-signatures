%----------------------------------------------------------------------------------------
%	SECTION 2
%----------------------------------------------------------------------------------------

\section{Non-negative Matrix Factorization}

Non-negative matrix Factorization reveals as a method to learn parts or components from high-dimensional data. In contrast to other matrix factorization methods such as Principal Component Analysis (PCA) or Vector Quantization (VQ), NMF does not learn holistic but part-based representations. In addition, NMF is constrained to have positive values as only additive combinations are considered \cite{Lee1999}. Basically, NMF consists on an approximation of an $m \times n$ $V$ matrix by the matrix multiplication of $W$ ($m \times k$) and $H$ ($k \times n$). This shall be accomplished by iteratively updating the rows and columns of $W$ and $H$ respectively. The constriction of the data to be positive also circumscribe the use of the method, though it has been applied in astronomy \cite{Ren2017}, language processing \cite{Bertin2010} or image processing \cite{Yang2007} among others. Gene expression, epigenetic modifications or mutation data have also been subject of analysis using NMF, considering the non-negativity of the data.

\medskip

Dealing with high dimensionality data using NMF implies taking certain decisions, starting from the choice of NMF dimensionality (number of $k$ signatures). An important measure used for this task is the reconstruction error, calculated by comparing the original data value with the result of the multiplication of the factorized matrices $W$ and $H$. Some of the studies compare the reconstruction error produced by each defined $k$ using the real data vs the one produced when taking random data, as it is the case in this analysis \cite{Frigyesi2008}. When examining the Residual Sum of Squares (RSS) as a function of $k$, the RSS decreases with increasing $k$ within the original dataset, however, this can be less the case for the random dataset. Ideally, for a given interval of $k$ values, there ought to be an optimal $k$ for which the slopes of the real and random data intersect. That is to say, we intend to find the highest $k$, reducing the reconstruction error and improving the accuracy of the model, before being within the noise. Other studies use the cophenetic correlation coefficient, where the similarity between results from several runs is compared for each $k$. In other words, the stability of the classification is studied choosing the highest $k$ before the robustness of the results drops \cite{Brunet2004}. 

\subsection{NMF and gene expression}

NMF was promptly introduced in bioinformatics as a method of dimensionality reduction of large-scale gene expression data \cite{Kim2003}. Here, functional relationships are yield from the analysis of 300 genome-wide expression measurements of yeast and compared with previous expertise. The 5346 genes analyzed in 300 samples were summarized in 50 patterns from which 12 were annotated with MIPS \cite{Mewes2002} functional categories, based on the frequency which genes from each category appeared in each of the signatures or patterns. In a similar way, NMF has also been applied to categorize tumor subtypes using 4651 human genes in 108 cases \cite{Frigyesi2008}, showing that NMF signatures may correspond specific disease gene expression patterns.

\subsection{NMF and mutations}

Besides finding shared patterns between diseases and gene expression, it has been shown that it is possible to characterize specific mutations produced in varied types of cancer \cite{Ramakrishna2012}. A catalogue of complete genomes for 21 primary breast cancer samples was sequenced and compared to the normal DNA of those same individuals. From these, 183,916 somatically acquired mutations were called and classified into 96 possible trinucleotide variations, composing the $V$ matrix of $183916 \times 21$ cells. NMF was then applied on this dataset and from the results they identified contribution of each of the mutational signatures in each of the patients, which yield similar arrangement for similar cancer types.  

\medskip

Although the number of identified signatures was five, further on this number got reduced to only four signatures \cite{Alexandrov2013}. In such a way, mutational processes operative in cancer genomes were modeled with great accuracy by a combination of this four mutational signatures. This allowed to make predictions based on the signatures and the data, as well as understand more deeply the biological processes involved on the different types of cancer.

\subsection{NMF and epigenetics}

The epigenetic modeling techniques explained before have shown the possibility to find genome-wide chromatin states based on the partition of the DNA. Nevertheless, they make the assumption that a small set of these chromatin states would be sufficient to describe the genomic expression, whereas several hundreds chromatin states have been estimated even with a small set of epigenetic marks used \cite{Ucar2011}. In addition, chromatin modifications tend to be highly correlated, hampering the task of assessing the importance of the chromatin marks and relating them to the biological mechanisms. NMF can be used in order to overcome these downsides by identifying combinatorial patterns of chromatin states. 

\medskip

This application was initially presented as a way to get combinatorial patterns of epigenetic marks from integrated epigenetic data sets \cite{Cieslik2014}. They characterized a small amount of combinatorial patterns, which could be displayed and interpreted, were statistically capable of regression and classification tasks. Each row of the $V$ matrix represents 2 kbp of regions flanking Transcription Start Sites (TSS) and columns represent the epigenetic marks used. In a case study for regression of \textit{Pol2} binding, ten epigenetic marks were used to identify seven quantitative epigenetic patterns. Using the 7 chromatin patterns in a regression model yield a performance of $r^2 = 0.85$, matching the one obtained by using 10 epigenetic marks but solving multicollinearity problems. 

\medskip

Later on, functional classification of the epigenome was performed by adding more epigenetic marks and in this case, segmenting the DNA into 200 bp regions as an attempt to resemble one nucleosome with each 200-bp bin \cite{Gandolfi2017}. NMF was here applied to 13 different epigenetic marks over 833,738 significant bins in human embryonic stem cells, finding 7 epigenetic signals or chromatin profiles. These seven epigenetic signals were then labeled based on the related biological process and then used for a wide range of tasks: study the genomic distribution of the signatures, investigate their recovery power on genomic features, association with the gene expression...

\medskip

Over the present analysis, aspects of these previous analysis have been reproduced, adopting the processing of the data scheme \{Gandolfi2017}, segmenting the genome into 200-bp bins and extending the analysis to three cancer cell types for which eleven common epigenetic marks were available. Results show that there are similar chromatin signatures found for these tissues, suggesting a possible generalization of the chromatin profiles in normal and cancer cells.



In order to achieve this, we need to map epigenetic modification reads (from Bisulfite-seq; DNA methylations, or ChIP-seq; histone modifications) into the human genome, specifically into genes and flanking/regulatory regions, which are the ones of interest. In this way we would obtain the counts of each epigenetic modification for different tissues. Non-negative Matrix Factorization (NMF) reveals as an ideal method for the task of finding combinatorial patterns of epigenetic modifications. We can then study the state of each epigenetic modification type in the defined loci of the tissue. From this information we would obtain the different epigenetic signatures which we will use for association and simulation analysis.

