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

Application of NMF to epigenetics,

Explain the use of Gandolfi, 2017 as guidance an explain differences

Gandolfi, Francesco, and Anna Tramontano. ``A computational approach for the functional classification of the epigenome.''