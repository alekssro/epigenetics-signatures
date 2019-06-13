%----------------------------------------------------------------------------------------
%	SECTION 4
%----------------------------------------------------------------------------------------

\section{Performance}

During the project, around 150 GB of sequencing data were downloaded, processed and analyzed. Therefore, it was crucial to use the high-performance computing facility provided by GenomeDK at Aarhus University. In particular, nodes with 64 GB of RAM were used. The large amount of data used made the task computationally intensive and some performance test were done.

\medskip

Pre-processing the data takes the major part of the time which could be avoided by saving the processed files into a database. It would also be possible to optimize the algorithm in the process of assigning the read counts for each bin. In any case, the whole process took approximately 8-10 hours of computation for the 87 ChIP-seq alignments used, using the mentioned computational cluster. Once the \(V\) matrices are generated, the NMF method lasted no more than 5 minutes to run.

\subsection{Pre-processing}

As mentioned, this is the most time-consuming step of the process. Figure~\ref{fig:PreProcess} shows the time it took for each read in the input to process the data. With changing numbers of reads (\textit{N}), the pre-processing time exhibits an asymptotic flat curve when the input size is in logarithmic. The time in seconds is divided by the input size (number of ChIP-seq reads), in order to get the time for each read. Consequently, the algorithm seems to perform in an \(\mathcal{O}(n \log{n})\) fashion, which can be caused by the indexing and sorting carried out in this step.

\begin{figure}[h]
    \centering
    \includegraphics[height=8cm]{Figures/performance/performance_per_read.png}
    \caption[Pre-processing performance by read]{\textbf{Pre-processing performance by read}. Figure description. Lots of things. Lots of things. Lots of things. Lots of things.Lots of things. Lots of things. Lots of things. Lots of things.Lots of things. Lots of things. Lots of things. Lots of things.}
    \label{fig:PreProcess}
\end{figure}

\medskip

\subsection{NMF analysis}

The NMF analysis is the cornerstone of the presented master thesis. Considering the high-dimensionality of the data, it is important that the method is able to find the signatures on the data with an acceptable performance. Additionally, we need to choose the \(k\) signatures used, and the way this number is chosen requires to perform NMF for every \(k\) value tested. Consequently, here we performed tests both to investigate the performance of the NMF with changing number of signatures \(k\), and with different number of bins input. Figure~\ref{fig:NMF}A shows how the different k values affect the running time and Figure~\ref{fig:NMF}B represents the performance dependent on the input size.

\begin{figure}[h]
    \centering
    \includegraphics[width=\textwidth]{Figures/performance/performance_nmf_join.png}
    \caption[NMF performance]{\textbf{NMF performance}. Figure description. Lots of things. Lots of things. Lots of things. Lots of things.Lots of things. Lots of things. Lots of things. Lots of things.Lots of things. Lots of things. Lots of things. Lots of things.Lots of things. Lots of things. Lots of things. Lots of things.}
    \label{fig:NMF}
\end{figure}

\medskip

Both plots suggest a linear time complexity behavior \(\mathcal{O}(n)\). It is important to take into account the stopping criteria used for the NMF, as it is based in convergence. The convergence of the algorithm can vary slightly from one run to another and it is revealed as some points lay out of the standard error boundary. However, the running times show that even with the full data \((\sim800,000)\) and using a personal computer the task run in short time.