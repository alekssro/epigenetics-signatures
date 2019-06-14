%----------------------------------------------------------------------------------------
%	SECTION 2
%----------------------------------------------------------------------------------------

\section{NMF Signatures}

The NMF produces two factorized matrices, \(H\) and \(W\). Among them, \(H\) matrix is the one that relates the epigenetic marks to the signatures, giving a \(7 \times 11\) array. From this matrix, we can study the combinatorial relations between the epigenetic modifications. One important aspect to bear in mind is that the signatures which yield from the method do not come off in the same order. This is due to the random initialization of the matrices. Therefore, referring to particular signature numbers in the following lines is just used for comparison purposes. A good way to compare the \(H\) matrix outcome is by analyzing the values in it using a heatmap plot. In order to relate the signatures from the different cell tissues, correlation between their values was used.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.8\textwidth]{Figures/NMF/vert_comb.png}
    \caption[Pre-processing performance by read]{\textbf{Heatmap of the obtained H matrices}. Figure description. Lots of things. Lots of things. Lots of things. Lots of things.Lots of things. Lots of things. Lots of things. Lots of things.Lots of things. Lots of things. Lots of things. Lots of things.}
    \label{fig:H_heatmaps}
\end{figure}

Figure~\ref{fig:H_heatmaps} shows the \(H\) matrix results for the three cell lines analyzed. The arrangement of the signatures was decided according in order to correspond to the one shown in \cite{Gandolfi2017}. As we can see, there are numerous coincidences between the heatmaps. Assuming we related the signatures correctly, we can adopt the same `genomic labels' for the signatures in the following order: (1) `Active Promoter’, (2) `Repressed Chromatin’,(3)  `Transcription Initiation’, (4) `Repressed Regulatory Regions', (5) `Gene Body Transcription’,  (6) `Enhancer Regions’ and (7) `Regulatory Elements’.