\section{Pipeline}

Explain the command line tool, the pipeline followed for the analysis and explain functioning of each script

You should follow the next pipeline for the analysis:

\begin{enumerate}
    \item Download the data (option \path{-d} or \path{--download})
    \item Prepare BED alignment files (option \path{-p} or \path{--process})
    \item Create V matrix of counts for marks (cols) by bins (rows) (option \path{-c} or \path{--counts})
    \item Filter V matrix bins to remove noise (option \path{-f} or \path{--filter})
    \item Choose the optimal `n' number of signatures for NMF (option \path{-k} or \path{--chooseN})
    \item NMF analysis (option \path{-n} or \path{--nmf})
\end{enumerate}