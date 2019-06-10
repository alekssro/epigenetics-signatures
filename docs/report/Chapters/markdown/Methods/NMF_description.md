%----------------------------------------------------------------------------------------
%	SECTION 1
%----------------------------------------------------------------------------------------

\section{NMF}

Non-Negative Matrix Factorization is the statistical framework in which this analysis is based on. Given a non-negative matrix \(V\), NMF is an unsupervised learning method which tries to find non-negative matrix factors \(W\) and \(H\) such that

\begin{equation}
     V \approx W \times H.
\end{equation}

The epigenetic data used, described in a later section, fulfills the precondition of non-negativity of the data. The \(V\) matrix used in the analysis is composed by 200-bp bins as rows and the epigenetic marks as columns. We are investigating the possibility of finding \(k\) signatures which will summarize combinatorial patterns of the data to the epigenetic marks by means of the \(H\) matrix and to the genome-wide bins by means of the \(W\) matrix.

\medskip

In order to perform the analysis, the \(V\) matrix was loaded into R environment and concretely the Rpackage \textit{NMF} was used \cite{Gaujoux2010}. The package used was select for consistency with previous similar studies \cite{Gandolfi2017} and seeing that attempts to implement the algorithm resulted in more time-consuming alternatives. Rpackage \textit{NMF} offers a framework with several NMF algorithms, including the ones explained in a previous section, from which \texttt{brunet} option was chosen. This algorithm is based on Kullback-Leibler divergence and matches the one described in \cite{Lee2001} and used in \cite{Brunet2004}, enhanced to avoid numerical overflow.

%-----------------------------------
%	SUBSECTION 1
%-----------------------------------
\subsection{Algorithm}

Passed on a non-negative matrix \(V\) of \(m \times n\) size and a chosen number of \(k\) signatures, \texttt{brunet} implementation will find the approximation \(V \approx WH\).

\begin{enumerate}
    \item First, both \(W\) and \(H\) matrices are randomly initialized.
    \item Then, every row in \(W\) is updated according to the correspondent multiplicative update rules mentioned in (1.7).
    \item Every column in \(H\) is updated according to (1.6).
    \item Repeat steps 2 and 3 until the default convergence criteria is reached.
\end{enumerate}

The stopping or convergence criteria for the NMF algorithm can be based on a fixed number of iterations or the invariability of the target values \[( \vert \vert WH \vert \vert  )_t = ( \vert \vert WH \vert \vert  )_{t+1}.\] Since no fixed number of iterations was assumed, the stopping criteria for the analysis was the invariability of the \(WH\) matrix multiplication, which means there were different number of iterations when NMF was applied to the alternative tissues, varying between 350-600 iterations. The implementation in \textit{NMF} package includes parallel computations to speed up the process.

%-----------------------------------
%	SUBSECTION 2
%-----------------------------------

\subsection{Choosing number of signatures}\label{chooseK}

In order to apply NMF to the data we foremost need to choose a suitable \(k\) number of signatures, also called rank. The number of clusters defined will largely influence the results and the explanation of them, therefore it is highly important to find the optimal number to produce ``meaningful'' clusters \cite{Brunet2004}. As explained in previous sections, there is no standard procedure to find the best \(k\) number. Therefore, an additional matrix was created with randomly permuted values from the original \(V\) matrix, which we will refer to as ``random'' \(V_R\) matrix. The \(V_R\) matrix is composed of the same values as \(V\) but column-wise permuted, meaning that the presumed chromatic patterns would not be able to be recognized.

\medskip

To compare the performance of NMF using different \(k\) values, the analysis is performed using values from 2 to 10 and then compared the reconstruction errors of the resultant models. As the NMF factorized matrix can vary from one run to another, the process is repeated 30 times and the results can be seen in [INSERT CROSS-REFERENCE]. The reconstruction error was calculated for comparison by residuals sum of squares (RSS),

\begin{equation}
     \vert \vert V - WH \vert \vert ^2 = \sum_{ij} (V_{ij} - (WH)_{ij})^2
\end{equation}

[INSERT CHOOSE N PLOT]

\medskip

As we can see from the results, for the same \(k\) value, the NMF model produces worst results for the \(V_R\) matrix as well as more variable results. Nevertheless, there is a \(k\) value for which the plots intersect due to the model being over-fitted. We chose \(k = 7\) because it is the value for which the reconstruction error is the lowest possible before modeling the noise in the data. The results found similar in all the three cell types.