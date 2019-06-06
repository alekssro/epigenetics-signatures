%----------------------------------------------------------------------------------------
%	SECTION 1
%----------------------------------------------------------------------------------------

\section{NMF}

Non-Negative Matrix Factorization is the statistical framework in which this analysis is based on. Given a non-negative matrix $V$, NMF is an unsupervised learning method which tries to find non-negative matrix factors $W$ and $H$ such that $$ V \approx W \times H $$. The epigenetic data used, described in a later section, fulfills the precondition of non-negativity of the data. The \(V\) matrix used in the analysis is composed by 200-bp bins as rows and the epigenetic marks as columns. We are investigating the possibility of finding \(k\) signatures which will summarize combinatorial patterns of the data to the epigenetic marks by means of the \(H\) matrix and to the genome-wide bins by means of the \(W\) matrix. 

\medskip

In order to perform the analysis, the \(V\) matrix was loaded into R environment and concretely the Rpackage \textit{NMF} was used \cite{Gaujoux2010}. The package used was select for consistency with previous similar studies \cite{Gandolfi2017} and seeing that attempts to implement the algorithm resulted in more time-consuming alternatives. Rpackage \textit{NMF} offers a framework with several NMF algorithms, including the ones explained in a previous section, from which \texttt{brunet} option was chosen. This algorithm is based on Kullback-Leibler divergence and matches the one described in \cite{Lee2001} and used in \cite{Brunet2004}, enhanced to avoid numerical overflow.

%-----------------------------------
%	SUBSECTION 1
%-----------------------------------
\subsection{Algorithm}

Passed on a non-negative matrix \(V\) of \(m \times n\) size and a chosen number of \(k\) signatures, \texttt{brunet} implementation will find the approximation \(V \approx WH\). 

\begin{enumerate}

​	\item First, both \(W\) and \(H\) matrices are randomly initialized.

​	\item Then, every row in \(W\) is updated according to the correspondent multiplicative update rules mentioned in (1.7).

​	\item Every column in \(H\) is updated according to (1.6).

​	\item Repeat steps 2 and 3 until the default convergence criteria is reached.

\end{enumerate}

The stopping or convergence criteria for the NMF algorithm can be based on a fixed number of iterations or the invariability of the target values (\(( \vert \vert WH \vert \vert  )_t = ( \vert \vert WH \vert \vert  )_{t+1}\)). Since no fixed number of iterations was assumed, the stopping criteria for the analysis was the invariability of the \(WH\) matrix multiplication, which means there were different number of iterations when NMF was applied to the alternative tissues, varying between 350-600 iterations. The implementation in \textit{NMF} package includes parallel computations to speed up the process. 

%-----------------------------------
%	SUBSECTION 2
%-----------------------------------

\subsection{Choosing number of signatures}

<https://www.academia.edu/238621/Non-negative_Matrix_Factorization_Assessing_Methods_for_Evaluating_the_Number_of_Components_and_the_Effect_of_Normalization_Thereon>

**Model Selection.** For any rank *k*, the NMF algorithm groups the samples into clusters. The key issue is to tell whether a given rank *k* decomposes the samples into “meaningful” clusters. For this purpose, we developed an approach to model selection that exploits the stochastic nature of the NMF algorithm. It is based on our group's previous work on consensus clustering ([11](https://www.pnas.org/content/101/12/4164.long#ref-11)) but adds a quantitative evaluation for robustness of the decomposition.

The NMF algorithm may or may not converge to the same solution on each run, depending on the random initial conditions. If a clustering into *k* classes is strong, we would expect that sample assignment to clusters would vary little from run to run. (Note that sample assignment depends only on the relative values in each column of *H*.)

For each run, the sample assignment can be defined by a connectivity matrix *C* of size *M* × *M*, with entry *cij* = 1 if samples *i*and *j* belong to the same cluster, and *cij* = 0 if they belong to different clusters. We can then compute the consensus matrix, C̄, defined as the average connectivity matrix over many clustering runs. (We select the number of runs by continuing until C̄ appears to stabilize; we typically find that 20–100 runs suffice in the applications below.) The entries of C̄ range from 0 to 1 and reflect the probability that samples *i* and *j* cluster together. If a clustering is stable, we would expect that *C* will tend not to vary among runs, and that the entries of C̄ will be close to 0 or 1. The dispersion between 0 and 1 thus measures the reproducibility of the class assignments with respect to random initial conditions. By using the off-diagonal entries of C̄ as a measure of similarity among samples, we can use average linkage HC to reorder the samples and thus the rows and columns of C̄.

We then evaluate the stability of clustering associated with a given rank *k*. Although visual inspection of the reordered matrix C̄ can provide substantial insight (see [Fig. 3](https://www.pnas.org/content/101/12/4164.long#F2)), it is important to have quantitative measure of stability for each value of *k*. We propose a measure based on the cophenetic correlation coefficient, ρ*k*(C̄), which indicates the dispersion of the consensus matrix C̄. ρ*k* is computed as the Pearson correlation of two distance matrices: the first, I-C̄, is the distance between samples induced by the consensus matrix, and the second is the distance between samples induced by the linkage used in the reordering of C̄. In a perfect consensus matrix (all entries = 0 or 1), the cophenetic correlation coefficient equals 1. When the entries are scattered between 0 and 1, the cophenetic correlation coefficient is <1. We observe how ρ*k* changes as *k* increases. W

Explain different ways to choose N, the way used, the N used and why.