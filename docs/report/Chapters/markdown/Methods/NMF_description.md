%----------------------------------------------------------------------------------------
%	SECTION 1
%----------------------------------------------------------------------------------------

\section{NMF}

Non-Negative Matrix Factorization is the statistical framework in which this analysis is based on. Given a non-negative matrix $V$, NMF is an unsupervised learning method which tries to find non-negative matrix factors $W$ and $H$ such that $$ V \approx W \times H $$. 

NMF can be applied to the statistical analysis of multivariate data in the following manner. Given a set of multivariate n-dimensional data vectors, the vectors are placed in the columns of an $n \times m$ matrix $V$ where $m$ is the number of examples in the data set. This matrix is the approximately factorized into an $n \times r$ matrix $W$ and an $r \times m$ matrix $H$. Usually $r$ is chosen to be smaller than $n$ or $m$ so that $W$ and $H$ dimensionality is smaller than the original matrix $V$. This results in a compressed version of the original data matrix.

What is the significance of the approximation in Eq. (1)? It can be rewritten column by column as $v \approx Wh$, where $v$ and $h$ are the corresponding columns of $V$ and $H$ respectively. In other words, each data vector $v$ is approximated by a linear combination of the columns of $W$, that is optimized for the linear approximation of the data in $V$. Since relatively few basis vectors are used to represent many data vectors, good approximation can only be achieved if the basis vectors discover structure that is latent in the data. 

Brunet 2004:

**Description of NMF Method.** We consider a data set consisting of the expression levels of *N* genes in *M* samples (which may represent distinct tissues, experiments, or time points). For gene expression studies, the number *N* of genes is typically in the thousands, and the number *M* of experiments is typically <100. The data are represented by an expression matrix *A* of size *N* × *M*, whose rows contain the expression levels of the *N* genes in the *M*samples.

Our goal is to find a small number of metagenes, each defined as a positive linear combination of the *N* genes. We can then approximate the gene expression pattern of samples as positive linear combinations of these metagenes.

Mathematically, this corresponds to factoring matrix *A* into two matrices with positive entries, *A* ∼ *WH*. Matrix *W* has size *N* × *k*, with each of the *k* columns defining a metagene; entry *wij* is the coefficient of gene *i* in metagene *j*. Matrix *H* has size *k* × *M*, with each of the *M* columns representing the metagene expression pattern of the corresponding sample; entry *hij* represents the expression level of metagene *i* in sample *j*. [Fig. 1](https://www.pnas.org/content/101/12/4164.long#F1) shows the simple case corresponding to *k* = 2.

Given a factorization *A* ∼ *WH*, we can use matrix *H* to group the *M*samples into *k* clusters. Each sample is placed into a cluster corresponding to the most highly expressed metagene in the sample; that is, sample *j* is placed in cluster *i* if the *hij* is the largest entry in column *j* ([Fig. 1](https://www.pnas.org/content/101/12/4164.long#F1)).

We note that there is a dual view of decomposition *A* ∼ *WH*, which defines metasamples (rather than metagenes) and clusters the genes (rather than the samples) according to the entries of *W*. We do not focus on this view here, but it is clearly of great interest.

NMF provides a natural way to cluster genes and samples, because it involves factorization into matrices with nonnegative entries. By contrast, principal component analysis provides a simple way to reduce dimensionality but requires that the matrices be orthogonal, which typically requires linear combination of components with arbitrary signs. NMF is more difficult algorithmically because of the nonnegativity requirement but provides a more intuitive decomposition of the data.

**NMF Algorithm.** Given a positive matrix *A* of size *N* × *M* and a desired rank *k*, the NMF algorithm iteratively computes an approximation *A* ∼ *WH*, where *W* and *H* are nonnegative matrices with respective sizes *N* × *k* and *k* × *M*. The method starts by randomly initializing matrices *W* and *H*, which are iteratively updated to minimize a divergence functional. The functional is related to the Poisson likelihood of generating *A* from *W* and *H, D*= Σ*i*,*j Ai*,*j*log(*Ai*,*j*/(*WH*)*i*,*j*) – *Ai*,*j* + (*WH*)*i*,*j*. At each step, *W* and *H* are updated by using the coupled divergence equations ([10](https://www.pnas.org/content/101/12/4164.long#ref-10)):![Math](https://www.pnas.org/sites/default/files/highwire/pnas/101/12/4164/embed/tex-math-1.gif)A simpler version of the NMF update equations that minimizes the norm of the residual ||*A*-*WH*||2 has also been derived in ref. [10](https://www.pnas.org/content/101/12/4164.long#ref-10). When applying the method to a medulloblastoma dataset (see *Results*), where we knew the underlying substructure, we observed that the divergence-based update equations were able to capture a subclass that the norm-based update equations did not. This is why our implementation of NMF uses the divergence form (see Data Sets and software).

**Model Selection.** For any rank *k*, the NMF algorithm groups the samples into clusters. The key issue is to tell whether a given rank *k* decomposes the samples into “meaningful” clusters. For this purpose, we developed an approach to model selection that exploits the stochastic nature of the NMF algorithm. It is based on our group's previous work on consensus clustering ([11](https://www.pnas.org/content/101/12/4164.long#ref-11)) but adds a quantitative evaluation for robustness of the decomposition.

The NMF algorithm may or may not converge to the same solution on each run, depending on the random initial conditions. If a clustering into *k* classes is strong, we would expect that sample assignment to clusters would vary little from run to run. (Note that sample assignment depends only on the relative values in each column of *H*.)

For each run, the sample assignment can be defined by a connectivity matrix *C* of size *M* × *M*, with entry *cij* = 1 if samples *i*and *j* belong to the same cluster, and *cij* = 0 if they belong to different clusters. We can then compute the consensus matrix, C̄, defined as the average connectivity matrix over many clustering runs. (We select the number of runs by continuing until C̄ appears to stabilize; we typically find that 20–100 runs suffice in the applications below.) The entries of C̄ range from 0 to 1 and reflect the probability that samples *i* and *j* cluster together. If a clustering is stable, we would expect that *C* will tend not to vary among runs, and that the entries of C̄ will be close to 0 or 1. The dispersion between 0 and 1 thus measures the reproducibility of the class assignments with respect to random initial conditions. By using the off-diagonal entries of C̄ as a measure of similarity among samples, we can use average linkage HC to reorder the samples and thus the rows and columns of C̄.

We then evaluate the stability of clustering associated with a given rank *k*. Although visual inspection of the reordered matrix C̄ can provide substantial insight (see [Fig. 3](https://www.pnas.org/content/101/12/4164.long#F2)), it is important to have quantitative measure of stability for each value of *k*. We propose a measure based on the cophenetic correlation coefficient, ρ*k*(C̄), which indicates the dispersion of the consensus matrix C̄. ρ*k* is computed as the Pearson correlation of two distance matrices: the first, I-C̄, is the distance between samples induced by the consensus matrix, and the second is the distance between samples induced by the linkage used in the reordering of C̄. In a perfect consensus matrix (all entries = 0 or 1), the cophenetic correlation coefficient equals 1. When the entries are scattered between 0 and 1, the cophenetic correlation coefficient is <1. We observe how ρ*k* changes as *k* increases. We select values of *k*where the magnitude of the cophenetic correlation coefficient begins to fall (see below).

%-----------------------------------
%	SUBSECTION 1
%-----------------------------------
\subsection{Algorithm}

We formally consider algorithms for solving the following problem:
Non-negative matrix factorization (NMF) Given a non-negative matrix
V, find non-negative matrix factors Wand H such that:
V~WH
(1)NMF can be applied to the statistical analysis of multivariate data in the following manner.
Given a set of of multivariate n-dimensional data vectors, the vectors are placed in the
columns of an n x m matrix V where m is the number of examples in the data set. This
matrix is then approximately factorized into an n x r matrix Wand an r x m matrix H.
Usually r is chosen to be smaller than nor m , so that Wand H are smaller than the original
matrix V. This results in a compressed version of the original data matrix.
What is the significance of the approximation in Eq. (1)? It can be rewritten column by
column as v ~ Wh, where v and h are the corresponding columns of V and H. In other
words, each data vector v is approximated by a linear combination of the columns of W,
weighted by the components of h. Therefore W can be regarded as containing a basis
that is optimized for the linear approximation of the data in V. Since relatively few basis
vectors are used to represent many data vectors, good approximation can only be achieved
if the basis vectors discover structure that is latent in the data.
The present submission is not about applications of NMF, but focuses instead on the tech-
nical aspects of finding non-negative matrix factorizations. Of course, other types of ma-
trix factorizations have been extensively studied in numerical linear algebra, but the non-
negativity constraint makes much of this previous work inapplicable to the present case
[8].
Here we discuss two algorithms for NMF based on iterative updates of Wand H. Because
these algorithms are easy to implement and their convergence properties are guaranteed,
we have found them very useful in practical applications. Other algorithms may possibly
be more efficient in overall computation time, but are more difficult to implement and may
not generalize to different cost functions. Algorithms similar to ours where only one of the
factors is adapted have previously been used for the deconvolution of emission tomography
and astronomical images [9, 10, 11, 12].
At each iteration of our algorithms, the new value of W or H is found by multiplying the
current value by some factor that depends on the quality ofthe approximation in Eq. (1). We
prove that the quality of the approximation improves monotonically with the application
of these multiplicative update rules. In practice, this means that repeated iteration of the
update rules is guaranteed to converge to a locally optimal matrix factorization.

3 Cost functions
To find an approximate factorization V ~ W H, we first need to define cost functions
that quantify the quality of the approximation. Such a cost function can be constructed
using some measure of distance between two non-negative matrices A and B . One useful
measure is simply the square of the Euclidean distance between A and B [13],
IIA - BI12 = L(Aij -
Bij)2
(2)
ij
This is lower bounded by zero, and clearly vanishes if and only if A = B .
Another useful measure is
D(AIIB)
=
2:
k·
( Aij log B:~
- Aij
+ Bij )
(3)
"J
Like the Euclidean distance this is also lower bounded by zero, and vanishes if and only
if A = B . But it cannot be called a "distance", because it is not symmetric in A and B,
so we will refer to it as the "divergence" of A from B. It reduces to the Kullback-Leibler
divergence, or relative entropy, when 2:ij Aij = 2:ij Bij = 1, so that A and B can be
regarded as normalized probability distributions.

We now consider two alternative formulations of NMF as optimization problems:
Problem 1 Minimize
IIV -
W
HI12 with
respect to Wand H, subject to the constraints
W,H~O.
Problem 2 Minimize D(VIIW H) with re.lpect to Wand H, subject to the constraints
W,H~O.
Although the functions IIV - W HI12 and D(VIIW H) are convex in W only or H only, they
are not convex in both variables together. Therefore it is unrealistic to expect an algorithm
to solve Problems 1 and 2 in the sense of finding global minima. However, there are many
techniques from numerical optimization that can be applied to find local minima.
Gradient descent is perhaps the simplest technique to implement, but convergence can be
slow. Other methods such as conjugate gradient have faster convergence, at least in the
vicinity of local minima, but are more complicated to implement than gradient descent
[8] . The convergence of gradient based methods also have the disadvantage of being very
sensitive to the choice of step size, which can be very inconvenient for large applications.
4 Multiplicative update rules
We have found that the following "multiplicative update rules" are a good compromise
between speed and ease of implementation for solving Problems 1 and 2.
Theorem 1 The Euclidean distance II V - W H II is non increasing under the update rules
(WTV)att
Hal' +- Hal' (WTWH)att
(V HT)ia
Wia +- Wia(WHHT)ia
(4)
The Euclidean distance is invariant under these updates if and only if Wand H are at a
stationary point of the distance.
Theorem 2 The divergence D(VIIW H) is nonincreasing under the update rules
H
att +-
H
att
2:i WiaVitt/(WH)itt
" W
L..Jk
ka
Wia +- Wia
2:1' HattVitt/(WH)itt
" H
L..Jv
av
(5)
The divergence is invariant under these updates if and only ifW and H are at a stationary
point of the divergence.
Proofs of these theorems are given in a later section. For now, we note that each update
consists of multiplication by a factor. In particular, it is straightforward to see that this
multiplicative factor is unity when V = W H, so that perfect reconstruction is necessarily
a fixed point of the update rules.
5 Multiplicative versus additive update rules
It is useful to contrast these multiplicative updates with those arising from gradient descent
[14]. In particular, a simple additive update for H that reduces the squared distance can be
written as
(6)
If 'flatt are all set equal to some small positive number, this is equivalent to conventional
gradient descent. As long as this number is sufficiently small, the update should reduce
IIV - WHII·

Now if we diagonally rescale the variables and set
Halt
"Ialt
(7)
= (WTW H)alt '
then we obtain the update rule for H that is given in Theorem 1. Note that this rescaling
results in a multiplicative factor with the positive component of the gradient in the denom-
inator and the absolute value of the negative component in the numerator of the factor.
For the divergence, diagonally rescaled gradient descent takes the form
Halt
f-
Halt
+ "Ialt
[~Wia (:;;)ilt - ~ Wia].
(8)
Again, if the "Ialt are small and positive, this update should reduce D (V II W H). If we now
set
Halt
"Ialt= ui
~ W. '
za
(9)
then we obtain the update rule for H that is given in Theorem 2. This rescaling can also
be interpretated as a multiplicative rule with the positive component of the gradient in the
denominator and negative component as the numerator of the multiplicative factor.
Since our choices for "Ialt are not small, it may seem that there is no guarantee that such a
rescaled gradient descent should cause the cost function to decrease. Surprisingly, this is
indeed the case as shown in the next section.
6 Proofs of convergence
To prove Theorems 1 and 2, we will make use of an auxiliary function similar to that used
in the Expectation-Maximization algorithm [15, 16].
Definition 1 G(h, h') is an auxiliary functionfor F(h)
G(h, h') ~ F(h),
G(h, h)
if the conditions
= F(h)
(10)
are satisfied.
The auxiliary function is a useful concept because of the following lemma, which is also
graphically illustrated in Fig. 1.
Lemma 1 IfG is an auxiliary junction, then F is nonincreasing under the update
ht+1 = argmlnG (h,h t )
Proof: F(ht+1) ~ G(ht+1, ht) ~ G(ht, ht)
(11)
= F(ht) •
Note that F(ht+1) = F(ht) only if ht is a local minimum of G(h, ht). If the derivatives
of F exist and are continuous in a small neighborhood of h t , this also implies that the
derivatives 'V F(ht) = O. Thus, by iterating the update in Eq. (11) we obtain a sequence
of estimates that converge to a local minimum h min = argminh F(h) of the objective
function:
We will show that by defining the appropriate auxiliary functions G(h, ht) for both IIV -
W HII and D(V, W H), the update rules in Theorems 1 and 2 easily follow from Eq. (11).

....

Explain the different types of algorithms to solve the problem and the one used.

%-----------------------------------
%	SUBSECTION 2
%-----------------------------------

\subsection{Choosing number of signatures}

<https://www.academia.edu/238621/Non-negative_Matrix_Factorization_Assessing_Methods_for_Evaluating_the_Number_of_Components_and_the_Effect_of_Normalization_Thereon>

Explain different ways to choose N, the way used, the N used and why.