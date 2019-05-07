# MSc. Thesis

## Title: Identifying signatures of human epigenetic modifications among tissues.

## Project description

The number of different epigenetic landscapes for a genome may be inestimable, but we can find correlations between specific epigenetic modifications which are typically associated in concrete functions and development states. In such a way, we reduce the dimensionality of the problem making it easier to draw conclusions from the analysis of the epigenetic modifications, as well as being able to use the smaller set of correlated modifications (or "signatures") as input for predictive modelling or supervised machine learning analysis.

Non-negative Matrix Factorization (NMF) reveals as an ideal method for the task of finding combinatorial patterns of epigenetic modifications. We can then study the state of each epigenetic modification type in the defined loci of the tissue. From this information we would obtain the different epigenetic signatures which we will use for association and simulation analysis.

## Preliminary study

### Considering different possibilities:

- Get raw/semi-process data from http://www.roadmapepigenomics.org/ and https://www.encodeproject.org/

  Bed files found in this webpages can be imported into R easily thanks to rtracklayer package. Then GenomicRanges & IRanges packages allow us to count how many epigenetic modifications are present in each gene.

- Get Arabidopsis data from https://datadryad.org/resource/doi:10.5061/dryad.80442

  Arabidopsis data includes Single Methylation Polymorfisms (SMPs) for Arabidopsis located in different locations, with following information: chr, position, strand, class, gclass, locations.

- Get alternative data

#### Plan A:

The one specified in the project description but it needs: define properly the elements in the V matrix and the final objective of the analysis.

Similar previous studies:

[A computational approach for the functional classification of the epigenome](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-017-0131-7) :

In this paper they try to highlight functional relationships between chromatin states and transcriptional activity. For that they applied NMF on multiple genomic datasets collecting 13 different epigenetic marks: H3K27me3, H3K9me3, CTCF (CCCTC binding factor), DNase-seq, H3K27ac, Pol2, H3K36me3, H3K79me1, H3K4me1, H2A.Z, H3K4me3, H3K9ac and H3K4me2 in human embryonic stem cells H1 (hESCs). They proceed in the following way:

1. Process BAM raw files to get BED files (BED files are already available).
2. They partitioned the genome in 200-bp non-overlapping bins (trying to approximate with each bin a single nucleosome along the DNA)
3. For each sample, each epigenetic mark read was assigned to one of the bins + correcting for the expected number of reads falling in a given bin.
4. To correct for variability in the signal ranges of the different epigenetic marks, they scaled the values from 0 to 1 using a sigmoid function.
5. At this point they obtained a V matrix with 833738 significant bins (rows) and 13 epigenetic marks (columns) in which they applied NMF (finding an optimal factorization rank of 7 when comparing to randomly generated data results).
6. From the results of the NMF, they could relate the 7 different profiles, in relation to the epigenetic marks (H matrix), to 7 different "genomic states": (1) Active promoter, (2) Repressed Chromatin, (3) Transcription Initiation, (4) Repressed Regulatory Regions, (5) Gene Body Transcription, (6) Enhancer Regions and (7) Regulatory Elements. For instance, in profile 1 marks Pol2 and H3K27ac are predominant, which are known to be related with active gene expression.
7. Moreover, they used these profiles in order to predict the function of already known regions and compared the performance to predictions based on individual epigenetic marks, showing that NMF are able to identify combinatorial interactions that are more informative than single mark contributions.
8. Furthermore, they studied the co-occurrence of chromatin profiles in bin assignment frequency, finding high co-occurrence between some profiles (such as Active promoter (1) and Transcription Initiation (3) or Repressed Chromatin (2) and Gene Body Transcription (5)).
9. Other analysis made in the paper include: spatial relationship between profiles, association between chromatin profiles and gene expression (finding a direct association between chromatin profiles and transcriptional status of the gene)...

[Combinatorial epigenetic patterns as quantitative predictors of chromatin biology](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-76 ) :

In this paper they applied NMF to 7 histone modifications and one histone variant in the A549 adenocarcinomic alveolar basal epithelial cell line (histone modifications = columns). The rows of the V matrix in this case are multiple loci; regions of TSS-proximal gene bodies since they contain epigenetic traces of transcription initiation and elongation. From this analysis they performed some case studies:

- Regression Pol2 binding: they try to predict levels of Pol2 binding at promoters of protein-coding genes in human embryonic stem cells (H1ESC). They used ridge regression comparing the performance of individual marks and combinatory patterns derived from NMF.
- Classification of Pol-2 bound enhancers vs promoters
- Gene set enrichment analysis.

[Intratumor heterogeneity in epigenetic patterns](https://www.sciencedirect.com/science/article/pii/S1044579X17302262)



Taking this into account I think it would be interesting to perform an NMF analysis in multiple cancer cells (A549, K562, MCF-7, HCT116, HeLa-S3, HepG2...), for which the 13 epigenetic marks (columns) are available, in the 833738 significant bins (rows). As they do in the first paper, relating the different profiles to biological processes and comparing the differences between the cancer types in terms of predominant metabolic pathways, by combining also gene expression data.

#### Plan B:

- Perform a SMPs analysis combined with the environmental and genetic features such as in https://datadryad.org/resource/doi:10.5061/dryad.80442





# Pipeline for processing files:

### create_conda_env.sh

Creates a conda environment with an input name and required packages for the analysis



### prepareBedAlignment.sh

Process raw BAM files into filtered BED files

1. Remove duplicates using `PicardCommandLine MarkDuplicates`
2. Filter by quality using `samtools view -q 20`
3. For TFs and HistMod signals apply tag extension of 200bp
   1. split by strand +/-
   2. apply tag extension of 200bp
   3. merge extended bed files

At the end we get "processed.bed" files as output.



### bedToNormalizedSignal.sh
Convert processed BED files into Normalized Signal tracks.

1. Calculate the total number of informative reads mapped across all the datasets: Take note of the number of informative reads in that sample to estimate total number of tags over all experiments in a cell line
2. Once the cell-type and the signal-type are defined, extract all samples from the cell-type X and signal-type Y from the ENCODE_project_datasheet_samples.tsv then loop row by row on the subset
3. Identify and retrieve processed.bed file and its path
4. Get the raw counts for each signal process file in each genomic bin (bins with no tags are discarted)
5. Add the rawCount full path name into a file (input_data_files.txt) for next steps.



- estimate both observed and expected read count distribution for each sample in a given epigenetic mark.
- Estimation of the normalized coverage signal for each chromatin mark.


#### TODO:
https://app.clickup.com/1701348/d/l?p=2259421&c=2661775&s=-1

## Selection of appropiate N:

Run 1: 

```R
[1] "N: 1 sum(V - W * H)^2 = 4576574736.30756"
[1] "N: 2 sum(V - W * H)^2 = 4756360089.77251"
[1] "N: 3 sum(V - W * H)^2 = 3894348833.50569"
[1] "N: 4 sum(V - W * H)^2 = 3409936797.94337"
[1] "N: 5 sum(V - W * H)^2 = 4201834121.94814"
[1] "N: 6 sum(V - W * H)^2 = 3630107685.789"
[1] "N: 7 sum(V - W * H)^2 = 2955505248.37936"
[1] "N: 8 sum(V - W * H)^2 = 2829057005.47625"
```

Run 2:

```R
[1] "N: 1 sum(V - W * H)^2 = 5629801823.51504"
[1] "N: 2 sum(V - W * H)^2 = 5064812567.71203"
[1] "N: 3 sum(V - W * H)^2 = 3913316199.03133"
[1] "N: 4 sum(V - W * H)^2 = 3753135654.59337"
[1] "N: 5 sum(V - W * H)^2 = 3455607394.97327"
[1] "N: 6 sum(V - W * H)^2 = 2496443481.38631"
[1] "N: 7 sum(V - W * H)^2 = 3327497061.56176"
[1] "N: 8 sum(V - W * H)^2 = 2885214901.65869"
```

Run 3:

```
[1] "N: 1 sum(V - W * H)^2 = 4672740729.33611"
[1] "N: 2 sum(V - W * H)^2 = 4906950231.04911"
[1] "N: 3 sum(V - W * H)^2 = 4263154041.06151"
[1] "N: 4 sum(V - W * H)^2 = 4027672185.15617"
[1] "N: 5 sum(V - W * H)^2 = 3410076391.85276"
[1] "N: 6 sum(V - W * H)^2 = 3879726657.22432"
[1] "N: 7 sum(V - W * H)^2 = 3740776225.74354"
[1] "N: 8 sum(V - W * H)^2 = 3266565144.27729"
```



**Filter the data**

https://app.clickup.com/t/vrd4r

![1557239241626](/home/alekssro/.config/Typora/typora-user-images/1557239241626.png)

<https://probabilityandstats.wordpress.com/tag/poisson-gamma-mixture/>

**Initiation of W and H**

A general rule of thumb for the stochastic initialization approach is to perform several runs of the NMF (i.e., several random initializations for matrices W and H) and keep the factorization that minimizes the reconstruction error across
multiple runs.
$$
min[δ = V − WH ]
$$

where δ is the difference between the real and the model output values of the epigenetic mark levels.



**Choosing of n**

A common strategy is to perform NMF in a pre-defined ranged of r values, estimate a quality measure of the results, and select the best value of r according to this quality criteria. 

Different strategies have been proposed to select the best factorization rank. The most common approach is based on the *cophenetic correlation coefficient*, which reflects the overall cluster stability obtained after the factorization process.

A more robust approach suggests to take the smallest value of r at which the decrease in the residual sum of squares (RSS) between V(j, k) and the NMF model is larger than the decrease observed in the random data.

