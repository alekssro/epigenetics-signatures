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

  Arabidopsis data includes Single Methylation Polymorfisms (SMPs) for Arabidopsis located in different locations, with following information: chr, position, strand, class, gclass, *locations.

- Get alternative data