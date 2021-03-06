% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{small_5050_mix}
\alias{small_5050_mix}
\title{A random subset of genes and cells from a 50:50 mixture of 293T:Jurkat cells}
\format{A Seurat object with 840 genes and 384 cells:
\describe{
  \item{1031 genes}{}
  \item{255 cells}{}
  ...
}}
\usage{
small_5050_mix
}
\description{
The dataset is a fully processed Seurat object of the mixture sub-sampled
data set by the standard Seurat procedure, therefore should contain all
normal fields expected from a Seurat object, sucha as identities, pca and
tsne dimensions.
}
\details{
Note that the subsampling was done such that variable genes were preserved
form the original dataset and 500 more were chosen randomly
}
\keyword{datasets}
