


#' A random subset of genes and cells from a 50:50 mixture of 293T:Jurkat cells
#'
#' The dataset is a fully processed seurat object of the mixture sub-sampled
#' data set by the standard seurat procedure, therefore should contain all
#' normal fields expected from a seurat object, sucha as identities, pca and
#' tsne dimensions.
#'
#' Note that the subsampling was done such that variable genes were preserved
#' form the original dataset and 500 more were chosen randomly
#'
#' @format A seurat object with 840 genes and 384 cells:
#' \describe{
#'   \item{1031 genes}{}
#'   \item{255 cells}{}
#'   ...
#' }
"small_5050_mix"


#' A random subset of genes and cells from a 99:1 mixture of 293T:Jurkat cells
#'
#' The dataset is a fully processed seurat object of the mixture sub-sampled
#' data set by the standard seurat procedure, therefore should contain all
#' normal fields expected from a seurat object, sucha as identities, pca and
#' tsne dimensions.
#'
#' Note that the subsampling was done such that variable genes were preserved
#' form the original dataset and 500 more were chosen randomly
#'
#' @format A seurat object with 840 genes and 384 cells:
#' \describe{
#'   \item{840 genes}{}
#'   \item{384 cells}{}
#'   ...
#' }
"small_9901_mix"
