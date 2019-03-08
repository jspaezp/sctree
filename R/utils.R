

#' TSNE plot generation
#'
#' slim wrapper arround Seurat::TSNEPlot that overwrites the return option and
#' by default uses a tenth of the data.
#'
#' @param object A seurat object
#' @param fraction fraction of the data to be used for plotting, defaults to 0.1
#'
#' @return a ggplot object with the tsne plot
#' @export
#'
#' @examples
#' > tsne_plot(Seurat::pbmc_small)
#' @importFrom Seurat TSNEPlot
tsne_plot <- function(object, fraction = 0.1, ...) {
    g <- Seurat::TSNEPlot(
        object = object,
        dim.1 = 1,
        do.return = TRUE,
        cells.use = sample(object@cell.names,
                           size = ceiling(length(object@cell.names)*fraction)),
        ...)
    g
}



#' Gets a data frame from a seurat object
#'
#' Provided a seurat object, returs a data frame of the count values, being
#' the columns each 'gene' and the rows each UMI/cell.
#'
#' It returns only the genes annotated as variable and the identity column.
#'
#' @param seurat A seurat object
#'
#' @return a data frame.
#' @export
#'
#' @examples
#' > as.data.frame.Seurat(pbmc_small, seurat@var.genes)[1:3,1:3]
#'                    LTB EAF2 CD19
#' ATGCCAGAACGACT 6.062788    0    0
#' CATGGCCTGTGCAT 6.714813    0    0
#' GAACCTGATGAACC 7.143118    0    0
#' @importFrom Seurat FetchData GetIdent
as.data.frame.Seurat <- function(seurat, genes) {
    # TODO possibly also a warning if it is not a variable gene and an error if it does not exist
    # also an argument to force though the error ...
    tmp <- Seurat::FetchData(seurat, vars.all = genes)
    tmp <-  as.data.frame(tmp, stringsAsFactors = FALSE)

    tmp$ident <- Seurat::GetIdent(seurat, uniq = FALSE, cells.use = rownames(tmp))
    return(tmp)
}

