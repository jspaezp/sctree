


#' plot_flowstyle
#'
#' Simulates flow cytometry plots for single cell rna seq data
#'
#' Note that it imputes a little bit of noise to values of 0 just for the
#' Sake of visualization
#'
#' @param object data, either a Seurat of a data frame object
#' @param markernames names of the markers to be used for the plot
#'
#' @return a ggplot grid with the plots
#' @export
#'
#' @examples
#' plot_flowstyle(Seurat::pbmc_small, c("ACRBP", "TSC22D1", "VDAC3"))
plot_flowstyle <- function(object, markernames) {
    # TODO add argument to change to natural scale ...
    UseMethod("plot_flowstyle", object)
}


#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 aes_string theme_bw
#' @describeIn plot_flowstyle Draw a flowstyle plot from a data.frame
#' @export
plot_flowstyle.data.frame <- function(df, markernames, classif_col = "ident") {

    tmp_ident <- df[[classif_col]]
    df <- df[,markernames]

    # TODO add argument decide when to add noise to the dataset ...
    df[df == 0] <- abs(rnorm(sum(df == 0),mean = 0, sd = 0.2))
    df$ident <- tmp_ident


    g <- GGally::ggpairs(
        as.data.frame(df),
        columns = 1:(ncol(df) - 1),
        ggplot2::aes_string(colour = "ident"),
        progress = FALSE,
        lower = list(
            continuous = GGally::wrap(
                "dot_no_facet",
                alpha = 0.2)),
        diag = list(
            continuous = GGally::wrap(
                'densityDiag',
                alpha = 0.3)),
        upper = list(
            continuous =  GGally::wrap(
                "density",
                alpha = 0.4))) +
        ggplot2::theme_bw()
    return(g)
}


#' @describeIn plot_flowstyle Draw a flowstyle plot from a seurat object
#' @export
plot_flowstyle.seurat <- function(Seurat, markernames, classif_col = "ident") {
    tmp <- as.data.frame.Seurat(Seurat, markernames)
    plot_flowstyle.data.frame(
        tmp, markernames = markernames,
        classif_col = classif_col)
}
