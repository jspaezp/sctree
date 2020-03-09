

#' include_roxygen_example
#'
#' Include the output of code in your documentation automatically
#'
#' This function takes a string that can be interpreted as r code and returns it
#' formatted as a roxygen exmaple, including tags and output commented out.
#'
#' It is meant to be used with the evalRd tag as
#' `#' @@evalRd include_roxygen_example("print('a')")`
#'
#' @param example_string the string to be interpreted as R code
#'
#' @return a string that can be interpreted literally as a roxygen comment
#' @importFrom utils capture.output
#' @keywords internal
include_roxygen_example <- function(example_string) {

    run_output <- capture.output(eval(parse(text = example_string)))
    run_output <- paste0("# ", run_output)

    return_string <- paste0(
        c("\\examples{", example_string, run_output, "}"),
        collapse = "\n")

    return(return_string)
}


#' Gets a data frame from a Seurat object
#'
#' Provided a Seurat object, returs a data frame of the count values, being
#' the columns each 'gene' and the rows each UMI/cell.
#'
#' It returns only the genes annotated as variable and the identity column.
#'
#' @param x A Seurat object
#' @param genes genes to extract to the data.frame
#' @param fix_names logical value indicating wether the gene names should be
#'     converted to R-compatible names. defaults to FALSE
#' @param ... additional arguments passed to `Seurat::FetchData`
#'
#' @return a data frame.
#' @export
#'
#' @evalRd include_roxygen_example({
#'     "as.data.frame(Seurat::pbmc_small,
#'     Seurat::VariableFeatures(Seurat::pbmc_small))[1:3,1:3]"
#'     })
#' @importFrom Seurat FetchData Idents VariableFeatures
as.data.frame.Seurat <- function(x, genes = Seurat::VariableFeatures(x), fix_names = TRUE, ...) {
    # TODO possibly also a warning if it is not a variable gene and an error
    # if it does not exist also an argument to force though the error ...
    tmp <- Seurat::FetchData(x, vars = genes, ...)
    tmp <-  as.data.frame(tmp, stringsAsFactors = FALSE)

    tmp$ident <- Seurat::Idents(x)[rownames(tmp)]

    if (fix_names) {
        colnames(tmp) <- make.names(colnames(tmp))
    }

    return(tmp)
}


#' Queries gene symbols and returns if they are annotated as membrane-localized
#'
#' returns a logical vector indicating wether there is a GO annotation
#' corresponding to *integral component of membrane*. the species to be queried
#' can be changed by using the corresponging database.
#'
#' Can also filter the annotations by the level of annotation desired.
#' For additional information on the accepted annotation codes please visit
#' http://geneontology.org/docs/guide-go-evidence-codes/
#'
#'
#' @param gene_symbols a character vector of gene symbols. zb: CD1A
#' @param db a database that inherits the select property.
#'     defaults to org.Hs.eg.db
#' @param evidence_codes the evidence codes that will be considered for the
#'     annotation. if set to null will return all annotations. defaults to
#'     c("EXP", "IDA", "IPI", "IMP", "IGI")
#'
#' @return logical vector indicating which gene_symbols match the criteria.
#' @export
#'
#' @examples
#' is_gene_membrane(c("CD4", "UPP1"))
#' # [1]  TRUE FALSE
#' @importFrom AnnotationDbi select
is_gene_membrane <- function(gene_symbols,
                             db = org.Hs.eg.db::org.Hs.eg.db,
                             evidence_codes = c("EXP", "IDA", "IPI",
                                                "IMP", "IGI")) {

    suppressMessages({
        go_annotated <- AnnotationDbi::select(
            db,
            keys = gene_symbols,
            columns = c("GOALL", "EVIDENCEALL"),
            keytype = "ALIAS")
    })


    if (!is.null(evidence_codes)) {
        go_annotated <- go_annotated[go_annotated[["EVIDENCEALL"]] %in%
                                         evidence_codes,]
    }

    # "GO:0005886" means plasma membrane
    # "GO:0044459" means plasma membrane part
    # "GO:0016021" integral component of membrane
    # "GO:0005887" integral component of plasma membrane
    surface_annotated <- go_annotated[
        go_annotated[["GOALL"]] %in%
            "GO:0005886",]

    return(gene_symbols %in% surface_annotated[["ALIAS"]])
}


# TODO: refactor this function to use standard evaluation to prevent the
# following error:
# Undefined global functions or variables:
# ALIAS SYMBOL plot.flowstyle

#' Gets the aliases for a gene names
#'
#' gets the aliases for a series of gene names
#'
#' @param gene_symbols a character vector of gene symbols. zb: CD1A
#' @param db a database that inherits the select property.
#'     defaults to org.Hs.eg.db
#'
#' @return a named list whose names are the provided symbols and the elements are
#'     the character vectors with the aliases
#' @export
#' @seealso get_genesymbols
#'
#' @examples
#' get_aliases(c("MAPK1", "CD4"))
#' # $MAPK1
#' # [1] "ERK"      "ERK-2"    "ERK2"     "ERT1"     "MAPK2"    "P42MAPK"
#' # [7] "PRKM1"    "PRKM2"    "p38"      "p40"      "p41"      "p41mapk"
#' # [13] "p42-MAPK" "MAPK1"
#' #
#' # $CD4
#' # [1] "CD4mut" "CD4"
#' @importFrom tidyr nest
#' @importFrom wrapr named_map_builder
#' @importFrom AnnotationDbi select
get_aliases <- function(gene_symbols,
                        db = org.Hs.eg.db::org.Hs.eg.db) {

    suppressMessages({
        alias_df <- AnnotationDbi::select(
            db,
            keys = gene_symbols,
            columns = c("ALIAS", "SYMBOL"),
            keytype = "SYMBOL")

        alias_df <- tidyr::nest(alias_df, data = "ALIAS")
        alias_vector <- wrapr::named_map_builder(
            names = alias_df[["SYMBOL"]],
            values = lapply(alias_df[["data"]], function(x) x[["ALIAS"]]))
    })

    return(alias_vector)
}


#' Gets the official gene symbol for protein or gene aliases
#'
#' gets the aliases for a series of common gene names
#'
#' @param gene_aliases a character vector of gene aliases zb: ERK
#' @param db a database that inherits the select property.
#'     defaults to org.Hs.eg.db
#'
#' @return a named list whose names are the provided aliases and the elements
#'     are the character vectors with the gene names
#' @export
#' @seealso get_aliases
#'
#' @examples
#' get_genesymbols("ERK")
#' # $ERK
#' # [1] "EPHB2" "MAPK1"
#' get_genesymbols(c("SUPERFAKEGENE", "ERK"))
#' # $SUPERFAKEGENE
#' # [1] NA
#' #
#' # $ERK
#' # [1] "EPHB2" "MAPK1"
#' @importFrom tidyr nest
#' @importFrom wrapr named_map_builder
#' @importFrom AnnotationDbi select
get_genesymbols <- function(gene_aliases,
                            db = org.Hs.eg.db::org.Hs.eg.db) {

    suppressMessages({
        alias_df <- AnnotationDbi::select(
            db,
            keys = gene_aliases,
            columns = c("ALIAS", "SYMBOL"),
            keytype = "ALIAS")

        alias_df <- tidyr::nest(alias_df, "SYMBOL")
        alias_vector <- wrapr::named_map_builder(
            names = alias_df[["ALIAS"]],
            values = lapply(alias_df[["data"]], function(x) x[["SYMBOL"]]))
    })

    return(alias_vector)
}


#' Check if the element has the structure of a confusion matrix
#'
#' Slim check for the element to assure it has 2 dimensions,
#' is a table and whose elements are exclusively integers
#'
#' @param x table object
#'
#' @return logical indicating wether the element matches the description
#' @export
#' @keywords internal
#'
#' @examples
#' is.confusion.matrix(table(1:3, 3:1, 1:3)) # FALSE
#' is.confusion.matrix(table(1:3, 3:1)) # TRUE
is.confusion.matrix <- function(x) {
    return(length(dim(x)) == 2 & is.table(x) & is.integer(x))
}

#' Convert confusion matrices and tables to frequency matrices
#'
#' @param confusion_matrix a matrix or table object that can be assumed to be a
#'     confusion matrix.
#'
#' @return a frequency matrix normalizing by the column total
#'     (expresed as percentages)
#' @export
#'
#' @examples
#' table(1:3, 3:1)
#' #   1 2 3
#' # 1 0 0 1
#' # 2 0 1 0
#' # 3 1 0 0
#' as.frequency.matrix(table(1:3, 3:1))
#' #
#' #     1   2   3
#' # 1   0   0 100
#' # 2   0 100   0
#' # 3 100   0   0
#' as.frequency.matrix(table(c(3,3,2), 3:1))
#' #
#' #     2   3
#' # 1 100   0
#' # 2   0 100
#' # 3   0 100
as.frequency.matrix <- function(confusion_matrix) {
    if (!is.confusion.matrix(confusion_matrix)) {
        warning(
            paste0(
                "Provided element is not configured ",
                "as a confusion matrix, attempting to follow though ",
                "but might result in unexpected results"
            )
        )
    }

    frequency_matrix <- t(
        100 * (t(as.matrix(confusion_matrix)) /
                   colSums(confusion_matrix)))

    class(frequency_matrix) <-
        c(class(frequency_matrix), "frequency.matrix")
    return(frequency_matrix)
}


#' Basic heatmap using ggplot
#'
#' @param df data.frame
#' @param x column name used in the x axis
#' @param y column name used in the y axis
#' @param value value of the color to be used
#' @param show_number logical indicating wether the value of z should be printed
#'     in the plot.
#'
#' @return a ggplot object
#' @keywords internal
#' @export
#'
#' @examples
#' my_df <- data.frame(x = 1:2, y = 1:2, z = 1:2)
#' ggheatmap_base(my_df, x = "x", y = "y", value = "z", show_number = TRUE)
#' ggheatmap_base(my_df, x = "x", y = "y", value = "z", show_number = FALSE)
#' @importFrom ggplot2 ggplot aes_string geom_tile theme_bw geom_text
#' @importFrom viridis viridis
ggheatmap_base <- function(df,
                           x = "cluster",
                           y = "predicted",
                           value = "Freq",
                           show_number = FALSE) {
    g <- ggplot2::ggplot(
        df,
        ggplot2::aes_string(
            x = x,
            y = y,
            fill = value,
            label = value)
    ) +
        ggplot2::geom_tile(colour = "#ffffff") +
        ggplot2::theme_bw() +
        viridis::scale_fill_viridis(direction = -1)


    if (show_number) {
        g <- g + ggplot2::geom_text()
    }

    return(g)
}


#' autoplot
#' wrapper arround ggplot2 autoplot generic
#' @param object an object ...
#' @param ... additional arguments passed to the correct method
#'
#' @export
#' @keywords internal
#' @importFrom ggplot2 autoplot
autoplot <- function(object, ...) {
    ggplot2::autoplot(object, ...)
}

#' autoplot.table
#'
#' Heatmap with ggplot from table objects
#'
#' @param object a table object
#' @param min_color minimum value to show colored, value sunder this will
#'     be shown as gray
#' @param show_number logical indicating wether the number should be printed
#'     in the box of the heatmap
#' @param ... additional arguments passed to downstream methods
#'
#'
#' @return ggpplot object
#' @method autoplot table
#' @export
#'
#' @examples
#' autoplot(table(1:10, 1:10), show_number = TRUE)
#' @importFrom ggplot2 geom_tile geom_point aes_string labs autoplot
autoplot.table <- function(object,
                           min_color = NULL,
                           show_number = FALSE,
                           ...) {
    df <- as.data.frame(object)

    df[["Freq"]] <- round(df[["Freq"]], 2)

    stopifnot(length(dim(object)) == 2)

    x_value <- colnames(df)[[1]]
    y_value <- colnames(df)[[2]]

    if (!is.null(min_color)) {
        df_under_min <- df[df[["Freq"]] < min_color, ]
        df[df[["Freq"]] < min_color, "Freq"] <- NA
    }

    g <- ggheatmap_base(
        df = df,
        x = x_value,
        y = y_value,
        value = "Freq",
        show_number = show_number
    )

    if (!is.null(min_color)) {
        g <- g +
            ggplot2::geom_tile(
                fill = "#666666",
                colour = "#ffffff",
                data = df_under_min,
                show.legend = TRUE
            ) +
            # adding the geom_point is the only hack I could find to make
            # the legend show up for the low frequency scores.
            ggplot2::geom_point(
                df_under_min,
                alpha = 0,
                mapping = ggplot2::aes_string(colour = paste0("Freq < ", min_color)),
                show.legend = TRUE,
                inherit.aes = TRUE
            ) +
            ggplot2::labs(colour = paste0("Freq < ", min_color))

    }

    return(g)
}


#' @describeIn autoplot.table plot a matrix object as a heatmap using ggplot2
#' @method autoplot matrix
#' @export
autoplot.matrix <- function(object, ...) {
    autoplot.table(as.table(object), ...)
}
