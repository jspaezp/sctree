

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
as.data.frame.Seurat <- function(seurat, genes, fix_names = FALSE) {
    # TODO possibly also a warning if it is not a variable gene and an error if it does not exist
    # also an argument to force though the error ...
    tmp <- Seurat::FetchData(seurat, vars.all = genes)
    tmp <-  as.data.frame(tmp, stringsAsFactors = FALSE)

    tmp$ident <- Seurat::GetIdent(seurat, uniq = FALSE, cells.use = rownames(tmp))

    if (fix_names) {
        colnames(tmp) <- make.names(colnames(tmp))
    }

    return(tmp)
}




#' is_gene_membrane
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
#' > is_gene_membrane(c("CD4", "GAPDH"))
#' [1]  TRUE FALSE
is_gene_membrane <- function(gene_symbols,
                             db = org.Hs.eg.db::org.Hs.eg.db,
                             evidence_codes = c("EXP", "IDA", "IPI",
                                                "IMP", "IGI")) {

    suppressMessages({
        go_annotated <- AnnotationDbi::select(
            # "GO:0005886" means plasma membrane
            # "GO:0044459" means plasma membrane part
            # "GO:0016021" integral component of membrane
            # "GO:0005887" integral component of plasma membrane
            db,
            keys = gene_symbols,
            columns = c("GOALL", "EVIDENCEALL"),
            keytype = "ALIAS")
    })


    if (!is.null(evidence_codes)) {
        go_annotated <- go_annotated[go_annotated[["EVIDENCEALL"]] %in%
                                         evidence_codes,]
    }

    surface_annotated <- go_annotated[
        go_annotated[["GOALL"]] %in%
            "GO:0005887",]

    return(gene_symbols %in% surface_annotated[["ALIAS"]])
}




#' get_aliases
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
#'
#' @examples
#' > get_aliases(c("MAPK1", "CD4"))
#' $MAPK1
#' [1] "ERK"      "ERK-2"    "ERK2"     "ERT1"     "MAPK2"    "P42MAPK"
#' [7] "PRKM1"    "PRKM2"    "p38"      "p40"      "p41"      "p41mapk"
#' [13] "p42-MAPK" "MAPK1"
#'
#' $CD4
#' [1] "CD4mut" "CD4"
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

        alias_df <- tidyr::nest(alias_df, ALIAS)
        alias_vector <- wrapr::named_map_builder(
            names = alias_df[["SYMBOL"]],
            values = lapply(alias_df[["data"]], function(x) x[["ALIAS"]]))
    })

    return(alias_vector)
}


#' get_genesymbols
#'
#' gets the aliases for a series of common gene names
#'
#' @param gene_symbols a character vector of gene aliases zb: ERK
#' @param db a database that inherits the select property.
#'     defaults to org.Hs.eg.db
#'
#' @return a named list whose names are the provided aliases and the elements
#'     are the character vectors with the gene names
#' @export
#'
#' @examples
#' > get_genesymbols("ERK")
#' $ERK
#' [1] "EPHB2" "MAPK1"
#' > get_genesymbols(c("SUPERFAKEGENE", "ERK"))
#' $SUPERFAKEGENE
#' [1] NA
#'
#' $ERK
#' [1] "EPHB2" "MAPK1"
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

        alias_df <- tidyr::nest(alias_df, SYMBOL)
        alias_vector <- wrapr::named_map_builder(
            names = alias_df[["ALIAS"]],
            values = lapply(alias_df[["data"]], function(x) x[["SYMBOL"]]))
    })

    return(alias_vector)
}

