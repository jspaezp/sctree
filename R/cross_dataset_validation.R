
# TODO refine what output is really needed
# TODO add option on how to handle missing values, either adding zeros or
#      filtering common names

#' Fits a decision tree in one data set and tests the performance in another
#'
#' @param train a Seurat object to be used for trainning.
#' @param test another Seurat object to be used for testing.
#' @param cluster the cluster whose equivalence needs to be found.
#' @param genes_use character vector specifying which genes to use for the
#'     classification, defaults to Seurat::VariableFeatures(train)
#' @param warn.gene.removal logical indicating wether to warn the user when
#'     genes are removed because they are missing in one of the datasets.
#'     defults to TRUE
#' @param ... additional arguments to be passed to ranger_importances.Seurat
#'
#' @return a list containing the
#'     (1) tree fit,
#'     (2) a summary_table
#'     (3) the concensus rules of the tree
#'     (4) ranger_significance_table
#'     (5) the suggested genes for the gating
#' @export
#'
#' @examples
#' cross_validate(small_5050_mix, small_9901_mix, cluster = "0")
#' cross_validate(small_5050_mix, small_9901_mix, cluster = "ALL")
#' @importFrom stats as.formula predict
cross_validate <- function(train, test,
                           cluster, genes_use = Seurat::VariableFeatures(train),
                           warn.gene.removal = TRUE,
                           ...) {

    rang_importances <- ranger_importances.Seurat(train,
                                                  cluster = cluster,
                                                  genes_use = genes_use,
                                                  ...)


    # TODO: implement a way to filter for memrane genes during the cross
    # validation
    #
    # imp_genes <- rang_importances$signif_importances_ranger$gene[
    #     is_gene_membrane(rang_importances$signif_importances_ranger$gene)]
    #

    ranger_significance_table <-
        rang_importances$signif_importances_ranger

    imp_genes <- rang_importances$signif_importances_ranger$gene

    comm_genes <- imp_genes[imp_genes %in% rownames(test@assays[[test@active.assay]]@data)]
    removed_genes <- imp_genes[!imp_genes %in% comm_genes]

    if (warn.gene.removal & length(removed_genes) > 0) {
        warning(paste0(
            "Some important genes were removed because they are not present in ",
            "the test dataset. \n",
            "Removed genes: ", paste(removed_genes, collapse = ", ")
        ))
    }

    partyfit <- fit_ctree(train, comm_genes, cluster = cluster)

    concensusrules <- get_concensus_rules(partyfit)

    # TODO add method to support data frames as inputs

    testset <-
        as.data.frame.Seurat(test, comm_genes, fix_names = FALSE)
    predicted <- predict(partyfit, testset)

    gating_genes <- names(partykit::varimp(partyfit[[1]]))

    confusion_matrix <- table(data.frame(predicted = predicted,
                                         cluster = testset$ident))

    return(
        list(
            party_fit = partyfit,
            confusion_matrix = confusion_matrix,
            concensus_rules = concensusrules,
            ranger_significance_table = ranger_significance_table,
            gating_genes = gating_genes
        )
    )
}

