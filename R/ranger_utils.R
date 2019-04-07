
# TODO: rename functions to make clear that they are fitting the forest
# TODO: add argument to whitelist, blacklist or filter the genes to use
# TODO: restructure the function to get method dispatchment (optional)

#' ranger_importances
#' @title ranger_importances
#'
#' Calculates ranger-based variable importances for data frames and
#' seurat objects
#'
#' @param object a seurat or data frame object
#' @param cluster a cluster name for which the markers will be found
#' @param pval_cutoff p value cutoff for the markers
#' @param imp_method importance method, either of "janitza" or "altmann"
#' @param num.trees number of trees to be build using ranger
#' @param genes_use a character vector indicating which genes to use in
#'     the classification. currently implemented only for seurat objects.
#'     (for data frames one can simply subset the input data frame)
#' @param warn.imp.method logical indicating wether warning should be issued
#'     when few negative importances are found to calculate the p.values in
#'     ranger.
#' @param ... additional arguments to be passed to ranger
#'
#' @return  list with 3 elements ranger_fit, importances_ranger, signif_importances_ranger
#' @export
#'
#' @examples
#' summary(ranger_importances.seurat(Seurat::pbmc_small, cluster = "ALL"))
#' # Length Class      Mode
#' # ranger_fit                15     ranger     list
#' # importances_ranger        62     -none-     numeric
#' # signif_importances_ranger  3     data.frame list
#' #' summary(ranger_importances.seurat(Seurat::pbmc_small, cluster = "0"))
#' @importFrom ranger ranger importance_pvalues
ranger_importances.df <- function(object, cluster = NULL,
                                  pval_cutoff = 0.05,
                                  imp_method = c("janitza", "altmann"),
                                  num.trees = 500, warn.imp.method = TRUE,
                                  ...) {

    base_ranger <- function(data, importance) {
        ranger_fit <- ranger::ranger(
            dependent.variable.name = "ident",
            data = data,
            num.trees = num.trees,
            mtry = floor(ncol(data)/5),
            importance = importance,
            classification = TRUE,
            ...)
        return(ranger_fit)
    }

    tmp <- object

    if (cluster == "ALL") {
        tmp$ident <- factor(tmp$ident)
    } else {
        if (!is.null(cluster)) {
            stopifnot(cluster %in% unique(tmp$ident))
            tmp$ident <- factor(make.names(tmp$ident == cluster))
        }
    }

    if (all(imp_method == c("janitza", "altmann"))) imp_method <- "janitza"
    stopifnot(imp_method[1] %in% c("janitza", "altmann"))

    if (imp_method == "altmann") {
        ranger_fit <- base_ranger(data = tmp, importance = "permutation", ...)
        importances_ranger <- ranger::importance_pvalues(
            ranger_fit, "altman",
            formula = ident ~ ., data = tmp)

    } else if (imp_method == "janitza") {
        ranger_fit <- base_ranger(
            data = tmp,
            importance = "impurity_corrected", ...)

        if (!warn.imp.method) {
            suppressWarnings(suppressMessages({{
                importances_ranger <- ranger::importance_pvalues(ranger_fit)
            }}))
        } else {
            importances_ranger <- ranger::importance_pvalues(ranger_fit)
        }


    }

    signif_importances_ranger <- importances_ranger[
        importances_ranger[,'pvalue'] < pval_cutoff,]

    signif_importances_ranger <- signif_importances_ranger[
        order(signif_importances_ranger[,"importance"], decreasing = TRUE), ]

    signif_importances_ranger <- as.data.frame(signif_importances_ranger)
    signif_importances_ranger[["gene"]] <- rownames(signif_importances_ranger)

    return(list(ranger_fit = ranger_fit,
                importances_ranger = importances_ranger,
                signif_importances_ranger = signif_importances_ranger))
}



#' @describeIn ranger_importances.df Calculate variable importances to calssify a seurat object
#' @export
ranger_importances.seurat <- function(object, cluster = NULL,
                                      pval_cutoff = 0.05,
                                      imp_method = c("janitza", "altmann"),
                                      num.trees = 500,
                                      genes_use = object@var.genes,
                                      warn.imp.method = TRUE,
                                      ...) {

    tmp <- as.data.frame.seurat(object, genes = genes_use, fix_names = FALSE)

    return(ranger_importances.df(tmp,
                                 cluster = cluster,
                                 pval_cutoff = pval_cutoff,
                                 imp_method = imp_method,
                                 num.trees = num.trees,
                                 warn.imp.method = warn.imp.method,
                                 ...))
}
