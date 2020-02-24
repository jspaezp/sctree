
# TODO: rename functions to make clear that they are fitting the forest
# TODO: add argument to whitelist, blacklist or filter the genes to use
# TODO: restructure the function to get method dispatchment (optional)


#' @title ranger_importances
#'
#' @description
#' Calculates ranger-based variable importances for data frames and
#' Seurat objects
#'
#' @param object a Seurat or data frame object
#' @param cluster a cluster name for which the markers will be found
#' @param pval_cutoff p value cutoff for the markers
#' @param imp_method importance method, either of "janitza" or "altmann"
#' @param num.trees number of trees to be build using ranger
#' @param return_what a subset of "ranger_fit", "importances_ranger",
#'     "signif_importances_ranger", defaults to signif_importances_ranger
#' @param genes_use a character vector indicating which genes to use in
#'     the classification. currently implemented only for Seurat objects.
#'     (for data frames one can simply subset the input data frame)
#' @param warn.imp.method logical indicating wether warning should be issued
#'     when few negative importances are found to calculate the p.values in
#'     ranger.
#' @param return ranger_fit, importances_ranger, signif_importances_ranger
#' @param ... additional arguments to be passed to ranger
#'
#' @return by default returns a data frame with the importances and p values
#'     but this behavior can be modified by the
#' @export
#'
#' @evalRd include_roxygen_example("
#'     head(ranger_importances(Seurat::pbmc_small, cluster = 'ALL', warn.imp.method = FALSE))
#'     ")
#'
#' @importFrom ranger ranger importance_pvalues
ranger_importances <- function(x, ...){
    UseMethod("ranger_importances")
}

#' @export
#' @method ranger_importances data.frame
ranger_importances.data.frame <- function(object, cluster = NULL,
                                  pval_cutoff = 0.05,
                                  imp_method = c("janitza", "altmann"),
                                  num.trees = 500,
                                  warn.imp.method = TRUE,
                                  identity_col_name = "ident",
                                  return_what = "signif_importances_ranger",
                                  sort_results = TRUE,
                                  ...) {

    # TODO find a better name for the `return_what` argument
    redirect_altmann_warn <- function(w) {
        warn_message1 <- paste0(
            "Only few negative importance values found, ",
            "inaccurate p-values. Consider the 'altmann' approach.")
        warn_message2 <- paste0(
            "No negative importance values found, ",
            "Consider the 'altmann' approach.")

        new_warning <- paste0(
            "Only few negative importance values found, ",
            "inaccurate p-values. Consider the 'altmann' approach.\n\n",
            "This can be done by setting the argument 'imp_method' to ",
            "'altmann', note that this method is extremely computationally ",
            "intensive.\n\n",
            "This warning can be disabled by setting the argument",
            " `warn.imp.method` to `FALSE`\n\n",
            "For more information please refer to ?ranger::ranger")

        if (grepl(paste(warn_message1, warn_message2, sep = "|"), w$message)) {
            warning(
                new_warning,
                call. = FALSE,
                immediate. = FALSE)
            invokeRestart("muffleWarning")
        }
    }

    base_ranger <- function(data, importance, ...) {
        ranger_fit <- ranger::ranger(
                dependent.variable.name = identity_col_name,
                data = data,
                num.trees = num.trees,
                mtry = floor(ncol(data)/5),
                importance = importance,
                classification = TRUE,
                ...)

        return(ranger_fit)
    }

    tmp <- object

    stopifnot(imp_method[1] %in%
                  c("janitza", "altmann"))
    stopifnot(return_what[1] %in%
                  c("ranger_fit", "importances_ranger",
                    "signif_importances_ranger"))

    if (is.null(cluster)) {
        stop("Please Specify a cluster to fit the forest (or assign cluster=\"ALL\")")
    }

    if (cluster == "ALL") {
        tmp$ident <- factor(tmp$ident)
    } else {
        if (!is.null(cluster)) {
            stopifnot(cluster %in% unique(tmp$ident))
            tmp$ident <- factor(make.names(tmp$ident == cluster))
        }
    }

    if (all(imp_method == c("janitza", "altmann"))) imp_method <- "janitza"

    if (imp_method == "altmann") {
        # This name shift seems to be necessary because the altman in ranger
        # method does not implement a dependent variable interface
        oldnames <- colnames(tmp)
        newnames <- make.names(oldnames)

        colnames(tmp) <- newnames
        ranger_fit <- base_ranger(data = tmp, importance = "permutation", ...)
        importances_ranger <- ranger::importance_pvalues(
            ranger_fit, "altmann",
            formula = "ident ~ .", data = tmp)

        stopifnot(all.equal(newnames[!newnames == "ident"],
                            rownames(importances_ranger)))

        rownames(importances_ranger) <- oldnames[!oldnames == "ident"]

    } else if (imp_method == "janitza") {
        ranger_fit <- base_ranger(
            data = tmp,
            importance = "impurity_corrected", ...)

        if (!warn.imp.method) {
            suppressWarnings(suppressMessages({{
                importances_ranger <- ranger::importance_pvalues(ranger_fit)
            }}))
        } else {
            withCallingHandlers({
                importances_ranger <- ranger::importance_pvalues(ranger_fit)
            }, warning = redirect_altmann_warn)
        }


    }

    importances_ranger <- as.data.frame(importances_ranger)
    importances_ranger[["gene"]] <- rownames(importances_ranger)

    signif_importances_ranger <- importances_ranger[
        importances_ranger[,'pvalue'] < pval_cutoff,]

    if (sort_results) {
        signif_importances_ranger <- signif_importances_ranger[
            order(signif_importances_ranger[,"importance"], decreasing = TRUE), ]
    }

    posible_returns <- list(
        ranger_fit = ranger_fit,
        importances_ranger = importances_ranger,
        signif_importances_ranger = signif_importances_ranger)

    ret <- posible_returns[return_what]
    if (length(ret) == 1 & is.list(ret) & !is.data.frame(ret)) ret <- ret[[1]]

    return(ret)
}


#' @export
#' @method ranger_importances Seurat
#' @importFrom Seurat VariableFeatures
ranger_importances.Seurat <- function(object, cluster = NULL,
                                      pval_cutoff = 0.05,
                                      imp_method = c("janitza", "altmann"),
                                      num.trees = 500,
                                      genes_use = Seurat::VariableFeatures(object),
                                      warn.imp.method = TRUE,
                                      ...) {

    tmp <- as.data.frame.Seurat(object, genes = genes_use, fix_names = FALSE)

    return(ranger_importances.data.frame(
        tmp,
        cluster = cluster,
        pval_cutoff = pval_cutoff,
        imp_method = imp_method,
        num.trees = num.trees,
        warn.imp.method = warn.imp.method,
        ...))
}


# TODO, make functions have the same interface as Seurat ...
#' @describeIn ranger_importances Calculate variable importances to each cluster in a Seurat object
#' @export
FindAllMarkers_ranger.Seurat <- function(object,
                                         genes_use = Seurat::VariableFeatures(object),
                                         ...) {

    tmp <- as.data.frame.Seurat(object, genes = genes_use, fix_names = FALSE)

    object_ranger_importances.Seurat <- function(cluster, ...) {
        ranger_importances.data.frame(
            tmp, cluster = cluster,
            return_what = "signif_importances_ranger", ...)
    }

    clusters <- sort(as.character(unique(tmp$ident)))

    results <- purrr::map(clusters,
                          .f = object_ranger_importances.Seurat,
                          ...)

    results <- purrr::map2(results, clusters, function(x, y) {
        x$cluster <- y
        return(x)})

    results <- do.call(what = rbind, args = results)

    return(results)
}


#' RangerDE
#'
#' A Helper function that interfaces with Seurat::FindMarkers
#'
#' @param data.use sparse matrix passed by seurat
#' @param cells.1 barcode names of the cells asigned as cluster 1
#' @param cells.2 barcode names of the cells asigned as cluster 2
#' @param verbose logical indicating wether the run should be verbose
#' @param ... additional arguments passed to sctree::ranger_importances.df
#'
#' @evalRd include_roxygen_example({
#'     "library(Seurat)
#'     library(sctree)
#'     head(FindMarkers(object = Seurat::pbmc_small,
#'     ident.1 = 0, test.use = 'RangerDE', warn.imp.method = FALSE))"
#'     })
#' @importFrom Matrix t
#' @keywords internal
RangerDE <- function(data.use, cells.1, cells.2, verbose, ...) {
    df <- as.data.frame(Matrix::t(data.use))
    df$ident <- as.character(rownames(df) %in% cells.1)
    ret <- ranger_importances.data.frame(
        df, cluster = "TRUE", verbose = verbose,
        return_what = "importances_ranger", ...)
    names(ret) <- gsub("pvalue", "p_val", names(ret))
    return(ret)
}
