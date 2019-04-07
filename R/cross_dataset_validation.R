
# TODO refine what output is really needed
# TODO add option on how to handle missing values, either adding zeros or
#      filtering common names

#' Fits a decision tree in one data set and tests the performance in another
#'
#' @param train a seurat object to be used for trainning.
#' @param test another seurat object to be used for testing.
#' @param cluster the cluster whose equivalence needs to be found.
#' @param ... additional arguments to be passed to ranger_importances.seurat
#'
#' @return a list containing the (1) tree fit,
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
cross_validate <- function(train, test, cluster, ...) {
    rang_importances <- ranger_importances.seurat(train,
                                                  cluster = cluster,
                                                  genes_use = train@var.genes,
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

    comm_genes <-
        train@var.genes[train@var.genes %in% test@var.genes]
    comm_imp_genes <-
        comm_genes[make.names(comm_genes) %in% imp_genes]

    # this make.names nonsense is only because ranger needs r-valid names
    # which is not allways the case for gene names

    treedata <- as.data.frame.seurat(
        train,
        genes = train@var.genes[train@var.genes %in% comm_imp_genes],
        fix_names = TRUE)

    # TODO make this its own function and place it in utils ...
    # something like "clean_ident_name"
    if (cluster == "ALL") {
        treedata$ident <- factor(treedata$ident)
    } else {
        classif_names <- c(
            # TODO add a way to bundle clusters ...
            # or to give a vector of potential clusters
            FALSE. = paste0("not clus ", cluster),
            TRUE. = paste0("indeed clus ", cluster)
        )

        treedata$ident <-
            factor(classif_names[make.names(treedata$ident == cluster)])
    }

    myformula <- as.formula("ident ~ .")
    # TODO add option to have another column name as the classifier
    partyfit <-
        partykit::ctree(formula = myformula, data = treedata)
    concensusrules <- get_concensus_rules(partyfit)

    # TODO add method to support data frames as inputs

    testset <-
        as.data.frame.seurat(test, comm_imp_genes, fix_names = TRUE)
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

#' Check if the element has the structure of a confusion matrix
#'
#' Slim check for the element to assure it has 2 dimensions,
#' is a table and whose elements are exclusively integers
#'
#' @param x table object
#'
#' @return logical indicating wether the element matches the description
#' @export
#'
#' @examples
#' is.confusion.matrix(table(1:3, 3:1, 1:3)) # FALSE
#' is.confusion.matrix(table(1:3, 3:1)) # TRUE
is.confusion.matrix <- function(x) {
    return(length(dim(x)) == 2 & is.table(x) & is.integer(x))
}

#' Convert confusion matrices and tables to frequency matrices
#'
#' @param confusion_matrix
#'
#' @return a frequency matrix ormalizing by the column total
#' @export
#'
#' @examples
#' table(1:3, 3:1)
#'
#' # 1 2 3
#' # 1 0 0 1
#' # 2 0 1 0
#' # 3 1 0 0
#' as.frequency.matrix(table(1:3, 3:1))
#' #
#' # 1   2   3
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


#' Heatmap with ggplot from table objects
#'
#' @param tbl a table object
#' @param min_color minimum value to show colored, value sunder this will
#'     be shown as gray
#' @param show_number logical indicating wether the number should be printed
#'     in the box of the heatmap
#'
#' @return ggpplot object
#' @export
#'
#' @examples
#' autoplot(table(1:10, 1:10), show_number = TRUE)
#' @importFrom ggplot2 geom_tile geom_point aes_string labs
autoplot.table <- function(tbl,
                           min_color = NULL,
                           show_number = FALSE) {
    df <- as.data.frame(tbl)

    df[["Freq"]] <- round(df[["Freq"]], 2)

    stopifnot(length(dim(tbl)) == 2)

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


