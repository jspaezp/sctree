

cross_validate <- function(train, test, cluster, ...) {
        rang_importances <- ranger_importances.seurat(
            train, cluster = cluster, genes_use = train@var.genes, ...)


        # TODO: implement a way to filter for memrane genes during the cross
        # validation
        #
        # imp_genes <- rang_importances$signif_importances_ranger$gene[
        #     is_gene_membrane(rang_importances$signif_importances_ranger$gene)]
        #

        ranger_significance_table <- rang_importances$signif_importances_ranger

        imp_genes <- rang_importances$signif_importances_ranger$gene

        comm_genes <- train@var.genes[train@var.genes %in% test@var.genes]
        comm_imp_genes <- comm_genes[make.names(comm_genes) %in% imp_genes]

        # this make.names nonsense is only because ranger needs r-valid names
        # which is not allways the case for gene names

        treedata <- as.data.frame.seurat(
            train, genes = train@var.genes[
                train@var.genes %in% comm_imp_genes],
            fix_names = TRUE)

        if (cluster == "ALL") {
            treedata$ident <- factor(treedata$ident)
        } else {
            classif_names <- c(
                FALSE. = paste0("not clus ", cluster),
                TRUE. = paste0("indeed clus ", cluster))

            treedata$ident <- factor(
                classif_names[make.names(treedata$ident == cluster)])
        }

        myformula <- as.formula("ident ~ .")
        partyfit <- partykit::ctree(formula = myformula, data = treedata)
        concensusrules <- get_concensus_rules(partyfit)
        plot(partyfit)

        testset <- as.data.frame.seurat(test, test@var.genes, fix_names = TRUE)
        predicted <- predict(partyfit, testset)

        gating_genes <- names(partykit::varimp(partyfit[[1]]))

        summary_table <- table(data.frame(predicted = predicted,
                                       cluster = testset$ident))
        summary_table <- 100 * (t(as.matrix(summary_table)) /
                                    colSums(summary_table))

        return(list(party_fit = partyfit, summary_table = summary_table,
                    concensus_rules = concensusrules,
                    ranger_significance_table = ranger_significance_table,
                    gating_genes = gating_genes))
}

# cross_validate(small_5050_mix, small_9901_mix, cluster = "0")
# cross_validate(small_9901_mix, small_5050_mix, cluster = "0")
# cross_validate(small_5050_mix, small_9901_mix, cluster = "ALL")
# cross_validate(small_9901_mix, small_5050_mix, cluster = "ALL")
