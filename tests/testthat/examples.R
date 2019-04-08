# File generated automatically, do not run by hand
### Name: as.data.frame.seurat
### Title: Gets a data frame from a seurat object
### Aliases: as.data.frame.seurat

### ** Examples

as.data.frame(Seurat::pbmc_small, Seurat::pbmc_small@var.genes)[1:3,1:3]



### Name: as.frequency.matrix
### Title: Convert confusion matrices and tables to frequency matrices
### Aliases: as.frequency.matrix

### ** Examples

table(1:3, 3:1)

as.frequency.matrix(table(1:3, 3:1))
#
as.frequency.matrix(table(c(3,3,2), 3:1))
#



### Name: autoplot.table
### Title: autoplot.table
### Aliases: autoplot.table

### ** Examples

require(ggplot2)
autoplot(table(1:10, 1:10), show_number = TRUE)



### Name: cross_validate
### Title: Fits a decision tree in one data set and tests the performance
###   in another
### Aliases: cross_validate

### ** Examples

cross_validate(small_5050_mix, small_9901_mix, cluster = "0")
cross_validate(small_5050_mix, small_9901_mix, cluster = "ALL")



### Name: featureplot_gadget
### Title: featureplot_gadget
### Aliases: featureplot_gadget

### ** Examples

## No test: 
##D featureplot_gadget(seurat_object = small_5050_mix,
##D                     genes.use = small_5050_mix@var.genes[
##D                         is_gene_membrane(small_5050_mix@var.genes)],
##D                     cache = "./.cache")
## End(No test)



### Name: fit_ctree
### Title: Fits a classification tree on a seurat object
### Aliases: fit_ctree

### ** Examples

fit_ctree(small_9901_mix, c("CCNB1", "PLK1", "AURKA"), cluster = "ALL")
fit_ctree(small_9901_mix, c("CCNB1", "PLK1", "AURKA"), cluster = "0")
#
#
fit_ctree(small_9901_mix, c("CCNB1", "PLK1", "AURKA"), cluster = "0", maxdepth = 2)
#
#
#



### Name: get_aliases
### Title: get_aliases
### Aliases: get_aliases

### ** Examples

get_aliases(c("MAPK1", "CD4"))
#



### Name: get_cluster_mapping
### Title: get_cluster_mapping
### Aliases: get_cluster_mapping

### ** Examples

iris_tree <- partykit::ctree(Species ~ ., data = iris)
get_cluster_mapping(iris_tree)



### Name: get_concensus_rules
### Title: get_concensus_rules
### Aliases: get_concensus_rules

### ** Examples

iris_tree <- partykit::ctree(Species ~ ., data = iris)
str(get_concensus_rules(iris_tree))

diamonds_tree <- partykit::ctree(cut ~ ., data = ggplot2::diamonds)
str(get_concensus_rules(diamonds_tree))



### Name: get_genesymbols
### Title: get_genesymbols
### Aliases: get_genesymbols

### ** Examples

get_genesymbols("ERK")
get_genesymbols(c("SUPERFAKEGENE", "ERK"))
#



### Name: ggheatmap_base
### Title: Basic heatmap using ggplot
### Aliases: ggheatmap_base

### ** Examples

my_df <- data.frame(x = 1:2, y = 1:2, z = 1:2)
ggheatmap_base(my_df, x = "x", y = "y", value = "z", show_number = TRUE)
ggheatmap_base(my_df, x = "x", y = "y", value = "z", show_number = FALSE)



### Name: is.confusion.matrix
### Title: Check if the element has the structure of a confusion matrix
### Aliases: is.confusion.matrix

### ** Examples

is.confusion.matrix(table(1:3, 3:1, 1:3)) # FALSE
is.confusion.matrix(table(1:3, 3:1)) # TRUE



### Name: is_gene_membrane
### Title: is_gene_membrane
### Aliases: is_gene_membrane

### ** Examples

is_gene_membrane(c("CD4", "GAPDH"))



### Name: plot_flowstyle
### Title: plot_flowstyle
### Aliases: plot_flowstyle plot_flowstyle.data.frame plot_flowstyle.seurat

### ** Examples

plot_flowstyle(Seurat::pbmc_small, c("ACRBP", "TSC22D1", "VDAC3"))



### Name: print.concensus.rules
### Title: Prints concensus rules
### Aliases: print.concensus.rules

### ** Examples

iris_tree <- partykit::ctree(Species ~ ., data = iris)
my_rules <- get_concensus_rules(iris_tree)
print(my_rules)



### Name: query_biocompare_antibodies
### Title: Query Biocompare for antibodies
### Aliases: query_biocompare_antibodies

### ** Examples

query_biocompare_antibodies("CD11bfakename")
head(query_biocompare_antibodies("CD11C"),3)



### Name: query_biolegend_antibodies
### Title: Query Biolegend for antibodies
### Aliases: query_biolegend_antibodies

### ** Examples

query_biolegend_antibodies("CD11bfakename")
head(query_biolegend_antibodies("CD11C"))



### Name: query_cc_antibodies
### Title: Query cell signaling for Single cell antibodies
### Aliases: query_cc_antibodies

### ** Examples

query_cc_antibodies("CD11bfakename")
query_cc_antibodies("CD11c")



### Name: query_sc_antibodies
### Title: Query Santa Cruz for antibodies
### Aliases: query_sc_antibodies

### ** Examples

query_sc_antibodies("CD11bfakename")
NULL
head(query_sc_antibodies("CD11C"))



### Name: ranger_importances.df
### Title: ranger_importances Calculates ranger-based variable importances
###   for data frames and seurat objects
### Aliases: ranger_importances.df ranger_importances.seurat

### ** Examples

summary(ranger_importances.seurat(Seurat::pbmc_small, cluster = "ALL"))
#' summary(ranger_importances.seurat(Seurat::pbmc_small, cluster = "0"))



### Name: top_n
### Title: top_n
### Aliases: top_n

### ** Examples

top_n(iris, 2, "Sepal.Width")
top_n(iris, -2, "Sepal.Width")



### Name: tsne_plot
### Title: TSNE plot generation
### Aliases: tsne_plot

### ** Examples

tsne_plot(Seurat::pbmc_small)



