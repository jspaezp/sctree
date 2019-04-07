
#' get_cluster_mapping
#'
#' returns a named vector, the names being the terminal node id in the tree
#' and the values being the classification it corresponds to
#'
#' @param tree a partykit::ctree fit
#'
#' @return a named vector
#' @export
#'
#' @examples
#' iris_tree <- partykit::ctree(Species ~ ., data = iris)
#' get_cluster_mapping(iris_tree)
#' # 2          5          6          7
#' # setosa versicolor versicolor  virginica
#' # Levels: setosa versicolor virginica
get_cluster_mapping <- function(tree) {
    # list of factor vectors, each element is a node and each elem in the
    # factor is the majority prediction
    predictions <- partykit::nodeapply(
      partykit::as.simpleparty(tree),
      ids = partykit::nodeids(tree, terminal = TRUE),
      function(x) partykit::info_node(x)$prediction)

    # Named vector, each name the node id and each element the predicted cluster
    cluster_mapping <- wrapr::named_map_builder(
      names(predictions),
      unlist(predictions))

    return(cluster_mapping)
}

# TODO, there is a known issue where suboptimas classification trees return
# 3: In max(elems_per_rule) :
# no non-missing arguments to max; returning -Inf

#' get_concensus_rules
#'
#' Returns a nested list of the concensus rules to define a cluster, given a
#' partykit decision tree.
#'
#' level 1 of the nested list is the group being predicted.
#' level 2 is either "all" or "majority".
#' each element is a character vector.
#'
#' @param tree partykit::ctree fit decision tree
#'
#' @return a nested list
#' @export
#'
#' @examples
#' iris_tree <- partykit::ctree(Species ~ ., data = iris)
#' str(get_concensus_rules(iris_tree))
#' # List of 3
#' # $ setosa    :List of 2
#' # ..$ all     : chr "Petal.Length <= 1.9"
#' # ..$ majority: chr(0)
#' # $ versicolor:List of 2
#' # ..$ all     : chr [1:2] "Petal.Length > 1.9" "Petal.Width <= 1.7"
#' # ..$ majority: chr "Petal.Length <= 4.8"
#' # $ virginica :List of 2
#' # ..$ all     : chr [1:2] "Petal.Length > 1.9" "Petal.Width > 1.7"
#' # ..$ majority: chr(0)
#'
#' diamonds_tree <- partykit::ctree(cut ~ ., data = ggplot2::diamonds)
#' str(get_concensus_rules(diamonds_tree))
#' # List of 5
#' # $ Fair     :List of 2
#' # ..$ all     : chr [1:2] "depth > 63" "depth > 64.3"
#' # ..$ majority: chr [1:3] "table <= 57" "table > 57" "table <= 62"
#' # $ Good     :List of 2
#' # ..$ all     : chr [1:2] "depth > 63" "depth <= 64.3"
#' # ..$ majority: chr [1:4] "table <= 57" "depth > 63.5" "table > 57" "table <= 62"
#' # $ Very Good:List of 2
#' # ..$ all     : chr [1:3] "depth > 63" "depth <= 64.3" "depth <= 63.5"
#' # ..$ majority: chr [1:3] "table <= 57" "table > 57" "table <= 62"
#' # $ Premium  :List of 2
#' # ..$ all     : chr [1:3] "table > 57" "table <= 62" "depth <= 63"
#' # ..$ majority: chr [1:3] "table <= 60" "depth > 58" "color %in% c(\"E\", \"G\", \"I\")"
#' # $ Ideal    :List of 2
#' # ..$ all     : chr [1:2] "table <= 57" "depth <= 63"
#' # ..$ majority: chr [1:2] "table <= 56.4" "z > 2.64"
#' @importFrom partykit as.simpleparty nodeids info_node
# @importFrom partykit as.simpleparty nodeids info_node .list.rules.party
#' @importFrom purrr map2_dbl
get_concensus_rules <-  function(tree) {
  # TODO the logic on this function is really hard to understand, we should
  # break it in several smaller ones ... and definitely find better names ...

  cluster_mapping <- get_cluster_mapping(tree)

  # List of lists of length 1, each sub-list being a named numeric vector
  number_per_node <- partykit::nodeapply(
    partykit::as.simpleparty(tree),
    ids = partykit::nodeids(tree, terminal = TRUE),
    function(x) partykit::info_node(x)$distribution)

  # Named vector, each name the node id and each element number of elements for
  # that predicted cluster
  number_per_node <- purrr::map2_dbl(
    number_per_node,
    cluster_mapping,
    .f = function(x,y) {
      x[[y]][[1]]} )

  # Character vector of the rules
  tree_rules <- partykit:::.list.rules.party(tree)


  # This is a list whose names are the cluster names and the elements the
  # name/number of the nodes
  nodes_per_clus <- lapply(split(cluster_mapping, cluster_mapping), names)

  # This returns the (overwhelming) rules for each cluster
  rules_per_clus <- lapply(nodes_per_clus, function(x) tree_rules[x])


  get_majority_rules <- function(rules_in_cluster) {
    split_rules <- strsplit(rules_in_cluster, "\\s&\\s")

    all_rules <- unique(unlist(split_rules))

    count_per_rule <- function(nodes_logical, number_per_node) {
      ### NOTE START
      ### basically counts how many cells are being classified
      ### in each terminal node
      ### Example imputs
      # > nodes_logical
      #   19   20   21   24   25   26
      # TRUE TRUE TRUE TRUE TRUE TRUE
      # > number_per_node
      #  19  20  21  24  25  26
      #  59 571  15  48   7   5
      ### NOTE END
      sum(number_per_node[names(nodes_logical)[which(nodes_logical)]])
      }

    .in <- function(x, y) {
      x %in% y
    }

    elems_per_rule <- purrr::map(
      all_rules, # TODO remove this annonimous function ... could use outer instead of a nested lapply ...
      function(x) purrr::map_lgl(split_rules, function(y) { x %in% y }))
    elems_per_rule <- purrr::map_dbl(
      elems_per_rule,
      count_per_rule,
      number_per_node)

    rules_in_majority <- all_rules[
      elems_per_rule > (max(elems_per_rule)/2) &
        elems_per_rule != max(elems_per_rule)]
    rules_in_all <- all_rules[elems_per_rule == max(elems_per_rule)]

    list(all = rules_in_all, majority = rules_in_majority)
    # TODO add a class to this and make a print method that prints it in
    # biologist-friendly form
  }

  purrr::map(rules_per_clus, get_majority_rules)
}


#' Fits a classification tree on a seurat object
#'
#' @param object a seurat object
#' @param genes_use a character vector indicating which genes to use in
#'     the classification. currently implemented only for seurat objects.
#'     (for data frames one can simply subset the input data frame)
#' @param cluster a cluster name for which the markers will be found
#' @param ... additional arguments to be passed to partykit::ctree_control
#'
#' @return a ctree fit
#' @export
#'
#' @examples
#' fit_ctree(small_9901_mix, c("CCNB1", "PLK1", "AURKA"), cluster = "ALL")
#' fit_ctree(small_9901_mix, c("CCNB1", "PLK1", "AURKA"), cluster = "0")
#' # Model formula:
#' #   ident ~ CCNB1 + PLK1 + AURKA
#' #
#' # Fitted party:
#' #   [1] root
#' # |   [2] PLK1 <= 2.31327
#' # |   |   [3] CCNB1 <= 3.01197
#' # |   |   |   [4] PLK1 <= 1.75363: indeed clus 0 (n = 202, err = 4.0%)
#' # |   |   |   [5] PLK1 > 1.75363: indeed clus 0 (n = 45, err = 20.0%)
#' # |   |   [6] CCNB1 > 3.01197
#' # |   |   |   [7] PLK1 <= 1.67634: indeed clus 0 (n = 26, err = 38.5%)
#' # |   |   |   [8] PLK1 > 1.67634: not clus 0 (n = 23, err = 17.4%)
#' # |   [9] PLK1 > 2.31327
#' # |   |   [10] AURKA <= 2.02766: not clus 0 (n = 21, err = 47.6%)
#' # |   |   [11] AURKA > 2.02766: not clus 0 (n = 67, err = 4.5%)
#' #
#' # Number of inner nodes:    5
#' # Number of terminal nodes: 6
#' fit_ctree(small_9901_mix, c("CCNB1", "PLK1", "AURKA"), cluster = "0", maxdepth = 2)
#' #
#' # Model formula:
#' #   ident ~ CCNB1 + PLK1 + AURKA
#' #
#' # Fitted party:
#' #   [1] root
#' # |   [2] PLK1 <= 2.31327
#' # |   |   [3] CCNB1 <= 3.01197: indeed clus 0 (n = 247, err = 6.9%)
#' # |   |   [4] CCNB1 > 3.01197: not clus 0 (n = 49, err = 40.8%)
#' # |   [5] PLK1 > 2.31327
#' # |   |   [6] AURKA <= 2.02766: not clus 0 (n = 21, err = 47.6%)
#' # |   |   [7] AURKA > 2.02766: not clus 0 (n = 67, err = 4.5%)
#' #
#' # Number of inner nodes:    3
#' # Number of terminal nodes: 4
fit_ctree <- function(object,
                      genes_use = object@var.genes,
                      cluster = "ALL", ...) {
    treedata <- as.data.frame.seurat(
      object, genes = genes_use, fix_names = FALSE)

    # TODO make this its own function and place it in utils ...
    # something like "clean_ident_name"
    # right now it is also hard coded in cross_dataset_validation.R
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
    # TODO add option to have another column name as the classifier
    partyfit <- partykit::ctree(
      formula = myformula,
      data = treedata,
      control = partykit::ctree_control(...))
    return(partyfit)
}
