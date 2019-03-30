
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
  }

  purrr::map(rules_per_clus, get_majority_rules)
}

