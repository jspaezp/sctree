
#' Gets the corresponding cluster for each terminl node in a classification tree
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
#' my_tree <- fit_ctree(small_5050_mix, c('CKB', 'ADA', 'ASNS', 'PRDX6', 'MZB1'))
#' get_cluster_mapping(my_tree)
#'
#' @evalRd paste("# ", capture.output(
#'     get_cluster_mapping(fit_ctree(small_5050_mix, c("CKB", "ADA", "ASNS", "PRDX6", "MZB1")))
#'     ))
get_cluster_mapping <- function(tree) {
    # list of factor vectors, each element is a node and each elem in the
    # factor is the majority prediction
    predictions <- partykit::nodeapply(
      partykit::as.simpleparty(tree),
      ids = partykit::nodeids(tree, terminal = TRUE),
      function(x) partykit::info_node(x)$prediction)

    # Named vector, each name the node id and each element the predicted cluster
    cluster_mapping <- wrapr::named_map_builder(
      names(predictions), # This is the node ids
      unlist(predictions)) # And this is the prediction (identity)

    return(cluster_mapping)
}


#' Converts a classification tree to a garnett output
#'
#' Converts a fitted ctree object to the classifier format supported by garnett
#' for more information on the specification please visit:
#' https://cole-trapnell-lab.github.io/garnett/docs/#submitting-a-classifier
#'
#' @param tree a fitted constparty/party object
#' @param digits an integer stating how many decimal places should be conserved
#' @param rules_keep a regular expression stating which rules should be kept,
#'     defaults to all rules (".*")
#'
#' @return a complex list that prints as a garnett file
#' @export
#'
#' @evalRd include_roxygen_example({
#'     "as.garnett(fit_ctree(Seurat::pbmc_small))"
#'     })
as.garnett <- function(tree, digits = 3, rules_keep = ".*") {
  cluster_mappings <- get_cluster_mapping(tree)
  node_ids <- names(cluster_mappings)
  cluster_mappings <- as.character(cluster_mappings)

  dup_names <- cluster_mappings %in%
    cluster_mappings[duplicated(cluster_mappings)]

  classif_in_many <- cluster_mappings[dup_names]
  replacement_names <- paste0(
    as.character(classif_in_many),
    "_node_",
    node_ids[dup_names])

  cluster_mappings[dup_names] <- replacement_names

  num_elements <- unlist(
    partykit::nodeapply(
      partykit::as.simpleparty(tree),
      partykit::nodeids(tree, terminal = TRUE),
      FUN = function(x) partykit::info_node(x)$n))

  cluster_mappings <- paste(
    cluster_mappings,
    paste0("\t(n = ", num_elements, ")"))

  rule_list <- partykit:::.list.rules.party(tree)
  split_rules <- strsplit(rule_list, "\\s&\\s")
  names(split_rules) <- cluster_mappings


  group_sub_rules <- function(split_rules) {

    rules_per_var <- split(
      split_rules,
      gsub(" (>|<).*$", "", split_rules))

    get_min_and_max_rule <- function(rule_group) {

      gt_lt_rules <- split(
        rule_group,
        gsub("(^.*)(>|<)(.*$)", "\\2", rule_group))

      get_nums <- function(x) {
        ret <- as.numeric(gsub(".* (-?[0-9.]+)\\s*$", "\\1", x, perl = TRUE))
        ret <- round(ret, digits = digits)
        return(ret)
      }

      if (length(gt_lt_rules[[">"]] > 0)) {
        max_gt <- max(get_nums(gt_lt_rules[[">"]]))
      } else {
        max_gt <- NA
      }

      if (length(gt_lt_rules[["<"]] > 0)) {
        min_lt <- min(get_nums(gt_lt_rules[["<"]]))
      } else {
        min_lt <- NA
      }

      ret <- c(min_lt, max_gt)

      ruletype <- c("TRUE_FALSE" = "gt",
                    "FALSE_FALSE" = "mixed",
                    "FALSE_TRUE" = "lt")[
        paste0(is.na(ret), collapse = "_")]

      ret <- ret[!is.na(ret)]

      attr(ret, "ruletype") <- ruletype

      return(ret)
    }

    max_n_mins <- lapply(rules_per_var, get_min_and_max_rule)
    names(max_n_mins) <- names(rules_per_var)

    ruletypes <- sapply(max_n_mins, function(x) attr(x, "ruletype"))

    ret_vars <- split(max_n_mins, ruletypes)

    list(`expressed above` = ret_vars[["gt"]],
         `expressed below` = ret_vars[["lt"]],
         `expressed between` =  ret_vars[["mixed"]])
  }

  garnett.list <- lapply(split_rules, group_sub_rules)

  garnett.list <- garnett.list[grepl(pattern = rules_keep, names(garnett.list))]

  class(garnett.list) <- "garnett.list"

  return(garnett.list)
}


#' @describeIn as.garnett Print a garnett classifier
#' @export
print.garnett.list <- function(x) {
  tmpfile <- tempfile()

  cat <- function(..., file = tmpfile) base::cat(..., file = file, append = TRUE)
  x <- x[sort(names(x))]
  for (entry in names(x)) {

    cat(paste0("> ", entry, "\n"))
    # For instance : "> clus 1 	(n = 54)\n"

    for (expression_type in names(x[[entry]])) {
      if (length(x[[entry]][[expression_type]]) == 0)  next()

      cat(paste0(expression_type, ": "))
      # for isntance: "expressed above: "

      pasted_expressions <- sapply(
        x[[entry]][[expression_type]],
        function(x) paste(x, collapse = " "))

      cat(paste(names(x[[entry]][[expression_type]]),
                pasted_expressions,
                collapse = ", "))
      # for instance:  "ASNS 1.432, ITGA4 0"

      cat("\n")
    }
    if (entry == names(x)[length(x)]) next()
    cat("\n")
  }

  ret <- readLines(tmpfile)
  unlink(tmpfile)

  cat(ret, file = "", sep = "\n")
  return(invisible(ret))
}


# TODO, there is a known issue where suboptimal classification trees return
# 3: In max(elems_per_rule) :
# no non-missing arguments to max; returning -Inf

#' Consolidates the rules in a decision tree per cluster
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
#' @evalRd include_roxygen_example("
#'     iris_tree <- partykit::ctree(Species ~ ., data = iris)
#'     str(get_concensus_rules(iris_tree))
#'     ")
#'
#' @importFrom partykit as.simpleparty nodeids info_node
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

  concensus_rules <- purrr::map(rules_per_clus, get_majority_rules)
  class(concensus_rules) <- c(class(concensus_rules), "concensus.rules")
  return(concensus_rules)
}

# TODO write `is.concensus.rules` function.
# ... a 3 level list with whared names in element 2 and only items interpretable
# as rules in level 3 ...

#' Prints concensus rules
#'
#' Prints in a human readable format a concensus.rules object, usually an output
#' of `sctree::get_concensus_rules`
#'
#' @param x a concensus.rules object output from `sctree::get_concensus_rules`
#'
#' @return returns silently the same object but prints the rules to the console.
#' @method print concensus.rules
#' @seealso get_concensus_rules
#' @export
#'
#' @evalRd include_roxygen_example("
#'     my_tree <- fit_ctree(small_5050_mix, c('CKB', 'ADA', 'ASNS', 'PRDX6', 'MZB1'))
#'     my_rules <- get_concensus_rules(my_tree)
#'     print(my_rules)
#'     ")
#'
print.concensus.rules <- function(x) {
  verbose_rules <- rapply(
    x,
    (function(y) gsub(" (>).*", " +", gsub(" (<=).*", " -", y))),
    classes = "character",
    how = "list")

  for (rule_name in names(verbose_rules)) {
    cat(paste0("Cluster-", rule_name, ": "), sep = "\n")
    rule <- verbose_rules[[rule_name]]
    for (ruletype in names(rule)) {
      sub_rules <- rule[[ruletype]]
      if (length(sub_rules) == 0) next()
      cat(paste0("\t", ruletype, " elements:\n"))
      cat(paste0("\t\t", sub_rules), sep = "\n")
    }
  }
  return(invisible(x))
}


#' Fits a classification tree on a Seurat object
#'
#' @param object a Seurat object
#' @param genes_use a character vector indicating which genes to use in
#'     the classification. currently implemented only for Seurat objects.
#'     (for data frames one can simply subset the input data frame)
#'     defaults to Seurat::VariableFeatures(object)
#' @param cluster a cluster name for which the markers will be found
#' @param ... additional arguments to be passed to partykit::ctree_control
#'
#' @return a ctree fit
#' @export
#'
#' @evalRd include_roxygen_example("
#'    fit_ctree(small_9901_mix, c('CCNB1', 'PLK1', 'AURKA'), cluster = 'ALL')
#'    fit_ctree(small_9901_mix, c('CCNB1', 'PLK1', 'AURKA'), cluster = '0')
#'    ")
#'
#' @importFrom Seurat VariableFeatures
#' @importFrom partykit ctree_control ctree
fit_ctree <- function(object,
                      genes_use = Seurat::VariableFeatures(object),
                      cluster = "ALL", ...) {
    treedata <- as.data.frame.Seurat(
      object, genes = genes_use, fix_names = FALSE)

    # TODO make this its own function and place it in utils ...
    # something like "clean_ident_name"
    # right now it is also hard coded in cross_dataset_validation.R
    if (cluster == "ALL") {
      treedata$ident <- factor(treedata$ident)
    } else {
      classif_names <- c(
        FALSE. = paste0("not clus ", cluster),
        TRUE. = paste0("clus ", cluster))

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


#' Plot Decision trees as gates
#'
#'
#' @param object a seurat object that will be used to draw the plots
#' @param tree a decision tree that will be used to construct the gates
#' @param terminal_node the name of the terminal node that will be used for the gating
#'
#' @return a ggplot object.
#' @export
#'
#' @examples
#' plot_gates(Seurat::pbmc_small, fit_ctree(Seurat::pbmc_small), '5')
#'
#' @importFrom ggplot2 geom_point ggplot aes_string
#' @importFrom partykit varimp
#' @importFrom cowplot plot_grid
#' @importFrom grDevices colorRampPalette
plot_gates <- function(object, tree, terminal_node) {

  plot_gate_subset <- function(data_subset,
                               next_expr,
                               cutoff_var,
                               next_cutoff_var) {
    data_subset$in_next <-
      with(data_subset, eval(parse(text = next_expr)))
    g <- ggplot2::ggplot(
      data_subset,
      ggplot2::aes_string(
        x = cutoff_var,
        y = next_cutoff_var,
        fill = "ident",
        colour = "ident",
        group = "ident"
      )
    )

    g <- g + ggplot2::geom_point(
      data = data_subset[data_subset$in_next,],
      alpha = 0.3,
      size = 4)
    g <- g + ggplot2::geom_point(
      data = data_subset[data_subset$in_next,],
      alpha = 1)
    g <- g + ggplot2::geom_point(
      data = data_subset[!data_subset$in_next,],
      alpha = 0.2)
    g <- g + ggplot2::theme_bw()
    return(g)
  }

  foo <- partykit:::.list.rules.party(tree)
  variable_names <- names(partykit::varimp(tree))

  split_rules <-
    unlist(strsplit(foo[[as.character(terminal_node)]], " & "))

  intermediates <- Reduce(function(x, y) {
    paste(x, y, sep = " & ")
  },
  split_rules, accumulate = TRUE)

  intermediates <-
    gsub("(?<=\\.\\d{2})(\\d+)", "", intermediates, perl = TRUE)

  last_var <- gsub(".*&\\s*", "", intermediates, perl = TRUE)

  # strips the conditional, for instance "FOO >= -0.001" becomes "FOO"
  last_var_names <- gsub("(.*)\\s([><=]+)\\s+(-?[0-9.]+)\\s*$",
                         "\\1",
                         last_var,
                         perl = TRUE)

  whole_data <- as.data.frame(object, genes = variable_names)

  intermediate_data <-
    lapply(intermediates, function(x, whole_data) {
      whole_data[with(whole_data, eval(parse(text = x))),]
    }, whole_data)

  g <- mapply(
    plot_gate_subset,
    data_subset = c(list(whole_data), intermediate_data[-length(intermediate_data)]),
    next_expr = c(intermediates),
    cutoff_var = c(last_var_names),
    next_cutoff_var = c(last_var_names[-1], last_var_names[length(last_var_names) - 1]),
    SIMPLIFY = FALSE
  )

  unique_idents <- unique(as.character(whole_data$ident))
  colourCount <- length(unique_idents)

  # This creates a function that will generate a character vector of colours
  # given a number of colour to generate
  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))

  Palette <- getPalette(colourCount)

  # This generates a named vector with values equal to the colours and names corresponding
  # To the identities that will have that colour
  names(Palette) <- unique_idents


  g <- lapply(g, function(g) {
    g + ggplot2::scale_colour_manual(values = Palette)} )

  g <- cowplot::plot_grid(
    plotlist = g,
    labels = paste0(seq_along(intermediates), ". ", intermediates),
    label_size = 10,
    align = "hv",
    scale = 0.9,
    label_x = 0.1,
    hjust = 0
  )
  return(g)
}


