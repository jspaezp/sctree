
#' featureplot_gadget
#'
#' Shiny gadget to explore clusters and features for a seurat object.
#'
#' @param seurat_object a seurat object to explore
#' @param starting_genes genes to start with
#'
#' @return TODO return something ... its a gadget after all ...
#' @export
#'
#'
#' @importFrom viridis viridis
#' @importFrom memoise cache_filesystem
#' @importFrom DT renderDataTable
featureplot_gadget <- function(seurat_object = Seurat::pbmc_small,
                               starting_genes = NULL,
                               cache = "./.cache") {

    plotting_pannel <- shiny::sidebarPanel(
        # TODO remove reactive button function to add genes or fix reactivity at startup
        shiny::h2("Cluster of interest"),
        shiny::uiOutput("clusters"),
        shiny::uiOutput("allgenes"),
        shiny::actionButton("add_gene", "Add to plottable genes"),
        shiny::uiOutput("genes"),
        shiny::actionButton("plot_command", "PLOT!"),
        shiny::checkboxInput(
            "hidelegends",
            "Hide Legends",
            value = FALSE),
        shiny::sliderInput(
            "pointsize",
            "Point Size:",
            min = 0.5,
            max = 10,
            value = 4
        )
    )


    ui <- shiny::fluidPage(
        shiny::titlePanel("Feature Plot for a Seurat Object"),

        shiny::sidebarLayout(
            plotting_pannel,

            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel("TSNE plot", shiny::plotOutput("generalTsne", height = '600px')),
                    shiny::tabPanel(
                        "Find Markers",
                        shiny::div(DT::dataTableOutput("rangerImpTable")),
                        shiny::plotOutput("rangerImpPlot")),
                    shiny::tabPanel("Feature plot", shiny::plotOutput("tsnePlot", height = '600px')),
                    shiny::tabPanel("Violin Plot", shiny::plotOutput("violinPlot", height = '600px')),
                    shiny::tabPanel(
                        "Pairs Plot",
                        shiny::plotOutput("generalTsne2", height = '600px', width = '600px'),
                        shiny::plotOutput("pairsPlot", height = '800px')))
            )

        )
    )

    server <- function(input, output) {
        object_df <- as.data.frame.Seurat(seurat_object, seurat_object@var.genes)
        filesys_cache <- memoise::cache_filesystem(cache)
        mem_ranger.df <- memoise::memoise(ranger_importances.df, cache = filesys_cache)


        colorscale <- viridis::viridis(2, direction = -1)
        # TODO add an option to select either all genes or only variables
        #total_genes <- rownames(seurat_object@raw.data)
        total_genes <- seurat_object@var.genes



        importance_df <- shiny::eventReactive(input$cluster, {
            importance_df <- mem_ranger.df(object = object_df,
                                           cluster = input$cluster,
                                           num.trees = 2000)[[3]]
            return(importance_df)
        })

        if (is.null(starting_genes)) {
            #plottable_genes <- rownames(head(importance_df))
            plottable_genes <- character()
        } else {
            plottable_genes <- starting_genes
        }

        get_plottable_genes <- shiny::eventReactive(input$add_gene, {

            plottable_genes <<- sort(unique(c(
                        plottable_genes, input$allgenes)))
            return(plottable_genes)
        })


        output$clusters <- shiny::renderUI({
            shiny::selectInput(
                inputId = "cluster",
                label = "Cluster to find markers for",
                choices = sort(unique(seurat_object@ident)))
        })


        output$generalTsne <- shiny::renderPlot({
            tsne_plot(seurat_object)
        })

        output$generalTsne2 <- shiny::renderPlot({
            tsne_plot(seurat_object)
        })

        output$rangerImpTable <- DT::renderDataTable({
            importance_df()
        })

        output$rangerImpPlot <- shiny::renderCachedPlot({
            importance_df <- importance_df()

            # TODO evaluate if reimplementing fct_reorder is worth one less dependency
            importance_df[["gene"]] <-  forcats::fct_reorder(
                importance_df[["gene"]],
                importance_df[["importance"]])

            g <- ggplot2::ggplot(
                dplyr::top_n(importance_df, 20, importance),
                ggplot2::aes_string(x = "importance", xend = 0,
                           y = "gene", yend = "gene")) +
                ggplot2::geom_point() +
                ggplot2::geom_segment() +
                ggplot2::ggtitle("Relative variable importance")
            g
        }, importance_df())

        output$genes <- shiny::renderUI({
            shiny::checkboxGroupInput(
                inputId = "genenames",
                label = "Genes to plot",
                choices = get_plottable_genes(),
                selected = get_plottable_genes())
        })

        output$allgenes <- shiny::renderUI({
            choices <- character()

            tryCatch({
                choices <- sort(total_genes[!total_genes %in%
                                                get_plottable_genes()])
            }, error = function(e) {
                warning("Recovering from error in `output$allgenes <- shiny::renderUI`")
                choices <<- sort(seurat_object@var.genes)
            })

            shiny::selectInput(
                'allgenes',
                "allgenes",
                choices = choices)
        })

        captured_input <- shiny::eventReactive(input$plot_command, {
            retlist <- list(
               genenames = input$genenames,
               pointsize = input$pointsize,
               hidelegends = input$hidelegends)
            return(retlist)
        })

       output$tsnePlot <- shiny::renderCachedPlot({
           capt_input <- captured_input()

           Seurat::FeaturePlot(
                object = seurat_object,
                features.plot = capt_input$genenames,
                cols.use = colorscale,
                reduction.use = "tsne",
                dim.1 = 1,
                no.legend = capt_input$hidelegends,
                do.return = FALSE,
                pt.size = capt_input$pointsize)
       }, cacheKeyExpr = {
           captured_input()
       })

       output$violinPlot <- shiny::renderCachedPlot({
           capt_input <- captured_input()

           Seurat::VlnPlot(
               object = seurat_object,
               features.plot = capt_input$genenames,
               do.return = FALSE)
       }, cacheKeyExpr = {
           captured_input()
       })

       output$pairsPlot <- shiny::renderCachedPlot({
           capt_input <- captured_input()

           plot.flowstyle(
               object = seurat_object,
               markernames = make.names(capt_input$genenames))
       }, cacheKeyExpr = {
           captured_input()
       })

    }

    # Run the application
    shiny::shinyApp(ui = ui, server = server)
}

##################### ================= NOTE ===================================
#####################  It seems that the slowest part of the process is to
#####################  serialize the seurat opbject, could try to run the
#####################  extraction of the data frame to cache the process ...
#####################
##################### ================= END NOTE ===============================
##################### library(shiny)
#library(Seurat)
#library(DT)
#library(memoise)

#filesys_cache <- memoise::cache_filesystem("./.cache")

#mem_ranger.seurat <- memoise::memoise(ranger_importances.seurat, cache = filesys_cache)

#system.time({mem_ranger.seurat(dataset_5050, cluster = "0", pval_cutoff = 0.05)})

#profvis({mem_ranger.seurat(dataset_5050, cluster = "1", pval_cutoff = 0.05)})

#profvis({shinyApp(ui = ui, server = server, options = "test.mode")} ,prof_input = '/path_to_save_output/random_name.Rprof')

#for (i in c("md5", "sha1", "crc32", "sha256", "sha512",
#            "xxhash32", "xxhash64", "murmur32")) {
#    print(i)
#    print(system.time({digest::digest(dataset_5050, i)}))
#}
#print(system.time({digest::digest(dataset_5050)}))
#require(profvis)
#profvis({
#    featureplot_gadget(dataset_5050,
#                       starting_genes = c("CD3D", "TMSB4X", "ARHGDIB", "CKB", "XIST", "CDKN2A"),
#                       cache = "./.cache")
#})
#for (i in c("md5", "sha1", "crc32", "sha256", "sha512",
#            "xxhash32", "xxhash64", "murmur32")) {
#    print(i)
#    print(system.time({digest::digest(object_df, i)}))
#}
#filesys_cache <- memoise::cache_filesystem("./cache")
#mem_ranger.df <- memoise::memoise(ranger_importances.df, cache = filesys_cache)
#ranger_importances.df(object_df, cluster = "0")
#system.time({mem_ranger.df(object_df, cluster = "1")})

#featureplot_gadget()


# dataset_5050 <- readRDS("C:/Users/Sebastian/Downloads/seurat_mix5050.RDS")


# object_df <- as.data.frame.Seurat(dataset_5050)
# filesys_cache <- memoise::cache_filesystem("./.cache")
# mem_ranger.df <- memoise::memoise(ranger_importances.df, cache = filesys_cache)
#
#
#
#
# clus_1 <- mem_ranger.df(object_df, cluster = "1")
#
# clus_0 <- mem_ranger.df(object_df, cluster = "0")
#
# clus_0_antibodies <- head(rownames(clus_0$signif_importances_ranger), 30) %>%
#     purrr::map(query_biolegend_antibodies)
#
# clus_1_antibodies <- head(rownames(clus_1$signif_importances_ranger), 30) %>%
#     purrr::map(query_biolegend_antibodies)
#
#
# featureplot_gadget(seurat_object = dataset_5050,
#                    starting_genes = c("CD3D", "TMSB4X", "ARHGDIB", "CKB", "XIST", "CDKN2A", "EIF5A"),
#                    cache = "./.cache")
#
#
#
# surface_annotated <- AnnotationDbi::select(
#     org.Hs.eg.db::org.Hs.eg.db,
#     keys = dataset_5050@var.genes,
#     columns = "GO",
#     keytype = "ALIAS") %>%
#     dplyr::filter(GO %in% c("GO:0009986", "GO:0005886", "GO:0044459"))
#
# head(surface_annotated)
#
#
# surface_annotated
#
# clus_1_surface <- mem_ranger.df(object_df[, c(surface_annotated$ALIAS, 'ident')], cluster = "1")
#
#
#
# #source("./R/antibodies.R")
#
#
# #query_biolegend_antibodies("CD3")
