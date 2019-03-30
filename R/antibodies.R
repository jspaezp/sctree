
# TODO: Fix biolegend function to return a data frame consistently with the rest
# of the other funtions
# TODO: implement caching the queries so they can quickly be recovered from the
# commonly searched antibodies
# TODO: write a "generic" function that unifies all companies.
# TODO: implement an automated way to search as well for the aliases of the gene
# names provided

#' Query cell signaling for Single cell antibodies
#'
#' Queries Cell Signaling Technologies and
#' filters for antibodies for IFC, IHC, and FLOW.
#'
#' Note that it is jsut a wrapper arround the web search of the site, so please
#' be nice on them and dont get your IP banned ...
#'
#' @param search_term The term used to query the cell signaling webpage
#'
#' @return NULL when no results are found or a DF with information
#' @export
#'
#' @examples
#' query_cc_antibodies("CD11bfakename")
#' # NULL
#' query_cc_antibodies("CD11c")
#' #   Cat_num                                       Product Name Reactivity
#' # 1   97585                           CD11c (D1V9Y) Rabbit mAb          M
#' # 2   97473 CD11c (3.9) Mouse mAb (violetFluor™ 450 Conjugate)          H
#' # 3   42847    CD11c (3.9) Mouse mAb (redFluor™ 710 Conjugate)          H
#' # 4   80342          CD11c (3.9) Mouse mAb (PE-Cy7® Conjugate)          H
#' # 5   56025               CD11c (3.9) Mouse mAb (PE Conjugate)          H
#' # 6   69627             CD11c (3.9) Mouse mAb (FITC Conjugate)          H
#' # 7   51562         CD11c (3.9) Mouse mAb (APC-Cy7® Conjugate)          H
#' # 8   36268              CD11c (3.9) Mouse mAb (APC Conjugate)          H
#' # Product_Type
#' # 1  Primary Antibodies
#' # 2 Antibody Conjugates
#' # 3 Antibody Conjugates
#' # 4 Antibody Conjugates
#' # 5 Antibody Conjugates
#' # 6 Antibody Conjugates
#' # 7 Antibody Conjugates
#' # 8 Antibody Conjugates
#' @importFrom rvest html_nodes html_table
#' @importFrom xml2 read_html
query_cc_antibodies <- function(search_term) {
    url <- paste0( "https://www.cellsignal.com/browse/?Ntt=",
                   search_term,
                   "&N=4294960091+4294960086+4294960087+4294964832+4294956287",
                   # ^^^ This filters for antibodies for IFC, IHC, and FLOW
                   "&No={offset}&Nrpp=1000") # and this one forces to show the top 1000

    antibody_df <- rvest::html_nodes(
        xml2::read_html(url),
        xpath='//*[@id="product-list"]')

    antibody_df <- rvest::html_table(antibody_df, fill = TRUE)
    # The former will be a list of length 0 if no results are found

    if (length(antibody_df) == 0) {return(NULL)}


    antibody_df <- antibody_df[[1]]
    colnames(antibody_df)[ncol(antibody_df)] <- "Product_Type"
    colnames(antibody_df)[1] <- "Cat_num"
    antibody_df <- antibody_df[grepl("Antibod", antibody_df$Product_Type),]
    antibody_df <- antibody_df[antibody_df[["Reactivity"]] != "",]

    # It does not contain the information of the usage so I will just delete it...
    # map(antibody_df[["Application"]], function(x) str_split(gsub("(\\t|\\n)+", ',', x), pattern = ","))
    antibody_df <- antibody_df[, ! colnames(antibody_df) %in% "Application" ]

    return(antibody_df)

}



#' Query Santa Cruz for antibodies
#'
#' Queries Santa Cruz for antibodies and
#' filters for antibodies for IFC, IHC, and FCM
#'
#' Note that it is jsut a wrapper arround the web search of the site, so please
#' be nice on them and dont get your IP banned ...
#'
#' @param search_term The term used to query the webpage
#'
#' @return NULL when no results are found or a DF with information
#' @export
#'
#' @examples
#' query_sc_antibodies("CD11bfakename")
#' NULL
#' head(query_sc_antibodies("CD11C"))
#' # Product_Name  Cat_num     Epitope    Species
#' # 1 Integrin alpha X Antibody (2Q862) sc-71454 Integrin aX  and human
#' # 2 Integrin alpha X Antibody (2Q865) sc-71455 Integrin aX      human
#' # 3   Integrin alpha X Antibody (3.9)  sc-1185 Integrin aX      human
#' # 4 Integrin alpha X Antibody (3H986) sc-71456 Integrin aX      mouse
#' # 5   Integrin alpha X Antibody (B-6) sc-46676 Integrin aX      human
#' # 6 Integrin alpha X Antibody (B-Iy6) sc-19989 Integrin aX      human
#' # Method
#' # 1                   IF and FCM
#' # 2                   IP and FCM
#' # 3               IP, IF and FCM
#' # 4                   IF and FCM
#' # 5 WB, IP, IF, IHC(P) and ELISA
#' # 6                   IF and FCM
#' @importFrom rvest html_nodes html_table
#' @importFrom xml2 read_html
query_sc_antibodies <- function(search_term) {
    url <- paste0("https://www.scbt.com/scbt/search?Ntt=",
                  search_term,
                  "&N=1354381666&Nrpp=50")

    antibody_df <- rvest::html_nodes(
        xml2::read_html(url),
        xpath='/html/body/div[1]/div/div/div/div[6]/table')

    antibody_df <- rvest::html_table(antibody_df, fill = TRUE)
    # The former will be a list of length 0 if no results are found

    if (length(antibody_df) == 0) {return(NULL)}


    antibody_df <- antibody_df[[1]]
    antibody_df <- antibody_df[,2:(ncol(antibody_df) - 1)]
    antibody_df <- antibody_df[grepl("detection", antibody_df$Description),]
    colnames(antibody_df)[2] <- "Cat_num"
    colnames(antibody_df)[1] <- "Product_Name"

    antibody_df <- antibody_df[grepl("IF|IHC|FCM", antibody_df$Description),]
    antibody_df <- antibody_df[!grepl("DISCONTINUED", antibody_df$Description),]


    descriptions <- strsplit(
        gsub(
            "recommended for detection of (.*) of (.*) origin by (.*)",
            "\\1::\\2::\\3",
            antibody_df$Description),
        split = "::")

    descriptions <- as.data.frame(do.call(rbind, descriptions))

    colnames(descriptions) <- c("Epitope", "Species", "Method")
    antibody_df <- cbind(
        antibody_df[, !colnames(antibody_df) %in%  "Description"],
        descriptions)

    return(antibody_df)
}


# ABCAM
query_ab_antibodies <- function() {
    stop("NOT IMPLEMENTED")


    # url <- paste0("https://www.scbt.com/scbt/search?Ntt=",
    #               search_term,
    #               "&N=1354381666&Nrpp=50")
    # url <- "https://www.abcam.com/products?sortOptions=Relevance&keywords=KLK3&selected.classification=Primary+antibodies"

    # antibody_df <- rvest::html_nodes(
    #     xml2::read_html(url),
    #     xpath='/html/body/div[1]/div/div/div/div[6]/table')

    # antibody_df <- rvest::html_table(antibody_df, fill = TRUE)
    # # The former will be a list of length 0 if no results are found

    # if (length(antibody_df) == 0) {return(NULL)}


    # antibody_df <- antibody_df[[1]]
    # antibody_df <- antibody_df[,2:(ncol(antibody_df)-1)]
    # antibody_df <- antibody_df[grepl("detection", antibody_df$Description),]
    # colnames(antibody_df)[2] <- "Cat_num"
    # colnames(antibody_df)[1] <- "Product_Name"

    # antibody_df <- antibody_df[grepl("IF|IHC|FCM", antibody_df$Description),]
    # antibody_df <- antibody_df[!grepl("DISCONTINUED", antibody_df$Description),]


    # descriptions <- strsplit(
    #     gsub(
    #         "recommended for detection of (.*) of (.*) origin by (.*)",
    #         "\\1::\\2::\\3",
    #         antibody_df$Description),
    #     split = "::")

    # descriptions <- as.data.frame(do.call(rbind, descriptions))

    # colnames(descriptions) <- c("Epitope", "Species", "Method")
    # antibody_df <- cbind(
    #     antibody_df[, !colnames(antibody_df) %in%  "Description"],
    #     descriptions)

    # return(antibody_df)
}


#' Query Biolegend for antibodies
#'
#' Queries Biolegend for antibodies and
#' filters for antibodies for Flow Cytometry
#'
#' Note that it is just a wrapper arround the web search of the site, so please
#' be nice on them and dont get your IP banned ...
#'
#' @param search_term The term used to query the webpage
#'
#' @return a length 0 character when no results are found or a character vector with the product names
#' @export
#'
#' @examples
#' query_biolegend_antibodies("CD11bfakename")
#' # NULL
#' head(query_biolegend_antibodies("CD11C"))
#' # [1] "MojoSortâ\u0084¢ Mouse CD11c Nanobeads"
#' # [2] "APC anti-human CD11c Antibody"
#' # [3] "Biotin anti-human CD11c Antibody"
#' # [4] "FITC anti-human CD11c Antibody"
#' # [5] "PE anti-human CD11c Antibody"
#' # [6] "PE/Cy5 anti-human CD11c Antibody"
#' @importFrom rvest html_nodes html_table
#' @importFrom xml2 read_html
query_biolegend_antibodies <- function(search_term) {
    url <- paste0("https://www.biolegend.com/en-us/search-results?Applications=FC&Keywords=",
                  search_term)

    antibody_list <- rvest::html_nodes(
        xml2::read_html(url),
        xpath='/html/body/div[1]/section/div/div/div/div/section/div/article/div[3]/ul[2]/li[*]/ul/li/h2/a')

    antibody_list <- rvest::html_text(antibody_list, trim = TRUE)
    # The former will be a list of length 0 if no results are found

    return(antibody_list)
}

