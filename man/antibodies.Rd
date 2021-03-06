% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antibodies.R
\name{antibodies}
\alias{antibodies}
\alias{query_cc_antibodies}
\alias{query_sc_antibodies}
\alias{query_biocompare_antibodies}
\alias{query_biolegend_antibodies}
\title{Query Antibodies from multiple vendors}
\usage{
query_cc_antibodies(search_term, sleep = 1)

query_sc_antibodies(search_term, sleep = 1)

query_biocompare_antibodies(search_term, sleep = 1)

query_biolegend_antibodies(search_term, sleep = 1)
}
\arguments{
\item{search_term}{The term used to query the antibodies}

\item{sleep}{an ammount of time to wait before submitting the query.
Check details for more information}
}
\value{
NULL when no results are found or a DF with information
}
\description{
Queries the respective antibody for each vendor and when posible
filters for antibodies for IFC, IHC, and FLOW.
}
\details{
Most of these functions have been implemented as queries to the
vendor website, therefore if many queries are done in quick
succession, it is posible that the vendor considers the queries as
an attack to their website. Therefore we have the sleep argument set
to 1 by default, the function will wait 1 second before doing the
query. If many queries will be done, feel free to modify the
parameter but be aware that some requests might be denied.
}
\examples{
query_cc_antibodies("CD11bfakename")
# NULL
query_cc_antibodies("CD11c")
#   Cat_num                                       Product Name Reactivity
# 1   97585                           CD11c (D1V9Y) Rabbit mAb          M
# 2   97473 CD11c (3.9) Mouse mAb (violetFluor™ 450 Conjugate)          H
# 3   42847    CD11c (3.9) Mouse mAb (redFluor™ 710 Conjugate)          H
# 4   80342          CD11c (3.9) Mouse mAb (PE-Cy7® Conjugate)          H
# 5   56025               CD11c (3.9) Mouse mAb (PE Conjugate)          H
# 6   69627             CD11c (3.9) Mouse mAb (FITC Conjugate)          H
# 7   51562         CD11c (3.9) Mouse mAb (APC-Cy7® Conjugate)          H
# 8   36268              CD11c (3.9) Mouse mAb (APC Conjugate)          H
# Product_Type
# 1  Primary Antibodies
# 2 Antibody Conjugates
# 3 Antibody Conjugates
# 4 Antibody Conjugates
# 5 Antibody Conjugates
# 6 Antibody Conjugates
# 7 Antibody Conjugates
# 8 Antibody Conjugates

query_sc_antibodies("CD11bfakename")
# NULL
head(query_sc_antibodies("CD11C"))
# Product_Name  Cat_num Citations Rating     Epitope
# 1   Integrin alpha X Antibody (B-6) sc-46676         3    (8) Integrin aX
# 2   Integrin alpha X Antibody (3.9)  sc-1185         2    (2) Integrin aX
# 3 Integrin alpha X Antibody (B-Iy6) sc-19989         1    (1) Integrin aX
# 4  Integrin alpha X Antibody (N418) sc-23951         3    (2) Integrin aX
# 6 Integrin alpha X Antibody (2Q865) sc-71455       NEW    (1) Integrin aX
# 7 Integrin alpha X Antibody (3H986) sc-71456       NEW    (1) Integrin aX
# Species                       Method
# 1   human WB, IP, IF, IHC(P) and ELISA
# 2   human               IP, IF and FCM
# 3   human                   IF and FCM
# 4   mouse                   IF and FCM
# 6   human                   IP and FCM
# 7   mouse                   IF and FCM

#' query_biocompare_antibodies("CD11bfakename")
# NULL
head(query_biocompare_antibodies("CD11C"),3)


query_biolegend_antibodies("CD11bfakename")
# NULL
head(query_biolegend_antibodies("CD11C"))
# [1] "MojoSortâ\u0084¢ Mouse CD11c Nanobeads"
# [2] "APC anti-human CD11c Antibody"
# [3] "Biotin anti-human CD11c Antibody"
# [4] "FITC anti-human CD11c Antibody"
# [5] "PE anti-human CD11c Antibody"
# [6] "PE/Cy5 anti-human CD11c Antibody"


}
