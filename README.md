sctree: a package to connect single cell rna-seq to biology using trees
================

  - [sctree](#sctree)
  - [Installation](#installation)
  - [Usage](#usage)
      - [Finding important variables to classify
        clusters](#finding-important-variables-to-classify-clusters)
      - [Visualizing the expected outcome of a flow cytometry
        experiment](#visualizing-the-expected-outcome-of-a-flow-cytometry-experiment)
      - [Suggesting a gating strategy for the
        markers](#suggesting-a-gating-strategy-for-the-markers)
      - [Finding equivalent clusters in two
        datasets](#finding-equivalent-clusters-in-two-datasets)
      - [Finding antibodies for the
        experiment](#finding-antibodies-for-the-experiment)
  - [Reproducing the runs in the purdue
    cluster](#reproducing-the-runs-in-the-purdue-cluster)
  - [Steps down the road](#steps-down-the-road)

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis build
status](https://travis-ci.org/jspaezp/sctree.svg?branch=master)](https://travis-ci.org/jspaezp/sctree)
[![Coverage
status](https://codecov.io/gh/jspaezp/sctree/branch/master/graph/badge.svg)](https://codecov.io/github/jspaezp/sctree?branch=master)

# sctree

The goal of sctree is to create a tool to accelerate the transition from
single cell rna-sequencing to calidation and new sub-population
discovery.

Features suggesting pseudo-gating strategies to purify found populations
via flow-cytometry, antibody querying and cross validations between
datasets.

Number of lines in roxygen comments: 652

Number of lines in R code: 1034

# Installation

    git clone https://github.rcac.purdue.edu/jpaezpae/sctree sctree
    cd sctree
    
    R -e "devtools::install('.')"

# Usage

I am assuming you have already done your clustering and dimensional
reduction using seurat and we have our seurat object.

For this examples we will use a dummy dataset that come bundled with the
package

`small_5050_mix`, this dataset comes originally from the 1:1 mixture of
Jurkat and 293T cells provided by 10x.

Original data can be found here:

1.  [1:1
    mixture](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/jurkat:293t_50:50)
2.  [99:1
    mixture](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/jurkat_293t_99_1)

<!-- end list -->

``` r
require(sctree)
#> Loading required package: sctree
#> Registered S3 methods overwritten by 'ggplot2':
#>   method         from 
#>   [.quosures     rlang
#>   c.quosures     rlang
#>   print.quosures rlang
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
#> Registered S3 method overwritten by 'R.oo':
#>   method        from       
#>   throw.default R.methodsS3
#> Registered S3 method overwritten by 'rvest':
#>   method            from
#>   read_xml.response xml2
require(Seurat)
#> Loading required package: Seurat

set.seed(6)

data(small_5050_mix)
small_5050_mix
#> An object of class Seurat 
#> 1031 features across 255 samples within 1 assay 
#> Active assay: RNA (1031 features)
#>  2 dimensional reductions calculated: pca, tsne

DimPlot(small_5050_mix, reduction = "tsne")
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Finding important variables to classify clusters

We base our importances on the “classification value” they give to a
random forest (using the implementation in the `ranger` package)

So lets fit the random forest … Here we are adding the `warn.imp.method`
to prevent a warning message sent by `ranger` when most of the variables
are correlated with the clustering.

Please reffer to the of the `importance_pvalues` section in the [ranger
documentation](https://cran.r-project.org/web/packages/ranger/ranger.pdf)
when addressing this issue and for more details.

``` r
rang_importances <- ranger_importances.Seurat(
    small_5050_mix,
    cluster = "ALL",
    warn.imp.method = FALSE)

names(rang_importances)
#> [1] "ranger_fit"                "importances_ranger"       
#> [3] "signif_importances_ranger"
```

This gives us a list with 3 elements.

1.  The ranger fit object itself (handy if you want to inspect its
    classification accuracy)
2.  The importance matrix deriven from ranger
3.  A data frame containing only importances with pvalues under 0.05
    (because biologists love p-values under 0.05)

<!-- end list -->

``` r
rang_importances[[1]]
#> Ranger result
#> 
#> Call:
#>  ranger::ranger(dependent.variable.name = "ident", data = data,      num.trees = num.trees, mtry = floor(ncol(data)/5), importance = importance,      classification = TRUE, ...) 
#> 
#> Type:                             Classification 
#> Number of trees:                  500 
#> Sample size:                      255 
#> Number of independent variables:  140 
#> Mtry:                             28 
#> Target node size:                 1 
#> Variable importance mode:         impurity_corrected 
#> Splitrule:                        gini 
#> OOB prediction error:             18.04 %
```

We can see that our classifier is 18.04 % accurate, as measured by its
*OOB prediction error*.

Your classifier USUALLY will be just as good as your initial clustering
is, therefore if your clustering is not all that consistent to start
with, it is unrealistic that predicions will be done accurately.

Now … taking a look at the importance matrix we can see that a `pvalue`
and a relative importance has been given to each gene.

``` r
head(rang_importances[[2]])
#>                importance     pvalue
#> TNFRSF4       -0.20663582 0.94318182
#> RP3-395M20.12  0.03391993 0.36363636
#> ID3            0.09088435 0.18181818
#> JUN            0.19175235 0.09090909
#> DEPDC1         0.15176841 0.09090909
#> CHI3L2         0.52477223 0.00000000
```

For simplicity we return by default a `data.frame` containing the ones
with `pvalue < 0.05`

``` r
head(rang_importances[[3]])
#>         importance pvalue    gene
#> ASNS      6.350801      0    ASNS
#> CD3D      4.477459      0    CD3D
#> ARHGDIB   4.367260      0 ARHGDIB
#> TMSB4X    4.344165      0  TMSB4X
#> ADA       3.731668      0     ADA
#> MZB1      3.222495      0    MZB1
```

Therefore, in this case we can say that the expression of the following
genes would be usefull to form the clusters.

As an analogous function to Seurat’s `FindAllMarkers`, we offer
`FindAllMarkers_ranger.Seurat`

``` r
markers <- FindAllMarkers_ranger.Seurat(
  small_5050_mix,
  warn.imp.method = FALSE)

head(markers)
#>         importance pvalue    gene cluster
#> ASNS      5.679208      0    ASNS       0
#> TMSB4X    5.281722      0  TMSB4X       0
#> ARHGDIB   5.002977      0 ARHGDIB       0
#> ADA       3.722655      0     ADA       0
#> CD3D      3.471757      0    CD3D       0
#> MZB1      3.033289      0    MZB1       0
```

## Visualizing the expected outcome of a flow cytometry experiment

Lets say we choose the top 5 markers from the former list and we did a
flow experiment … HYPOTHETICALLY the marker distribution would resemble
the rna expression profile for which we have the function
`plot_flowstyle`

``` r
top_markers <- head(rang_importances[[3]]$gene)
top_markers
#> [1] "ASNS"    "CD3D"    "ARHGDIB" "TMSB4X"  "ADA"     "MZB1"
g <- plot_flowstyle(small_5050_mix, markernames = top_markers)
g
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Based on this, we can see that the red cluster in this plot is
predominantly CD3+ ADA+ and ARHGDIB+, as well as ASNS-

We can also focus in one of the pannels (and check the color
conventions)

``` r
g[1,2]
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Suggesting a gating strategy for the markers

A general strategy to get separate all clusters

``` r

top_markers <- head(rang_importances[[3]]$gene)

tree_fit <- fit_ctree(small_5050_mix,
                      genes_use = top_markers, 
                      cluster = "ALL")
```

Visualizing the tree as … a tree … we can see how our model is a simple
series of yes/no questions.

If we wanted to classifiy a random cell: in the first `node`, we check
if the expression of that gene is higher or lower than a given value, if
it is lower, we proceed to the left, if not we go right. We keep doing
that until we have no more `branches`. This final node will have a
predicted cluster, in this plot we can also see how pure can we expect
this group to be and how many of the cells in our training set clasify
as part of it.

``` r
plot(tree_fit)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

When inspecting the tree\_fit, we can see a more detailed text
representation of this tree.

``` r
print(tree_fit)
#> 
#> Model formula:
#> ident ~ ASNS + CD3D + ARHGDIB + TMSB4X + ADA + MZB1
#> 
#> Fitted party:
#> [1] root
#> |   [2] ADA <= 2.83724
#> |   |   [3] ASNS <= 1.43241: 0 (n = 40, err = 35.0%)
#> |   |   [4] ASNS > 1.43241
#> |   |   |   [5] ARHGDIB <= 1.90619: 1 (n = 55, err = 10.9%)
#> |   |   |   [6] ARHGDIB > 1.90619: 0 (n = 7, err = 28.6%)
#> |   [7] ADA > 2.83724
#> |   |   [8] ASNS <= 2.10406: 0 (n = 146, err = 3.4%)
#> |   |   [9] ASNS > 2.10406: 1 (n = 7, err = 28.6%)
#> 
#> Number of inner nodes:    4
#> Number of terminal nodes: 5
```

Sometimes one might think that the proposed strategy is too complicated
or not implementable in the experimental settings, in order to add
constrians to the fit one can give additional arguments that will be
passed to `partykit::ctree_control`, such as `maxdepth = 2` (maximum 2
questions per cell)

``` r
tree_fit <- fit_ctree(
  small_5050_mix, genes_use = top_markers, 
  cluster = "ALL", maxdepth = 2)
print(tree_fit)
#> 
#> Model formula:
#> ident ~ ASNS + CD3D + ARHGDIB + TMSB4X + ADA + MZB1
#> 
#> Fitted party:
#> [1] root
#> |   [2] ADA <= 2.83724
#> |   |   [3] ASNS <= 1.43241: 0 (n = 40, err = 35.0%)
#> |   |   [4] ASNS > 1.43241: 1 (n = 62, err = 17.7%)
#> |   [5] ADA > 2.83724
#> |   |   [6] ASNS <= 2.10406: 0 (n = 146, err = 3.4%)
#> |   |   [7] ASNS > 2.10406: 1 (n = 7, err = 28.6%)
#> 
#> Number of inner nodes:    3
#> Number of terminal nodes: 4
plot(tree_fit)
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Since not all variables are ultimately used in our classifier, one can
acces the ones that were by using `varimp(tree_fit)`

``` r
partykit::varimp(tree_fit)
#>       ADA      ASNS 
#> 0.2231765 0.2185317
plot_flowstyle(small_5050_mix, names(partykit::varimp(tree_fit)))
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

One can also request the package to suggest a specific strategy only for
a given cluster. This function is not expected to give drastically
different results in datasets with few clusters, but it can definitely
come usefull when many clusters are present and one is interested in a
specific one.

``` r
tree_fit <- fit_ctree(small_5050_mix, genes_use = top_markers, cluster = "0")
print(tree_fit)
#> 
#> Model formula:
#> ident ~ ASNS + CD3D + ARHGDIB + TMSB4X + ADA + MZB1
#> 
#> Fitted party:
#> [1] root
#> |   [2] ADA <= 2.83724
#> |   |   [3] ASNS <= 1.43241: clus 0 (n = 40, err = 35.0%)
#> |   |   [4] ASNS > 1.43241
#> |   |   |   [5] ARHGDIB <= 1.90619: not clus 0 (n = 55, err = 10.9%)
#> |   |   |   [6] ARHGDIB > 1.90619: clus 0 (n = 7, err = 28.6%)
#> |   [7] ADA > 2.83724
#> |   |   [8] ASNS <= 2.10406: clus 0 (n = 146, err = 3.4%)
#> |   |   [9] ASNS > 2.10406: not clus 0 (n = 7, err = 28.6%)
#> 
#> Number of inner nodes:    4
#> Number of terminal nodes: 5
```

## Finding equivalent clusters in two datasets

``` r
data(small_9901_mix)
small_9901_mix
#> An object of class Seurat 
#> 840 features across 384 samples within 1 assay 
#> Active assay: RNA (840 features)
#>  2 dimensional reductions calculated: pca, tsne
```

``` r
validation_results <- cross_validate(
    small_5050_mix, small_9901_mix, 
    cluster = "ALL",
    warn.imp.method = FALSE)
#> Warning in cross_validate(small_5050_mix, small_9901_mix, cluster = "ALL", : Some important genes were removed because they are not present in the test dataset. 
#> Removed genes: ASNS, CD3D, ADA, HEY1, RPL26, CA2, XIST, CDKN2A, TSC22D3, AIF1, GAL, PSMB9, FAM127B, MAP1A, PSMB8, TSTD1, MDK, CSRP2, DMKN, HOXA9, RNF138, ZNF503, CTC1, C21orf90, HSPA1B, PYGL, SPRED2, ZFAT

validation_results[[1]]
#> 
#> Model formula:
#> ident ~ ARHGDIB + TMSB4X + MZB1 + SOX4 + FYB + UBE2C + HIST1H1E + 
#>     CD1E + ITGA4 + ITM2A + HIST1H4C + CXCR4 + CDK1 + C12orf57 + 
#>     HIST1H2BK + ARPP21 + HIST1H1C + CCNB1 + IGLL1 + CHI3L2 + 
#>     JUN
#> 
#> Fitted party:
#> [1] root
#> |   [2] ARHGDIB <= 2.29381
#> |   |   [3] SOX4 <= 2.74241: 1 (n = 90, err = 27.8%)
#> |   |   [4] SOX4 > 2.74241: 0 (n = 8, err = 0.0%)
#> |   [5] ARHGDIB > 2.29381: 0 (n = 157, err = 6.4%)
#> 
#> Number of inner nodes:    2
#> Number of terminal nodes: 3
```

``` r
validation_results$confusion_matrix
#>          cluster
#> predicted   0   1   2
#>         0 263  98   8
#>         1   0   0  15
```

``` r
freq_matrix <- as.frequency.matrix(validation_results$confusion_matrix)
freq_matrix
#>          cluster
#> predicted         0         1         2
#>         0 100.00000 100.00000  34.78261
#>         1   0.00000   0.00000  65.21739
```

This would mean that cells in cluster 0 in the `small_9901_mix` dataset
are classified 100 % of the time as bleonging to cluster 0 of the
`small_5050_mix`

These frequencies can be visualized by plotting them in a heatmap

``` r

autoplot(freq_matrix,
         min_color = 50,
         show_number = TRUE)
#> Warning: Removed 3 rows containing missing values (geom_text).
```

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Here we can see that in the *9901* dataset (predicted), both clusters 0
and 1 are classified mostrly as cluster 0 in the *5050* dataset, while
the cluster 2 is mainly classified as 1.

(remember that the numbers are arbitrary and only mean something within
each dataset)

``` r
print(validation_results[[3]])
#> Cluster-0: 
#>  all elements:
#>      ARHGDIB +
#> Cluster-1: 
#>  all elements:
#>      ARHGDIB -
#>      SOX4 -
```

``` r
gating_genes <- validation_results$gating_genes
gating_genes
#> [1] "ARHGDIB" "SOX4"
```

``` r
g1 <- plot_flowstyle(small_5050_mix, markernames = gating_genes)
g2 <- plot_flowstyle(small_9901_mix, markernames = gating_genes)
 
g1
```

![](README_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
g2
```

![](README_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

``` r
g2[1,2]
```

![](README_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## Finding antibodies for the experiment

Since we acknowledge most experimental workflows need antibodies. We
have implemented several functions to look for antibodies in vendor
websites, as well as some helper functions to find the other posible
aliases a gene might have.

Here is a simple example for a gene widely know to have an antibody
available

``` r
head(query_biocompare_antibodies("CD11b"))
#>                                  title            vendor
#> 1        Anti-CD11b antibody [EPR1344]             Abcam
#> 2         Anti-CD11b/c antibody [OX42]             Abcam
#> 3     InVivoMab anti-mouse/human CD11b        Bio X Cell
#> 4     InVivoMab anti-mouse/human CD11b        Bio X Cell
#> 5 Monoclonal Antibody to CD11b (human)   MyBioSource.com
#> 6               Anti-CD11b (Mouse) mAb MBL International
#>                                                                                                               specification
#> 1    Applications: WB, IHC-p; Reactivity: Hu, Ms, Rt, Pg, RhMk; Conjugate/Tag: Unconjugated; Quantity: 10 µl, 40 µl, 100 µl
#> 2 Applications: ELISA, FCM, ICC, IF, IHC-fr, IP; Reactivity: Rat (Rt); Conjugate/Tag: Unconjugated; Quantity: 10 µg, 100 µg
#> 3                             Applications: FCM, AfP, Neut; Reactivity: Hu, Ms; Conjugate/Tag: Unconjugated; Quantity: 1 mg
#> 4                            Applications: FCM, AfP, Neut; Reactivity: Hu, Ms; Conjugate/Tag: Unconjugated; Quantity: 50 mg
#> 5                 Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 0.1 mg
#> 6                 Applications: Flow Cytometry (FCM); Reactivity: Mouse (Ms); Conjugate/Tag: Unconjugated; Quantity: 100 ug
```

``` r
sessionInfo()
#> R version 3.6.0 (2019-04-26)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 18356)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.1252 
#> [2] LC_CTYPE=English_United States.1252   
#> [3] LC_MONETARY=English_United States.1252
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.1252    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Seurat_3.0.0      sctree_0.0.2.0001
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rtsne_0.15           colorspace_1.4-1     selectr_0.4-1       
#>   [4] ggridges_0.5.1       listenv_0.7.0        npsurv_0.4-0        
#>   [7] ggrepel_0.8.1        DT_0.6               bit64_0.9-7         
#>  [10] AnnotationDbi_1.46.0 mvtnorm_1.0-10       ranger_0.11.2       
#>  [13] xml2_1.2.0           codetools_0.2-16     splines_3.6.0       
#>  [16] R.methodsS3_1.7.1    lsei_1.2-0           libcoin_1.0-4       
#>  [19] knitr_1.22           Formula_1.2-3        jsonlite_1.6        
#>  [22] ica_1.0-2            cluster_2.0.8        png_0.1-7           
#>  [25] R.oo_1.22.0          shiny_1.3.2          sctransform_0.2.0   
#>  [28] compiler_3.6.0       httr_1.4.0           assertthat_0.2.1    
#>  [31] Matrix_1.2-17        lazyeval_0.2.2       later_0.8.0         
#>  [34] htmltools_0.3.6      tools_3.6.0          rsvd_1.0.0          
#>  [37] igraph_1.2.4.1       partykit_1.2-3       gtable_0.3.0        
#>  [40] glue_1.3.1           RANN_2.6.1           reshape2_1.4.3      
#>  [43] dplyr_0.8.0.1        Rcpp_1.0.1           Biobase_2.44.0      
#>  [46] gdata_2.18.0         ape_5.3              nlme_3.1-139        
#>  [49] gbRd_0.4-11          lmtest_0.9-37        inum_1.0-1          
#>  [52] xfun_0.6             stringr_1.4.0        globals_0.12.4      
#>  [55] rvest_0.3.3          mime_0.6             irlba_2.3.3         
#>  [58] gtools_3.8.1         future_1.13.0        MASS_7.3-51.4       
#>  [61] zoo_1.8-5            scales_1.0.0         promises_1.0.1      
#>  [64] parallel_3.6.0       RColorBrewer_1.1-2   curl_3.3            
#>  [67] yaml_2.2.0           memoise_1.1.0        reticulate_1.12     
#>  [70] pbapply_1.4-0        gridExtra_2.3        ggplot2_3.1.1       
#>  [73] rpart_4.1-15         reshape_0.8.8        stringi_1.4.3       
#>  [76] RSQLite_2.1.1        S4Vectors_0.22.0     caTools_1.17.1.2    
#>  [79] BiocGenerics_0.30.0  bibtex_0.4.2         Rdpack_0.11-0       
#>  [82] SDMTools_1.1-221.1   rlang_0.3.4          pkgconfig_2.0.2     
#>  [85] bitops_1.0-6         evaluate_0.13        lattice_0.20-38     
#>  [88] ROCR_1.0-7           purrr_0.3.2          labeling_0.3        
#>  [91] htmlwidgets_1.3      cowplot_0.9.4        bit_1.1-14          
#>  [94] tidyselect_0.2.5     GGally_1.4.0         wrapr_1.8.6         
#>  [97] plyr_1.8.4           magrittr_1.5         R6_2.4.0            
#> [100] IRanges_2.18.0       gplots_3.0.1.1       DBI_1.0.0           
#> [103] pillar_1.4.0         fitdistrplus_1.0-14  survival_2.44-1.1   
#> [106] tibble_2.1.1         future.apply_1.2.0   tsne_0.1-3          
#> [109] crayon_1.3.4         KernSmooth_2.23-15   plotly_4.9.0        
#> [112] rmarkdown_1.12       viridis_0.5.1        grid_3.6.0          
#> [115] data.table_1.12.2    blob_1.1.1           metap_1.1           
#> [118] digest_0.6.18        xtable_1.8-4         httpuv_1.5.1        
#> [121] tidyr_0.8.3          R.utils_2.8.0        stats4_3.6.0        
#> [124] munsell_0.5.0        viridisLite_0.3.0
```

# Reproducing the runs in the purdue cluster

To reproduce the runs in the purdue cluster run as follows …

1.  We get the data from the temporary directory

<!-- end list -->

    git clone https://github.rcac.purdue.edu/jpaezpae/data_sctree data
    cd data
    bash untar_data.bash

2.  We run the standard seurat workflow.

This will output a report and generate an .RDS file for each of the
final seurat objects

    bash ./bash/build_jobs_seurat_workflow.bash ./data/filtered_matrices_mex_5050/hg19/ mix5050
    bash ./bash/build_jobs_seurat_workflow.bash ./data/filtered_matrices_mex_9901/hg19/ mix9901 

3.  Whenever those are done, run this …

This will run the benchmarks for the datasets. Will also generate 2 .RDS
files containing a list with a lot of stuff in it.

    for i in seurat*.RDS ; do bash ./bash/build_jobs_acc_benchmark.bash $i ; done

# Steps down the road

3.  Address some of the TODO’s in this repository
4.  Reduce dependecies by replacing functions to base equivalents.
5.  Add links to the documentation to make nicer to explore the package
    from inside R
6.  Implement plot that actually illustrates the progressive gating in
    the decision tree
7.  Implement a way to find markers for clusters exclusively upregulated
