sctree: a package to connect single cell rna-seq to biology using trees
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# sctree

The goal of sctree is to create a tool to accelerate the transition from
single cell rna-sequencing to calidation and new sub-population
discovery.

Features suggesting pseudo-gating strategies to purify found populations
via flow-cytometry, antibody querying and cross validations between
datasets.

Number of lines in roxygen comments: 0

Number of lines in R code: 0

# Installation

    git clone https://github.rcac.purdue.edu/jpaezpae/sctree sctree
    cd sctree
    
    R -e "devtools::install('.')"

# Usage

I am assuming you have already done your clustering and dimensional
reduction using seurat and we have our seurat object.

For this examples we will use a dummy dataset that come bundled with the
package

`seurat_mix5050`, this dataset comes originally from the 1:1 mixture of
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
require(Seurat)
#> Loading required package: Seurat
#> Loading required package: ggplot2
#> Loading required package: cowplot
#> 
#> Attaching package: 'cowplot'
#> The following object is masked from 'package:ggplot2':
#> 
#>     ggsave
#> Loading required package: Matrix

set.seed(6)

seurat_mix5050 <- readRDS("seurat_mix5050.RDS")
seurat_mix9901 <- readRDS("seurat_mix9901.RDS")

TSNEPlot(seurat_mix9901)
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
TSNEPlot(seurat_mix5050)
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

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
rang_importances <- ranger_importances.seurat(
    seurat_mix5050,
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
#> Sample size:                      3388 
#> Number of independent variables:  538 
#> Mtry:                             107 
#> Target node size:                 1 
#> Variable importance mode:         impurity_corrected 
#> Splitrule:                        gini 
#> OOB prediction error:             0.24 %
```

We can see that our classifier is 0.24 % accurate, as measured by its
*OOB prediction error*.

Your classifier USUALLY will be just as good as your initial clustering
is, therefore if your clustering is not all that consistent to start
with, it is unrealistic that predicions will be done accurately.

Now … taking a look at the importance matrix we can see that a `pvalue`
and a relative importance has been given to each gene.

``` r
head(rang_importances[[2]])
#>                importance      pvalue
#> HES4           1.40955457 0.000000000
#> ISG15          3.92726194 0.000000000
#> TNFRSF4        0.11080444 0.005780347
#> RP3-395M20.12 -0.02833254 0.965317919
#> KCNAB2         0.02682665 0.052023121
#> CTNNBIP1      -0.03562707 0.988439306
```

For simplicity we return by default a `data.frame` containing the ones
with `pvalue < 0.05`

``` r
head(rang_importances[[3]])
#>         importance pvalue    gene
#> TMSB4X    141.5518      0  TMSB4X
#> ARHGDIB   139.4114      0 ARHGDIB
#> CKB       133.3990      0     CKB
#> XIST      132.2081      0    XIST
#> CDKN2A    119.2037      0  CDKN2A
#> ADA       111.5039      0     ADA
```

Therefore, in this case we can say that the expression of the following
genes would be usefull to form the clusters.

## Visualizing the expected outcome of a flow cytometry experiment

Lets say we choose the top 5 markers from the former list and we did a
flow experiment … HYPOTHETICALLY the marker distribution would resemble
the rna expression profile for which we have the function
`plot_flowstyle`

``` r
top_markers <- head(rang_importances[[3]]$gene)
top_markers
#> [1] "TMSB4X"  "ARHGDIB" "CKB"     "XIST"    "CDKN2A"  "ADA"
g <- plot_flowstyle(seurat_mix5050, markernames = top_markers)
g
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Based on this, we can see that the red cluster in this plot is
predominantly CD3+ ADA+ and ARHGDIB+, as well as ASNS-

We can also focus in one of the pannels (and check the color
conventions)

``` r
g[1,2]
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Suggesting a gating strategy for the markers

A general strategy to get separate all clusters

``` r

top_markers <- head(rang_importances[[3]]$gene)

tree_fit <- fit_ctree(seurat_mix5050, genes_use = top_markers, cluster = "ALL")
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

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

When inspecting the tree\_fit, we can see a more detailed text
representation of this tree.

``` r
print(tree_fit)
#> 
#> Model formula:
#> ident ~ TMSB4X + ARHGDIB + CKB + XIST + CDKN2A + ADA
#> 
#> Fitted party:
#> [1] root
#> |   [2] CKB <= 1.32554
#> |   |   [3] XIST <= 1.12625
#> |   |   |   [4] ARHGDIB <= 1.20526: 0 (n = 18, err = 5.6%)
#> |   |   |   [5] ARHGDIB > 1.20526: 0 (n = 1734, err = 0.0%)
#> |   |   [6] XIST > 1.12625: 1 (n = 18, err = 50.0%)
#> |   [7] CKB > 1.32554
#> |   |   [8] ARHGDIB <= 1.29467
#> |   |   |   [9] XIST <= 0
#> |   |   |   |   [10] TMSB4X <= 0.7876: 3 (n = 19, err = 31.6%)
#> |   |   |   |   [11] TMSB4X > 0.7876: 3 (n = 7, err = 28.6%)
#> |   |   |   [12] XIST > 0
#> |   |   |   |   [13] TMSB4X <= 1.77781
#> |   |   |   |   |   [14] ARHGDIB <= 0.87204
#> |   |   |   |   |   |   [15] XIST <= 1.23346
#> |   |   |   |   |   |   |   [16] CDKN2A <= 1.92909: 1 (n = 25, err = 4.0%)
#> |   |   |   |   |   |   |   [17] CDKN2A > 1.92909: 1 (n = 7, err = 28.6%)
#> |   |   |   |   |   |   [18] XIST > 1.23346
#> |   |   |   |   |   |   |   [19] ARHGDIB <= 0.69905: 1 (n = 1437, err = 0.0%)
#> |   |   |   |   |   |   |   [20] ARHGDIB > 0.69905: 1 (n = 30, err = 3.3%)
#> |   |   |   |   |   [21] ARHGDIB > 0.87204
#> |   |   |   |   |   |   [22] TMSB4X <= 0.8868: 1 (n = 19, err = 0.0%)
#> |   |   |   |   |   |   [23] TMSB4X > 0.8868: 1 (n = 7, err = 42.9%)
#> |   |   |   |   [24] TMSB4X > 1.77781: 2 (n = 16, err = 50.0%)
#> |   |   [25] ARHGDIB > 1.29467
#> |   |   |   [26] CDKN2A <= 0.73944: 0 (n = 7, err = 28.6%)
#> |   |   |   [27] CDKN2A > 0.73944: 2 (n = 44, err = 0.0%)
#> 
#> Number of inner nodes:    13
#> Number of terminal nodes: 14
```

Sometimes one might think that the proposed strategy is too complicated
or not implementable in the experimental settings, in order to add
constrians to the fit one can give additional arguments that will be
passed to `partykit::ctree_control`, such as `maxdepth = 2` (maximum 2
questions per cell)

``` r
tree_fit <- fit_ctree(
  seurat_mix5050, genes_use = top_markers, 
  cluster = "ALL", maxdepth = 2)
print(tree_fit)
#> 
#> Model formula:
#> ident ~ TMSB4X + ARHGDIB + CKB + XIST + CDKN2A + ADA
#> 
#> Fitted party:
#> [1] root
#> |   [2] CKB <= 1.32554
#> |   |   [3] XIST <= 1.12625: 0 (n = 1752, err = 0.1%)
#> |   |   [4] XIST > 1.12625: 1 (n = 18, err = 50.0%)
#> |   [5] CKB > 1.32554
#> |   |   [6] ARHGDIB <= 1.29467: 1 (n = 1567, err = 2.2%)
#> |   |   [7] ARHGDIB > 1.29467: 2 (n = 51, err = 9.8%)
#> 
#> Number of inner nodes:    3
#> Number of terminal nodes: 4
plot(tree_fit)
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Since not all variables are ultimately used in our classifier, one can
acces the ones that were by using `varimp(tree_fit)`

``` r
partykit::varimp(tree_fit)
#>        CKB       XIST    ARHGDIB 
#> 0.82866521 0.03625885 0.33081234
plot_flowstyle(seurat_mix5050, names(partykit::varimp(tree_fit)))
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

One can also request the package to suggest a specific strategy only for
a given cluster. This function is not expected to give drastically
different results in datasets with few clusters, but it can definitely
come usefull when many clusters are present and one is interested in a
specific
one.

``` r
tree_fit <- fit_ctree(seurat_mix5050, genes_use = top_markers, cluster = "0")
print(tree_fit)
#> 
#> Model formula:
#> ident ~ TMSB4X + ARHGDIB + CKB + XIST + CDKN2A + ADA
#> 
#> Fitted party:
#> [1] root
#> |   [2] CKB <= 1.32554
#> |   |   [3] XIST <= 1.12625
#> |   |   |   [4] ARHGDIB <= 1.20526: indeed clus 0 (n = 18, err = 5.6%)
#> |   |   |   [5] ARHGDIB > 1.20526: indeed clus 0 (n = 1734, err = 0.0%)
#> |   |   [6] XIST > 1.12625: not clus 0 (n = 18, err = 33.3%)
#> |   [7] CKB > 1.32554
#> |   |   [8] ARHGDIB <= 1.9863
#> |   |   |   [9] TMSB4X <= 2.72634: not clus 0 (n = 1584, err = 0.0%)
#> |   |   |   [10] TMSB4X > 2.72634: not clus 0 (n = 17, err = 5.9%)
#> |   |   [11] ARHGDIB > 1.9863: not clus 0 (n = 17, err = 29.4%)
#> 
#> Number of inner nodes:    5
#> Number of terminal nodes: 6
```

## Finding equivalent clusters in two datasets

``` r
data(seurat_mix9901)
#> Warning in data(seurat_mix9901): data set 'seurat_mix9901' not found
seurat_mix9901
#> An object of class seurat in project SeuratProject 
#>  32738 genes across 4185 samples.
```

``` r
validation_results <- cross_validate(
    seurat_mix5050, seurat_mix9901, 
    cluster = "ALL")
#> Warning in ranger::importance_pvalues(ranger_fit): Only few negative
#> importance values found, inaccurate p-values. Consider the 'altmann'
#> approach.
#> Warning in max(elems_per_rule): no non-missing arguments to max; returning
#> -Inf

#> Warning in max(elems_per_rule): no non-missing arguments to max; returning
#> -Inf

#> Warning in max(elems_per_rule): no non-missing arguments to max; returning
#> -Inf

validation_results[[1]]
#> 
#> Model formula:
#> ident ~ CKB + TMSB4X + CD3D + EIF5A + XIST + ARHGDIB + CDKN2A + 
#>     ADA + MZB1 + OCIAD2 + LCK + NUCB2 + RPS4X + GAPDH + ITM2A + 
#>     CA2 + RPL29 + CD3G + BEX1 + PSMB8 + GMFG + FTH1 + LEF1 + 
#>     CBR1 + RNASEH2B + GAL + GTF3A + TSPO + `MT-CO3` + `MT-CO1` + 
#>     CD1E + NT5C3B + CHCHD10 + FYB + ISG15 + RAB13 + DDT + CALD1 + 
#>     RPL12 + DMKN + MDK + RPL26 + TUBB + MALAT1 + TYMS + SH3BGRL3 + 
#>     CDK6 + MAP1A + PTPN7 + PCSK1N + HEY1 + CHCHD5 + AIF1 + UCHL1 + 
#>     EIF4EBP1 + C3orf14 + PRDX4 + CALR + RAB34 + HES4 + GYPC + 
#>     PCNA + SELM + HOXA10 + IGLL1 + RBBP7 + UBE2L6 + CDC42EP1 + 
#>     PGK1 + HSP90AA1 + RPL36A + CHI3L2 + FKBP10 + HOXA9 + TUBA1B + 
#>     ZNF703 + GLO1 + DNAJA1 + UBC + RAC3 + `HHIP-AS1` + ARPP21 + 
#>     NME3 + CD63 + B2M + TSTD1 + CD3E + GLRX + TCF7 + NEURL1 + 
#>     CORO1A + RHOH + ACTB + SPN + PLS3 + ID2 + RPP25 + S100A11 + 
#>     STMN1 + TFDP2 + RAB40B + ECI1 + MORF4L2 + RPLP1 + CCNB2 + 
#>     MGMT + CD151 + PGRMC1 + GSDMD + AK1 + S100A13 + TCEAL4 + 
#>     PTTG1 + AC002454.1 + TFAP2A + KRT18 + TUBB4B + MYO1G + SLIT2 + 
#>     STRA13 + HSP90B1 + MYB + FTL + ENAH + ARF1 + PODXL2 + CNN3 + 
#>     ACTG1 + VAMP5 + SHQ1 + CDK1 + HSD17B12 + NAA10 + CD1B + DCAF13 + 
#>     PRDX6 + HSPA1B + ZNF503 + SH2D1A + PERP + MRP63 + LCP1 + 
#>     CXCR4 + SIVA1 + VIM + PPP1R18 + CKLF + ASF1A + ZNF22 + TPTEP1 + 
#>     CENPF + MT2A + USP14 + ANXA2 + PSMB9 + PYGL + PSME1 + CAPN2 + 
#>     SYNGR1 + PPIC + MCM5 + EPM2AIP1 + HSPA8 + CSNK1E + IFI16 + 
#>     SCML1 + ME2 + ACAT2 + BTG1 + MAP1B + SMC2 + IL2RG + CST3 + 
#>     FUCA2 + PFN2 + LINC00632 + FABP5 + LPXN + DEPDC1 + ACSF3 + 
#>     CEP70 + NDFIP2 + BIRC5 + CARHSP1 + DNAJC15 + COA3 + NMRAL1 + 
#>     DEGS1 + KPNA2 + AEBP1 + MLF1 + CBX5 + RAC2 + HOTAIRM1 + EIF5AL1 + 
#>     HNRNPF + TNFRSF4 + TIMP1 + NMI + C14orf64 + RFC4 + LDLRAD4 + 
#>     KIAA1430 + CBS + KNSTRN + MTAP + C12orf65 + HOXB6 + SEPW1 + 
#>     CDC20 + DPP7 + CSRP2 + PPP2CA + KCNAB2 + NEAT1 + TPX2 + MFNG + 
#>     PTPRCAP + RRM2 + ARHGAP30 + PINLYP + CKS2 + FAM127B + TM7SF3 + 
#>     FBXO25 + HIST1H4C + SKP1 + NME4 + LAPTM5 + XBP1 + FHL1 + 
#>     GNA15 + TRAF3IP3 + CCNB1 + TSC22D3 + TMSB10 + CD1D + ANXA5 + 
#>     CCND3 + SRSF6 + MSN + SELL + HIST1H2BK + C12orf57 + FNBP1 + 
#>     TCTEX1D2 + MYO7B + IDH2 + CENPE + LIN28B + NAPRT1 + ITGA4 + 
#>     FDFT1 + NUF2 + EPDR1 + ANKRD12 + HCLS1 + MYC + SPTBN1 + KIAA0226L + 
#>     EMP3 + HMGB3 + BID + TSC22D1 + FSCN1 + AES + JUN + KIAA0020 + 
#>     GPSM3 + SLA + PLP2 + CDCA8 + ATP6V1H + OXCT1 + MT1X + TEX30 + 
#>     TK1 + ICAM2 + S100A10 + PAK2 + SOX4 + ARID5B + NUDT8 + PCID2 + 
#>     ALDH1A2 + BAG2 + ID4 + RASGRP2 + `RP11-973H7.4` + CITED1 + 
#>     ZAP70 + HENMT1 + DNTT + `PRKAG2-AS1` + CD247 + DHRS4 + ATP6AP1L + 
#>     SEPT1 + FES + ELOVL5 + NDN + MAL + SGOL2 + CDKN2AIP + NUSAP1 + 
#>     AURKA + SLC43A3 + RAP1GDS1 + CELF2 + SELPLG + PSME2 + CFL2 + 
#>     CD7 + IER2 + EIF4A3 + PBK + JAKMIP1 + GRAP2 + VAMP8 + CTSG + 
#>     BEX2 + HIST2H2AC + SMC1B + TBC1D10C + RCSD1 + CTNNBIP1 + 
#>     TUBA1A + EVL + SIT1 + FXYD5 + HMGB2 + RHOB + FERMT3 + CPNE1 + 
#>     APLP2 + CENPA + CD47 + DUSP1 + C5orf56 + SPAG5 + ICAM3 + 
#>     ETHE1 + CMPK2 + PTPRC + CASP8 + PSAT1 + CCZ1B + RAG1 + RBM3 + 
#>     `RP11-1094M14.11`
#> 
#> Fitted party:
#> [1] root
#> |   [2] CD3D <= 1.23381
#> |   |   [3] OCIAD2 <= 0.90609: 1 (n = 1565, err = 1.8%)
#> |   |   [4] OCIAD2 > 0.90609: 0 (n = 16, err = 18.8%)
#> |   [5] CD3D > 1.23381
#> |   |   [6] XIST <= 1.20355
#> |   |   |   [7] GAL <= 0.76111: 0 (n = 1734, err = 0.1%)
#> |   |   |   [8] GAL > 0.76111: 0 (n = 18, err = 33.3%)
#> |   |   [9] XIST > 1.20355
#> |   |   |   [10] ARHGAP30 <= 0.38927
#> |   |   |   |   [11] EIF5A <= 2.09015: 2 (n = 11, err = 27.3%)
#> |   |   |   |   [12] EIF5A > 2.09015: 2 (n = 37, err = 0.0%)
#> |   |   |   [13] ARHGAP30 > 0.38927: 0 (n = 7, err = 42.9%)
#> 
#> Number of inner nodes:    6
#> Number of terminal nodes: 7
```

``` r
validation_results$confusion_matrix
#>          cluster
#> predicted    0    1    2    3    4
#>         0 1814 1068  945  327    1
#>         1    0    0    1    0   27
#>         2    0    2    0    0    0
#>         3    0    0    0    0    0
```

``` r
freq_matrix <- as.frequency.matrix(validation_results$confusion_matrix)
freq_matrix
#>          cluster
#> predicted           0           1           2           3           4
#>         0 100.0000000  99.8130841  99.8942918 100.0000000   3.5714286
#>         1   0.0000000   0.0000000   0.1057082   0.0000000  96.4285714
#>         2   0.0000000   0.1869159   0.0000000   0.0000000   0.0000000
#>         3   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
```

This would mean that cells in cluster 0 in the `seurat_mix9901` dataset
are classified 99.23 % of the time as bleonging to cluster 0 of the
`seurat_mix5050`

These frequencies can be visualized by plotting them in a heatmap

``` r

autoplot.table(freq_matrix,
               min_color = 50,
               show_number = TRUE)
#> Warning: Removed 15 rows containing missing values (geom_text).
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Here we can see that in the *9901* dataset (predicted), both clusters 0
and 1 are classified mostrly as cluster 0 in the *5050* dataset, while
the cluster 2 is mainly classified as 1.

(remember that the numbers are arbitrary and only mean something within
each dataset)

``` r
validation_results[[3]]
#> $`0`
#> $`0`$all
#> [1] "CD3D > 1.23380572030749"
#> 
#> $`0`$majority
#> [1] "XIST <= 1.20354519639701" "GAL <= 0.761105335702687"
#> 
#> 
#> $`1`
#> $`1`$all
#> [1] "CD3D <= 1.23380572030749"    "OCIAD2 <= 0.906093418104814"
#> 
#> $`1`$majority
#> character(0)
#> 
#> 
#> $`2`
#> $`2`$all
#> [1] "CD3D > 1.23380572030749"       "XIST > 1.20354519639701"      
#> [3] "ARHGAP30 <= 0.389265177607656"
#> 
#> $`2`$majority
#> [1] "EIF5A > 2.09014605183273"
#> 
#> 
#> $`3`
#> $`3`$all
#> NULL
#> 
#> $`3`$majority
#> NULL
```

``` r
gating_genes <- validation_results$gating_genes
gating_genes
#> [1] "CD3D"     "OCIAD2"   "XIST"     "GAL"      "ARHGAP30" "EIF5A"
```

``` r
g1 <- plot_flowstyle(seurat_mix5050, markernames = gating_genes)
g2 <- plot_flowstyle(seurat_mix9901, markernames = gating_genes)
 
g1
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
g2
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

``` r
g2[1,2]
```

![](sctree_readme_full_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

## Finding antibodies for the experiment

Since we acknowledge most experimental workflows need antibodies. We
have implemented several functions to look for antibodies in vendor
websites, as well as some helper functions to find the other posible
aliases a gene might have.

Here is a simple example for a gene widely know to have an antibody
available

``` r
head(query_biocompare_antibodies("CD11b"))
#>                                                                       title
#> 1                                                 Anti-CD11b/ITGAM Antibody
#> 2                                        Anti-CD11b/ITGAM Picoband Antibody
#> 3 Anti-Human CD11b Monoclonal Antibody mFluor450 Conjugated, Flow Validated
#> 4     Anti-CD11b (integrin alpha-M) Rabbit Monoclonal Antibody, Clone#RM290
#> 5                                      Monoclonal Antibody to CD11b (human)
#> 6                                                      Mouse Anti-Rat CD11b
#>            vendor
#> 1       BosterBio
#> 2       BosterBio
#> 3       BosterBio
#> 4       BosterBio
#> 5 MyBioSource.com
#> 6      RayBiotech
#>                                                                                                          specification
#> 1           Applications: Western Blot (WB); Reactivity: Hu, Ms, Rt; Conjugate/Tag: Unconjugated; Quantity: 100ug/vial
#> 2 Applications: WB, FCM, ICC, IHC-fr, IHC-p; Reactivity: Hu, Ms, Rt; Conjugate/Tag: Unconjugated; Quantity: 100ug/vial
#> 3  Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: mFluor450; Quantity: 25 Tests, 100 Tests
#> 4                          Applications: WB, IHC; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100uL
#> 5            Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 0.1 mg
#> 6                                    Applications: Flow Cytometry (FCM); Reactivity: Rat (Rt); Quantity: 500 <U+00B5>g
```

``` r
aliases <- get_aliases(gating_genes[[1]])

print(aliases)
#> $CD3D
#> [1] "CD3-DELTA" "IMD19"     "T3D"       "CD3D"

lapply(aliases[[1]], function(alias) {
     head(query_biocompare_antibodies(alias))
})
#> [[1]]
#>                                                                title
#> 1  Anti-Human CD3 Monoclonal Antibody APC Conjugated, Flow Validated
#> 2 Anti-Human CD3 Monoclonal Antibody FITC Conjugated, Flow Validated
#> 3   Anti-Human CD3 Monoclonal Antibody PE Conjugated, Flow Validated
#> 4                               CD3d / CD3 delta Antibody, Mouse MAb
#> 5                                                 CD3 delta antibody
#> 6                                      CD3-delta Polyclonal Antibody
#>                  vendor
#> 1             BosterBio
#> 2             BosterBio
#> 3             BosterBio
#> 4 Sino Biological, Inc.
#> 5               Biorbyt
#> 6       MyBioSource.com
#>                                                                                                                      specification
#> 1         Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: APC; Quantity: 25 tests, 100 tests, 500 tests
#> 2                   Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: FITC; Quantity: 25 tests, 100 tests
#> 3          Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: PE; Quantity: 25 tests, 100 tests, 500 tests
#> 4 Applications: Immunohistochemistry-paraffin (IHC-p); Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100<U+00B5>L
#> 5                                   Applications: WB, ELISA; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100 ug
#> 6                                  Applications: WB, ELISA; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 0.05 ml
#> 
#> [[2]]
#>                                                                title
#> 1                               Anti-CD3D Rabbit Monoclonal Antibody
#> 2  Anti-Human CD3 Monoclonal Antibody APC Conjugated, Flow Validated
#> 3 Anti-Human CD3 Monoclonal Antibody FITC Conjugated, Flow Validated
#> 4   Anti-Human CD3 Monoclonal Antibody PE Conjugated, Flow Validated
#> 5                                           CD3D Monoclonal Antibody
#> 6                                           CD3D Monoclonal Antibody
#>                 vendor
#> 1            BosterBio
#> 2            BosterBio
#> 3            BosterBio
#> 4            BosterBio
#> 5 LifeSpan BioSciences
#> 6 LifeSpan BioSciences
#>                                                                                                              specification
#> 1                      Applications: WB, IF, IP; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100ug/vial
#> 2 Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: APC; Quantity: 25 tests, 100 tests, 500 tests
#> 3           Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: FITC; Quantity: 25 tests, 100 tests
#> 4  Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: PE; Quantity: 25 tests, 100 tests, 500 tests
#> 5                         Applications: ELISA; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 50 <U+00B5>l
#> 6                    Applications: WB, ELISA; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100 <U+00B5>g
#> 
#> [[3]]
#>                                                                title
#> 1                               Anti-CD3D Rabbit Monoclonal Antibody
#> 2  Anti-Human CD3 Monoclonal Antibody APC Conjugated, Flow Validated
#> 3 Anti-Human CD3 Monoclonal Antibody FITC Conjugated, Flow Validated
#> 4   Anti-Human CD3 Monoclonal Antibody PE Conjugated, Flow Validated
#> 5                                                 anti-CD3D Antibody
#> 6                                                 anti-CD3D Antibody
#>              vendor
#> 1         BosterBio
#> 2         BosterBio
#> 3         BosterBio
#> 4         BosterBio
#> 5 antibodies-online
#> 6 antibodies-online
#>                                                                                                              specification
#> 1                      Applications: WB, IF, IP; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100ug/vial
#> 2 Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: APC; Quantity: 25 tests, 100 tests, 500 tests
#> 3           Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: FITC; Quantity: 25 tests, 100 tests
#> 4  Applications: Flow Cytometry (FCM); Reactivity: Human (Hu); Conjugate/Tag: PE; Quantity: 25 tests, 100 tests, 500 tests
#> 5                  Applications: WB, IHC; Reactivity: Hu, Ms, Rt, Mu; Conjugate/Tag: Unconjugated; Quantity: 200 <U+00B5>L
#> 6                   Applications: WB, IHC; Reactivity: Hu, Ms, Rt, Mu; Conjugate/Tag: Unconjugated; Quantity: 50 <U+00B5>L
#> 
#> [[4]]
#>                                  title                vendor
#> 1 Anti-CD3D Rabbit Monoclonal Antibody             BosterBio
#> 2     CD3D & CD3E Antibody, Rabbit PAb Sino Biological, Inc.
#> 3 CD3d / CD3 delta Antibody, Mouse MAb Sino Biological, Inc.
#> 4 CD3d / CD3 delta Antibody, Mouse MAb Sino Biological, Inc.
#> 5             CD3D Monoclonal Antibody  LifeSpan BioSciences
#> 6             CD3D Monoclonal Antibody  LifeSpan BioSciences
#>                                                                                                                      specification
#> 1                              Applications: WB, IF, IP; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100ug/vial
#> 2                                 Applications: ELISA; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100<U+00B5>L
#> 3 Applications: Immunohistochemistry-paraffin (IHC-p); Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100<U+00B5>L
#> 4                                 Applications: ELISA; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100<U+00B5>L
#> 5                                 Applications: ELISA; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 50 <U+00B5>l
#> 6                            Applications: WB, ELISA; Reactivity: Human (Hu); Conjugate/Tag: Unconjugated; Quantity: 100 <U+00B5>g
```

``` r
sessionInfo()
#> R version 3.5.2 (2018-12-20)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: CentOS Linux 7 (Core)
#> 
#> Matrix products: default
#> BLAS: /home/jpaezpae/opt/R-3.5.2/lib/libRblas.so
#> LAPACK: /home/jpaezpae/opt/R-3.5.2/lib/libRlapack.so
#> 
#> locale:
#> [1] C
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] Seurat_2.3.4      Matrix_1.2-15     cowplot_0.9.4     ggplot2_3.1.0    
#> [5] sctree_0.0.1.9002
#> 
#> loaded via a namespace (and not attached):
#>   [1] snow_0.4-3           backports_1.1.3      Hmisc_4.2-0         
#>   [4] selectr_0.4-1        wrapr_1.8.6          plyr_1.8.4          
#>   [7] igraph_1.2.4         lazyeval_0.2.2       splines_3.5.2       
#>  [10] digest_0.6.18        foreach_1.4.4        htmltools_0.3.6     
#>  [13] viridis_0.5.1        lars_1.2             gdata_2.18.0        
#>  [16] magrittr_1.5         checkmate_1.9.1      memoise_1.1.0       
#>  [19] cluster_2.0.7-1      mixtools_1.1.0       ROCR_1.0-7          
#>  [22] R.utils_2.8.0        colorspace_1.4-1     blob_1.1.1          
#>  [25] rvest_0.3.2          xfun_0.6             dplyr_0.8.0.1       
#>  [28] crayon_1.3.4         jsonlite_1.6         libcoin_1.0-4       
#>  [31] survival_2.43-3      zoo_1.8-5            iterators_1.0.10    
#>  [34] ape_5.3              glue_1.3.1           gtable_0.3.0        
#>  [37] kernlab_0.9-27       prabclus_2.2-7       BiocGenerics_0.26.0 
#>  [40] DEoptimR_1.0-8       scales_1.0.0         mvtnorm_1.0-10      
#>  [43] DBI_1.0.0            GGally_1.4.0         bibtex_0.4.2        
#>  [46] Rcpp_1.0.1           metap_1.1            dtw_1.20-1          
#>  [49] viridisLite_0.3.0    xtable_1.8-3         htmlTable_1.13.1    
#>  [52] reticulate_1.11.1    foreign_0.8-71       bit_1.1-14          
#>  [55] proxy_0.4-23         mclust_5.4.3         SDMTools_1.1-221    
#>  [58] Formula_1.2-3        stats4_3.5.2         tsne_0.1-3          
#>  [61] DT_0.5               htmlwidgets_1.3      httr_1.4.0          
#>  [64] gplots_3.0.1.1       RColorBrewer_1.1-2   fpc_2.1-11.1        
#>  [67] acepack_1.4.1        modeltools_0.2-22    ica_1.0-2           
#>  [70] pkgconfig_2.0.2      reshape_0.8.8        R.methodsS3_1.7.1   
#>  [73] flexmix_2.3-15       nnet_7.3-12          labeling_0.3        
#>  [76] tidyselect_0.2.5     rlang_0.3.3          reshape2_1.4.3      
#>  [79] later_0.8.0          AnnotationDbi_1.42.1 munsell_0.5.0       
#>  [82] tools_3.5.2          RSQLite_2.1.1        ranger_0.11.2       
#>  [85] ggridges_0.5.1       evaluate_0.13        stringr_1.4.0       
#>  [88] yaml_2.2.0           npsurv_0.4-0         org.Hs.eg.db_3.4.0  
#>  [91] knitr_1.22           bit64_0.9-7          fitdistrplus_1.0-14 
#>  [94] robustbase_0.93-4    caTools_1.17.1.2     purrr_0.3.2         
#>  [97] RANN_2.6.1           pbapply_1.4-0        nlme_3.1-137        
#> [100] mime_0.6             R.oo_1.22.0          xml2_1.2.0          
#> [103] hdf5r_1.1.1          compiler_3.5.2       rstudioapi_0.10     
#> [106] curl_3.3             png_0.1-7            lsei_1.2-0          
#> [109] tibble_2.1.1         stringi_1.4.3        lattice_0.20-38     
#> [112] trimcluster_0.1-2.1  pillar_1.3.1         Rdpack_0.10-1       
#> [115] lmtest_0.9-36        data.table_1.12.0    bitops_1.0-6        
#> [118] irlba_2.3.3          gbRd_0.4-11          httpuv_1.5.1        
#> [121] R6_2.4.0             latticeExtra_0.6-28  promises_1.0.1      
#> [124] KernSmooth_2.23-15   gridExtra_2.3        IRanges_2.14.12     
#> [127] codetools_0.2-16     MASS_7.3-51.1        gtools_3.8.1        
#> [130] assertthat_0.2.1     withr_2.1.2          S4Vectors_0.18.3    
#> [133] diptest_0.75-7       parallel_3.5.2       doSNOW_1.0.16       
#> [136] grid_3.5.2           rpart_4.1-13         tidyr_0.8.3         
#> [139] class_7.3-14         rmarkdown_1.11       inum_1.0-0          
#> [142] segmented_0.5-3.0    Rtsne_0.15           partykit_1.2-3      
#> [145] Biobase_2.40.0       shiny_1.2.0          base64enc_0.1-3
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
final seurat
    objects

    bash ./bash/build_jobs_seurat_workflow.bash ./data/filtered_matrices_mex_5050/hg19/ mix5050
    bash ./bash/build_jobs_seurat_workflow.bash ./data/filtered_matrices_mex_9901/hg19/ mix9901 

3.  Whenever those are done, run this …

This will run the benchmarks for the datasets. Will also generate 2 .RDS
files containing a list with a lot of stuff in
    it.

    for i in seurat*.RDS ; do bash ./bash/build_jobs_acc_benchmark.bash $i ; done

# Steps down the road

1.  Make figure list.
2.  Address some of the TODO’s in this repository
3.  Reduce dependecies by replacing functions to base equivalents.
4.  Add links to the documentation to make nicer to explore the package
    from inside R
5.  Implement plot that actually illustrates the progressive gating in
    the decision tree
6.  DRASTICALLY increase code coverage (right now everything runs
    without errors due to R CMD check) but expectations are not tested
7.  Implement a way to find markers for clusters exclusively upregulated
