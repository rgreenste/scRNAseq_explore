Reading in the data and removing empty droplets
================

Overview
--------

An exploratory analysis of public single-cell RNA-seq data from [Guiu et al](https://www.nature.com/articles/s41586-019-1212-5#Sec2) and retrieved from the Single Cell Expression Atlas

Dependencies
------------

Load the major packages needed for this analysis first. Additional packages will be loaded as needed.

``` r
# dependencies needed
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)
```

Obtain in the data from ArrayExpress
------------------------------------

``` r
## -- download and unzip data  -- ##

library(BiocFileCache)
bfc <- BiocFileCache("../raw_data", ask = FALSE)
guiu.zip <- bfcrpath(bfc, 
                    file.path("https://www.ebi.ac.uk/arrayexpress/files",
                              "E-MTAB-7660/E-MTAB-7660.processed.1.zip"))
unzip(guiu.zip, exdir = "../raw_data")
tarPath <- file.path("../raw_data/guiu2019_10xProcessed.tar.gz")
untar(tarPath, exdir = "../raw_data")
```

Read the 10x data into a SingleCellExperiment container
-------------------------------------------------------

Using a dedicated function for data collected from 10X Genomics Chromiun platform, read the raw files into R as a SingleCellExperiment (sce). Take a look at the sce to see what's in there.

``` r
## -- reading in the 10x data  -- ##

library(DropletUtils)
library(Matrix)
files <- file.path(getwd(), "../raw_data/guiu2019_10xProcessed/raw_gene_bc_matrices/mm10")
sce.guiu <- read10xCounts(files, col.names=TRUE)
sce.guiu
```

    ## class: SingleCellExperiment 
    ## dim: 27998 737280 
    ## metadata(1): Samples
    ## assays(1): counts
    ## rownames(27998): ENSMUSG00000051951 ENSMUSG00000089699 ...
    ##   ENSMUSG00000096730 ENSMUSG00000095742
    ## rowData names(2): ID Symbol
    ## colnames(737280): AAACCTGAGAAACCAT-1 AAACCTGAGAAACCGC-1 ...
    ##   TTTGTCATCTTTAGTC-1 TTTGTCATCTTTCCTC-1
    ## colData names(2): Sample Barcode
    ## reducedDimNames(0):
    ## spikeNames(0):
    ## altExpNames(0):

Additional exploration of the sce.

1.  How many droplets are present?

2.  How are the genes annotated?

3.  Are there spike in transcripts?

``` r
# do some exploring of this sce 
dim(sce.guiu) # how many rows and columns
```

    ## [1]  27998 737280

``` r
assays(sce.guiu) # what assays are present
```

    ## List of length 1
    ## names(1): counts

``` r
head(rownames(sce.guiu))# confirmed annotation with Ensembl notation
```

    ## [1] "ENSMUSG00000051951" "ENSMUSG00000089699" "ENSMUSG00000102343"
    ## [4] "ENSMUSG00000025900" "ENSMUSG00000109048" "ENSMUSG00000025902"

``` r
head(rowData(sce.guiu)) # annotated with gene ID and gene Symbol
```

    ## DataFrame with 6 rows and 2 columns
    ##                                    ID      Symbol
    ##                           <character> <character>
    ## ENSMUSG00000051951 ENSMUSG00000051951        Xkr4
    ## ENSMUSG00000089699 ENSMUSG00000089699      Gm1992
    ## ENSMUSG00000102343 ENSMUSG00000102343     Gm37381
    ## ENSMUSG00000025900 ENSMUSG00000025900         Rp1
    ## ENSMUSG00000109048 ENSMUSG00000109048         Rp1
    ## ENSMUSG00000025902 ENSMUSG00000025902       Sox17

``` r
names(colData(sce.guiu)) # sample names and cell barcode
```

    ## [1] "Sample"  "Barcode"

``` r
table(grepl("^ERCC", rownames(sce.guiu))) # there are no ERCC spike ins
```

    ## 
    ## FALSE 
    ## 27998

``` r
table(grepl("^SIRV", rownames(sce.guiu))) # there are no SIRV spike ins
```

    ## 
    ## FALSE 
    ## 27998

Removing empty droplets
-----------------------

First plot the UMI's per barcode by barcode rank. Allows us to determine "knee" and "inflection point" in the curve. Beyond the inflection point, there is a sharp decline in UMI's per barcode.

``` r
# -- compare barcode rank to total UMI count  -- #

barcodeRank <- barcodeRanks(counts(sce.guiu))
# Only showing unique points to speed up the plotting.

uniqBc <- !duplicated(barcodeRank$rank)
plot(barcodeRank$rank[uniqBc], barcodeRank$total[uniqBc], log="xy",
     xlab="Barcode Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(barcodeRank)$inflection, col="red", lty=2)
abline(h=metadata(barcodeRank)$knee, col="blue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("red", "blue"), lty=2, cex=1.2)
```

![](1_readData_countCells_files/figure-markdown_github/unnamed-chunk-5-1.png)

Next use the emptyDrops function to test whether each barcode expression profile is significantly different from the pool of ambient RNA in flowcell.Compute p-values and apply an FDR threshold to retain cells and remove empty drops.

``` r
# use emptyDrops function to remove empty drops (uses monte carlo sim so have to set.seed)
set.seed(101)

# calculates pValues to determine cells from empty drops with FDR 0.001
e.out <- emptyDrops(counts(sce.guiu))
summary(e.out$FDR <= 0.001)
```

    ##    Mode   FALSE    TRUE    NA's 
    ## logical    6639    5289  725352

``` r
# subset the sce based on the columns which pass our FDR cutoff of 0.1%
sce.guiu <- sce.guiu[,which(e.out$FDR <= 0.001)]
dim(sce.guiu) # now there are only 5289 cells remaining
```

    ## [1] 27998  5289

Output
------

After removing empty drops, output the resulting sce as unfiltered and save the file to the processed\_data folder

``` r
unfiltered <- sce.guiu # save the sce post removing empty's but before additional QC
saveRDS(unfiltered, file="../processed_data/guiu_unfiltered.rds")
```

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] DropletUtils_1.6.1          BiocFileCache_1.10.2       
    ##  [3] dbplyr_1.4.2                Rtsne_0.15                 
    ##  [5] uwot_0.1.5                  Matrix_1.2-18              
    ##  [7] scran_1.14.5                scater_1.14.6              
    ##  [9] ggplot2_3.2.1               SingleCellExperiment_1.8.0 
    ## [11] SummarizedExperiment_1.16.1 DelayedArray_0.12.2        
    ## [13] BiocParallel_1.20.1         matrixStats_0.55.0         
    ## [15] Biobase_2.46.0              GenomicRanges_1.38.0       
    ## [17] GenomeInfoDb_1.22.0         IRanges_2.20.2             
    ## [19] S4Vectors_0.24.2            BiocGenerics_0.32.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-6             bit64_0.9-7              httr_1.4.1              
    ##  [4] tools_3.6.2              backports_1.1.5          R6_2.4.1                
    ##  [7] irlba_2.3.3              HDF5Array_1.14.1         vipor_0.4.5             
    ## [10] DBI_1.1.0                lazyeval_0.2.2           colorspace_1.4-1        
    ## [13] withr_2.1.2              tidyselect_0.2.5         gridExtra_2.3           
    ## [16] bit_1.1-14               curl_4.3                 compiler_3.6.2          
    ## [19] BiocNeighbors_1.4.1      scales_1.1.0             rappdirs_0.3.1          
    ## [22] stringr_1.4.0            digest_0.6.23            R.utils_2.9.2           
    ## [25] rmarkdown_2.0            XVector_0.26.0           pkgconfig_2.0.3         
    ## [28] htmltools_0.4.0          limma_3.42.0             rlang_0.4.2             
    ## [31] RSQLite_2.2.0            DelayedMatrixStats_1.8.0 R.oo_1.23.0             
    ## [34] dplyr_0.8.3              RCurl_1.95-4.12          magrittr_1.5            
    ## [37] BiocSingular_1.2.1       GenomeInfoDbData_1.2.2   Rhdf5lib_1.8.0          
    ## [40] Rcpp_1.0.3               ggbeeswarm_0.6.0         munsell_0.5.0           
    ## [43] viridis_0.5.1            R.methodsS3_1.7.1        lifecycle_0.1.0         
    ## [46] stringi_1.4.5            yaml_2.2.0               edgeR_3.28.0            
    ## [49] zlibbioc_1.32.0          rhdf5_2.30.1             grid_3.6.2              
    ## [52] blob_1.2.0               dqrng_0.2.1              crayon_1.3.4            
    ## [55] lattice_0.20-38          locfit_1.5-9.1           zeallot_0.1.0           
    ## [58] knitr_1.26               pillar_1.4.3             igraph_1.2.4.2          
    ## [61] glue_1.3.1               evaluate_0.14            RcppParallel_4.4.4      
    ## [64] vctrs_0.2.1              gtable_0.3.0             purrr_0.3.3             
    ## [67] assertthat_0.2.1         xfun_0.12                rsvd_1.0.2              
    ## [70] viridisLite_0.3.0        tibble_2.1.3             beeswarm_0.2.3          
    ## [73] memoise_1.1.0            statmod_1.4.33
