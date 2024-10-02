# UKBB-Diverse-pops
GWAS across all analyzable phenotypes and diverse ancestry groups in the UK Biobank.

A full description of how these summary statistics were generated and how to work with them can be found here: http://pan.ukbb.broadinstitute.org/. 

# Code overview
The code in this repo can be used to generate all the figures in the manuscript. They look for data files in ./data, and if they are not found, downloads them as needed from the public data repository.

Each figure panel can be generated individually, or figures as a whole. For instance, in R/fig3_spectrum.R, we provide a function called `figure3()` which can generate the entirety of figure 3. Alternatively, running the code inside the function can generate each figure panel separately. Note that for some figures, on some R setups, attempting to generate the full figure by calling the function directly can crash R: we are uncertain of the cause of the issue, but it can be resolved by running the code inside the function step-wise.

# System Requirements
## Hardware requirements
There are no specific hardware requirements, except enough RAM to fit the datasets in memory (~16Gb).

## Software requirements
### OS Requirements
This package has been tested on macOS Sonoma (14.5).

### R Dependencies
The code depends on R packages, including those listed at the top of `constants.R`, using R 4.3.3 (`sessionInfo()` used below).

# Installation Guide:

```
git clone https://github.com/atgu/ukbb_pan_ancestry
cd ukbb_pan_ancestry
R
```
- Install any packages needed, and run `constants.R` - if this runs without error, the remaining files will load.

## Citation

To cite the Pan-UKB resource, please use the citation: [Karczewski, Gupta, Kanai et al., 2024](https://medrxiv.org/content/10.1101/2024.03.13.24303864).

### Session info

```
> sessionInfo()
R version 4.3.3 (2024-02-29)
Platform: aarch64-apple-darwin23.2.0 (64-bit)
Running under: macOS Sonoma 14.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /opt/local/Library/Frameworks/R.framework/Versions/4.3/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] googlesheets4_1.1.1  colorspace_2.1-0     reticulate_1.35.0    magick_2.8.1         png_0.1-8           
 [6] IsoplotR_5.6         ggpmisc_0.6.0        ggpp_0.5.8-1         ggbreak_0.1.2        ggthemes_5.0.0      
[11] ggwordcloud_0.6.1    cowplot_1.1.1        RMySQL_0.10.27       DBI_1.1.3            ggrepel_0.9.4       
[16] pbapply_1.7-2        rlang_1.1.4          tidygraph_1.3.1      meta_6.5-0           ggrastr_1.0.2       
[21] ggpubr_0.6.0         ggridges_0.5.4       readxl_1.4.3         corrr_0.4.4          corrplot_0.92       
[26] patchwork_1.1.3.9000 naniar_1.0.0         plotROC_2.3.1        gghighlight_0.4.0    skimr_2.1.5         
[31] gapminder_1.0.0      trelliscopejs_0.2.10 scales_1.3.0         magrittr_2.0.3       slackr_3.3.1        
[36] plotly_4.10.3        broom_1.0.5          lubridate_1.9.3      forcats_1.0.0        stringr_1.5.1       
[41] dplyr_1.1.4          purrr_1.0.2          readr_2.1.4          tidyr_1.3.0          tibble_3.2.1        
[46] ggplot2_3.5.1        tidyverse_2.0.0      plyr_1.8.9           RCurl_1.98-1.16      Hmisc_5.1-1         

loaded via a namespace (and not attached):
  [1] splines_4.3.3           bitops_1.0-7            ggplotify_0.1.2         locusviz_0.1.0         
  [5] cellranger_1.1.0        rpart_4.1.21            lifecycle_1.0.4         rstatix_0.7.2          
  [9] doParallel_1.0.17       vroom_1.6.4             lattice_0.22-5          MASS_7.3-60            
 [13] backports_1.4.1         metafor_4.4-0           rmarkdown_2.25          minqa_1.2.6            
 [17] RColorBrewer_1.1-3      abind_1.4-5             fidelius_0.0.2          BiocGenerics_0.48.1    
 [21] yulab.utils_0.1.0       nnet_7.3-19             circlize_0.4.15         IRanges_2.36.0         
 [25] S4Vectors_0.40.2        MatrixModels_0.5-3      codetools_0.2-19        xml2_1.3.6             
 [29] tidyselect_1.2.0        shape_1.4.6             aplot_0.2.2             farver_2.1.1           
 [33] lme4_1.1-35.1           matrixStats_1.1.0       stats4_4.3.3            base64enc_0.1-3        
 [37] googledrive_2.1.1       webshot_0.5.5           mathjaxr_1.6-0          jsonlite_1.8.8         
 [41] GetoptLong_1.0.5        Formula_1.2-5           survival_3.5-7          iterators_1.0.14       
 [45] foreach_1.5.2           tools_4.3.3             progress_1.2.2          Rcpp_1.0.11            
 [49] glue_1.6.2              gridExtra_2.3           xfun_0.41               withr_2.5.2            
 [53] numDeriv_2016.8-1.1     fastmap_1.1.1           boot_1.3-28.1           fansi_1.0.5            
 [57] SparseM_1.81            digest_0.6.33           timechange_0.2.0        R6_2.5.1               
 [61] gridGraphics_0.5-1      visdat_0.6.0            scattermore_1.2         utf8_1.2.4             
 [65] generics_0.1.3          data.table_1.14.8       prettyunits_1.2.0       httr_1.4.7             
 [69] htmlwidgets_1.6.3       pkgconfig_2.0.3         gtable_0.3.4            ComplexHeatmap_2.18.0  
 [73] htmltools_0.5.7         carData_3.0-5           autocogs_0.1.4          clue_0.3-65            
 [77] ggfun_0.1.3             knitr_1.45              rstudioapi_0.15.0       tzdb_0.4.0             
 [81] rjson_0.2.21            curl_5.1.0              checkmate_2.3.1         nlme_3.1-164           
 [85] nloptr_2.0.3            repr_1.1.6              cachem_1.0.8            GlobalOptions_0.1.2    
 [89] parallel_4.3.3          vipor_0.4.5             metadat_1.2-0           foreign_0.8-86         
 [93] pillar_1.9.0            vctrs_0.6.5             car_3.1-2               cluster_2.1.6          
 [97] beeswarm_0.4.0          htmlTable_2.4.2         evaluate_0.23           cli_3.6.1              
[101] compiler_4.3.3          crayon_1.5.2            ggsignif_0.6.4          labeling_0.4.3         
[105] mclust_6.0.1            fs_1.6.3                ggbeeswarm_0.7.2        stringi_1.8.4          
[109] viridisLite_0.4.2       munsell_0.5.0           lazyeval_0.2.2          CompQuadForm_1.4.3     
[113] quantreg_5.97           Matrix_1.6-4            hms_1.1.3               bit64_4.0.5            
[117] gridtext_0.1.5          gargle_1.5.2            igraph_2.0.3            memoise_2.0.1          
[121] bit_4.0.5               DistributionUtils_0.6-1 polynom_1.4-1
```
