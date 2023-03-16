# Abiotic and biotic regulators of species distribution

## Overview
This repository provides the data and analysis code for the manuscript "Experimental heatwaves and warming induce distinctive community responses through their interactions with a novel species"

As the climate warms, species shift their distributions at different rates, re-organising ecological communities. The resulting novel interactions will shape the local community’s response to ongoing climate change. The distinction between extreme events and a rising mean temperature in driving range expansion of the neighbouring species has not been examined empirically, nor has the resulting ecological impact propagating through multi-trophic networks been addressed. 

In this study, we recreated a high-elevation host-parasitoid community comprising Drosophila species and their associated parasitoid species from the Australian Wet Tropics, and subjected them to either heatwaves or warming in combination with the introduction of a low-elevation-specific Drosophila species. Specifically, we 1) analysed how population size and single-generation reproduction success changed by the combination of temperature and invasion treatments; 2) used NMDs to reveal the differences in community compositions associated with different treatments; 3) conducted path analysis (by piecewiseSEM) to understand the direct and indirect effects that different treatments have on reproduction of every species within a host-parasitoid community. 


## Authors information
Jinlin Chen (Department of Biology, University of Oxford)
Owen T. Lewis (Department of Biology, University of Oxford)

Authors contribution: JC and OTL both contributed to the development of ideas. JC designed and conducted the experimental work. JC analyzed the results and led the writing of the manuscript. OTL contributed to the writing.

Correspondence: Jinlin Chen (chen.jin.lin**`AT`**hotmail.com)


## Layout
The repository is split into 3 main directories. Each section is described below. 

### **`Analysis`** 
Where all of the *executed* analysis lives. This includes one scripts and Data. 
 * **`commCage_analysis.R`**: all code to analyze the data and generate tables, figures, statisticla results in the main text and the SI. 
 * **`commCageCensus.csv`**: data of end population size stored in csv for analysis.
 * **`commCageSample.csv`**: data of reproductive success stored in csv for analysis.
 * **`hmSample.csv`**: data of reproductive success before, during and after heatwave event stored in csv for analysis.

### **`Raw data`** 
The "JinlinChen_invasion experiment_entered on 0316.xlsx" file contains the full set of data from all three different group of measurement. Metadata is also provided in one of the spreadsheet in this document.
The "longtermClimate_formated.csv" file contains air temperature data during the field survey.

### **`Results`** 
Figures and Tables of the main text and supplementary materials. 

## Version control
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] piecewiseSEM_2.1.2 multcompView_0.1-8 vegan_2.5-7        lattice_0.20-41   
 [5] permute_0.9-7      brms_2.14.4        Rcpp_1.0.6         lme4_1.1-26       
 [9] Matrix_1.2-18      emmeans_1.7.2      cowplot_1.1.1      forcats_0.5.1     
[13] stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4        readr_2.1.2       
[17] tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.3      tidyverse_1.3.1   

loaded via a namespace (and not attached):
  [1] readxl_1.3.1         backports_1.2.1      plyr_1.8.6           igraph_1.3.1        
  [5] splines_4.0.3        crosstalk_1.1.1      TH.data_1.1-1        rstantools_2.1.1    
  [9] inline_0.3.17        digest_0.6.27        htmltools_0.5.2      rsconnect_0.8.16    
 [13] fansi_0.4.2          magrittr_2.0.1       cluster_2.1.0        openxlsx_4.2.3      
 [17] tzdb_0.3.0           modelr_0.1.8         RcppParallel_5.0.2   matrixStats_0.58.0  
 [21] xts_0.12.1           sandwich_3.0-1       prettyunits_1.1.1    colorspace_2.0-0    
 [25] rvest_1.0.2          haven_2.3.1          callr_3.7.0          crayon_1.4.1        
 [29] jsonlite_1.7.2       survival_3.2-7       zoo_1.8-8            glue_1.6.2          
 [33] gtable_0.3.0         V8_3.4.0             car_3.0-10           pkgbuild_1.2.0      
 [37] rstan_2.21.1         abind_1.4-5          scales_1.1.1         mvtnorm_1.1-1       
 [41] DBI_1.1.1            miniUI_0.1.1.1       xtable_1.8-4         foreign_0.8-80      
 [45] stats4_4.0.3         StanHeaders_2.21.0-7 DT_0.17              htmlwidgets_1.5.3   
 [49] httr_1.4.2           threejs_0.3.3        DiagrammeR_1.0.9     RColorBrewer_1.1-2  
 [53] ellipsis_0.3.2       pkgconfig_2.0.3      loo_2.4.1            dbplyr_2.1.1        
 [57] utf8_1.1.4           tidyselect_1.1.2     rlang_1.0.2          reshape2_1.4.4      
 [61] later_1.1.0.1        munsell_0.5.0        cellranger_1.1.0     tools_4.0.3         
 [65] visNetwork_2.1.0     cli_3.3.0            generics_0.1.0       broom_0.8.0         
 [69] ggridges_0.5.3       fastmap_1.1.0        yaml_2.2.1           processx_3.5.2      
 [73] fs_1.5.0             zip_2.1.1            nlme_3.1-149         mime_0.9            
 [77] projpred_2.0.2       xml2_1.3.2           compiler_4.0.3       bayesplot_1.8.0     
 [81] shinythemes_1.2.0    rstudioapi_0.13      gamm4_0.2-6          curl_4.3            
 [85] reprex_2.0.1         statmod_1.4.35       stringi_1.5.3        ps_1.5.0            
 [89] Brobdingnag_1.2-6    nloptr_1.2.2.2       markdown_1.1         shinyjs_2.0.0       
 [93] vctrs_0.4.1          pillar_1.7.0         lifecycle_1.0.1      bridgesampling_1.0-0
 [97] estimability_1.3     data.table_1.14.0    httpuv_1.5.5         R6_2.5.0            
[101] promises_1.1.1       gridExtra_2.3        rio_0.5.26           codetools_0.2-16    
[105] boot_1.3-25          colourpicker_1.1.0   MASS_7.3-53          gtools_3.8.2        
[109] assertthat_0.2.1     withr_2.4.3          shinystan_2.5.0      multcomp_1.4-19     
[113] mgcv_1.8-33          parallel_4.0.3       hms_1.0.0            grid_4.0.3          
[117] coda_0.19-4          minqa_1.2.4          carData_3.0-4        shiny_1.6.0         
[121] lubridate_1.8.0      base64enc_0.1-3      dygraphs_1.1.1.6 


# License Information
To the extent possible under law, *Jinlin Chen* has waived all copyright to use the included code as a template for related analysis. Please cite the archived or published version of the paper if adapting the code for your research. Please contact *Jinlin Chen* for potential collaboration if you want to use any of the datasets involved in this project. 

Copyright (c) 2023 JINLIN CHEN
