CAMDA 2017: Neuroblastoma challenge

For a summary, see `doc/CAMDA_presentation`. Session info from the analysis is below.

All output is saved to the output/ directory. To recreate the output, delete the output/ folder and run `src/main.R`. All models will be refitted (using the processed data in the data/ folder), and the output figures and tables will be created.
 
The raw datasets can be obtained from the GEO database with accession number GSE49711. Processed data is available and stored in data/processed/

Project organization:
  bin/ 
    | - currently empty (compiled code/external scripts go here).
  data/ 
    | raw/ 
    | - | ... [Not uploaded due to large file sizes.]
    | processed/
    | - | ... [The raw data is processed by `src/rna_seq_data_scripts/run_processing.R`, and the results are saved here.]
    | patient_data : extracted patient data from raw data. See `src/patient_data_scripts/patient_data.R`.
    | gene_cards : list of genes related to neuroblastoma, obtained from the GeneCards database.
  docs/ 
    | CAMDA_presentation : Powerpoint presentation given at the 2017 CAMDA conference.
  logs/
    | ... [When scripts are run, output is directed to text files in this directory.]
  neuroblastoma.Rproj 
  output/
    | ensemble/
    | - | ... [Models fit via ensemble method are saved here.]
    | figures/
    | - | km/
    | - | - | ... [Kaplan-Meier curves for each model.]
    | - | measures/
    | - | - | ... [Plots for performance measures for each model.]
    | results/
    | - | ... [R objects containing list of all fitted models.]
    | tables/
    | - | ... [Table summarizing Kaplan-Meier estimates for 2 year survival for each model, formatted using LaTex.]
  README.txt
  src/
    | ensemble.R
    | fit.R
    | imputation.R
    | init_required_packages_and_files.R
    | loss_functions.R
    | main.R
    | miscellaneous.R
    | models/
    | - | ... [Wrappers for each model used: elastic net, lasso, PLS, linear regression, and sparse PLS.]
    | output_scripts/
    | - | ... [scripts for generating output figures and tables.]
    | patient_data_scripts/
    | - | ... [scripts for extracting patient clinical data.]
    | plots.R
    | rna-seq_data_scripts/
    | - | ... [scripts for extracting patient RNA-seq data.]

    
    
sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.4

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] survcomp_1.28.5   prodlim_1.6.1     gdata_2.18.0      Hmisc_4.1-1      
 [5] Formula_1.2-2     lattice_0.20-35   km.ci_0.5-2       boot_1.3-20      
 [9] survminer_0.4.2   ggpubr_0.1.6      magrittr_1.5      ggplot2_2.2.1    
[13] RankAggreg_0.6.4  doParallel_1.0.11 iterators_1.0.9   glmnet_2.0-16    
[17] foreach_1.4.4     Matrix_1.2-14     spls_2.2-2        truncnorm_1.0-8  
[21] maxLik_1.3-4      miscTools_0.6-22  survival_2.41-3   dplyr_0.7.4      

loaded via a namespace (and not attached):
 [1] tidyr_0.8.0          splines_3.4.4        gtools_3.5.0        
 [4] microbenchmark_1.4-4 rmeta_3.0            assertthat_0.2.0    
 [7] latticeExtra_0.6-28  bootstrap_2017.2     yaml_2.1.18         
[10] survivalROC_1.0.3    pillar_1.2.1         backports_1.1.2     
[13] glue_1.2.0           digest_0.6.15        RColorBrewer_1.1-2  
[16] checkmate_1.8.5      colorspace_1.3-2     sandwich_2.4-0      
[19] cmprsk_2.2-7         htmltools_0.3.6      plyr_1.8.4          
[22] psych_1.8.3.3        pkgconfig_2.0.1      broom_0.4.4         
[25] purrr_0.2.4          xtable_1.8-2         scales_0.5.0        
[28] lava_1.6.1           KMsurv_0.1-5         tibble_1.4.2        
[31] htmlTable_1.11.2     nnet_7.3-12          lazyeval_0.2.1      
[34] mnormt_1.5-5         nlme_3.1-137         foreign_0.8-69      
[37] SuppDists_1.1-9.4    tools_3.4.4          data.table_1.10.4-3 
[40] stringr_1.3.0        munsell_0.4.3        cluster_2.0.7-1     
[43] bindrcpp_0.2.2       pls_2.6-0            compiler_3.4.4      
[46] rlang_0.2.0          grid_3.4.4           rstudioapi_0.7      
[49] htmlwidgets_1.0      base64enc_0.1-3      gtable_0.2.0        
[52] codetools_0.2-15     reshape2_1.4.3       R6_2.2.2            
[55] gridExtra_2.3        zoo_1.8-1            knitr_1.20          
[58] survMisc_0.5.4       bindr_0.1.1          KernSmooth_2.23-15  
[61] stringi_1.1.7        Rcpp_0.12.16         rpart_4.1-13        
[64] acepack_1.4.1    

