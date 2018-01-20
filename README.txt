CAMDA 2017: Neuroblastoma challenge

For a summary, see `doc/CAMDA_presentation`. 

All output is saved to the output/ directory. To recreate the output, simply delete the output/ folder and run `src/main.R`. All models will be refitted (using the processed data in the data/ folder), and the output figures and tables will be created.

Note that the current output is from fitting the models without using the GeneCards dataset to subset genes. Not all figures and tables have been generated. 
 
The raw datasets can be obtained from the GEO database with accession number GSE49711.



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
