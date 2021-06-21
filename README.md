# Application of an in vitro assay to identify chemicals that increase estradiol and progesterone synthesis and are potential breast cancer risk factors
## Authors
Bethsaida Cardona, Ruthann A. Rudel<br/>
Silent Spring Institute, Newton, MA<br/>
Corresponding author: Ruthann Rudel, Silent Spring Institute, Newton, MA 02460, USA.<br/>

## Sources version
May 2021<br/>
R version 3.6.3 (2020-02-29)<br/>

## libraries 
- tidyverse
- tidyr
- janitor
- readxl
- lubridate
- stringi
- ggforce
- ggrepel
- ggmosaic
- ggstatsplot
- ggpubr
- ggplot2
- ggpattern (#remotes::install_github("coolbutuseless/ggpattern"))


## Input files
Input files should be in the input folder in the "../" directory 
(Script 1)
- Supp7_Global_H295R_ANOVA_OECD_strings_filtered_output_2017-08-09.txt: binary code with the significance associated with each chemical dose (Haggard 2018 supplemental https://academic.oup.com/toxsci/article/162/2/509/4690670)
- Supp4_OECD_GLOBAL_ANOVA_output_pValues2017-08-09.txt: the chemical concentrations and their associated p-values for all hormones (Haggard 2018 supplemental https://academic.oup.com/toxsci/article/162/2/509/4690670)
- Supp3_H295R_master_table_2017-08-08.txt: hormone data concentrations for each chemical dose (Haggard 2018 supplemental https://academic.oup.com/toxsci/article/162/2/509/4690670)
- Supp11_Global_OECD_Mahalanobis_distances_2017-08-09.txt:  Supplemental file from Derik's paper containing adjusted max Mahalanobis distance (Haggard 2018 supplemental https://academic.oup.com/toxsci/article/162/2/509/4690670)
- Supp2_Global_H295R_ANOVA_OECD_directionality_filtered_output_2018-08-09.txt: the overall direction of chemical effect on hormone concentrations in H295R cells (Haggard 2019 supplemental https://www.sciencedirect.com/science/article/abs/pii/S0273230019302740#appsec1)
- ac50_Matrix_190226.csv: File from Toxcast database containing the ac50s and chemical code. Also used to find the number of active assays and number of total assays. NA means not tested, "1000000" for not active and numerical values less than 1000000 indicating the AC50 in microMolar; downloaded from the EPA Comptox chemical dashboard: https://comptox.epa.gov/dashboard/downloads (last update 2-26-2019)
- Chemical_Summary_190226.csv: File to match chemical code and CAS/chemical name from the EPA chemical dashboard https://comptox.epa.gov/dashboard/downloads (last update 2-26-2019)
- SupTable-all.chem.preds-2018-11-28.txt: Supplemental file from Ring 2019 paper containing HT exposure predictions (https://pubs.acs.org/doi/10.1021/acs.est.8b04056)
(Script 2 - exposure)
- cpdat_source_substances.xlsx: CPDAT- links CASN to the associated CPDat record id (last updated 2018, data from https://comptox.epa.gov/dashboard/downloads)
- chemical_and_product_categories.xlsx: CPDAT- product type as well as the CPDat record id (last updated 2018, data from https://comptox.epa.gov/dashboard/downloads)
- products.xlsx: CPDAT- list of products with no ingredients (last updated 2018, data from https://comptox.epa.gov/dashboard/downloads)
- ingredients.xlsx: CPDAT- list of ingredients with the CPDat record id but no CASN (last updated 2018, data from https://comptox.epa.gov/dashboard/downloads)
- SupTable-all.chem.preds-2018-11-28.txt: Supplemental file from Ring 2019 paper containing HT exposure predictions (https://pubs.acs.org/doi/10.1021/acs.est.8b04056)
(Script 2 - invivo)
- toxval_pod_summary_human health_2019-08-20.csv: ToxValDB data Downloaded ftp://newftp.epa.gov/COMPTOX/STAFF/rjudson/datasets/ToxValDB/ (last updated 08-20-19)
- toxval_cancer_summary_2019-08-20.xlsx: ToxValDB data Downloaded ftp://newftp.epa.gov/COMPTOX/STAFF/rjudson/datasets/ToxValDB/ (last updated 08-20-19)
- dev list.csv: mammary developmental toxicants from Rudel 2011 (data compiled by authors using paper results https://pubmed.ncbi.nlm.nih.gov/21697028/
- MC list.csv: mammary carcinogens from Rudel 2007 (data compiled by authors using paper results https://pubmed.ncbi.nlm.nih.gov/17503434/
- list_nhaneschemicals-2020-01-17-15-52-55.csv: NHANES chemicals (last updated 01-17-20; data download Comptox Dashboard https://comptox.epa.gov/dashboard/chemical_lists/NHANES2019)
- p65chemicalslist2021p.xlsx: Prop65 chemicals (last updated 03-19-21; data download https://oehha.ca.gov/proposition-65)
(Script 3)
- modl_ac10_Matrix_190226.csv: ac10 and ac50 values for the H295R assay (and other Toxcast assays) downloaded from Comptox website: https://comptox.epa.gov/dashboard/downloads) (last update 2-26-2019)
- hitc_Matrix_190226.csv:  ac10 and ac50 values for the H295R assay (and other Toxcast assays) downloaded from Comptox website: https://comptox.epa.gov/dashboard/downloads (last update 2-26-2019)
- pesticide product status2.csv: dataset with number of pesticide products that are approved by EPA and contain a certain chemical ingredient (casn shown), also contrains number of products not currently approved but which have been approved in the past data compiled by authors using data files downloaded from https://www.epa.gov/ingredients-used-pesticide-products/ppis-download-product-information-data (last update 5-5-20)
- products.txt: file listing FDA approved drug products. downloaded from https://www.fda.gov/drugs/drug-approvals-and-databases/approved-drug-products-therapeutic-equivalence-evaluations-orange-book last update 6-30-2020
- supp_kfv168_toxsci-15-0258-File002_adj.xlsx: ER AUC data from Judson 2015 supplemental (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4635633/) some CAS were converted to dates in original file. These were manually adjusted in new column 
- pesticides_weffects.csv: list of pesticides with observed mammary effects, from Cardona and Rudel 2020 (data compiled by authors using paper results https://pubmed.ncbi.nlm.nih.gov/32645345/)

## Scripts
Scripts need to be run sequentially following this order
- R Script 1: prep data. Code identifying E2/P4-up data
- R Script 2: exposures. Code focusing on exposure sources
- R Script 2: in vivo. Code focusing only on invivo effects of E2/P4-up chemicals
- R Script 3: E2/P4. Code puts together all tables relevant to E2/P4 up chemicals 


