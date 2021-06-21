
# Application of an in vitro assay to identify chemicals that increase estradiol 
# and progesterone synthesis and are potential breast cancer risk factors
# Bethsaida Cardona1, Ruthann A. Rudel1
# 1Silent Spring Institute, Newton, MA 
# Corresponding author: Ruthann Rudel, Silent Spring Institute, Newton, MA 02460, USA. 

# May 2021
# R version 3.6.3 (2020-02-29)

#R Script 1: prep data. Code identifying E2/P4-up data


# Set working directory where the R script lives to workingdir
workingdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workingdir)

# Set the working directory one directory back of workingdir
setwd("..") 


#load packages
library(tidyverse)
library(tidyr)
library(janitor)
library(readxl)
library(lubridate)
library(stringi)


#effectively disables scientific notation 
options(scipen = 999)


my_signif = function(x, digits) floor(x) + signif(x %% 1, digits)


#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read in data  -----
#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#binary code with the significance associated with each chemical dose (Haggard 2018 supplemental)
hits <- read_tsv("./input/Supp7_Global_H295R_ANOVA_OECD_strings_filtered_output_2017-08-09.txt")

#the chemical concentrations and their associated p-values for all hormones (Haggard 2018 supplemental)
concentrations <- read_tsv("./input/Supp4_OECD_GLOBAL_ANOVA_output_pValues2017-08-09.txt")

#hormone data concentrations for each chemical dose (Haggard 2018 supplemental)
hormd <- read_tsv("./input/Supp3_H295R_master_table_2017-08-08.txt")

# Supplemental file from Derik's paper containing adjusted max Mahalanobis 
# distance (Haggard 2018 supplemental)
mmd <- read_tsv("./input/Supp11_Global_OECD_Mahalanobis_distances_2017-08-09.txt")

#the overall direction of chemical effect on hormone concentrations in H295R cells (Haggard 2019 supplemental)
direction <- read_tsv("./input/Supp2_Global_H295R_ANOVA_OECD_directionality_filtered_output_2018-08-09.txt")


#File from Toxcast database containing the ac50s and chemical code. Also used to find
#the number of active assays and number of total assays. NA means not tested, "1000000"
#for not active and numerical values less than 1000000 indicating the AC50 in microMolar
#downloaded from Comptox website:  (last update 2-26-2019)
ac50_and_assaycount <- read_csv("./input/ac50_Matrix_190226.csv")

#File to match chemical code and CAS/chemical name
CASN <- read_csv("./input/Chemical_Summary_190226.csv")


# Supplemental file from Ring 2019 paper containing HT exposure predictions
htp_pred <- read_csv("./input/SupTable-all.chem.preds-2018-11-28.txt")

# Warning about parsing failures for htp_pred. Investigate with problems()
problems(htp_pred) %>% pull(col) %>% unique

# parsing failure only applies to "Pred.SHEDS.Indirect", "Pred.Food.Contact",
# "Pred.REDS", "Pred.FINE", "Pred.RAIDAR.ICE", "Pred.SHEDS.Direct"
# read_tsv is incorrectly predicting that these columns contain logical
# variables when in fact they are numeric. See help for read_tsv() for more
# explanation. Currently we don't care about these columns so it is okay
# if they are parsed incorrectly. If need to use these columns, see read_tsv
# for help about how to specify column types rather than let read_tsv guess



#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Potency Measures  -----
#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1) Direction of chemical effect (Increase v. Decrease) ------------------------

# Want to gather this data to identify chems that increase estradiol and/or 
# progesterone production

# This info is summarized in the direction file.

# a. Transpose 

#pulled out the steroid column and transposed the data set 
direction_t <- 
  direction %>% 
  select(-steroid) %>%  
  t() %>% 
  as.data.frame() 

#gave a name to the rows by pulling the name of the variables in the steroid column 
names(direction_t) <- direction %>% pull(steroid)

#gave a name (data_chnm_plate) to the first columm 
direction_t <-
  direction_t %>%
  rownames_to_column %>%
  rename(date_chnm_plate = rowname)

#b.For cases where chemicals were run multiple times, pick first plate for each 
#chemical for this exercise. Will consider reproducibility across plates later. 
# To do this, need to split date/chem/plate into separate columns and then pick
# first plate
direction_t2 <-
  direction_t %>%
  separate(date_chnm_plate, into = c("date", "chnm", "plate"), sep = "_", remove = FALSE) %>%
  # remove the "Plate" string from plate to only leave the number
  mutate(plate = gsub("Plate", "", plate) %>% as.numeric()) %>%
  group_by(chnm) %>%
  # top_n: If n is negative, selects the bottom n rows (where "bottom" = lowest
  # number rank)
  top_n(wt = plate, n = -1) %>%
  ungroup()

# Check to make sure we ended up with one row for each chem
stopifnot(nrow(direction_t2) == length(unique(direction_t2$chnm)))


# c.Make direction long -- that is collapse individual steroid columns and direction
# info into one column with steroid name, one column with direction info. 
# This will allow us to match up with hits file
direction_long <-
  direction_t2 %>%
  gather(steroid, direction, `OH-Pregnenolone`:`Estradiol`)


#2) Hit (Significance) associated with each chemical dose ----------------------


# a. Separate out the hit info at each concentration -- currently all contained
# in single cell, we want one row with the hit info for 
# each chemical concentration and steroid

hits_long <-
  hits %>%
  gather(steroid, hit, `OH-Pregnenolone`:`Estradiol`)

hits_long2 <-
  hits_long %>%
  mutate(hit = gsub("c\\(|\\)", "", hit)) %>%
  separate(hit, into = c("1", "2", "3", "4", "5", "6"), sep = ", ") %>%
  gather(dose_level, hit, `1`:`6`) %>%
  mutate(dose_level = as.numeric(dose_level))

# b. Merge direction data in, using left_join to keep only plates that appear 
# in the filtered direction data (to avoid having multiple plates for chemicals that 
# were run more than once). 
hits_and_direction <-
  left_join(direction_long %>% select(date_chnm_plate, steroid, direction, plate), 
            hits_long2)


# c. Need to link up concentrations file to hits file. To do this, 
# need to add dose level indicator (1 for lowest, 6 for highest) 

concentrations2 <- 
  concentrations %>%
  select(date_chnm_plate, conc) %>%
  group_by(date_chnm_plate) %>%
  arrange(date_chnm_plate, conc) %>% # I don't think this arrange step is 
  # necessary, rank would do the re-arranging. Just in case. 
  mutate(dose_level = rank(conc)) %>%
  ungroup()

# d. Bring concentrations and hits/direction together.

hits_dir_conc <-
  left_join(hits_and_direction, concentrations2)

stopifnot(nrow(hits_dir_conc) == nrow(hits_and_direction))


#3) Fold Change Calculation --------------------------------------------------- 

# Need to compile info about the hormone concentration measured at each chemical
#  dose to calculate fold change  

# Info summarized in the hormd file 

hormd_tidy <- hormd %>%
  separate(apid, c("date_text","trash","plate"), "\\.") %>%
  mutate(date_text_fixed = case_when(date_text=="04112017" ~ "11Apr2011",
                                     date_text=="20130320" ~ "20Mar2013",
                                     date_text=="20130321" ~ "21Mar2013",
                                     TRUE ~ date_text)) %>%
  mutate(date = dmy(date_text_fixed)) %>%
  mutate(flag = case_when(ceetox_raw=="ND" ~ 1,    # make flag =1 for non-detects
                        ceetox_raw=="NQ" ~ 1,
                        TRUE ~ 0)) %>%
  mutate(plate = as.numeric(plate)) %>% 
  select(spid, chnm, casn, plate, steroid, conc_index, conc, rowi, coli, uM, date, flag)



# We have two observations at each concentration level - these are the replicates - 
# just going to take the mean of these two
hormd2 <- hormd_tidy %>%
  group_by(steroid, casn, chnm, date, plate, conc_index, conc) %>%
  summarize(uM_bar=mean(uM), flag_bar=sum(flag))



# We also need to get the DMSO (baseline) values for each steroid
DMSO <- hormd_tidy %>%
  filter(spid=="DMSO") %>%
  select(date, plate, casn, rowi, coli, steroid, uM, flag) %>%
  rename(DMSO_uM=uM) %>%
  rename(DMSO_flag=flag) %>% 
  #For a given hormone-date-plate, need to average across the duplicate samples
  group_by(steroid, date, plate) %>%
  summarize(DMSO_uM_bar=mean(DMSO_uM), DMSO_flag_bar=sum(DMSO_flag))

# Merge in the DMSO's
hormd3 <- hormd2 %>%
  left_join(DMSO, by=c("steroid","date", "plate")) %>%
  # Now we need to compute the fold change of the hormone (comparing DMSO to treatments) - 
  # up to 2 decimal points
  mutate(foldchange = (uM_bar/DMSO_uM_bar),
         foldchange = my_signif(foldchange, 2), 
         foldchange = as.double(foldchange)) %>% 
  #for merging with hits_dir_conc dataset 
  group_by(casn, steroid, plate, date) %>%
  mutate(dose_level = rank(conc))


stopifnot(sum(is.na(hormd3$chnm)) ==  sum(is.na(hormd2$chnm))) #rows with chnm is NA are colums with info about DMSO 



#4) Combine fold change info with that of hits and direction ------------------------------------

# also want to indicate whether fold change went up or down or had no change (fc_dir) and 
# add the associated direction to "hits" (hitd)


data_long <- hits_dir_conc %>% 
  # merging by rank of conc (dose_level) rather than conc because these vary by decimals 
  left_join(hormd3, by = c("casn", "steroid", "plate", "dose_level")) %>% 
  # chmn.y varies from chmn.x due to how they are written but they refer to the same chemical 
  select(-c(date, conc_index, flag_bar, DMSO_flag_bar, conc.y, chnm.y )) %>% 
  #want to indicate direction of fold change  
  mutate(fc_dir = case_when(foldchange > 1 ~ 1, 
                            foldchange < 1 ~ -1, 
                            foldchange == 1 ~ 0,
                            TRUE ~ foldchange)) %>% 
  mutate(hit = as.numeric(hit)) %>% 
  #currently the hit tells us whether a dose signifcantly changed the hormone 
  #level but does not tell us in what direction, we can look at the fold change 
  #to answer this 
  mutate(hitd = case_when(hit != 0 ~ fc_dir, 
                          TRUE ~ hit)) 

stopifnot(nrow(hits_dir_conc) == nrow(data_long))



#5) Calculate the max fold change and LEC (Lowest effective concentration)

#to make it easier, we want to remove any data associated with those chemical 
#doses which are not significant following OECD logic we only want data from 
#chemical doses which had consecutive significant increases or which were 
#significant at the highest chemical dose administered

data_wide <- data_long %>% 
  select(-c(hit, uM_bar, DMSO_uM_bar, fc_dir)) %>% 
  #only interested in finding the LEC and fold change for chemicals that 
  #increased the hormone concentration 
  filter(direction == 1) %>% 
  rename(fc = foldchange, conc = conc.x, chnm = chnm.x) %>% 
  #each chemical dose level has a different associated chemical concentration, 
  #fold change, and direction of hit
  pivot_wider(names_from = dose_level, values_from = c(hitd, conc, fc)) %>% 
  #are not interested in data were significance of fold change was not 
  #consecutive or at highest level 
  #remove hits we are not interested, check that chems all follow OECD logic - they do 
  mutate(hitd_1.a = case_when(hitd_1 == 1 & hitd_2 == 1 ~ hitd_1, 
                              TRUE ~ NA_real_),
         hitd_2.a = case_when(hitd_2 == 1 & (hitd_1 == 1 | hitd_3 %in% c(NA_real_, 1))  ~ hitd_2, 
                              TRUE ~ NA_real_), 
         hitd_3.a = case_when(hitd_3 == 1 & (hitd_2 == 1 | hitd_4 %in% c(NA_real_, 1)) ~ hitd_3, 
                              TRUE ~ NA_real_), 
         hitd_4.a = case_when(hitd_4 == 1 & (hitd_3 == 1 | hitd_5 %in% c(NA_real_, 1)) ~ hitd_4, 
                              TRUE ~ NA_real_), 
         hitd_5.a = case_when(hitd_5 == 1 & (hitd_4 == 1 | hitd_6 %in% c(NA_real_, 1)) ~ hitd_5, 
                              TRUE ~ NA_real_), 
         hitd_6.a = case_when(hitd_6 == 1  ~ hitd_6, 
                              TRUE ~ NA_real_)) %>% 
  #remove concentration data if hit is not 1 
  mutate(conc_1 = case_when(hitd_1.a == 1 ~ conc_1,
                            TRUE ~ NA_real_),
         conc_2 = case_when(hitd_2.a == 1 ~ conc_2,
                            TRUE ~ NA_real_),
         conc_3 = case_when(hitd_3.a == 1 ~ conc_3,
                            TRUE ~ NA_real_),
         conc_4 = case_when(hitd_4.a == 1 ~ conc_4,
                            TRUE ~ NA_real_),
         conc_5 = case_when(hitd_5.a == 1 ~ conc_5,
                            TRUE ~ NA_real_),
         conc_6 = case_when(hitd_6.a == 1 ~ conc_6,
                            TRUE ~ NA_real_)) %>% 
  #remove fc data if hit is not 1 
  mutate(fc_1 = case_when(hitd_1.a == 1 ~ fc_1,
                          TRUE ~ NA_real_),
         fc_2 = case_when(hitd_2.a == 1 ~ fc_2,
                          TRUE ~ NA_real_),
         fc_3 = case_when(hitd_3.a == 1 ~ fc_3,
                          TRUE ~ NA_real_),
         fc_4 = case_when(hitd_4.a == 1 ~ fc_4,
                          TRUE ~ NA_real_),
         fc_5 = case_when(hitd_5.a == 1 ~ fc_5,
                          TRUE ~ NA_real_),
         fc_6 = case_when(hitd_6.a == 1 ~ fc_6,
                          TRUE ~ NA_real_)) 


#identify lec and MFC
#add indicator for fc > 1.5
maxfc_lec <-  data_wide %>% 
  rowwise  %>% 
  mutate(MFC = max(c(fc_1, fc_2, fc_3, fc_4, fc_5, fc_6), na.rm = T)) %>% 
  mutate(LEC = min(c(conc_1, conc_2, conc_3, conc_4, conc_5, conc_6), na.rm = T)) %>% 
  ungroup() %>% 
  select(date_chnm_plate:casn, MFC:LEC) 



#5) Identify the max and min chem dose and bring back the hits vector ----------

#to get column with the max dose tested: bring back in maximum dose tested 
#from the "concentration" dataset

maxd_hits_maxfc_lec <- maxfc_lec %>% 
  #max dose
  left_join(., 
            concentrations2) %>%
  group_by(date_chnm_plate, steroid) %>%
  mutate(max_dose = max(dose_level)) %>%
  filter(dose_level %in% max_dose) %>%
  rename(max_tested_conc = conc) %>%
  select(-max_dose, -dose_level) %>%
  ungroup() %>%
  #min dose
  left_join(., 
            concentrations2) %>%
  group_by(date_chnm_plate, steroid) %>%
  mutate(min_dose = min(dose_level)) %>%
  filter(dose_level %in% min_dose) %>%
  rename(min_tested_conc = conc) %>%
  select(-min_dose, -dose_level) %>%
  ungroup() %>%
  # Bring back the hits vector
  left_join(., hits_long) %>%
  rename(hits = hit) 



#6) Identify chemicals that increase E2 and/or P

E2_P_chems <- maxd_hits_maxfc_lec %>%
  #filter for chemicals of interest
  filter(steroid %in% c("Progesterone", "Estradiol") & direction == 1) %>%
  #direction should be 1 for all the chemicals on here so we can remove it for now. 
  #Will read it in later on... 
  select(-direction) %>%
  group_by(date_chnm_plate) %>%
  pivot_wider(names_from = steroid, values_from = c(LEC, MFC, hits)) %>% 
  #NA under Estradiol or Progesterone indicates that for that 
  #chemical, the 'direction' for that hormone was not 1 
  ungroup() %>%
  #adding in the direction of progesterone and estradiol for each chemical 
  left_join(., direction_long %>%
              select(-c(date_chnm_plate, date, plate)) %>%
              filter(steroid %in% c("Progesterone", "Estradiol")) %>%
              spread(key = steroid, value = direction)) %>% 
  filter(Progesterone == 1 | Estradiol == 1)




#7) Adjusted maxMd ------------------------------------------------------------------------------------


# Adding the adjusted max Mahalanobis distance from Haggard paper


chems_with_maxMd  <-
  left_join(E2_P_chems, mmd %>% select(date_chnm_plate:casn, adj_maxmMD)) %>% 
  #round to two decimal pleces 
  mutate(adj_maxmMD = sprintf("%0.2f", adj_maxmMD))


# 8) AC50 ---------------------------------------------------------------------------------------------------

#creating new databases to manipulate; contains columns of interest 

# selected only the H295R assays which increase progesterone and estradiol 
# this dataset contains the code name but does not link to the casn 
ac50_data <- ac50_and_assaycount %>%
  select(1, contains ("H295R")) %>%
  select(1, contains("PROG"), contains("ESTRADIOL")) %>%
  select(1, contains("ESTRADIOL_up"), contains("_PROG_up")) %>%
  rename(code = X1)

#contains the code name and the casn 
ac50_key <- CASN %>%
  select(chnm:code)

#making sure all of the codes and casn are unique and that each casn/code 
#combination is unique if they are we will match based on code 
ac50_key$casn[duplicated(ac50_key$casn)]
ac50_key$code[duplicated(ac50_key$code)]  
nrow(distinct(ac50_key, casn, code))

#combines the casn from ac50_key and the code name from ac50_data
ac50 <- ac50_data %>%
  left_join(ac50_key, by = "code") %>%
  rename("AC50_Estradiol" = "CEETOX_H295R_ESTRADIOL_up",
         "AC50_Progesterone" = "CEETOX_H295R_PROG_up") %>%
  select(-code) #no longer need the code  

stopifnot(nrow(ac50) == nrow(ac50_data))


#combining the ac50 data with the dataset of interest
chems_with_ac50 <- chems_with_maxMd  %>%
  left_join(ac50, by = c("casn")) %>% 
  #only want the AC50 up to 2 decimal pleces
  mutate(AC50_Progesterone = sprintf("%0.2f", AC50_Progesterone), 
         AC50_Estradiol = sprintf("%0.2f", AC50_Estradiol))


#checks post-join()
stopifnot(nrow(chems_with_maxMd) == nrow(chems_with_ac50))



#############################################################################################################
## Creating tables for E2 and P up 
#############################################################################################################

E2_P_data <- chems_with_ac50 %>% 
  rename("chnm" = "chnm.x") %>% 
  select(-chnm.y) 


# #############################################################################################################
# ## Exposure Predictions 
# #############################################################################################################


# HT Exposure Predictions  ------------------------------------------------ 


#Adding high throughput exposure predictions from Ring 2019

# Note that there are some chemicals that don't have high throughput 
# exposure pedictions - 9 chemicals 
(cas_no_htp <- setdiff(E2_P_data$casn, htp_pred$CAS %>% unique()))

chems_no_htp <-
  E2_P_data %>% 
  filter(casn %in% cas_no_htp) %>% 
  select(chnm) %>% 
  unique

# clean up htp_pred chemical names -- just a few small changes required
# to match up.
setdiff(unique(E2_P_data$chnm), unique(htp_pred$Substance_Name))

htp_pred <- 
  htp_pred %>%
  # htp_pred has double quotes around some chemical names...remove
  mutate(chnm = gsub('"', "", Substance_Name))

setdiff(unique(E2_P_data$chnm), unique(htp_pred$chnm))

htp_pred <-
  htp_pred %>%
  mutate(chnm = 
           # also drop "_" in PharmaGSID_48509 and similar
           case_when(grepl("PharmaGSID_", chnm) ~  gsub("_", "", chnm), 
                     # add a "1-" in front of "Bromo-3-chloro-5,5,-dimethylhydantoin"
                     grepl("Bromo-3-chloro-5,5-dimethylhydantoin", chnm) ~ "1-Bromo-3-chloro-5,5-dimethylhydantoin",
                     # lastly, drop "h" from Terbuthylazine (I matched the CAS 
                     # numbers so am sure it is same substance)
                     grepl("Terbuthylazine", chnm) ~ "Terbutylazine",
                     TRUE ~ chnm))

setdiff(unique(E2_P_data$chnm), unique(htp_pred$chnm))


exposure_potency <-
  left_join(E2_P_data, 
            htp_pred %>% 
              select(CAS, chnm, seem3:Pathway), 
            by = c("casn" = "CAS", "chnm" = "chnm")) %>%
  #want to trunctuate to significant digits for seem: 
  mutate(seem3.u95 = my_signif(seem3.u95, 2), 
         seem3 = my_signif(seem3, 2)) %>% 
  ##for cases when seem3.up95 is NA we want to replece the NA with seem3 which is the median; none of the chemicals for which 
  ##the seem3.u95 was NA had a seem3 < or = to 0.01
  mutate(seem = case_when(is.na(seem3.u95) ~ sub("$", " (m)", seem3), 
                          TRUE ~ sub("$", " (u)", seem3.u95))) %>%#create indicator for whether value used was the median or upper 95
  mutate(casn = paste0("'", casn))


# check that only chems in 'chems_no_htp' have NAs for seem

stopifnot(exposure_potency %>% filter(is.na(seem3)) %>% pull(chnm) %>% unique() 
          %in% chems_no_htp$chnm)


#output for use with other R scripts
write_csv(exposure_potency, "./output/exposure_potency.csv") 


