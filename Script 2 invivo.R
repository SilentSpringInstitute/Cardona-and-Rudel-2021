
# Application of an in vitro assay to identify chemicals that increase estradiol 
# and progesterone synthesis and are potential breast cancer risk factors
# Bethsaida Cardona1, Ruthann A. Rudel1
# 1Silent Spring Institute, Newton, MA 
# Corresponding author: Ruthann Rudel, Silent Spring Institute, Newton, MA 02460, USA. 

# May 2021
# R version 3.6.3 (2020-02-29)

#R Script 2: in vivo. Code focusing only on invivo effects of E2/P4-up data
#-Needs Script 1 to run. Can run concurrently with Script 2: exposures


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
library(ggrepel)
library(lubridate)
library(ggforce)


#effectively disables scientific notation
options(scipen = 999)


#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read in data  -----
#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#data with potency measures and exposure predictions from Script 1: data prep
E2_P_data <- read_csv("./output/exposure_potency.csv") %>%
  mutate(casn = gsub("'", "", casn))

#ToxvalDB for invivo tests 
#Downloaded ftp://newftp.epa.gov/COMPTOX/STAFF/rjudson/datasets/ToxValDB/ (last updated 08-20-19)
invivo_ToxVal <- read_csv("./input/toxval_pod_summary_human health_2019-08-20.csv")
cancer_ToxVal <- read_excel("./input/toxval_cancer_summary_2019-08-20.xlsx")

#mammary developmental toxicants from Rudel 2011
dev_tox<- read_csv("./input/dev list.csv")

#mammary carcinogens from Rudel 2007 
mamm_car <- read_csv("./input/MC list.csv")

#NHANES chemicals (last updated 01-17-20; data download Comptox Dashboard)
nhanes <- read_csv("./input/list_nhaneschemicals-2020-01-17-15-52-55.csv")

#Prop65 chemicals (last updated 03-19-21; data download https://oehha.ca.gov/proposition-65)
prop65 <- read_excel("./input/p65chemicalslist2021p (1).xlsx", skip = 10)

#####################################################################################################################
##In-vivo
#####################################################################################################################

#using ToxValDb,identify chemicals that are likely/unlikely/not enough data to be developmental toxicants 
#based on toxicity values reported in invivo studies. label each study as a dev/repro study conducted 
#in mammals, as reporting an effect or a no effect level and if value reported was above or under 
#100 mg/kg-day 

filtered_invivo_ToxVal <- invivo_ToxVal %>%
  #subset to chemicals of interest
  filter(casrn %in% E2_P_data$casn) %>%  
  #want our columns of interest in the front
  select(risk_assessment_class, species_supercategory, toxval_type, toxval_numeric_qualifier, 
         toxval_numeric, toxval_units, everything()) %>% 
  #identify if study is a developmental or reproductive study conducted in mammals, 
  #we are interested in these
  mutate(risk_assessment_class = tolower(risk_assessment_class), 
         `Repro or dev study in mammals` = case_when (grepl("repro|dev", risk_assessment_class) == TRUE & 
                                                        risk_assessment_class != "developmental neurotoxicity" &
                                                        grepl("mammals", species_supercategory) == TRUE ~ 1, 
                                                      TRUE ~ 0)) %>% 
  #identify the type of toxicity value reported, is it reported an effect level or a no effect level
  mutate(`Effect value type` = case_when(`Repro or dev study in mammals` == 1 & 
                                           toxval_type %in% c("BMDL", "ED10", "LEL", "LOAEC", "LOAEL", "LOEC",  "LOEL", "BMC05") ~ "effect", 
                                         `Repro or dev study in mammals` == 1 & 
                                           toxval_type %in% c("HNEL", "NEL", "NOAEC", "NOAEL", "NOEL", "NOEC") ~ "no effect",
                                         TRUE ~ NA_character_)) %>% 
  #specify if reproductive or developmental study
  mutate(`Study type` = case_when(`Repro or dev study in mammals` == 1 & 
                                    risk_assessment_class == "reproductive developmental" ~ "repro/dev" ,
                                    grepl("repro", risk_assessment_class) ~ "reproductive", 
                                  `Repro or dev study in mammals` == 1 &
                                    grepl("dev", risk_assessment_class) ~ "developmental", 
                                  TRUE ~ "")) %>% 
  select(`Repro or dev study in mammals`, `Effect value type`, `Study type`, everything()) %>% 
  #categorize each study as not reporting an effect, reporting an effect level (greater or less than 100?), or 
  #reporting a no effect level (greater or less than 100?)
  mutate(`Effect type` = case_when(`Repro or dev study in mammals` == 0 
                                   ~ "not a developmental or reproductive risk assessment", 
                                   #effect level under 100 
                                   `Repro or dev study in mammals` == 1 &
                                     `Effect value type` == "effect" &
                                     toxval_numeric_qualifier %in% c("=", "<", "<=", "~") &
                                     toxval_numeric < 100 &
                                     toxval_units == "mg/kg-day" 
                                   ~ "toxicant with effect level < 100",   
                                   #effect level over 100
                                   `Repro or dev study in mammals` == 1 &
                                     `Effect value type` == "effect" &
                                     toxval_numeric_qualifier %in% c("=", ">", ">=", "~") &
                                     toxval_numeric >= 100 &
                                     toxval_units == "mg/kg-day" 
                                   ~ "toxicant with effect level >= 100", 
                                   #any with no effect reported? 
                                   #under 100 NOEL 
                                   `Repro or dev study in mammals` == 1 &
                                     `Effect value type` == "no effect" &
                                     toxval_numeric_qualifier %in% c("=", "<", "<=", "~") &
                                     toxval_numeric < 100 &
                                     toxval_units == "mg/kg-day" 
                                   ~ "NOEL less than 100 mg/kg-day", 
                                   #over 100 noel 
                                   `Repro or dev study in mammals` == 1 &
                                     `Effect value type` == "no effect" &
                                     toxval_numeric_qualifier %in% c("=", ">", ">=", "~") &
                                     toxval_numeric >= 100 &
                                     toxval_units == "mg/kg-day" 
                                   ~ "unlikely (NOEL over 100 mg/kg-day)", 
                                   #everything else will not be categorizd e.g. not units of mg/kg-day, 
                                   #indicates ">" or ">=" some value greater than a number under 100 but does not specify
                                   #or indicates a no effect level < some number but doesn't specify
                                   TRUE ~ "uncategorized (vague toxicity value or different units)")) %>% 
  mutate(summary = case_when(`Repro or dev study in mammals` == 1 ~ paste(`Study type`, `Effect type`, sep = " - "), 
                             TRUE ~ paste(`Effect type`))) %>% 
  select(-c(`Study type`, `Effect type`)) %>% 
  select(summary, everything())


#checking if there are any chemicals which increased estradiol or progesterone 
#but which are not in the filtered in-vivo dataset 
not_in_ToxvalDb <- E2_P_data%>%
  filter(!casn %in% filtered_invivo_ToxVal$casrn) #71 chemicals 

#chems in database but with no reproductive or developmental tests in mammals
chem_no_test <- filtered_invivo_ToxVal %>% 
  group_by(casrn, `Repro or dev study in mammals`) %>% count(name="count") %>% 
  pivot_wider(names_from = "Repro or dev study in mammals", values_from = "count") %>% 
  filter(is.na(`1`))  #159 chemicals


#we can see which studies are in the categories we defined above to double check again 
#use them in tables below

repro_unlikely_toxval <- filtered_invivo_ToxVal %>% 
  filter(summary %in% c("reproductive - unlikely (NOEL over 100 mg/kg-day)", 
                        "repro/dev - unlikely (NOEL over 100 mg/kg-day)")) 


devo_unlikely_toxval <- filtered_invivo_ToxVal %>% 
  filter(summary %in% c("developmental - unlikely (NOEL over 100 mg/kg-day)", 
                        "repro/dev - unlikely (NOEL over 100 mg/kg-day)")) 

repro_under100 <- filtered_invivo_ToxVal %>% 
  filter(summary %in% c("reproductive - toxicant with effect level < 100", 
                        "repro/dev - toxicant with effect level < 100")) 

devo_under100 <- filtered_invivo_ToxVal %>%
  filter(summary %in% c("developmental - toxicant with effect level < 100", 
                        "repro/dev - toxicant with effect level < 100")) 


#every other study in the filtered_invivo_ToxVal dataset will be "imprecise" because using study alone we can't determine 
#if LOEL is under 100 or NOEL above 100, we would need additional data
#examples of why
# 1) if only a NOEL below 100 mg/kg-day is reported but no LOEL, what was the LOEL? is it also under 100 mg/kg-day or above? we can't say
# 2)  we can't interpret a study reporting a LOEL that is less than 100 mg/kg-day and using ">" or a ">=" sign because 
#     this can indicate a LOEL that is indeed less than 100 or greater than 100
# 3) if only a LOEL above 100 mg/kg-day is reported but no NOEL, does this mean the NOEL is above 100 mg/kg-day as well or under 100 mg/kg-day
# 4) if only a NOEL under 100 mg/kg-day is reported and no LOEL mg/kg-day, does this mean the LOEL is also under 100 mg/kg-day or above? 
# 5) units are not in mg/kg-day

imprecise_ToxVal_repro <- filtered_invivo_ToxVal %>% 
  filter(`Repro or dev study in mammals` == 1, #select only studies with reproductive or developmental tests in mammals 
         ! summary %in% c("reproductive - unlikely (NOEL over 100 mg/kg-day)", 
                          "developmental - unlikely (NOEL over 100 mg/kg-day)", 
                          "repro/dev - unlikely (NOEL over 100 mg/kg-day)", 
                          "reproductive - toxicant with effect level < 100", 
                          "developmental - toxicant with effect level < 100", 
                          "repro/dev - toxicant with effect level < 100")) %>% #remove the test summaries that have enough information on their own
  filter(grepl("repro", summary)) #only reproductive related studies


imprecise_ToxVal_dev <- filtered_invivo_ToxVal %>% 
  filter(`Repro or dev study in mammals` == 1, #select only studies with reproductive or developmental tests in mammals 
         ! summary %in% c("reproductive - unlikely (NOEL over 100 mg/kg-day)", 
                          "developmental - unlikely (NOEL over 100 mg/kg-day)", 
                          "repro/dev - unlikely (NOEL over 100 mg/kg-day)", 
                          "reproductive - toxicant with effect level < 100", 
                          "developmental - toxicant with effect level < 100", 
                          "repro/dev - toxicant with effect level < 100"))  %>% #remove the test summaries that have enough information on their own
  filter(grepl("dev", summary)) #only developmental related studies

#Carcinogenicity in Toxvaldb ----------------------------------------------------

#use the carcinogenicity assessment file in ToxVal DB to find the the cancer 
#classification of chemicals 

#want to filter to only our chemicals of interest (E2 and P4 up chemicals)
cancer_ToxVal_filt <- cancer_ToxVal %>%
  filter(casrn %in% E2_P_data$casn)

unique(cancer_ToxVal_filt$cancer_call) #33 unique cancer calls to be assigned 

#assign each cancer call as indicating "likely" carcinogencity, "unlikely" or 
#"not enough info to assess"

#calls indicating likely
likely_call <- c("A (Human carcinogen)",
                 "Carcinogenic to humans",
                 "Known Human Carcinogen",
                 "Group 1 - Carcinogenic to humans",
                 "potential occupational carcinogen",
                 "B2 (Probable human carcinogen - based on sufficient evidence of carcinogenicity in animals)",
                 "C (Possible human carcinogen)",
                 "Likely to be carcinogenic to humans",
                 "Likely to be carcinogenic to humans (oral route)"  ,
                 "Suggestive evidence of carcinogenic potential"  ,
                 "Suggestive evidence of carcinogenicity in humans" ,
                 "Reasonably Anticipated To Be Human Carcinogen",
                 "Group 2A - Probably carcinogenic to humans",
                 "Group 2B - Possibly carcinogenic to humans",
                 "Group IIIB: (possibly carcinogenic to humans)",
                 "Group II: CEPA (probably carcinogenic to humans)",
                 "Group B2 Probable Human Carcinogen",
                 "Group C Possible Human Carcinogen",
                 "Likely to be Carcinogenic in Humans at High Doses; Not Likely to be Carcinogenic to Humans at Low Doses",
                 "Likely to be Carcinogenic to Humans",
                 "Suggestive Evidence of Carcinogenicity but Not Sufficient to Assess Human Carcinogenic Potential",
                 "Suggestive Evidence of Carcinogenicity to Humans",
                 "Likely Human Carcinogen")

#chemicals indicating likely 
likely <- cancer_ToxVal_filt %>% 
  filter(cancer_call %in% likely_call)


#cancer calls indicating not likely
notcancer_call <- c("Group E Evidence of Non-carcinogenicity for Humans", 
                    "Not Likely to Be Carcinogenic in Humans")

#chemicals indicating not likely
notcancer <- cancer_ToxVal_filt %>% 
  filter(cancer_call %in% notcancer_call)


#cancer calls indicating insufficient info
insufficient_call <- 
  c("D (Not classifiable as to human carcinogenicity)", 
    "Inadequate for an assessment of carcinogenic potential",     
    "Group 3 - Not classifiable as to its carcinogenicity to humans", 
    "Group VA: CEPA (inadequate data for evaluation)",                                                        
    "Group D: IRIS (not classifiable as to human carcinogenicity)", 
    "Data are Inadequate for an Assessment of Human Carcinogenic Potential",
    "Group D Not Classifiable as to Human Carcinogenicity",   
    "Not Yet  Determined")

#chemicals indicating insufficient detail
insufficient <- cancer_ToxVal_filt %>% 
  filter(cancer_call %in% insufficient_call)

#Prop65 ------------------------------------------------------------------------

#for prop65 want to characterize chemicals for their designation as "carcinogens" 
#or "developmental toxicants"

prop65_tidy <- prop65 %>% 
  #lets clean up data a bit, delet chemicals that have been delisted
  filter(!grepl("Delisted", Chemical)) %>% 
  #some chemicals contain multiple casn, lets separate them (gives warning message 
  #some cells can't be separated but we can just ignore)
  separate(`CAS No.`, c("A","B"), sep = "([/;])") %>%
  #recombine separated casn into one column, remove white space, and remove any chemicals without a casn
  pivot_longer(cols = c(A, B), values_to = "CAS No.") %>% 
  filter(!(is.na(`CAS No.`) | `CAS No.` == "---")) %>% 
  mutate(across(where(is.character), str_trim))


#chemicals indicated potential carcinogens
prop65_canc <- prop65_tidy %>% 
  group_by(Chemical, `CAS No.`) %>% 
  filter(grepl("cancer", `Type of Toxicity`))

#chemicals indicated developmental toxicants
prop65_dev <- prop65_tidy %>% 
  group_by(Chemical, `CAS No.`) %>% 
  filter(grepl("development", `Type of Toxicity`))



#table of all the E2-up and P4-up chemicals ------------------------------------


#gather all of the invivo data into one file 

#matching by casn or chemical name 

invivo_comb <- E2_P_data %>%
  select(chnm, casn) %>%
  #is chemical in NHANES? 1 = yes, 0 = no
  mutate(Nhanes = case_when (casn %in% nhanes$CASRN  ~ 1, 
                             TRUE ~ 0)) %>% 
  #is chemical identified as mammary carcinogen, developmental toxicant, or both 
  #in Rudel 2007 or Rudel 2011
  mutate(Rudel = case_when(casn %in% dev_tox$CAS_No & casn %in% mamm_car$CAS_No
                           ~ "Mammary carcinogen and developmental toxicant",
                           casn %in% dev_tox$CAS_No & !casn %in% mamm_car$CAS_No
                           ~ "Mammary developmental toxicant",
                           !casn %in% dev_tox$CAS_No & casn %in% mamm_car$CAS_No
                           ~ "Mammary carcinogen",
                           !casn %in% dev_tox$name & !casn %in% mamm_car$name
                           ~ "",
                           TRUE ~ "fix")) %>%
  #what is the carcinogenicity assessment in ToxValDb? Giving priority to "Likely", followed by "Unlikely" and then
  #"inadequate evidence" if there is a conflict between authoritative sources
  mutate(Carcinogencity = case_when(casn %in% likely$casrn ~ "Likely",
                                    casn %in% notcancer$casrn ~ "Unlikely",
                                    casn %in% insufficient$casrn  ~ "Inadequate evidence",
                                    !casn %in% cancer_ToxVal_filt$casrn ~ "No assessment",
                                    TRUE ~ "FIX")) %>% 
  #are there any reported invivo effects at chemical concentrations 
  #less than 100 mg/kg-day, as reported in ToxValDb. Are they repro, dev or both?
  mutate(`Invivo effect < 100 mg/kg-day` = 
           case_when((casn %in% devo_under100$casrn & casn %in% repro_under100$casrn) 
                     ~ "Reproductive and developmental",
                     casn %in% repro_under100$casrn
                     ~ "Reproductive", 
                     casn %in% devo_under100$casrn
                     ~ "Developmental",
                     TRUE ~ "")) %>% 
  #are there any reported no effect levels at chemical concentrations greater than 
  #100 mg/kg-day.? Repro, dev or both? Only writing if there was no in vivo effect < 100 mg/kg-day based on previous column
  mutate(`No invivo effect >= 100 mg/kg-day` = 
           case_when((casn %in% devo_unlikely_toxval$casrn & casn %in% repro_unlikely_toxval$casrn) & 
                       `Invivo effect < 100 mg/kg-day` == "" 
                     ~ "Reproductive and developmental", 
                     casn %in% repro_unlikely_toxval$casrn &
                       !grepl("Repro", `Invivo effect < 100 mg/kg-day`)
                     ~ "Reproductive",
                     casn %in% devo_unlikely_toxval$casrn  &
                       !grepl("mental", `Invivo effect < 100 mg/kg-day`)
                     ~ "Developmental",
                     TRUE ~ "")) %>% 
  #if a chemical is unaccounted for in previous ToxVal re columns, is it because a study was reporting inprecise units 
  #or because chemical was not tested in a reproductive/developmental test? 
  mutate(`Other (ToxValDb)` = 
           
           #repro/dev study available in ToxValDb for chemical but we couldn't assign it to LOEL under 100 or NOEL over 100
           case_when((casn %in% imprecise_ToxVal_repro$casrn & casn %in% imprecise_ToxVal_dev$casrn) & 
                       `Invivo effect < 100 mg/kg-day` == "" & 
                       `No invivo effect >= 100 mg/kg-day` == "" 
                     ~ "repro/dev - imprecise toxicity value or different units", 
                     
                      casn %in% imprecise_ToxVal_repro$casrn &
                       !grepl("(R|r)epro", `Invivo effect < 100 mg/kg-day`) & 
                       !grepl("(R|r)epro", `No invivo effect >= 100 mg/kg-day`)
                     ~"repro - imprecise toxicity value or different units",
                     
                     casn %in% imprecise_ToxVal_dev$casrn &
                       !grepl("(D|d)ev", `Invivo effect < 100 mg/kg-day`) & 
                       !grepl("(D|d)ev", `No invivo effect >= 100 mg/kg-day`)
                     ~"dev - imprecise toxicity value or different units",
                     
                     #chemicals with no repro/dev tests found in ToxValDb
                     
                     casn %in% not_in_ToxvalDb$casn | casn %in% chem_no_test$casrn ~
                       "no repro nor dev tests in ToxValDb", 
                     TRUE ~ "")) %>% 
  
  #Is the chemical in Prop65, and is it in there as a carcinogen, developmental toxicant or both? 
  mutate(Prop65= case_when(casn %in% prop65_canc$`CAS No.` & casn %in% prop65_dev$`CAS No.` 
                           ~ "Cancer and developmental", 
                           casn %in% prop65_canc$`CAS No.` & !casn %in% prop65_dev$`CAS No.`
                           ~ "Cancer", 
                           !casn %in% prop65_canc$`CAS No.` & casn %in% prop65_dev$`CAS No.` 
                           ~ "Developmental", 
                           !casn %in% prop65_canc$`CAS No.` & !casn %in% prop65_dev$`CAS No.`
                           ~ "", 
                           TRUE ~ "fix")) 

#Lets use the dataset we created to assign a simplified carcinogenicity and dev/repro
#toxicity assignment for each chemical  

invivo_arrange <- invivo_comb %>% 
  #if chemical was reported as carcinogen in ToxValDb, Rudel 2007 or Prop65 assign as likely carcinogenic
  #if chemical was reported as unlikely in ToxValDb, assign as unlikely 
  #if chemical not assigned a carcinogenicity assessment or reported as having insufficient evidence in 
  #ToxValDB, assign "inadequate evidence" 
  mutate(`Carcinogenicity Assessment (Summary)` = case_when(Carcinogencity == "Likely" | grepl("carcinogen", Rudel) == TRUE |
                                                              grepl("Cancer", Prop65) == TRUE  
                                                            ~ "Likely",
                                                            Carcinogencity == "Unlikely" ~ "Unlikely",
                                                            Carcinogencity %in% c("Inadequate evidence", "No assessment") 
                                                            ~ "Inadequate evidence", 
                                                            TRUE ~ "fix")) %>% 
  
  #If chemical reported in ToxValDB as having repro and/or dev effect under 100mg-kg/day, 
  #in Rudel 2011 as mammary gland developmental toxicant or in Prop65 as developmental toxicant,
  #assign "likely" dev/repro tox. 
  #
  #If a study reports a NOEL above 100 mg/kg-day 
  #for both repro and dev study, assign unlikely to be dev/repro tox. 
  #
  #if chemical was neither likely or unlikely, and missing a repro or dev study, or had
  #vague toxicity value or different units than mg/kg-day, assign inadequate evidence
  #
#We are prioritizing likely over unlikely over inadequate evidence (okay to overwrite)
mutate (`Developmental or reproducive toxicity (Summary)` = case_when(grepl("(D|d)evelopment|Reproductive", `Invivo effect < 100 mg/kg-day`) == TRUE |
                                                                        grepl("(D|d)evelopment", Rudel) == TRUE |
                                                                        grepl("(D|d)evelopmental", Prop65) == TRUE
                                                                      ~ "Likely", 
                                                                      
                                                                        `No invivo effect >= 100 mg/kg-day` == "Reproductive and developmental" 
                                                                      ~ "Unlikely",
                                                                      
                                                                      `Other (ToxValDb)` != "" |
                                                                        `No invivo effect >= 100 mg/kg-day` %in% c("Reproductive", "Developmental")
                                                                      ~ "Inadequate evidence",
                                                                      TRUE ~ "FIX")) %>% 
  #adjust column names for export
  rename("Invivo effect < 100 mg/kg-day (ToxValDb)" = "Invivo effect < 100 mg/kg-day",
         "No invivo effect >= 100 mg/kg-day (ToxValDb)" = "No invivo effect >= 100 mg/kg-day",
         "Carcinogencity (ToxValDb)" = "Carcinogencity")


write_csv(invivo_arrange, "./output/rearrange_invivo.csv")



#supplemental ------------------------------------------------------------------

#for supplemental, gather all the carcinogenicity calls given to a chemical by 
#each authoritative source (as listed in ToxValDB). 
#will indicate which sources contradict with each other 

cancer_ToxVal_rearrange <- cancer_ToxVal %>% 
  filter(casrn %in%  E2_P_data$casn) %>% 
  mutate(cancer_call2 = case_when(exposure_route != "-" ~  paste(cancer_call, " (", exposure_route, ")", sep = ""), 
                                  TRUE ~ cancer_call)) %>% 
  select(casrn, name, source, cancer_call2) %>% 
  group_by(casrn, name, source) %>%
  summarise(cancer_call2, sep = ", ") %>% 
  pivot_wider(names_from = "source", values_from = "cancer_call2") %>% 
  select(-sep)

write_csv(cancer_ToxVal_rearrange, "./output/ToxValDb_cancercalls.csv")
