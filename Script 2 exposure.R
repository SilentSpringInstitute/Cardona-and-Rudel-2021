
# Application of an in vitro assay to identify chemicals that increase estradiol 
# and progesterone synthesis and are potential breast cancer risk factors
# Bethsaida Cardona1, Ruthann A. Rudel1
# 1Silent Spring Institute, Newton, MA 
# Corresponding author: Ruthann Rudel, Silent Spring Institute, Newton, MA 02460, USA. 

# May 2021
# R version 3.6.3 (2020-02-29)

#R Script 2: exposures. Code focusing on exposure sources
#-Needs Script 1 to run. Can run concurrently with Script 2: invivo

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
library(stringi)


#effectively disables scientific notation
options(scipen = 999)

#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read in data  -----
#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#data with potency measures and exposure predictions from Script 1
exposure_potency <- read_csv("./output/exposure_potency.csv") %>% 
  mutate(casn = gsub("'", "", casn))


#Files for data on chemical exposure sources ----------------------------

#Following files from CPDat/CPCat, last updated in 2018 (data download: https://comptox.epa.gov/dashboard/downloads)

#CPCat
#links CASN to the associated CPDat record id
casn_expo <- read_excel("./input/cpdat_source_substances.xlsx")
#has the product type as well as the CPDat record id 
product_type <- read_excel("./input/chemical_and_product_categories.xlsx")

#CPDAt 
#list of products with no ingredients
products <- read_excel("./input/products.xlsx")
#list of ingredients with the CPDat record id but no CASN
ingredients <- read_excel("./input/ingredients.xlsx")


# Supplemental file from Ring 2019 paper containing HT-exposure predictions
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


################################################################################
# Chemical exposure sources - broad ---------------------------------------------
################################################################################

# Using data from CPCat 

#this will produce the dataset used for both the broad and specific categories
chem_sourcesa <- product_type %>% 
  #some chemical/source/cassette pairs repeat, so we only want to count 1 instance,
  #we don't care about web_url or website download info
  group_by(source, cassette, cpdat_substance_record_id) %>% count(name = "repeats") %>% 
  #remove any cassette with term "discontinued", "prohibited", 'restricted"
  filter(!grepl("discontinued|prohibited|restricted", cassette)) %>% 
  left_join (casn_expo) %>% 
  left_join (., exposure_potency %>% select(casn,chnm), by = c("casrn" = "casn")) %>% 
  #filter to chemicals that increase E2 or P4
  filter(casrn %in% exposure_potency$casn) %>% 
  #for chemicals with a sector specific source, we want to assign it into its general category 
  mutate(fromsource = gsub("_ACToRUseDB", "", cassette))  %>% 
  mutate(fromsource = case_when(
    grepl("Industrial", source) ~ "ind.",
    grepl("NACE", source) ~ "ind.",
    grepl("Retail", source) & !grepl("food|beverage", cassette) ~ "cons.",
    grepl("Drug", source) ~ "pharma.", 
    grepl("DfE", source) ~ "functional", 
    grepl("Dow", source) ~ "functional", 
    grepl("Kem", source) ~ "functional",
    grepl("DfE", source) ~ "functional", 
    grepl("Toxome", source) ~ "functional", 
    grepl("Functional", source) ~ "functional", 
    TRUE ~ cassette)) %>% 
  #for the following we are mainly interested in the first term... 
  mutate(split = sapply(str_split(fromsource, pattern = "\\|", n = Inf), "[", 1)) %>% 
  #unless the product has an industrial or functional use, if it has a pesticide associated term we are
  #interested in extracting it (does not matter if it is categorized as consumer_use as well)
  mutate(cassette3 = case_when(grepl("pesticide", split) ~ "pest.",
                               grepl("herbicide", split) ~ "pest.",
                               grepl("extermination", split) ~ "pest.",
                               grepl("fungicide", split) ~ "pest.",
                               grepl("antimicrobial", split) ~ "pest.",
                               grepl("drug", split) ~ "pharma.",
                               grepl("food", split) ~ "diet", #we also want to include anything that has contact with food 
                               grepl("beverage", split) ~ "diet",
                               grepl("drinking_water", split) ~ "diet",
                               TRUE ~ split)) %>% 
  #if any of the above are also categorized as "consumer_use" we want to take note of that...
  mutate(cassette3 = case_when(cassette3 %in% c("pest.", "pharma", "diet.") & grepl("consumer_use", cassette) ~ paste(cassette3, "cons.", sep = " ,"), 
                               grepl("consumer_use", fromsource) ~ "cons.",
                               TRUE ~ cassette3)) %>% 
  #for the remaining chemicals in the "general uses" databses, we will manually 
  #assigned them into a category based on terms in cassette
  mutate(cassette4 = case_when(grepl("apparel", cassette) &!grepl("ind.", fromsource) ~ "cons.",
                               grepl("personal_care", cassette3) ~ "cons.",
                               grepl("arts_crafts", cassette3) ~ "cons.",
                               grepl("furniture", cassette) & !grepl("ind.", fromsource) ~ "cons.",
                               grepl("child_use", cassette3) ~ "cons.",
                               grepl("decor", cassette3) ~ "cons.", 
                               grepl("toy", cassette3) ~ "cons.",
                               grepl("electronics", cassette3) ~ "cons.",
                               grepl("lawn_garden", cassette3) ~ "cons.",
                               grepl("sports_equipment", cassette3) ~ "cons.",
                               grepl("baby_use", cassette3) ~ "cons.",
                               grepl("pet", cassette3) ~ "cons.",
                               grepl("tools|personal_care|dental|toothbrush|detected", cassette3) ~ "cons.",
                               grepl("cleaning_washing", cassette3) ~ "cons.",
                               grepl("industrial", cassette3) ~ "ind.",
                               grepl("raw_material", cassette3) ~ "ind.",
                               grepl("industrial_manufacturing", cassette3) ~ "ind.",
                               grepl("industrial_fluid", cassette3) ~ "ind.",
                               grepl("mining", cassette3) ~ "ind.",
                               grepl("manufacturing", cassette3) ~ "ind.",
                               grepl("resource_extraction", cassette3) ~ "ind.",
                               grepl("rubber_processing", cassette3) ~ "ind.",
                               grepl("soap", cassette3) ~ "cons.",
                               grepl("disinfectant", cassette3) ~ "pest.",
                               grepl("automotive_care", cassette) &!grepl("ind.", fromsource)  ~ "cons.",
                               TRUE ~ cassette3), 
         #sometimes something labeled as pesticide or consumer use is also found in diet through food contact
         cassette4 = case_when(grepl("food|beverage", cassette) & !cassette4 %in% c("diet", "ind.") & !cassette3 %in% c("ind.") ~ paste(cassette4, "diet", sep = ", "),
                               grepl("extermination", cassette) ~ paste(cassette4, "pest.", sep = ", "),
                               grepl("disinfectant", cassette) & !cassette4 %in% c("disinfectant", "pest.", "feed") ~ paste(cassette4, "pest.", sep = ", "),
                               TRUE ~ cassette4))

                            

#bring in CPDat... everything is a consumer product, some may also be pesticides
chem_sourcesb <- casn_expo %>% 
  right_join(ingredients) %>% 
  select(3:6) %>% 
  right_join(products) %>% 
  filter(casrn %in% exposure_potency$casn, 
         !is.na(product_use_category)) %>% 
  group_by(casrn, product_use_category, source) %>% count(name = "repeats") %>% ungroup() %>% 
  mutate(cassette4 = case_when(grepl("pesticide", product_use_category) ~ "pest., cons.", 
                               grepl("herbicide", product_use_category) ~ "pest., cons.", 
                               TRUE ~ "cons.")) 



#broad categories for CPCat
broad_cpcat <- chem_sourcesa %>% 
  group_by(cassette4, casrn) %>% slice(1) %>% ungroup() 

#broad categories for CPDat
broad_cpdat <- chem_sourcesb %>% 
  group_by(cassette4, casrn) %>% slice(1) %>% ungroup() 

#combine two "broad" datasets 
chem_sources_broad <- plyr::rbind.fill(broad_cpcat, broad_cpdat) %>% 
  #capitalize the sources, remove "-"
  mutate(cassette4 = stri_trans_totitle(cassette4), 
         cassette4 = case_when(grepl("_", cassette4) ~  gsub("_", " ", cassette4),
                               TRUE ~ cassette4)) %>% 
  group_by(cassette4, casrn) %>% slice(1) %>% ungroup() %>% 
  #we want to list all of the broad sources into one column based on the chemical
  group_by(casrn) %>% 
  summarize(broad = paste(cassette4,collapse=", ")) %>% 
  ungroup() %>% 
  #for each of the "main" exposure sources (diet, cons. ind. pharma. pest.) we want to separate into their own columns 
  mutate(Cons. = case_when(grepl("Cons", broad) == TRUE ~ "1", 
                           TRUE ~ "0"), 
         Diet = case_when(grepl("Diet", broad)  == TRUE | grepl("Drinking", broad)  ~ "1", 
                          TRUE ~ "0"),
         Ind. = case_when(grepl("Ind.", broad)  == TRUE ~ "1", 
                          TRUE ~ "0"),
         Pharma. = case_when(grepl("Pharma.", broad)  == TRUE ~ "1", 
                             TRUE ~ "0"),
         Pest. = case_when(grepl("Pest.", broad)  == TRUE ~ "1", 
                           TRUE ~ "0"),
         'Uncategorized Use' = case_when(broad != " " & Diet == "0" & Ind. == "0" & Cons. == "0" & Pharma. == "0" & Pest. == "0" ~ "1",
                                         TRUE ~ "0"))  


################################################################################
# Chemical exposure sources - Specific -----------------------------------------
################################################################################

#we are also interested in gathering more specific information about the product uses of these chemicals 
#we want to see if any of these chemicals are in the following: fragrance, microbial, 
#hair dyes, personal care, drinking water contaminant, cigarettes, human metabolites, 
#flame-retardant and plastics. For some of these, the purpose may be industrial or functional

#using CPCAT data first 
specific_cpcat <- chem_sourcesa %>% 
  mutate(term1 = sapply(str_split(cassette, pattern = "\\|", n = Inf), "[", 1)) %>% 
  mutate(Other = case_when (grepl("personal_care", term1) ~ "pers.",
                            grepl("fragrance", cassette) ~ "fragrance", #some fragrance is also personal care, may refer to "fragrance
                            #as a functional use
                            grepl("air_freshener", cassette) ~ "fragrance", 
                            grepl("antimicrobial|disinfectant", cassette) ~ "antimicrobial", 
                            grepl("cigarettes", term1) ~ "cigarettes",
                            grepl("human_metabolite", cassette) ~ "human_metabolite", 
                            grepl("plastic", cassette) ~ "plastic", 
                            grepl("drinking_water", cassette) ~ "drinking_water_contaminant",
                            grepl("flame_retardant", cassette) ~ "flame_retardant", #functional use included
                            grepl("food_contact", cassette) ~ "food_contact", #includes peticides used on food,
                            grepl("food_additive", cassette) ~ "food_additive",
                            grepl("food_residue", cassette) ~ "food_residue",
                            grepl("textile", term1) ~ "textile",
                            TRUE ~ " ")) %>% 
  #the others have some interferences with the above so we will paste additional material next to them 
  mutate(Other = case_when(grepl("hair_dye", cassette) & Other != "hair_dye" ~ paste(Other, "hair_dye", sep = ", "),
                           grepl("fragrance", cassette) & Other != "fragrance" ~ paste(Other, "fragrance", sep = ", "),
                           grepl("drinking_water", cassette) & Other != "drinking_water_contaminant" ~ paste(Other, "drinking_water_contaminant", sep = ", "),
                           grepl("food_contact", cassette) & Other != "food_contact" ~ paste(Other, "food_contact", sep = ", "),
                           TRUE ~ Other)) 



#now lets gather the specific data for cpdat
specific_cpdat <- casn_expo %>% 
  right_join(ingredients) %>% 
  select(3:6) %>% 
  right_join(products) %>% 
  filter(casrn %in% exposure_potency$casn, 
         !is.na(product_use_category)) %>% 
  group_by(casrn, product_use_category, source) %>% slice(1) %>% ungroup() %>% 
  mutate(Other = case_when(grepl("fragrance", product_use_category) ~ "fragrance", 
                           #since personal products may apply to things such as fragrance, we want to separate it 
                           grepl("air freshener", product_use_category) ~ "fragrance", 
                           grepl("hair color", product_use_category) ~ "hair_dye", 
                           grepl("disinfectant", product_use_category) ~ "antimicrobial", 
                           TRUE ~ " ")) %>% 
  #since personal products may apply to things such as fragrance, we want to separate it 
  mutate(Other = case_when(grepl("personal care", product_use_category) & Other != " " ~ paste("pers.", Other, sep = ", "),
                           grepl("personal care", product_use_category) & Other == " " ~ "pers.", 
                           TRUE ~ Other)) %>% 
  group_by(casrn, Other) %>% slice(1) %>% ungroup() 




#bind the specific data from cpdat and cpcat
chem_sources_specific <- plyr::rbind.fill(specific_cpcat, specific_cpdat) %>% 
  select(casrn, chnm, cassette, Other) %>% 
  #capitalize the sources, remove "-"
  mutate(Other = stri_trans_totitle(Other), 
         Other = case_when(grepl("_", Other) ~  gsub("_", " ", Other),
                           TRUE ~ Other)) %>% 
  group_by(Other, casrn) %>% slice(1) %>% ungroup() %>% 
  #we want to list all of the other sources into one column based on the chemical
  group_by(casrn) %>% 
  summarize(Other = paste(Other, collapse= ", ")) %>% 
  ungroup() %>% 
  #we want to remove any commas at the start or any repeats.
  mutate(Other= gsub("^ , ", " ", Other)) 



################################################################################
## Combine the broad and specific data by the casn -----------------------------
################################################################################


#we also want to add source predictions based on Ring data - only adding if pathway assigned using CPDat/CPCat is 0
chem_sources <- chem_sources_broad %>% 
  left_join(chem_sources_specific) %>%
  right_join(exposure_potency %>% select(casn), by = c("casrn" = "casn")) %>% 
  left_join(htp_pred %>% select(CAS, Pathway), by = c("casrn" = "CAS")) %>% 
  mutate(Pathway = case_when(Pathway == "All Four" ~ "Cons., Ind., Pest., Diet.", 
                             TRUE ~ Pathway)) %>% 
  mutate_at(c("Cons.", "Diet", "Pest.", "Ind."), as.character) %>%
  mutate(Cons. = case_when(grepl("Cons", Pathway) == TRUE & Cons. %in% c("0", NA_real_) ~ "1*", 
                           TRUE ~ Cons.), 
         Diet = case_when(grepl("Diet", Pathway) == TRUE & Diet %in% c("0", NA_real_) ~ "1*", 
                          TRUE ~ Diet), 
         Pest. = case_when(grepl("Pest", Pathway) == TRUE & Pest. %in% c("0", NA_real_) ~ "1*", 
                           TRUE ~ Pest.),
         Ind. = case_when(grepl("Ind.", Pathway) == TRUE & Ind. %in% c("0", NA_real_) ~ "1*", 
                          TRUE ~ Ind.)) %>% 
  mutate(`No Toxcast Info.` = case_when(Pathway == "Unknown" & is.na(Cons.) == TRUE & is.na(Diet) == TRUE 
                                        & is.na(Pest.) == TRUE & is.na(Ind.) == TRUE ~ "1",
                                        TRUE ~ " ")) %>% 
  mutate(row = 1:n())


#to graph the specific sources we need to seperate them from the one column they are all listed in and 
#pivot_longer

chem_sources_graph_specific <- as_tibble(str_split(chem_sources$Other, pattern = ",", n= Inf, simplify = TRUE)) %>% 
  mutate(row = 1:n()) %>%  
  pivot_longer(cols = contains("V", ignore.case = FALSE), names_to = "source") %>%  
  filter(!value %in% c("", " ", NA_real_)) %>% 
  group_by(value, row) %>% slice(1) %>% ungroup() %>%  
  mutate(value =  gsub("^\\s+", "", value), 
         value = case_when(grepl("Pers.", value) ~ "Personal care",
                           TRUE ~ value)) %>% 
  right_join(chem_sources) %>% 
  select(-c(broad, Pathway, `Uncategorized Use`)) 

write_csv(chem_sources_graph_specific %>% mutate(casrn = paste0("'", casrn)), "./output/chem_sources_graph_specific.csv")



chem_sources_tidy <- chem_sources_graph_specific %>% 
  pivot_wider(names_from = "source", values_from = "value") %>% 
  unite(col = "Sources of Interest", contains("V", ignore.case = FALSE), sep = ", ", remove = TRUE, na.rm = TRUE) %>%  
  mutate(`Sources of Interest` = tolower(`Sources of Interest`)) %>% 
  select(-c(row, Other, `NA`))
  
write_csv(chem_sources_tidy %>% mutate(casrn = paste0("'", casrn)), "./output/chem_sources.csv")






