
# Application of an in vitro assay to identify chemicals that increase estradiol 
# and progesterone synthesis and are potential breast cancer risk factors
# Bethsaida Cardona1, Ruthann A. Rudel1
# 1Silent Spring Institute, Newton, MA 
# Corresponding author: Ruthann Rudel, Silent Spring Institute, Newton, MA 02460, USA. 

# May 2021
# R version 3.6.3 (2020-02-29)

#R Script 3: E2/P4. Code puts together all tables relevant to E2/P4 up chemicals 
#run in previous scripts 
#-Needs data output from Script 1:prep data, Script 2:invivo and Script 2:exposure 


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
library(ggmosaic)
library(ggstatsplot)
library(ggpubr)
library(ggplot2)
library(ggpattern)

#effectively disables scientific notation
options(scipen = 999)

#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read in data  -----
#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#data with potency measures and median exposure predictions. Output from R Script 1
exposure_potency <- read_csv("./output/exposure_potency.csv")

#data with chemical exposure sources. Output from R script 2:exposures
chem_sources <- read_csv("output/chem_sources.csv") %>% 
  mutate(casrn = gsub("'", "", casrn))

#data for graphing the specific exposure sources. Output from R script 2: exposures
chem_sources_graph_specific <- read_csv("output/chem_sources_graph_specific.csv") %>% 
  mutate(casrn = gsub("'", "", casrn))

#chemical invivo data. From Script 2: invivo
in_vivo <- read_csv("output/rearrange_invivo.csv")

#ac10 and ac50 values for the H295R assay (and other Toxcast assays) 
#downloaded from Comptox website:  (last update 2-26-2019)
ac10 <- read_csv("input/modl_ac10_Matrix_190226.csv")
CASN <- read_csv("input/Chemical_Summary_190226.csv")
hitc <- read_csv(("input/hitc_Matrix_190226.csv"))

#dataset with number of pesticide products that are approved by EPA and contain 
#a certain chemical ingredient (casn shown), also contrains number of products 
#not currently approved but which have been approved in the past
#data compiled by authors using data files downloaded from 
#https://www.epa.gov/ingredients-used-pesticide-products/ppis-download-product-information-data (last update 5-5-20)
status <- read_csv("input/pesticide product status2.csv") %>% 
  mutate(CASN = gsub("'", "", CASN))


#file listing FDA approved drug products. 
#downloaded from https://www.fda.gov/drugs/drug-approvals-and-databases/approved-drug-products-therapeutic-equivalence-evaluations-orange-book
#last update 6-30-2020
fdadrugs <-  read.table("input/products.txt", sep = "~" , header = T , stringsAsFactors= F, fill = TRUE, quote = "")

#ER AUC data from Judson 2015 supplemental
#some CAS were converted to dates in original file. 
#these were manually adjusted in new column 
auc <- read_xlsx("input/supp_kfv168_toxsci-15-0258-File002_adj.xlsx")

#list of pesticides with observed mammary effects, from Cardona and Rudel 2020
mammary_pesticides <- read_csv("input/pesticides_weffects.csv")


#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Join data for E2 and P4 with with CPdat and in vivo data -----
#### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data_long <- exposure_potency %>% 
  mutate(casn = gsub("'", "", casn)) %>% 
  pivot_longer(cols = c(contains("_Estradiol"), contains("_Progesterone")), 
               names_to = c(".value", "hormone"), 
               names_sep = "_", 
               values_drop_na = TRUE) %>% 
  select(-c(date_chnm_plate, plate)) %>% 
  filter(!(is.na(LEC))) %>% 
  rename("Chemical name" = "chnm", 
         "CASN" = "casn",
        "Adj.maxmMd" = "adj_maxmMD",
        "Max. tested conc." = "max_tested_conc", 
        "Min. tested conc." = "min_tested_conc") %>% 
  #recoding the AC50
  mutate(AC50 = case_when(AC50 == 1000000 ~ "-", 
                                    is.na(AC50) == TRUE ~ " ", 
                                    TRUE ~ as.character(AC50))) %>% 
  #bringing in ER AUC information
  left_join(auc %>% select(`CASRN_Adj (4-30-21)`, AUC.Agonist, AUC.Antagonist) %>% 
              group_by(`CASRN_Adj (4-30-21)`) %>% 
              mutate(AUC.Bioactivity = max(AUC.Agonist, AUC.Antagonist)) %>% 
              pivot_longer(contains("AUC")) %>%
              mutate(value = case_when(value >= 0.1 ~ paste(as.character(value), "(active)", sep = " "),
                                       value >= 0.01 ~ paste(as.character(value), "(ambiguous)", sep = " "),
                                       TRUE ~ as.character(value))) %>%
              pivot_wider(names_from = "name", values_from = "value")  , by = c("CASN" = "CASRN_Adj (4-30-21)")) %>% 

  #chem_sources has info on chemical exposure 
  left_join(chem_sources, by = c("CASN" = "casrn")) %>% 
  rename("Predicted median intake rate (mg/kg BW/day)" = "seem",
         "Consumer" = "Cons.",
         "Industrial" = "Ind.", 
         "No exposure source data" = "No Toxcast Info.", 
         "Additional exposure sources" = "Sources of Interest") %>% 
  #mammary pesticides info
  left_join(mammary_pesticides %>% select (CASRN, status), by = c("CASN" = "CASRN")) %>% 
  mutate(`Pesticide effects on mammary gland` = case_when(status == "considered - tumor" ~ "Tumor", 
                                                          status == "not considered - both" ~ "Tumor and other effect (potential)", 
                                                          status == "not considered - tumor" ~ "Tumor (potential)", 
                                                          status == "considered - other" ~ "Other effect",
                                                          TRUE ~ "")) %>% 
  mutate(`Pesticide effects on mammary gland` = case_when(`Pesticide effects on mammary gland` != "" ~ 
                                                            paste(`Pesticide effects on mammary gland`, " (Cardona and Rudel 2020)", sep = ""), 
                                                          TRUE ~ NA_character_))  %>% 
    
  #joining with invivo sources 
  left_join(in_vivo %>% select(-"chnm"), by = c("CASN" = "casn")) %>% 
  mutate(Rudel = gsub("\\bMammary carcinogen\\b", "Tumor (Rudel 2007)", Rudel),
         Rudel = gsub("and developmental toxicant", "and alters mammary gland development (Rudel 2011)", Rudel), 
         Rudel = gsub("Mammary developmental toxicant", "Alters mammary gland development (Rudel 2011)", Rudel)) %>% 
  mutate(`Effects on mammary gland` = case_when(!is.na(Rudel) ~ Rudel, 
                                                TRUE ~ `Pesticide effects on mammary gland`)) %>% 
  select(-c(Rudel, status, `Pesticide effects on mammary gland`)) %>% 
  #ordering by fc 
  arrange(desc(MFC))
           
length(unique(data_long$CASN)) #418 unique E2/P4-up chemicals

length(unique(data_long[which((data_long$hormone == "Estradiol")),]$CASN)) #274 chemicals increase estradiol
length(unique(data_long[which((data_long$hormone == "Progesterone")),]$CASN)) #283 chemicals increase progesterone


################################################################################
##Figures for paper 
################################################################################


##Filter chemicals by efficacy and potency -------------------------------------

ds_mfc_lac <- data_long %>% 
  #removing chemicals based on exposure from the data set we made at the beginning
  filter(`MFC` >= 1.5) %>% 
  filter(Adj.maxmMd > 0) %>% 
  filter(LEC <= 33) %>% 
  #remove hormones
  filter(!CASN %in% c("53-16-7", "474-86-2", "50-27-1", "57-85-2", "57-83-0", "57-63-6", "58-18-4", 	
                      "521-18-6", "53-41-8", "68-96-2", "63-05-8", "53-43-0", "50-28-2", "57-91-0")) %>% 
  mutate(logLEC = log10(LEC), 
         seem = `Predicted median intake rate (mg/kg BW/day)`, 
         seem = parse_number(seem)) %>% 
  #for purposes of graphing, we are going to shorten the names of a few common chemicals 
  mutate(`Chemical name` = case_when(`Chemical name` == "2,2-Bis(4-hydroxyphenyl)-1,1,1-trichloroethane" ~ "HPTE", 
                                     `Chemical name` == "3,3,5-Trimethylcyclohexyl salicylate" ~ "Homosalate", 
                                     `Chemical name` == "3,3',5,5'-Tetrabromobisphenol A" ~ "3,3',5,5'-TBBPA", 
                                     `Chemical name` == "Pentachlorophenol" ~ "PCP", 
                                     `Chemical name` == "Triphenyl phosphate" ~ "TPhP", 
                                     `Chemical name` == "Di(2-ethylhexyl) phthalate" ~ "DEHP", 
                                     `Chemical name` == "7,12-Dimethylbenz(a)anthracene" ~ "DMBA", 
                                     `Chemical name` == "2-Hydroxy-4-methoxybenzophenone" ~ "BP-3", 
                                     `Chemical name` == "3,3'-Dimethoxybenzidine" ~ "o-Dianisidine", 
                                     `Chemical name` == "3,3'-Dimethylbenzidine" ~ "2-Tolidine", 
                                     `Chemical name` == "4,4'-Sulfonylbis[2-(prop-2-en-1-yl)phenol]" ~ "TGSA",
                                     `Chemical name` == "1,4-Bis(butylamino)anthracene-9,10-dione" ~ "Oil Blue 35", 
                                     `Chemical name` == "(+/-)-cis-Permethrin" ~ "cis-Permethrin", 
                                     `Chemical name` == "(E)-beta-Damascone" ~ "Damascone beta", 
                                     `Chemical name` == "17-Methyltestosterone" ~ "Methyltestosterone", 
                                     `Chemical name` == "2-Amino-5-azotoluene" ~ "Solvent Yellow 3", 
                                     `Chemical name` == "2-Hydroxyethyl acrylate" ~ "2-HEA", 
                                     `Chemical name` == "2-Methoxy-5-methylaniline" ~ "p-cresidine", 
                                     `Chemical name` == "2,2-Dibromo-3-nitrilopropionamide" ~ "DBNPA", 
                                     `Chemical name` == "4-Chloro-2-methylaniline" ~ "4-COT",
                                     `Chemical name` == "4-Chloro-2-methylphenol" ~ "4-Chloro-o-cresol",
                                     `Chemical name` == "4-Methoxyaniline hydrochloride" ~ "p-Anisidinium chloride",
                                     `Chemical name` == "4,4'-Sulfonyldiphenol" ~ "Bisphenol S", 
                                     `Chemical name` == "all-trans-Retinoic acid" ~ "Tretinoin", 
                                     `Chemical name` == "Tri-o-cresyl phosphate" ~ "TOCP", 
                                     `Chemical name` == "3,3'-Dimethoxybenzidine dihydrochloride" ~ "o-Dianisidine dihydrochloride", 
                                     `Chemical name` == "4-Chloro-3-methylphenol" ~ "4-Chloro-m-cresol", 
                                     `Chemical name` == "Nordihydroguaiaretic acid" ~ "NDGA", 
                                     `Chemical name` == "4-(Butan-2-yl)phenol" ~ "4-sec-Butylphenol", 
                                     `Chemical name` == "2-Ethoxy-5-(1-propenyl)phenol" ~ "Vanitrope", 
                                     `Chemical name` == "5-Chloro-2-methyl-3(2H)-isothiazolone" ~ "MCI", 
                                     `Chemical name` == "N-Phenyl-2-naphthylamine" ~ "Neozone",
                                     `Chemical name` == "Azinphos-methyl" ~ "Gusathion", 
                                     `Chemical name` == "Methyl 2,4-dihydroxy-3,6-dimethylbenzoate" ~ "Methyl 3-methylorsellinate", 
                                     `Chemical name` == "Bisphenol A" ~ "BPA", 
                                     `Chemical name` == "10-Undecenoic acid" ~ "Undecylenic acid", 
                                     `Chemical name` == "17alpha-Ethinylestradiol"	~ "Ethinyl estradiol",
                                     `Chemical name` == "2-Naphthalenol" ~ "2-Naphthol", 
                                     `Chemical name` == "2,4-D 1-butyl ester" ~ "2,4-D Butyl ester", 
                                     `Chemical name` == "3-Phenyl-2-propen-1-ol" ~ "Cinnamyl alcohol", 
                                     `Chemical name` == "3,3'-Dimethylbenzidine dihydrochloride" ~	"o-Tolidine dihydrochloride", 
                                     `Chemical name` == "4-(2-Methylbutan-2-yl)phenol" ~ "p-tert-Pentylphenol", 
                                     `Chemical name` == "5-Amino-2-methylphenol" ~ "5-Amino-o-cresol", 
                                     `Chemical name` == "5,7-Dimethoxy-2H-chromen-2-one" ~ "5,7-Dimethoxycoumarin",
                                     `Chemical name` == "Diethylstilbestrol" ~ "DES", 
                                     `Chemical name` == "4-(2-Phenylpropan-2-yl)-N-[4-(2-phenylpropan-2-yl)phenyl]aniline" ~ 
                                       "4-(2-Phenylpropan-2-yl)-N-\n[4-(2-phenylpropan-2-yl)phenyl]aniline",
                                     TRUE ~ `Chemical name`), 
         `Chemical name` = gsub("C.I. ", "", `Chemical name`), 
         `Chemical name` = gsub("Dichlorophenol", "DCP", `Chemical name`), 
         `Chemical name` = gsub("Trichlorophenol", "TCP", `Chemical name`), 
         `Chemical name` = gsub("Tribromophenol", "TBP", `Chemical name`), 
         `Chemical name` = gsub("Dinitrotoluene", "DNT", `Chemical name`), 
         `Chemical name` = gsub("Dimethylphenol", "xylenol", `Chemical name`)) 



length(unique(ds_mfc_lac[which(ds_mfc_lac$hormone == "Estradiol"),]$`Chemical name`)) #182 E2-up chemicals 
length(unique(ds_mfc_lac[which(ds_mfc_lac$hormone == "Progesterone"),]$`Chemical name`)) #185 P4-up chemicals


length(unique(ds_mfc_lac$CASN)) #296 unique chems after filtering based on our efficacy/potency criteria


#group_by potency and efficacy

pot_eff <- ds_mfc_lac %>% 
  group_by(hormone) %>% 
  nest() %>% 
  mutate(data_mut = purrr::map(data, function(x) x %>% select(`Chemical name`, CASN, MFC, LEC, logLEC) %>% 
                                 mutate(rank_MFC = percent_rank(MFC)) %>% 
                                 mutate(rank_LEC = percent_rank(-LEC)) %>% 
                                 mutate(cumulative_rank = (rank_MFC+rank_LEC)/2) %>% 
                                 arrange(-cumulative_rank) %>% 
                                 mutate(rank = 1:nrow(.)) %>% 
                                 mutate(rank_percent = rank/nrow(.)) %>% 
                                 #clustering, 25% as "higher", middle 25-75% of chemicals as "intermediate". 
                                 #bottom 25% as "lower"
                                 mutate(cluster = case_when(rank_percent < .25 ~ "higher",
                                                            rank_percent >= .25 & rank_percent < .75 ~ "intermediate",
                                                            rank_percent >= .75 ~ "lower", 
                                                            TRUE ~ "0")))) %>% 
  
  unnest(data_mut) %>% 
  select(-data)
  

#separate e2-up and p4-up chemicals 
e2_cluster <- pot_eff %>% 
  filter(hormone == "Estradiol")


p4_cluster <- pot_eff %>% 
  filter(hormone == "Progesterone")



## graph ac50 and ac10 values --------------------------------------------------

#lets gather ac50 information for those chemicals that were prioritized by 
#efficacy and potency

ac10info <- ac10 %>%  
  left_join(CASN %>% select(code, casn), by = c("X1" = "code" )) %>%  
  #select the columns reporting ac10 values for estradiol up and progesterone up
  select(X1, casn, contains("ESTRADIOL_up"), `CEETOX_H295R_PROG_up`) %>% 
  pivot_longer(cols = contains("up"), values_to = "AC10") %>% 
  #joing with hit call data. hitcall = 1 = increase hormone conc; 0 = no effect; 
  #-1 = decrease conc; NA = chemical not tested
  left_join(hitc %>% select (X1, contains("ESTRADIOL_up"), `CEETOX_H295R_PROG_up`) %>% 
            pivot_longer(cols = contains("up"), values_to = "hitcall")) %>% 
  #simplify column name
  mutate(hormone = case_when(name == "CEETOX_H295R_ESTRADIOL_up" ~ "Estradiol", 
                             name == "CEETOX_H295R_PROG_up" ~ "Progesterone", 
                             TRUE ~ "1000000")) %>% 
  rename("tcpl direction hit call" = "name" )
  



# combine ac10 data with E2-up/P4-up chemicals data
tcplval_sub <- data_long %>% 
  # log the AC50 and LEC values (AC10s are already logged) to compare with each other
  mutate(AC50 = log10(as.double(AC50)), 
         LEC = log10(as.double(LEC))) %>% 
  select(`Chemical name`, CASN, AC50, LEC, hormone) %>% 
  left_join(ac10info, by = c("CASN" =  "casn", "hormone" = "hormone")) %>% 
  #want to indicate ac10 values for chemicals with positive hits 
  #will be the same ones that have an ac50 value
  mutate(AC10_active = case_when(hitcall == 1 ~ AC10, 
                                 TRUE ~ NA_real_)) %>% 
  rename("AC50_active" = "AC50", 
         "AC10_all" = "AC10") 


#separate by E2 and p4 to graph

#E2 ---------------------------------------------

e2_tcpl <- tcplval_sub %>% 
  filter(hormone == "Estradiol") %>% 
  filter(CASN %in% e2_cluster$CASN)
  

long_e2_tcpl <- e2_tcpl %>% 
  pivot_longer(c(AC50_active, AC10_all, AC10_active)) %>% 
  rename("POD" = "name")


ggscatter(long_e2_tcpl, x = "LEC", y = "value",
          add = "reg.line",                         # Add regression line
          conf.int = TRUE,                          # Add confidence interval
          color = "POD", palette = "jco",           # Color by groups "name"
          shape = "POD"                             # Change point shape by groups "name"
)+
  stat_cor(aes(color = POD)) + 
  theme(panel.grid.major = element_line(colour="black", size = (1.0)),
        panel.grid.minor = element_line(size = (0.2), colour="grey"), 
        legend.title = element_blank()) + 
  labs(x = "log(LEC)", 
       y = "log(value)", 
       title = "E2-up chemicals") +
  geom_abline(intercept = 0, slope = 1, color="red", size = 1)

ggsave("./output/correlations_subchms_E2.png", height = 6, width = 7.5, units = "in")


cor.test(e2_tcpl$AC50_active, e2_tcpl$LEC, 
         method = "pearson")

cor.test(e2_tcpl$AC10_active, e2_tcpl$LEC, 
         method = "pearson")

cor.test(e2_tcpl$AC10_all, e2_tcpl$LEC, 
         method = "pearson")


#P4 ----------------------------------

p4_tcpl <- tcplval_sub %>% 
  filter(hormone == "Progesterone") %>% 
  filter(CASN %in% p4_cluster$CASN)


long_p4_tcpl <- p4_tcpl %>% 
  pivot_longer(c(AC50_active, AC10_all, AC10_active)) %>% 
  rename("POD" = "name")


ggscatter(long_p4_tcpl, x = "LEC", y = "value",
          add = "reg.line",                         # Add regression line
          conf.int = TRUE,                          # Add confidence interval
          color = "POD", palette = "jco",           # Color by groups "name"
          shape = "POD"                             # Change point shape by groups "name"
)+
  stat_cor(aes(color = POD)) + 
  theme(panel.grid.major = element_line(colour="black", size = (1.0)),
        panel.grid.minor = element_line(size = (0.2), colour="grey"), 
        legend.title = element_blank()) + 
  labs(x = "log(LEC)", 
       y = "log(value)", 
       title = "p4-up chemicals") +
  geom_abline(intercept = 0, slope = 1, color="red", size = 1)

ggsave("./output/correlations_subchms_P4_update.png", height = 6, width = 7.5, units = "in")


cor.test(p4_tcpl$AC50_active, p4_tcpl$LEC, 
         method = "pearson")

cor.test(p4_tcpl$AC10_active, p4_tcpl$LEC, 
         method = "pearson")

cor.test(p4_tcpl$AC10_all, p4_tcpl$LEC, 
         method = "pearson")



#-------------------------------------------------------------------------------
#update our table with the clusters 
#and ac50 values/chemical abbreviations

data_long_update <- data_long %>%  
  left_join(ac10info, by = c("CASN" = "casn", "hormone")) %>% 
  mutate(AC10_notlog = as.character(round((10^AC10), digits = 3))) %>% 
  mutate(AC10 = case_when(hitcall == 0 & !is.na(AC10_notlog) ~ paste(AC10_notlog, "(not active)", sep = " "), 
                          hitcall == 1 ~ AC10_notlog, 
                          TRUE ~ NA_character_)) %>% 
  select(-c(AC10_notlog, hitcall, `tcpl direction hit call`, X1)) %>% 
  left_join(pot_eff %>% rename ("Synonym" = "Chemical name") %>%  select(Synonym, CASN, cluster, rank, rank_percent, hormone)) %>% 
  mutate(Synonym = case_when(`Chemical name` != Synonym ~ Synonym,
                             TRUE ~ "")) %>% 
  select(hormone, `Chemical name`, Synonym, CASN, cluster, rank, rank_percent, MFC,
         AC50, AC10, hits, everything()) %>%  
  rename("Efficacy/potency" = "cluster") %>% 
  mutate(`Failed drug candidate` = case_when(grepl("UK-|CP-|Pharma|GSK|MK-|CI-|SSR|PD-|FR|CJ-|SAR|PD |AVE|PK|HMR" , `Chemical name`) ~ 1, 
                                             TRUE ~ 0)) %>% 
  rename("ER_Bioactivity" = "AUC.Bioactivity", 
         "ER_Agonist_AUC" = "AUC.Agonist", 
         "ER_Antagonist_AUC" = "AUC.Antagonist")


#Plot 1 ------------------------------------------------------------------------
#Create a plot with three variables showing potency, efficacy and exposure rate 
#
#looking at chemicals with a mahalanobis > 0, an mfc >= 1.5, a lac<= 33 
#also removing hormones in the pathway 

#want to remove chemicals with pharmaceutical prefixes as these are pharmaceutical drug candidates

#e2 ----------------------------------------------------------------------------


e2_plot <- e2_cluster %>% 
  #remove 19 failed drug candidates from graph 
  mutate(Delete = case_when(grepl("UK-|CP-|Pharma|GSK|MK-|CI-|SSR|PD-|FR|CJ-|SAR|PD |AVE|PK|HMR" , `Chemical name`) ~ "yes", 
                            TRUE ~ "no")) %>% 
  filter(Delete == "no") %>% 
  mutate(cluster = factor(cluster, levels = (c("higher", "intermediate", "lower"))))



e2_cluster_num <- e2_plot %>% 
  group_by(cluster) %>% 
  count()

#this code below will allow us to label the scatterplot with the names 
#of the chemicals; 
e2_text_zoom <- e2_plot %>% 
  mutate(parameters = case_when(MFC <= 3.84 & logLEC > -.7 ~ TRUE, 
                                TRUE ~ FALSE)) 


ggplot(e2_plot, aes(y = MFC, x = logLEC, shape = cluster, fill = cluster), color = "black") + 
  geom_point(position = "jitter") +
  theme_bw(base_size = 12) +
  #to show up in the plot that is not zoomed
  geom_label_repel(data= e2_text_zoom, aes(x= logLEC, y=MFC, label=`Chemical name`),
                   segment.color = "gray66", size = 2.3, label.padding = unit(0.1, "lines")) + 
  facet_zoom(ylim = c(1.3, 3.75), xlim = c(-.7, 2.3), zoom.data = parameters, 
             zoom.size = 3) + 
  labs(x = expression(paste("log (Lowest effective concentration (", mu, "M))")), 
       y = "Maximum E2 fold change", 
       title = "      E2-up chemicals") + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 10, vjust = -0.7),
        axis.title.y = element_text(size = 10, vjust = 1), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        plot.title = element_text(size = 16, face = 4), 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10)) +
  scale_shape_manual(labels = c("higher", "intermediate", "lower"), values = c(21, 22, 24), name = "Efficacy/\npotency") +
  scale_fill_manual(labels = c("higher", "intermediate", "lower"), values = c("#6baed6", "#c6dbef","#f7fbff"), name = "Efficacy/\npotency") + 
  guides(fill = FALSE) +
  guides(shape = guide_legend(override.aes = list(fill = c(a = "#6baed6", b="#c6dbef", c="#f7fbff"), size = 3)))


ggsave("./output/Figure2_E2.pdf", height = 8.5, width = 11, units = "in")




#p4 ------------------------------------------

p4_plot <- p4_cluster %>% 
  #remove 20 failed drug candidates from plot
  mutate(Delete = case_when(grepl("UK-|CP-|Pharma|GSK|MK-|CI-|SSR|PD-|FR|CJ-|SAR|PD |AVE|PK|HMR" , `Chemical name`) ~ "yes", 
                            TRUE ~ "no")) %>% 
  filter(Delete == "no") %>% 
  mutate(cluster = factor(cluster, levels = (c("higher", "intermediate", "lower")))) 

p4_cluster_num <- p4_plot  %>% 
  group_by(cluster) %>% 
  count()


p4_text_zoom <- p4_plot %>% 
  mutate(parameters = case_when(MFC <= 11.72 & logLEC > -.65 ~ TRUE, 
                                TRUE ~ FALSE)) 

ggplot(p4_plot, aes(y = MFC, x = logLEC, shape = cluster, fill = cluster), color = "black") + 
  geom_point(position = "jitter") +
  theme_bw(base_size = 12) +
  #to show up in the plot that is not zoomed
  geom_label_repel(data= p4_text_zoom, aes(x= logLEC, y=MFC, label=`Chemical name`),
                   segment.color = "gray66", size = 2.3, label.padding = unit(0.1, "lines")) + 
  facet_zoom(ylim = c(0, 11.72), xlim = c(-.5, 1.7), zoom.data = parameters, 
             zoom.size = 3) + 
  labs(x = expression(paste("log (Lowest effective concentration (", mu, "M))")), 
       y = "Maximum P4 fold change", 
       title = "      P4-up chemicals") + 
  theme(axis.title.x = element_text(size = 10, vjust = -0.7),
        axis.title.y = element_text(size = 10, vjust = 1), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        plot.title = element_text(size = 16, face = 4), 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10)) +
  scale_shape_manual(labels = c("higher", "intermediate", "lower"), values = c(21, 22, 24), name = "Efficacy/\npotency") +
  scale_fill_manual(labels = c("higher", "intermediate", "lower"), values = c("#6baed6", "#c6dbef","#f7fbff"), name = "Efficacy/\npotency") + 
  guides(fill = FALSE) +
  guides(shape = guide_legend(override.aes = list(fill = c(a = "#6baed6", b="#c6dbef", c="#f7fbff"), size = 3)))

ggsave("./output/Figure2_P4.pdf", height = 8.5, width = 11, units = "in")


  
################################################################################  
#exposure ----------------------------------------------------------------------
################################################################################

# median exposure rates ------------------------

exposure <- data_long_update %>%
  filter(!is.na(`Efficacy/potency`)) %>%
  select(`Chemical name`, CASN, hormone, `Efficacy/potency`, `Predicted median intake rate (mg/kg BW/day)`, 
         Consumer:`No exposure source data`) %>%
  pivot_wider(names_from = hormone, values_from = `Efficacy/potency`) %>%
  mutate(hormone = case_when (!is.na(Progesterone) & !is.na(Estradiol) == 1 ~ "Both",
                              !is.na(Progesterone) ~ "Progesterone",
                              !is.na(Estradiol) ~ "Estradiol",
                              TRUE ~ "fix")) 

medexp <- exposure %>%
  filter(`Predicted median intake rate (mg/kg BW/day)` >= 0.01) %>%
  mutate(`Efficacy/potency` = case_when(Estradiol == "higher" | Progesterone == "higher" ~ "Higher",
                                        Estradiol == "intermediate" | Progesterone == "intermediate" ~ "Intermediate",
                                        Estradiol == "lower" | Progesterone == "lower" ~ "Lower",
                                        TRUE ~ "Fix")) %>%
  mutate(`Efficacy/potency` = factor(`Efficacy/potency`, levels = c("Lower", "Intermediate", "Higher"))) %>% 
  mutate(`Predicted median intake rate (mg/kg BW/day)` = parse_number(`Predicted median intake rate (mg/kg BW/day)`)) %>% 
  mutate(facet = case_when(`Predicted median intake rate (mg/kg BW/day)` > 2.5 ~ "high", 
                           `Predicted median intake rate (mg/kg BW/day)` > .3 ~ "intermediate",
                           TRUE ~ "low"), 
         facet = factor(facet, levels = c("low", "intermediate", "high"))) %>% 
  rename("Hormone" = "hormone")

#functions to adjust plot scales
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

#plot
ggplot(medexp, aes(x = `Predicted median intake rate (mg/kg BW/day)`, 
                   y = reorder(`Chemical name`, `Predicted median intake rate (mg/kg BW/day)`), size = `Efficacy/potency`)) + 
  geom_point(aes(fill = Hormone, colour = Hormone), shape = 21, colour = "black") +
  facet_wrap_custom(~facet, scales = "free_x", scale_overrides = list(
    scale_override(1, scale_x_continuous(breaks=seq(0, .35, by = .05), expand = c(.1, 0))),
    scale_override(2, scale_x_continuous(breaks=seq(.5, 4, by = .5), expand = c(.1, 0))),
    scale_override(3, scale_x_continuous(breaks = seq(5, 60, by = 10), expand = c(.1, 0)))
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(y = "Chemical name", x = "Median exposure rate (mg/kg BW/day)") + 
  scale_fill_manual(values = c("#252525", "#969696", "#f7f7f7")) +
  scale_size_manual(values = c(2.5, 4, 6.5)) + 
  theme(strip.background = element_blank(), strip.text = element_blank(), 
        text = element_text(size = 16), 
        panel.grid.major = element_line(colour="grey90"), 
        panel.grid.minor = element_line(colour="grey100")) +
  guides(fill = guide_legend(override.aes = list(size=5)))


ggsave("./output/exposure_rate.png", height = 11, width = 10, units = "in")


#exposure sources -----------------------------------------------------

exposure <- data_long_update %>% 
  filter(!is.na(`Efficacy/potency`)) %>%
  select(`Chemical name`, CASN, hormone, `Efficacy/potency`, `Predicted median intake rate (mg/kg BW/day)`, 
         Consumer:`No exposure source data`, LEC) %>% 
  rename( "chnm" = "Chemical name") %>% 
  mutate(`Pharma.` = as.character(`Pharma.`), 
         `No exposure source data` = as.character(`No exposure source data`)) %>% 
  pivot_longer(Consumer:`No exposure source data`, names_to = "category") %>% 
  filter(value %in% c("1", "1*"))



plot_sources <- exposure %>% 
  select(chnm, CASN, category, hormone) %>% 
  arrange(hormone) %>% 
  group_by(chnm, CASN, category) %>% 
  summarise(hormones = paste(hormone, collapse = ", ")) %>% 
  group_by(category, hormones) %>% 
  count() %>% ungroup() %>% 
  mutate(categories = factor(category, levels = c("Consumer", "Diet", "Industrial", "Pest.", "Pharma.", "No exposure source data"), 
                             labels = c("Consumer \nproducts", "Diet", "Industrial", "Pesticides", "Pharmaceuticals", "No data")),
         Hormone = gsub("Estradiol, Progesterone", "Both", hormones)) 


#plot graph
ggplot(plot_sources, aes(x= reorder(categories, n), fill = Hormone, y = n)) +
  geom_bar_pattern(aes(fill= Hormone,
                       pattern_angle = Hormone,
                       pattern_density = Hormone,
                       pattern_fill = Hormone),
                   pattern_spacing = .01,
                   pattern_alpha = .3,
                   position = "stack", stat = "identity",
                   pattern_colour  = "white", #color of pattern border
                   color = "black", #border of bars color
                   size = .01) +
  xlab("Source of exposure") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0.5)) +
  theme_bw(base_size = 11 ) +
  scale_fill_manual(values = c("#252525", "#969696", "#f7f7f7")) +
  scale_pattern_fill_manual(values = c('blue', "#252525", '#969696')) +
  labs(x = "", y = "Number of chemicals") + 
  theme(strip.text.x = element_text(
    size = 12, color = "black"), 
    axis.title.x = element_text(margin = margin(t = 15 , r = 0, b = 0, l = 0)))
ggsave("./output/Figure3_broad.png", width = 7.25, height = 3, units = "in")



plot_specific <- chem_sources_graph_specific %>%  
  select(-source) %>% 
  right_join(exposure %>% group_by(CASN, hormone) %>% slice(1), by = c("casrn" ="CASN")) %>%  
  select(chnm, casrn, value.x, hormone) %>% 
  group_by(chnm, casrn, value.x) %>% 
  summarise(hormones = paste(hormone, collapse = ", ")) %>% 
  filter(!is.na(value.x)) %>% 
  group_by(hormones, value.x) %>% 
  count() %>% 
  mutate(Hormone = gsub("Estradiol, Progesterone", "Both", hormones)) 



ggplot(plot_specific, aes(x= reorder(value.x, n), fill = Hormone, y = n)) +
  geom_bar_pattern(aes(fill= Hormone,
                       pattern_angle = Hormone,
                       pattern_density = Hormone,
                       pattern_fill = Hormone),
                   pattern_spacing = .01,
                   pattern_alpha = .3,
                   position = "stack", stat = "identity",
                   pattern_colour  = "white", #color of pattern border
                   color = "black", #border of bars color
                   size = .01) +
  xlab("Source of exposure") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0.5)) +
  theme_bw(base_size = 11 ) +
  scale_fill_manual(values = c("#252525", "#969696", "#f7f7f7")) +
  scale_pattern_fill_manual(values = c('blue', "#252525", '#969696')) +
  labs(x = "", y = "Number of chemicals") +
  theme(
    strip.text.x = element_text(
      size = 11, color = "black"),
    axis.title.x = element_text(margin = margin(t = 15 , r = 0, b = 0, l = 0)))
ggsave("./output/Figure3_specific.png", width = 7.25, height = 5, units = "in")



################################################################################
## Prioritization --------------------------------------------------------------
################################################################################

#what drugs are currently ingredients in fda approved drug products? -----------

drugs <- fdadrugs %>% 
  mutate(Ingredient = tolower(Ingredient)) %>% 
  filter(!Type == "DISCN") %>%  
  mutate(row = 1:n())
  
#some drug products contain multiple drug ingredients, we are going to separate them 
drug_split <- as_tibble(str_split(drugs$Ingredient, pattern = "; ", n= Inf, simplify = TRUE)) %>% 
  mutate(row = 1:n()) %>% 
  pivot_longer(cols = contains("V", ignore.case = FALSE), values_to = "ingredients") %>%
  filter(!ingredients == "") %>% 
  group_by(ingredients) %>%  count()


drug_numbers <- data_long_update %>% 
  #remove chem duplicate
  group_by(CASN) %>%
  slice(1) %>% 
  select(`Chemical name`, CASN, Pharma.) %>%
  mutate(`Chemical name2` = tolower(`Chemical name`)) %>% 
  left_join(drug_split, by = c("Chemical name2" = "ingredients")) 

#what pesticide active ingredients are currently in epa approved products? -----

pest_prods <- data_long_update %>% 
  group_by(CASN) %>%
  slice(1) %>% 
  select(`Chemical name`, CASN) %>% 
  left_join(status) %>% 
  mutate(products = case_when(is.na(`Active products`) & is.na(`Not active products`) ~ "",
                              !is.na(`Active products`) ~ "yes", 
                              !is.na(`Not active products`) ~ "no", 
                              TRUE ~ "fix")) %>% 
  filter(`products` == "yes")
  

#Prioritization ----------------------------------------------------------------

#what are the pesticides which have active products in the US? 

data_long_update2 <- data_long_update %>% 
  mutate(`Chemical name2` = tolower(`Chemical name`)) %>% 
  mutate(`FDA approved drug products` = case_when(`Chemical name2` %in% c(drug_split$ingredients) ~ "yes", 
                           TRUE ~ "no")) %>% 
  select(-`Chemical name2`) %>% 
  mutate(`EPA approved pesticide products` = case_when(CASN %in% pest_prods$CASN ~ "yes", 
                                                       TRUE ~ "no")) %>% 
  # select(`Chemical name`:LEC, Consumer:Pest., Nhanes, contains("Summary")) %>% \
  mutate(Prioritize = case_when(!`Efficacy/potency` == "" &
                                  ((Consumer == 1 | Nhanes == 1 | 
                                     `EPA approved pesticide products` == "yes" |
                                     `FDA approved drug products` == "yes")) ~ "yes", 
                                TRUE ~ "no")) %>% 
  rename ("Prioritize (based on efficacy, potency, and exposure potential)" = "Prioritize") %>% 
  select(`Chemical name`, CASN, everything()) 



#Finalizing the table for publication; separating E2 and P4  ----------------------------------------------


data_long_final <- data_long_update2 %>% 
  rename("Efficacy/potency rank" = "rank", 
         "Efficacy/potency percentile" = "rank_percent", 
         "Hitcall per chemical concentration" = "hits",
         "NHANES" = "Nhanes", 
         "Effect on mammary gland" = "Effects on mammary gland", 
         "Dev/repro in vivo effect level < 100 mg/kg-day (ToxValDb)" = "Invivo effect < 100 mg/kg-day (ToxValDb)",
         "Dev/repro in vivo no effect level >= 100 mg/kg-day (ToxValDb)" = "No invivo effect >= 100 mg/kg-day (ToxValDb)",
         "Other (ToxValDb re dev/repro in vivo effect)" = "Other (ToxValDb)",
         "Carcinogenicity assessment (Summary)" = "Carcinogenicity Assessment (Summary)",
         "Carcinogenicity (ToxValDb)" = "Carcinogencity (ToxValDb)", 
         "Chemical of concern (based on hazard and exposure potential)" = "Prioritize (based on efficacy, potency, and exposure potential)", 
         "Developmental or reproductive toxicity (Summary)" = "Developmental or reproducive toxicity (Summary)", 
         "ER Bioactivity" = "ER_Bioactivity", 
         "ER Agonist AUC" = "ER_Agonist_AUC", 
         "ER Antagonist AUC" = "ER_Antagonist_AUC") %>% 
  select(hormone, `Chemical name`, Synonym, CASN, `MFC`, `LEC`, AC50, AC10, `Adj.maxmMd`, 
         `Efficacy/potency`, 
         `Efficacy/potency percentile`,
         `Efficacy/potency rank`, 
         `Min. tested conc.`, 
         `Max. tested conc.`, 
         `Hitcall per chemical concentration`,
          Estradiol, Progesterone,
         `Predicted median intake rate (mg/kg BW/day)`, 
          Consumer, Diet, Industrial, Pharma., Pest., `No exposure source data`, 
         `Failed drug candidate`,
         `Additional exposure sources`, NHANES, `Effect on mammary gland`, 
         `Carcinogenicity (ToxValDb)`, 
         `Dev/repro in vivo effect level < 100 mg/kg-day (ToxValDb)`, 
         `Dev/repro in vivo no effect level >= 100 mg/kg-day (ToxValDb)`,
         `Other (ToxValDb re dev/repro in vivo effect)`,
         `Prop65`, `Carcinogenicity assessment (Summary)`, 
         `Developmental or reproductive toxicity (Summary)`, 
         `FDA approved drug products`, `EPA approved pesticide products`, 
         `Chemical of concern (based on hazard and exposure potential)`, 
         `ER Bioactivity`, `ER Agonist AUC`, `ER Antagonist AUC`) %>% 
  mutate(`Efficacy/potency` = case_when(CASN %in% c("53-16-7", "474-86-2", "50-27-1", "57-85-2", "57-83-0", "57-63-6", "58-18-4", 	
                                         "521-18-6", "53-41-8", "68-96-2", "63-05-8", "53-43-0", "50-28-2", "57-91-0") ~ "hormone substrate", 
                                        TRUE ~ `Efficacy/potency`), 
         `Efficacy/potency` = case_when(is.na(`Efficacy/potency`) ~ "borderline active",
                                        TRUE ~ `Efficacy/potency`))

#E2 table -------------------------------------------------------


final_E2 <- data_long_final %>% 
  filter(hormone == "Estradiol") %>% 
  select(-c(hormone, Estradiol)) %>% 
  rename("Direction of P4 synthesis" = "Progesterone")

write_csv(final_E2 %>% mutate(CASN = paste0("'", CASN)), "./output/E2_supp.csv") 

#P4 table -------------------------------------------------------

final_P4 <- data_long_final %>% 
  filter(hormone == "Progesterone") %>% 
  select(-c(hormone, Progesterone)) %>% 
  rename("Direction of E2 synthesis" = "Estradiol")

write_csv(final_P4 %>% mutate(CASN = paste0("'", CASN)), "./output/P4_supp.csv") 
