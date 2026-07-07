####R code for 'Integration of metabolomics into early phase clinical trials from 2010 to 2025'####

#### 1.) PREPARING WORKSPACE####
### 1a. Install required tidyverse and janitor packages
install.packages("tidyverse")
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
install.packages("janitor") 
library(janitor)
library(lubridate)
library(purrr)
library(tidyr)
library(tibble)
library(ggrepel)

### 1b. Set appropriate working directory using setwd(xxxxxxx)

### 1c. Upload ICTRP_Full: the csv file of all studies in the ICTRP database. This includes 1173018 obs and 63 variables
ICTRP1 <- read_csv(
  "ICTRP_Full_R.csv",
  na = c("", "NA"),
  quote = "\"",   
  col_types = cols(.default = "c"),  
  locale = locale(encoding = "UTF-8"),  
  trim_ws = TRUE, col_names = FALSE)

#### 2.) FORMATTING DATAFRAME  
### 2a. Add appropriate column names 
colnames(ICTRP1) <- c(
  "TrialID","last_update","SecondaryIDs","public_title","Scientific_title","url",
  "Public_Contact_Firstname","Public_Contact_Lastname","Public_Contact_Address",
  "Public_Contact_Email","Public_Contact_Tel","Public_Contact_Affiliation",
  "Scientific_Contact_Firstname","Scientific_Contact_Lastname",
  "Scientific_Contact_Address","Scientific_Contact_Email",
  "Scientific_Contact_Tel","Scientific_Contact_Affiliation",
  "study_type","study_design","phase","Date_registration","Date_enrolment",
  "Target_size","Recruitment_status","Primary_sponsor","Secondary_sponsors",
  "Source_Support","Countries","Conditions","Interventions","Agemin","Agemax",
  "Gender","Inclusion_Criteria","Exclusion_Criteria","Primary_Outcome",
  "Secondary_Outcomes","Bridging_flag","Bridged_type","Childs","type_enrolment",
  "Prospective_registration","results_actual_enrolment","results_url_link",
  "results_summary","results_date_posted","results_date_first_publication",
  "results_baseline_char","results_participant_flow","results_adverse_events",
  "results_outcome_measures","results_url_protocol","results_IPD_plan",
  "results_IPD_description","results_date_completed","results_yes_no",
  "Ethics_Status","Ethics_Approval_Date","Ethics_Contact_Name",
  "Ethics_Contact_Address","Ethics_Contact_Phone","Ethics_Contact_Email"
)
### 2b. Clean column names so I only have lower case
ICTRP1 <- ICTRP1 %>% clean_names()

### 2c. Select only the columns I'm interested in
ICTRP2 <- ICTRP1 %>% select(trial_id, public_title, scientific_title, secondary_i_ds, study_type, primary_sponsor, secondary_sponsors, primary_outcome, secondary_outcomes, date_enrolment, date_registration, study_design, study_type, study_design, phase, countries, conditions, interventions)

### 2d. Renaming column name to make neater 
ICTRP2 <- ICTRP2 %>%
  rename(secondary_ids= secondary_i_ds)

## 2e. Adding row numbers to make it easier to relate back to original csv file for cross referencing (and then putting this as the first column with relocate()) 
ICTRP2 <- ICTRP2 %>% mutate (row_number= row_number()) %>% relocate(row_number, .before=1)
dim(ICTRP2)

### 2f. Create new column called 'year' that picks the appropriate year from the unstructured 'date_enrolment' column. 

ICTRP3 <- ICTRP2 %>%
  mutate(
    date_enrolment_clean = str_trim(str_replace_all(date_enrolment, '"', '')),
    year = coalesce(
      year(parse_date_time(
        date_enrolment_clean,
        orders = c( "Y-m-d","d-m-Y","m-d-Y","d/m/Y","m/d/Y",
                    "d-m-y","m-d-y","d/m/y","m/d/y",
                    "B d, Y","b d, Y","d B Y","d b Y","B Y","b Y","Y"),
        quiet = TRUE)),
      as.numeric(str_extract(date_enrolment_clean, "\\b(19|20)\\d{2}\\b")))) %>% select(-date_enrolment_clean)

## 2g. To check this is done correctly 
ICTRP3 %>%filter(is.na(year)) %>% select(date_enrolment)
##This tells me 8,810 rows (about 0.6% of all) have missing data for year 

ICTRP3 %>% filter(is.na(year)) %>% distinct(date_enrolment) %>% head(20)
##This confirms missing data is the problem rather than misreading/inappropriate code.

#### 3.) FILTERING ICTRP_Full for early phase clinical trials from 2010-2025 (ICTRP-All)

### 3a. Select only INTERVENTIONAL studies (i.e clinical trials)-> ICTRP4
ICTRP4 <- ICTRP3 %>% filter(str_detect(study_type, regex("intervention", ignore_case=TRUE)))
str(ICTRP4)
#This produces a dataframe with 911,889 rows of all interventional studies (across all years and all phases) 


### 3b. Select only early PHASE studies (Phase 0, 1, 2 or a combination of these)
## First the code looks for if there is a 'Phase X: Yes or no' pattern to ensure no false positives. 
## Second the code searches for any mention of 0,1/I,2/II and categorises accordingly. 

ICTRP5 <- ICTRP4 %>% mutate(phase_lower = str_to_lower(phase),
has_phase_yesno = str_detect(phase_lower, "phase\\s*(i{1,3}|iv|\\d+).*?\\b(yes|no)\\b"),
phase_clean = case_when(has_phase_yesno ~ map_chr(phase_lower, function(x) {yes_matches <- str_match_all(x, "phase\\s*(i{1,3}|iv|\\d+)\\)?\\s*:?\\s*yes" )[[1]][, 2]
 if (length(yes_matches) == 0) { NA_character_ } else {paste(yes_matches, collapse = "/")} }),
           TRUE ~ phase_lower) %>%
           str_replace_all("early phase 1", "phase 1") %>%
           str_replace_all("\\biii\\b", "3") %>%
           str_replace_all("\\bii\\b", "2") %>%
           str_replace_all("\\biv\\b", "4") %>%
           str_replace_all("\\bi\\b", "1") %>%
           str_replace_all("phase", " ") %>%
           str_replace_all("([0-9])[a-z]", "\\1") %>%
           str_replace_all("[,;\\-]", "/") %>%
           str_replace_all("\\s*/\\s*", "/") %>%
           str_replace_all("[^0-9/ ]", " ") %>%
           str_squish(),
         
         phase_group = case_when(
           str_detect(phase_clean, "\\b0/1\\b|\\b1/0\\b") ~ "Phase 0/1",
           str_detect(phase_clean, "\\b1/2\\b|\\b2/1\\b") ~ "Phase 1/2",
           str_detect(phase_clean, "\\b2/3\\b|\\b3/2\\b") ~ "Phase 2/3",
           
           str_detect(phase_clean, "^0$|\\b0\\b") ~ "Phase 0",
           str_detect(phase_clean, "^1$|\\b1\\b") ~ "Phase 1",
           str_detect(phase_clean, "^2$|\\b2\\b") ~ "Phase 2",
           
           TRUE ~ NA_character_) )

ICTRP6 <- ICTRP5 %>%
  filter(phase_group %in% c(
    "Phase 0",
    "Phase 0/1",
    "Phase 1",
    "Phase 1/2",
    "Phase 2"))

#3c. To check this is done correctly 
ICTRP5 %>%filter(has_phase_yesno) %>% select(phase, phase_clean, phase_group) %>% head(20)
#This provided ICTRP6 240,463 rows for phase 0,1 or 2 studies (+ combos)

###3d. Select only YEARS 2010-2025 making ICTRP7 (this leaves 196,244 early phase interventional studies)
ICTRP7 <- ICTRP6 %>% filter(year>=2010, year <=2025) 

#### 4.) FURTHER FORMATTING: Creating a new column called registry_name using each trial_id

##4a. Full registry name lookup
registry_lookup <- tribble(
  ~registry_code, ~registry_name,
  "ANZCTR",   "Australian New Zealand Clinical Trials Registry",
  "REBEC",    "Brazilian Clinical Trials Registry",
  "ChiCTR",   "Chinese Clinical Trial Register",
  "NCT",      "Clinical Trials.gov (USA)",
  "CTRI",     "Clinical Trials Registry India",
  "CRIS",     "Clinical research information centre (Republic of Korea)",
  "RPCEC",    "Cuban Public Registry of Clinical Trials",
  "EUCTR",    "EU Clinical Trials Register",
  "CTIS",     "Clinical Trials Information System",
  "DRKS",     "German Clinical Trials Register",
  "IRCT",     "Iranian Registry of Clinical Trials",
  "ISRCTN",   "ISRCTN (UK)",
  "JPRN",     "Japan Primary Registries Network",
  "PACTR",    "Pan African Clinical Trial Register",
  "REPEC",    "Peruvian Clinical Trials Register",
  "SLCTR",    "Sri Lanka Clinical Trials Register",
  "TCTR",     "Thai Clinical Trials Register",
  "OMON",     "The Netherlands National Trial Register",
  "LBCTR",    "Lebanese Clinical Trials Registry",
  "ITMCTR",   "International Traditional Medicine Clinical Trial Registry")

## 4b. Regex lookup table for TrialID prefixes. Note PER- has been included as part of peruvian registry
registry_patterns <- tribble(
  ~pattern,                   ~registry_code,
  "^ACTRN",                   "ANZCTR",
  "^ANZCTR",                  "ANZCTR",
  "^RBR[-]",                  "REBEC",
  "^REBEC",                   "REBEC",
  "^CHICTR[-](TRC|IIR|IPR)",  "ChiCTR",
  "^CHICTR[-]",               "ChiCTR",
  "^CHICTR",                  "ChiCTR",
  "^NCT",                     "NCT",
  "^CTRI",                    "CTRI",
  "^CRIS",                    "CRIS",
  "^KCT",                     "CRIS",
  "^RPCEC",                   "RPCEC",
  "^EUCTR",                   "EUCTR",
  "^CTIS",                    "CTIS",
  "^DRKS",                    "DRKS",
  "^IRCT",                    "IRCT",
  "^ISRCTN",                  "ISRCTN",
  "^JPRN[-](JRCT|JRCTS|UMIN)","JPRN",
  "^JPRN",                    "JPRN",
  "^PACTR",                   "PACTR",
  "^PER[-]",                  "REPEC",
  "^REPEC",                   "REPEC",
  "^SLCTR",                   "SLCTR",
  "^TCTR",                    "TCTR",
  "^NL[-]OMON",               "OMON",
  "^OMON",                    "OMON",
  "^LBCTR",                   "LBCTR",
  "^ITMCTR",                  "ITMCTR")


##4c. Function to assign registry code from TrialID
get_registry_code <- function(trial_id, registry_patterns) {trial_id <- str_to_upper(trial_id)

map_chr(trial_id, function(x) {hit <- registry_patterns %>% filter(str_detect(x, pattern))

if (nrow(hit) == 0) {NA_character_} else {hit$registry_code[1] } })}

##4d. Apply to the ICTRP dataframe
ICTRP8 <- ICTRP7 %>% mutate(registry_code = get_registry_code(trial_id, registry_patterns)) %>% left_join(registry_lookup, by = "registry_code")


## 4e. Checking how many studies for each registry 
ICTRP8 %>% count(registry_code, registry_name, sort = TRUE)

## 4f.To check there are no trials left without an associated registry_code
which(is.na(ICTRP8$registry_code)) #correct, no studies left 


#### 5.) REMOVING REPLICATED STUDIES: to remove studies registered with multiple national regitries and keep the one with the earliest registration date

### 5a. Take ICTRP8 and make a new column for a standardised/readable registration date label (known as date_regitration_std) -> label the new dataframe as ICTRP9

ICTRP9 <- ICTRP8 %>% mutate(date_registration_clean = str_trim(str_replace_all(date_registration, '"', '')), date_registration_std = case_when(
                              
                              # Full dates (many formats)
                              !is.na(parse_date_time(date_registration_clean,
                                                     orders = c("Y-m-d","d-m-Y","m-d-Y","d/m/Y","m/d/Y",
                                                                "d-m-y","m-d-y","d/m/y","m/d/y",
                                                                "B d, Y","b d, Y","d B Y","d b Y"),
                                                     quiet = TRUE)) ~ as.Date(parse_date_time(
                                                       date_registration_clean,
                                                       orders = c("Y-m-d","d-m-Y","m-d-Y","d/m/Y","m/d/Y",
                                                                  "d-m-y","m-d-y","d/m/y","m/d/y",
                                                                  "B d, Y","b d, Y","d B Y","d b Y"),
                                                       quiet = TRUE))))



### 5b.  identify all the replicated rows and label as dup_rows.
#This means all the rows where BOTH public title and scientific title are identical while ignoring those with empty/missing public and scientific titles (as well as .,-,/, 'n/a' etc.) 
clean_title <- function(x) {
  x %>%
    str_replace_all('"', '') %>%
    str_to_lower() %>%
    str_squish() %>%
    na_if("") %>%
    if_else(str_detect(., "^(n/?a|none|unknown)$"), NA_character_, .) %>% 
    if_else(str_detect(., "^[ .-]*$"), NA_character_, .) %>%
    str_replace_all("\\s+", "")}

dup_rows <- ICTRP9 %>%
  mutate(public_title_norm = clean_title(public_title),scientific_title_norm = clean_title(scientific_title)) %>%
  filter(!is.na(public_title_norm), !is.na(scientific_title_norm)) %>%
  group_by(public_title_norm, scientific_title_norm) %>%
  filter(n() > 1) %>%
  arrange(public_title_norm, scientific_title_norm, year) %>%
  ungroup()
#dup_rows -> there are 19988 replicated rows 

## 5c. Find out how many replicate groups in total (because one study can be in multiple registries) in dup_rows  
dup_rows %>% distinct(public_title_norm, scientific_title_norm) %>% nrow() 
#6785 discrete groups
# If there are 6785 discrete groups and 19988 replicate rows...
## 19988-6785= 13203 rows should be removed. 


### 5d. Select the replicated rows with the earliest registration date (as per column date_registration_std)
## I have used cleaned title names and ignore the rows where the public or scientific title are missing. 
dup_rows_refined <- ICTRP9 %>%
  mutate(
    public_title_norm = clean_title(public_title),
    scientific_title_norm = clean_title(scientific_title)
  ) %>%
  mutate(
    public_title_norm = str_replace_all(public_title_norm, "\\s+", ""),
    scientific_title_norm = str_replace_all(scientific_title_norm, "\\s+", "")
  ) %>%
  filter(!is.na(public_title_norm), !is.na(scientific_title_norm)) %>%
  group_by(public_title_norm, scientific_title_norm) %>%
  filter(n() > 1) %>%
  arrange(is.na(date_registration_std), date_registration_std, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

### 5e. Remove dup_rows from ICTRP9 and then add back dup_rows_refined
ICTRP9_nodups <- ICTRP9 %>% filter(!row_number %in% dup_rows$row_number)
ICTRP10 <- bind_rows(ICTRP9_nodups, dup_rows_refined)
#ICTRP10: 183,041 early phase interventional studies with replicates removed and the row with the earliest registration date kept.


#### 6.) FINDING METABOLOMIC STUDIES

###6a. Build my search terms and look for these character strings in select columns. (These two pieces of code take a couple of hours to run.)

search_terms <- c("metabolom\\w*",
                  "metabonom\\w*",
                  "metabolite profil\\w*",
                  "metabolite fingerprint\\w*",
                  "small molecule profiling",
                  "pharmacometabolom\\w*",
                  "lipidom\\w*",
                  "glycom\\w*",
                  "steroidom\\w*",
                  "metaboli\\w* phenotyping")

search_cols <- c("public_title",
                 "scientific_title",
                 "primary_outcome",
                 "secondary_outcomes",
                 "conditions",
                 "interventions")

ICTRP_metabolomics <- ICTRP10 %>%
  rowwise() %>%
  mutate(
    found_in = paste(
      search_cols[
        sapply(search_cols, function(col) {
          str_detect(
            coalesce(get(col), ""),
            regex(paste(search_terms, collapse = "|"), ignore_case = TRUE)) })],collapse = "; "),
    
    matched_terms = paste(
      search_terms[
        sapply(search_terms, function(term) {
          any(
            sapply(search_cols, function(col) {
              str_detect(
                coalesce(get(col), ""),
                regex(term, ignore_case = TRUE))}) )})], collapse = "; ")) %>%
  ungroup() %>%
  filter(found_in != "")  


###6b. To see how often each field matches 
ICTRP_metabolomics %>% count(found_in, sort = TRUE) %>% print(n=Inf)
ICTRP_metabolomics %>% count(matched_terms, sort = TRUE)

#### 7.) WRITE ICTRP_Met and ICTRP_All into files to save work for subsequent analysis
### 7a. remove unnecessary columns 
ICTRP_Met <- ICTRP_metabolomics %>%
  select(
    -public_title_norm,
    -scientific_title_norm,
    -phase_clean,
    -has_phase_yesno,
    -phase_lower,
    -secondary_sponsors)

write_csv(ICTRP_Met, "ICTRP_Met.csv")


ICTRP_All <- ICTRP10 %>%
  select(
    -public_title_norm,
    -scientific_title_norm,
    -phase_clean,
    -has_phase_yesno,
    -phase_lower,
    -secondary_sponsors)

write_csv(ICTRP_All, "ICTRP_All.csv")



#####ANALYSIS#####

#### 8. Breakdown of studies by country (ICTRP_ALl)

### 8a. The country names are not structured and inputted differently. Consequently, I have produced a separate CSV file as a look-up table for countries with inconsistent naming to make them readable.  
#This helps group the ICTRP_All into 187 countries and 1 Missing data categories. Missing data includes vague denominations for example entire continents rather than the specific country or unrecognisable country names.

country_lookup <- tribble(
  ~raw_country, ~country_standard,
  
  "Czechia", "Czech Republic",
  "Czech", "Czech Republic",
  
  "Russian Federation", "Russia",
  "Russian Fed.", "Russia",
  
  "Iran, Islamic Republic of", "Iran",
  "Iran (Islamic Republic Of)", "Iran",
  
  "Hong Kong", "China",
  
  "United States of America", "United States",
  "Us", "United States",
  "United States Minor Outlying Islands", "United States",
  "The United States of America", "United States",
  "The United States", "United States",
  "The Us", "United States",
  "United State", "United States",
  "Unitedstates Of America", "United States",
  "Perto Rico", "United States",
  "United Status", "United States",
  "U.s,A", "United States",
  "United States,", "United States",
  "Unitedstates Ofamerica", "United States",
  "United States Of Americ", "United States",
  "Usa (Including Puerto Rico)", "United States",
  "Puerto Rico", "United States",
  "Guam", "United States",
  "United Staes", "United States",
  "Usa", "United States",
  "U.s.a", "United States",
  "U.s.", "United States",
  "U.s.a.", "United States",
  "American Samoa", "United States",
  
  "Holland, The Netherlands", "Holland",
  "Netherland", "Holland",
  "Nederland", "Holland",
  "Republic Of Netherlands", "Holland",
  "Netherlands", "Holland",
  "Benelux", "Holland",
  "Netherlands Kingdom Of The", "Holland",
  "Nertherlands", "Holland",
  "Thenetherlands", "Holland",
  "The Netherlands", "Holland",
  "Netherlands Antilles", "Holland",
  "Aruba", "Holland",
  
  "Moldova, Republic of", "Moldova",
  "Republic Of Moldova", "Moldova",
  "Moldova Republic Of", "Moldova",
  "Modalvia", "Moldova",
  
  "Western Sahara", "Morocco", 
  
  "Turkey (Türkiye)", "Turkey",
  "Turkey (Turkiye)", "Turkey",
  "Türkiye", "Turkey",
  
  "United Republic of Tanzania", "Tanzania",
  "Tanzania, United Republic Of", "Tanzania",
  
  "Korea, Democratic People’s Republic of", "North Korea",
  "Democratic People's Republic Of Korea", "North Korea",
  "Korea, Democratic People's Republic Of", "North Korea",
  "Korea North", "North Korea",
  "Democratic People Republic Of Korea", "North Korea",
  "Democratic People's Republic Of", "North Korea",
  
  "Syrian Arab Republic", "Syria",
  "Phillipines", "Philippines", 
  
  "None", "Missing data",
  "??", "Missing data",
  "À¹„À¸—À¸¢", "Missing data",
  "Not Applicable", "Missing data",
  "No", "Missing data",
  "Other", "Missing data",
  "-", "Missing data",
  "Non", "Missing data",
  "Not Stated", "Missing data",
  "The People's Republic Of", "Missing data",
  "Many", "Missing data",
  "Noe", "Missing data",
  "None.", "Missing data",
  "Nore", "Missing data",
  "Nothing", "Missing data",
  "Refer To The Section, Other.", "Missing data",
  "A", "Missing data",
  "N", "Missing data",
  "North America, Europe, Asia (Except Japan), Oceania, South America, Africa", "Missing data",
  "Apac", "Missing data",
  "Part A: Usa Part B C: Asia, North America, Eu", "Missing data",
  "East Asia(Taiwan, Hong Kong, Singapore And Other)", "Missing data",
  "Asia Exept Japan", "Missing data",
  "Asia", "Missing data",
  "Emea", "Missing data",
  "Expansion Phase :", "Missing data",
  "Null", "Missing data",
  "Asia/ Europe/ North America/ Oceania/ European Union / South America", "Missing data",
  "Japan,Asia(Except Japan)", "Missing data",
  "Japan,Asia(Except Japan),North America", "Missing data",
  "Dose Escalation Phase :", "Missing data",
  "	Japan,Asia(Except Japan),Europe", "Missing data",
  "Japan,Europe", "Missing data",
  "Japan,North America", "Missing data",
  "Japan,North America,South America", "Missing data",
  "Japan,South America", "Missing data",
  " America", "Missing data",
  "European Union", "Missing data",
  "Africa", "Missing data",
  "Asia(Except Japan),Africa", "Missing data",
  "South America", "Missing data",
  "North America", "Missing data",
  "Europe", "Missing data",
  "Oceania", "Missing data",
  "Nate", "Missing data",
  "Asia(Except Japan)", "Missing data",
  "Asia Except Japan", "Missing data",
  "Japan, Asia(Exept Japan)", "Missing data",
  
  "Congo, The Democratic Republic Of The", "Congo",
  "Democratic Republic Of The Congo", "Congo",
  "Congo - Kinshasa", "Congo",
  "Republic Of The Congo", "Congo",
  
  "Lao People’s Democratic Republic", "Laos",
  
  "Chnia", "China",
  "P.r. China", "China",
  "P.r.china", "China",
  "Beijing", "China",
  "Cnina", "China",
  "Chain", "China",
  "Chine", "China",
  "Guangdong", "China",
  "Hina", "China",
  "People’s Republic Of China", "China",
  "People's Republic Of China", "China",
  "Shanghai", "China",
  "Republic of China", "China",
  "Anhui", "China",
  "Chaina", "China",
  "Chengdu", "China",
  "Chiina", "China",
  "Chin", "China",
  "Chin’a", "China",
  "China Hong Kong", "China",
  "China Mainland (Prc)", "China",
  "China·", "China",
  "Ching", "China",
  "Chinia", "China",
  "Chins", "China",
  "Chna?Cna", "China",
  "Chia", "China",
  "Cina", "China",
  "Chna", "China",
  "Can", "China",
  "Guangxi Zhuang Autonomous Region", "China",
  "Hebei", "China",
  "Hksar", "China",
  "Hong Kong Sar, China", "China",
  "Hong Kong, China", "China",
  "Hong Kong,China", "China",
  "Maternal and Child Health Hospital of Changxing Couty", "China",
  "Pepole's Republic Of China", "China",
  "Rpc", "China",
  "Tianjin", "China",
  "Zhong Guo", "China",
  "Zhongguo", "China",
  "Hongkong", "China",
  "Macau", "China",
  "Chian", "China",
  "Chinese", "China",
  "The First Affiliated Hospital Of Hunan University Of Traditional Chinese Medicine", "China",
  "The People's Republic Of China", "China",
  "Pr China", "China",
  "Cna", "China",
  
  "Cote Divoire", "Côte D'Ivoire",
  "Cote D’Ivoire", "Côte D'Ivoire",
  "Côte D'ivoire", "Côte D'Ivoire",
  
  "Macedonia, The Former Yugoslav Republic of", "Macedonia",
  "North Macedonia", "Macedonia",
  "The Former Yugoslav Rep Of Macedonia", "Macedonia",
  "The Former Yugoslav Republic Of Macedonia", "Macedonia",
  
  "Viet Nam", "Vietnam",
  
  "Korea, South", "South Korea",
  "Korea South", "South Korea",
  "Korea Republic Of", "South Korea",
  "Korea (The Republic Of)", "South Korea",
  "Sounth Korea", "South Korea",
  "Korea (Republic Of)", "South Korea",
  "Sourth Korea", "South Korea",
  "Korea,Republic Of", "South Korea",
  "Suoth Korea", "South Korea",
  "S Korea", "South Korea",
  "Sk", "South Korea",
  "Sout Korea", "South Korea",
  "Southkorea", "South Korea",
  "S. Korea", "South Korea",
  "Rep. Of Korea", "South Korea",
  "Korea", "South Korea",
  "Koria", "South Korea",
  "Korea, Republic Of", "South Korea",
  "Republic Of Korea,", "South Korea",
  "Republic Of Korea", "South Korea",
  "Sount Korea", "South Korea",
  "Rupublic Of Korea", "South Korea",
  
  "United Kindgdom", "United Kingdom",
  "Isle of Man", "United Kingdom",
  "United Kingdom Of Great Britain And Northern Ireland", "United Kingdom",
  "United Kingdom Of Great Britain", "United Kingdom",
  "The United Kingdom", "United Kingdom",
  "U.k.", "United Kingdom",
  "United Kingdum", "United Kingdom",
  "Unitedkingdom Of Great Britain", "United Kingdom",
  "United Kingdom Of Great Britain And Northern Irel", "United Kingdom",
  "United Of Kingdom", "United Kingdom",
  "Unitred Kingdom", "United Kingdom",
  "Scotland", "United Kingdom",
  "Wales", "United Kingdom",
  "Northern Ireland", "United Kingdom",
  "England", "United Kingdom",
  "Gibraltar", "United Kingdom",
  "Uk", "United Kingdom",
  "Cayman Islands", "United Kingdom",
  
  "Japan Only", "Japan",
  "None (Japan Only)", "Japan",
  
  "Palestinian Territory, Occupied", "Palestine",
  "Palestinian Territories", "Palestine",
  "Palestinian Territories, Occupied", "Palestine",
  
  "Thai", "Thailand",
  "Maetha", "Thailand",
  "Thailand.", "Thailand",
  "Thailan", "Thailand",
  "Thaiand", "Thailand",
  
  "Aland Islands", "Finland",
  "Republic Of Finland", "Finland",
  
  "Mone", "Myanmar",
  
  "México", "Mexico",
  "Republic Of Mexico", "Mexico",
  
  "Taiwan R.o.c", "Taiwan",
  "Taiwan Province Of China", "Taiwan",
  "Taiwan,Province Of China", "Taiwan",
  "Taiwian", "Taiwan",
  "Taiwan, Province Of China", "Taiwan",
  "Taiwan (Province Of China)", "Taiwan",
  
  "Venezuela, Bolivarian Republic Of", "Venezuela",
  "Venezuela (Bolivarian Republic Of)", "Venezuela",
  
  "Malaysa", "Malaysia",
  "Malasya", "Malaysia",
  
  "Serbia And Montenegro", "Serbia",
  "Former Serbia and Montenegro", "Serbia",
  "Republic Of Serbia", "Serbia",
  
  "Isreal", "Israel",
  "State of Irael", "Israel",
  "Israelc", "Israel",
  "State Of Israel", "Israel",
  
  "Italia", "Italy",
  "Republic Of Italy", "Italy",
  "Itary", "Italy",
  
  "Spein", "Spain",
  "Spain,", "Spain",
  "Kingdom Of Spain", "Spain",
  
  "Sweden Etc", "Sweden",
  "Sweeden", "Sweden",
  "Swizerland", "Switzerland",
  
  "The Gambia", "Gambia",
  
  "Canda", "Canada",
  "Columbia", "Colombia",
  
  "Paraguar", "Paraguay",
  
  "Hungaria", "Hungary",
  
  "Republic Of Portugal", "Portugal",
  "Portugal,", "Portugal",
  
  "Bulugaria", "Bulgaria",
  
  "Commonwealth Of Australia", "Australia",
  "Australia", "Australia",
  "Christmas Island", "Australia",
  "Australia," , "Australia",
  
  "Reunion", "France",
  "Rèunion", "France",
  "Réunion", "France",
  
  
  "Federal Republic Of Germany", "Germany",
  "German", "Germany",
  
  "Kingdom Of Belgium", "Belgium",
  "Beugium", "Belgium",
  "Berugium", "Belgium",
  
  "Republic Of Singapore", "Singapore",
  "Shingapore", "Singapore",
  
  "Latovia", "Latvia",
  
  "Ukrine", "Ukraine",
  
  "Ausria", "Austria",
  "Austraria", "Austria",
  
  "Republic Of South Africa", "South Africa",
  
  "Bolivia, Plurinational State Of", "Bolivia",
  
  "Libyan Arab Jamahiriya", "Libya",
  
  "Arabia Saudi", "Saudi Arabia",
  
  "Bosnial And Herzegovina", "Bosnia And Herzegovina",
  
  "Lebano", "Lebanon",
  
  "Roumania", "Romania",
  
  "Newzealand", "New Zealand",
  "Niue", "New Zealand", 
  
  "Slovakia,", "Slovakia",
  
  "Brunei Darussalam", "Brunei",
  
  "Faroe Islands", "Denmark",
  "Greenland", "Denmark",
  "French Guiana", "France", 
  "Martinique", "France", 
  "Guadeloupe", "France",
  "Swaziland", "Eswatini")

write_csv(country_lookup, "country_lookup.csv")


### 8b. Fixing 7 problematic rows where the countries are separated by , rather than ;. 
problem_rows <- ICTRP_All %>%
  filter(
    str_detect(countries, "Argentina,Belarus") |
      str_detect(countries, "Australia,Austria") |
      str_detect(countries, "Australia,China") |
      str_detect(countries, "Australia,India") | 
      str_detect(countries, "United States, China") |
      str_detect(countries, "South Korea, Taiwan") |
      str_detect(countries, "India,Italy,Malaysia,Sri Lanka"))


problem_rows_fixed <- problem_rows %>%
  mutate(
    countries = str_replace_all(countries, ",", ";"))

ICTRP_All_country <- ICTRP_All %>% filter(!countries %in% problem_rows$countries) %>% bind_rows(problem_rows_fixed)


### 8c. Run code to make function using country_lookup file while ensuring spaces, capital letters, symbols and certain punctuation do not prevent matches between countries 

norm_country <- function(x) {
  x %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_to_lower() %>%
    str_replace_all("[’`]", "'") %>%
    str_replace_all('"', "") %>%
    str_squish()}

country_lookup <- read_csv("country_lookup.csv", show_col_types = FALSE) %>%
  mutate(raw_key = norm_country(raw_country))

###8d. Breakdown by country for ICTRP_All. 
#Note that ICTRP_ALl_countries_fractional makes a new row for each additional hosting country under the same study (which explains why there are 327,339 rows)

ICTRP_All_countries_fractional <- ICTRP_All_country %>%
  mutate(
    row_id = row_number(),
    countries_clean = countries %>%
      str_replace_all('"', "") %>%
      str_replace_all("<br\\s*/?>", ";") %>%
      str_replace_all("\\s*;\\s*", ";") %>%
      str_squish(),
    countries_clean = if_else(is.na(countries_clean) | countries_clean == "", "Missing data", countries_clean))  %>%
  separate_rows(countries_clean, sep = ";") %>%
  mutate(
    countries_clean = str_trim(countries_clean),
    countries_clean = if_else(is.na(countries_clean) | countries_clean == "", "Missing data", countries_clean),
    raw_key = norm_country(countries_clean)
  ) %>%
  left_join(
    country_lookup %>% select(raw_key, country_standard),
    by = "raw_key"
  ) %>%
  mutate(
    countries_final = case_when(
      !is.na(country_standard) ~ country_standard,
      raw_key == "missing data" ~ "Missing data",
      TRUE ~ str_to_title(countries_clean)
    )
  ) %>%
  select(-countries_clean, -raw_key, -country_standard) %>%
  rename(countries_clean = countries_final) %>%
  distinct(trial_id, countries_clean, .keep_all = TRUE) %>%
  group_by(trial_id) %>%
  mutate(
    n_countries = n_distinct(countries_clean),
    country_weight = 1 / n_countries
  ) %>%
  ungroup()


###8e. Create the final counts table including missing data

All_country_counts_fractional <- ICTRP_All_countries_fractional %>%
  group_by(countries_clean) %>%
  summarise(
    weighted_studies = sum(country_weight),
    n_study_country_records = n(),
    .groups = "drop"
  ) %>%
  mutate(
    percent_of_studies = 100 * weighted_studies / n_distinct(ICTRP_All_countries_fractional$trial_id)
  ) %>%
  arrange(desc(weighted_studies))

### 8f. CHECK: To make sure i've not lost studies in the process 
sum(All_country_counts_fractional$weighted_studies)
n_distinct(ICTRP_All$trial_id) #Correct, 183041 studies remain in both.

#### 9. Breakdown by COUNTRY for ICTRP_Met
### 9a. Repeat the same code from section 8 for ICTRP_Met

ICTRP_Met_countries_fractional <- ICTRP_Met %>%
  mutate(
    row_id = row_number(),
    countries_clean = countries %>%
      str_replace_all('"', "") %>%
      str_replace_all("<br\\s*/?>", ";") %>%
      str_replace_all("\\s*;\\s*", ";") %>%
      str_squish(),
    countries_clean = if_else(is.na(countries_clean) | countries_clean == "", "Missing data", countries_clean))  %>%
  separate_rows(countries_clean, sep = ";") %>%
  mutate(
    countries_clean = str_trim(countries_clean),
    countries_clean = if_else(is.na(countries_clean) | countries_clean == "", "Missing data", countries_clean),
    raw_key = norm_country(countries_clean)
  ) %>%
  left_join(
    country_lookup %>% select(raw_key, country_standard),
    by = "raw_key"
  ) %>%
  mutate(
    countries_final = case_when(
      !is.na(country_standard) ~ country_standard,
      raw_key == "missing data" ~ "Missing data",
      TRUE ~ str_to_title(countries_clean)
    )
  ) %>%
  select(-countries_clean, -raw_key, -country_standard) %>%
  rename(countries_clean = countries_final) %>%
  distinct(trial_id, countries_clean, .keep_all = TRUE) %>%
  group_by(trial_id) %>%
  mutate(
    n_countries = n_distinct(countries_clean),
    country_weight = 1 / n_countries
  ) %>%
  ungroup()


### 9b. Create the final counts table including missing data

Met_country_counts_fractional <- ICTRP_Met_countries_fractional %>%
  group_by(countries_clean) %>%
  summarise(
    weighted_studies = sum(country_weight),
    n_study_country_records = n(),
    .groups = "drop"
  ) %>%
  mutate(
    percent_of_studies = 100 * weighted_studies / n_distinct(ICTRP_Met_countries_fractional$trial_id)
  ) %>%
  arrange(desc(weighted_studies))

### 9c. CHECK: To make sure i've not lost studies in the process 
sum(Met_country_counts_fractional$weighted_studies)
n_distinct(ICTRP_Met$trial_id) #Correct, 640 studies remain in both.


#### 10. Analysis of country breakdown
### 10a. The percent of metabolomic studies as a fraction of all studies per hosting country 
country_met_fraction <- All_country_counts_fractional %>%
  select(countries_clean, n_all = weighted_studies) %>%
  left_join(
    Met_country_counts_fractional %>%
      select(countries_clean, n_met = weighted_studies),
    by = "countries_clean"
  ) %>%
  mutate(
    n_met = if_else(is.na(n_met), 0, n_met),  # countries with no metabolomics studies
    fraction_metabolomics = n_met / n_all,
    percent_metabolomics = 100 * fraction_metabolomics
  ) %>%
  arrange(desc(percent_metabolomics))

write_csv(country_met_fraction, "country_met_fraction.csv")

### 10b. Create a table (for supplementary data) of the percent metabolomic studies for each country (which also includes both weighted and total counts) 
country_met_fraction2 <- country_met_fraction %>%
  left_join(
    All_country_counts_fractional %>%
      select(countries_clean, n_study_country_records) %>%
      rename(All_total_counts = n_study_country_records),
    by = "countries_clean"
  ) %>%
  left_join(
    Met_country_counts_fractional %>%
      select(countries_clean, n_study_country_records) %>%
      rename(Met_total_counts = n_study_country_records),
    by = "countries_clean")

### 10c. Select the 8 countries with the highest percent of metabolomic studies (with at least 1000 hosted trials)
top8_fraction <- country_met_fraction2 %>%
  filter(All_total_counts > 1000) %>%          # only countries with >1000 trials
  arrange(desc(fraction_metabolomics)) %>%    # highest fraction first
  slice_head(n = 8)     

### 10d. Set standardised colours for each country for representation on the graphs and ensure the registry colours align with the country colours where possible. 

country_colours <- c(
  "United States" = "#4C78A8",   # blue
  "China" = "#E45756",           # red
  "India" = "#F58518",           # orange
  "Germany" = "#54A24B",         # green
  "United Kingdom" = "#B279A2",  # purple
  "Spain" = "#FF9DA6",           # pink
  "France" = "#9D755D",          # brown
  "Italy" = "#72B7B2",           # teal
  "Canada" = "#B5BD61",          # olive
  "Holland" = "#86BCB6",         # soft cyan-green
  "Australia" = "#F2CF5B",       # muted yellow
  "Denmark" = "#E17C05",         # burnt orange
  "Missing data" = "#BAB0AC")   # grey

registry_colours <- c(
  "Clinical Trials.gov (USA)" = "#4C78A8",   # United States
  "Australian New Zealand Clinical Trials Registry" = "#F2CF5B",   # Australia
  "Clinical Trials Registry India" = "#F58518",   # India
  "Pan African Clinical Trial Register" = "#7f7f7f",   # Africa → neutral grey
  "EU Clinical Trials Register/Clinical Trials Information System" = "#54A24B",   # Europe → Germany green
  "ISRCTN (UK)" = "#B279A2",   # United Kingdom
  "Thai Clinical Trials Register" = "#FF9DA6",   # Thailand → pink (closest available)
  "Japan Primary Registries Network" = "#72B7B2",   # Japan → teal
  "German Clinical Trials Register" = "#54A24B",   # Germany
  "Brazilian Clinical Trials Registry" = "#17becf",   # Brazil → teal-blue
  "Cuban Public Registry of Clinical Trials" = "#8c564b",   # Cuba → brown
  "Peruvian Clinical Trials Register" = "#bcbd22",   # Peru → olive
  "Lebanese Clinical Trials Registry" = "#9D755D",   # Lebanon → brown variant
  "Chinese Clinical Trial Register" = "#E45756",   # China
  "Iranian Registry of Clinical Trials" = "#e377c2",   # Iran → pink-purple
  "Clinical research information centre (Republic of Korea)" = "#1f77b4",   # South Korea → blue variant
  "International Traditional Medicine Clinical Trial Registry" = "#B22222",   # Close to chinese red 
  "Sri Lanka Clinical Trials Register" = "#ff7f0e",   # Sri Lanka → orange variant
  "The Netherlands National Trial Register" = "#86BCB6",  # Holland
  "Other" = "#BAB0AC")

#### 11.) GRAPHS for stratified analysis by country

### 11a. Select top countries by metabolomics count (excluding Missing data) with horizontal ranked bar chart with colour for percentage metabolomic 
plot_data <- country_met_fraction %>%
  filter(countries_clean != "Missing data") %>%
  arrange(desc(n_met)) %>%
  slice_head(n = 20)

ggplot(plot_data, aes(x = reorder(countries_clean, n_met), y = n_met)) +
  geom_col(aes(fill = percent_metabolomics)) +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    x = "Country",
    y = "Number of metabolomics studies",
    fill = "% metabolomics studies",
    title = "Metabolomics studies by country") +
  theme_minimal()


### 11c. Percent Metabolomic studies by country: Same again but instead a scatter plot (using weighted scores)

plot_data <- country_met_fraction %>% 
  filter(countries_clean != "Missing data", n_all >= 4000)

ggplot(plot_data, aes(x = n_all, y = n_met)) +
  geom_point(
    aes(colour = percent_metabolomics),  # 👈 keep colour mapping
    size = 4,                            # 👈 fixed size for all points
    alpha = 0.7
  ) +
  geom_text_repel(aes(label = countries_clean), size = 3) + scale_colour_gradient(low = "#deebf7", high = "#08519c") +
  labs(
    x = "Total studies",
    y = "Metabolomics studies",
    colour = "% metabolomics"
  ) +
  theme_minimal() 

### 11d. Number of metabolomic studies by COUNTRY: simple bar chart inlcuding only top 8 countries by absolute number of metabolomic studies 

plot_data <- country_met_fraction2 %>%
  filter(countries_clean != "Missing data") %>%
  arrange(desc(Met_total_counts)) %>%
  slice_head(n = 8)

ggplot(plot_data, aes(x = reorder(countries_clean, Met_total_counts), y = Met_total_counts, fill = countries_clean)) +
  geom_col() +
  geom_text(aes(label = Met_total_counts), hjust = -0.2) +
  coord_flip() + scale_fill_manual(values = country_colours) +
  labs(
    x = "Country",
    y = "Number of metabolomic studies",
    fill = "Country",
    title = "Top 8 countries by number of metabolomic studies"
  ) +
  theme_minimal() +
  guides(fill = "none") +
  expand_limits(y = max(plot_data$Met_total_counts) * 1.1)

### 11e. A bubble plot for countries with >1000 trials (absolute count), where percent metabolomic is on weighted numbers (top 12 countries have same colours as subsequent national registries)

plot_data <- country_met_fraction2 %>%
  filter(All_total_counts > 1000, countries_clean != "Missing data") %>%
  arrange(percent_metabolomics) %>%
  mutate(countries_clean = factor(countries_clean, levels = countries_clean),
         
         # Assign colours only to countries in my palette, others become "Other"
         colour_group = ifelse(countries_clean %in% names(country_colours), as.character(countries_clean), "Other" ))

ggplot(plot_data, aes(x = percent_metabolomics, y = countries_clean)) +
  geom_point(aes(size = All_total_counts, colour = colour_group), alpha = 0.7) +
  
  scale_colour_manual(
    values = c(country_colours, "Other" = "grey70"),
    guide = "none") + scale_size_continuous(name = "All studies") + labs (
      x = "Percentage",
      y = "Hosting country",
      title = "Adoption of metabolomics by country") +
  
  theme_minimal() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold"))


#### 12.) GRAPHS for stratified analysis by PHASE

### 12a. ICTRP_All Phases 

plot_data_all <- ICTRP_All %>%
  filter(!is.na(phase_group)) %>%
  count(phase_group, name = "n") %>%
  arrange(desc(n))

ggplot(plot_data_all, aes(x = reorder(phase_group, n), y = n, fill = phase_group)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2") +
  geom_text(aes(label = n), hjust = -0.2, fontface = "bold") +
  coord_flip() +
  expand_limits(y = max(plot_data$n) * 1.1) +   # 👈 FIX
  labs(
    x = "Phase of clinical study",
    y = "Number of studies",
    title = "Distribution of all clinical studies by phase",
    fill = "Phase"
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"))

### 12b. ICTRP_Met phases
plot_data_met <- ICTRP_Met %>%
  filter(!is.na(phase_group)) %>%
  count(phase_group, name = "n") %>%
  mutate(
    phase_group = factor(
      phase_group,
      levels = c("Phase 0", "Phase 1/2", "Phase 1", "Phase 2")))


ggplot(plot_data_met, aes(x = phase_group, y = n, fill = phase_group)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2") +
  geom_text(aes(label = n), hjust = -0.2, fontface = "bold") +
  coord_flip() +
  expand_limits(y = max(plot_data$n) * 1.1) +   # 👈 FIX
  labs(
    x = "Phase of clinical study",
    y = "Number of studies",
    title = "Distribution of Metabolomic clinical studies by phase",
    fill = "Phase"
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"))

### 12c. Proportion of metabolomic studies by phase 
phase_levels <- c("Phase 0", "Phase 1/2", "Phase 1", "Phase 2")

plot_data_phase_percent <- plot_data_all %>%
  rename(n_all = n) %>%
  left_join(
    plot_data_met %>% rename(n_met = n),
    by = "phase_group"
  ) %>%
  mutate(
    n_met = if_else(is.na(n_met), 0, n_met),
    percent_metabolomic = 100 * n_met / n_all,
    phase_group = factor(phase_group, levels = phase_levels))

ggplot(plot_data_phase_percent, aes(x = phase_group, y = percent_metabolomic, fill = phase_group)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2") +
  geom_text(aes(label = sprintf("%.1f%%", percent_metabolomic)),
            hjust = -0.2, fontface = "bold") +
  coord_flip() +
  expand_limits(y = max(plot_data_phase_percent$percent_metabolomic) * 1.1) +
  labs(
    x = "Phase of clinical study",
    y = "Metabolomic studies (%)",
    title = "Increasing corporation of metabolomics at earlier phases",
    fill = "Phase"
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"))

### 12d Dot plot for percent metabolomic by phase 
ggplot(plot_data_phase_percent, aes(x = percent_metabolomic, y = reorder(phase_group, percent_metabolomic))) +
  geom_point(size = 4) +
  labs(
    x = "Metabolomic studies (%)",
    y = "Phase"
  ) +
  theme_minimal()

### 12e. Number (and percentage of total number) of Metabolomic studies by phase 

plot_data_met2 <- plot_data_met %>% mutate (percentage= (n / 640) *100)

### 12f. Pie chart of metabolomic studies by phase
ggplot(plot_data_met2, aes(x = "", y = percentage, fill = phase_group)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Trial Phase")

#### 13. GRAPHS for stratified analysis by REGISTRY
### 13a. Breakdown by registry
all_registry_counts <- ICTRP_All %>% count(registry_name, name = "n_studies", sort = TRUE)
met_registry_counts <- ICTRP_Met %>% count(registry_name, name= "n_studies", sort=TRUE)

registry_met_fraction <- all_registry_counts %>%
  rename(n_all = n_studies) %>%
  left_join(
    met_registry_counts %>% rename(n_met = n_studies),
    by = "registry_name"
  ) %>%
  mutate(
    n_met = if_else(is.na(n_met), 0, n_met),  # registries with no metabolomics studies
    percent_metabolomics = 100 * n_met / n_all
  ) %>%
  arrange(desc(percent_metabolomics))

### 13b. Breakdown by registry with EUCTR and CTIS combined 
#Rewriting registry_name to combine EUCTR and CTIS 
ICTRP_All_Registry <- ICTRP_All %>%
  mutate(
    registry_name = recode(
      registry_name,
      "EU Clinical Trials Register" = "EU Clinical Trials Register/Clinical Trials Information System",
      "Clinical Trials Information System" = "EU Clinical Trials Register/Clinical Trials Information System"))

ICTRP_Met_Registry <- ICTRP_Met %>%
  mutate(
    registry_name = recode(
      registry_name,
      "EU Clinical Trials Register" = "EU Clinical Trials Register/Clinical Trials Information System",
      "Clinical Trials Information System" = "EU Clinical Trials Register/Clinical Trials Information System"))


all_registry_counts_combined <- ICTRP_All_Registry %>% count(registry_name, name = "n_studies", sort = TRUE)
met_registry_counts_combined  <- ICTRP_Met_Registry %>% count(registry_name, name= "n_studies", sort=TRUE)

### 13c. To produce dataframe for supplementary information on the fraction of metabolomic studies per registry
combined_registry <- all_registry_counts_combined %>%
  rename(n_all = n_studies) %>% 
  left_join(
    met_registry_counts_combined %>%
      rename(n_met = n_studies),
    by = "registry_name")

combined_registry <- combined_registry %>%
  mutate(
    n_met = coalesce(n_met, 0),   # replace NA with 0
    fraction_metabolomics = n_met / n_all)

combined_registry <- combined_registry %>%
  mutate(
    percent_metabolomics = 100 * fraction_metabolomics)

### 13d. To make a piechart of the breakdown by REGISTRY of metabolomic studies (excluding registries with 5 or less metabolomic studies)
pie_data <- combined_registry %>%
  mutate(
    registry_grouped = if_else(n_met <= 5, "Other registries", registry_name)
  ) %>%
  group_by(registry_grouped) %>%
  summarise(n_met = sum(n_met), .groups = "drop") %>%
  arrange(desc(n_met))

### 13d (continued) Add 'Other registries" colour to my predefined palette'
registry_colours <- c(registry_colours,"Other registries" = "#BAB0AC")  # grey

### 13d (continued) Plot pie chart   
ggplot(pie_data, aes(x = "", y = n_met, fill = registry_grouped)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = registry_colours) +
  labs(fill = "Registry") +
  theme_void()

### 13e. To make a ranked bar chart of metabolomic studies by registry
bar_data <- combined_registry %>%
  mutate(
    registry_grouped = if_else(n_met <= 5, "Other registries", registry_name)) %>%
  group_by(registry_grouped) %>%
  summarise(n_met = sum(n_met), .groups = "drop") %>%
  arrange(desc(n_met)) %>%
  mutate(
    registry_grouped = factor(registry_grouped, levels = rev(registry_grouped)))  

registry_colours <- c(
  registry_colours,
  "Other registries" = "#BAB0AC")

ggplot(bar_data, aes(x = registry_grouped, y = n_met, fill = registry_grouped)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = registry_colours) +
  geom_text(aes(label = n_met), hjust = -0.2, size = 4) +
  labs(
    x = "Registry",
    y = "Number of metabolomic studies",
    title = "Distribution of metabolomic studies by registry"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")) +
  expand_limits(y = max(bar_data$n_met) * 1.1)

### 13f. Graph: Bar chart of top 8 registries by number of metabolomic studies 

plot_data <- combined_registry %>%arrange(desc(n_met)) %>% slice_head(n = 8)


ggplot(plot_data, aes(x = reorder(registry_name, n_met), y = n_met, fill = registry_name)) +
  geom_col() +
  geom_text(aes(label = n_met), hjust = -0.2) +   
  coord_flip() + scale_fill_manual(values = registry_colours) +
  labs(
    x = "Registry",
    y = "Number of metabolomic studies",
    title = "Top 8 registries by number of metabolomic studies",
    fill = "Registry"
  ) +
  theme_minimal() +
  guides(fill = "none") +   
  expand_limits(y = max(plot_data$n_met) * 1.1)

#### 14.) Stratified analysis of by year and registry

### 14a. Start with ICTRP_Met:produce a stacked bar chart of metabolomic studies per year with breakdown of the contributions of the top 6 registries (all "Other" registries are lumped together for visual simplicity)

## 14a (continued) Identify top 6 registries (excluding Multiple registries if desired)
top_registries <- ICTRP_Met_Registry %>%
  filter(!is.na(registry_name)) %>%
  count(registry_name, sort = TRUE) %>%
  slice_head(n = 6) %>%
  pull(registry_name)

## 14a (continued) Create new dataframe with "Other" category (so that all registries not in the top 6 are lumped into other category)
ICTRP_Met_Registry2 <- ICTRP_Met_Registry %>%
  mutate(
    registry_name = case_when(
      is.na(registry_name) ~ "Other",
      registry_name %in% top_registries ~ registry_name,
      TRUE ~ "Other"))

## 14a (continued) Summarise for plotting
Met_plot_data <- ICTRP_Met_Registry2 %>%
  filter(!is.na(year)) %>%
  count(year, registry_name, name = "n")

## 14a (continued) Plot stacked bar chart
ggplot(Met_plot_data, aes(x = year, y = n, fill = registry_name)) +
  geom_col() + scale_fill_manual(values = registry_colours) +
  labs(
    x = "Year",
    y = "Metabolomic studies",
    fill = "Registry",
    title = "Metabolomic studies over time by registry"
  ) + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"))


#### 14b. Repeating the process now for ICTRP_All 

top_registries <- ICTRP_All_Registry %>%
  filter(!is.na(registry_name)) %>%
  count(registry_name, sort = TRUE) %>%
  slice_head(n = 6) %>%
  pull(registry_name)

## 14b(continued) Create new dataframe with "Other" category
ICTRP_All_Registry2 <- ICTRP_All_Registry %>%
  mutate(
    registry_name = case_when(
      is.na(registry_name) ~ "Other",   
      registry_name %in% top_registries ~ registry_name,
      TRUE ~ "Other"))

## 14b (continued) Summarise for plotting
All_plot_data <- ICTRP_All_Registry2 %>%
  filter(!is.na(year)) %>%
  count(year, registry_name, name = "n")

## 14b (continued) Plot stacked bar chart
ggplot(All_plot_data, aes(x = year, y = n, fill = registry_name)) +
  geom_col() + scale_fill_manual(values = registry_colours) +
  labs(
    x = "Year",
    y = "All early phase clinical studies ",
    fill = "Registry",
    title = "All early phase studies over time"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"))


### 15. Producing a Combination stacked bar chart for total studies per year and registry with overlying metabolomic percentage on the second axis 

#15a ICTRP_All_registry breakdown by year for subsequent trend of metabolomic percent per year 
All_studies_per_year <- ICTRP_All_Registry %>% count(year, name = "n_all_studies") %>% arrange(year)

Met_studies_per_year <- ICTRP_Met_Registry %>%
  count(year, name = "n_met_studies") %>%
  arrange(year)


studies_per_year_combined <- All_studies_per_year %>%
  rename(n_all = n_all_studies) %>%   # adjust if your column name differs
  left_join(
    Met_studies_per_year %>% rename(n_met = n_met_studies),
    by = "year"
  ) %>%
  mutate(
    n_met = if_else(is.na(n_met), 0, n_met),
    percent_metabolomic = 100 * n_met / n_all
  ) %>%
  arrange(year)


### 15b. Top registries + Other for all early phase studies
top_registries <- ICTRP_All_Registry %>%
  filter(!is.na(registry_name)) %>%
  count(registry_name, sort = TRUE) %>%
  slice_head(n = 6) %>%
  pull(registry_name)

ICTRP_All_Registry_plot <- ICTRP_All_Registry %>%
  mutate(
    registry_name = case_when(
      is.na(registry_name) ~ "Other",
      registry_name %in% top_registries ~ registry_name,
      TRUE ~ "Other"),
    registry_name = factor(registry_name, levels = c(top_registries, "Other")))

### 15c. Stacked bar data
plot_data <- ICTRP_All_Registry_plot %>%
  filter(!is.na(year)) %>%
  count(year, registry_name, name = "n")

### 15d. Secondary-axis line data
line_data <- studies_per_year_combined %>%
  filter(!is.na(year), !is.na(percent_metabolomic)) %>%
  arrange(year)

### 15e. Scale factor to map % onto left y-axis
scale_factor <- max(plot_data$n %>% tapply(plot_data$year, sum), na.rm = TRUE) /
  max(line_data$percent_metabolomic, na.rm = TRUE)

ggplot(plot_data, aes(x = year, y = n, fill = registry_name)) +
  geom_col() +
  geom_line(
    data = line_data,
    aes(x = year, y = percent_metabolomic * scale_factor, group = 1),
    inherit.aes = FALSE,
    linewidth = 1.1
  ) + scale_fill_manual(values = registry_colours) +
  geom_point(
    data = line_data,
    aes(x = year, y = percent_metabolomic * scale_factor),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_y_continuous(
    name = "All early phase trials",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Metabolomic (%)"
    )
  ) +
  labs(
    x = "Year",
    fill = "Registry",
    title = "Growing inclusion of metabolomics in early phase clinical trials"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"))


#### 16.) Stratified analysis by OUTCOME

ICTRP_metabolomics %>%
  summarise(
    n_primary_outcome = sum(str_detect(found_in, "primary_outcome"), na.rm = TRUE),
    n_secondary_outcomes = sum(str_detect(found_in, "secondary_outcomes"), na.rm = TRUE))


### 16a. Horizontal bar plot for breakdown by outcome 

outcome_data <- tibble(
  category = c("Primary", "Secondary"),
  n = c(198, 640 - 198)) %>%
  mutate(
    percent = n / sum(n),
    label = paste0(round(percent * 100, 1), "%"))

ggplot(outcome_data, aes(x = 1, y = percent, fill = category)) + geom_col(width = 0.3) + coord_flip() +
  geom_text( aes(label = label), position = position_stack(vjust = 0.5), colour = "white",fontface = "bold",
             size = 5) +
  
  scale_fill_manual(values = c(
    "Primary" = "#4C78A8",      
    "Secondary" = "#D3D3D3")) + labs(
      x = NULL,
      y = NULL,
      fill = NULL,
      title = "Use of metabolomics in early phase trials by outcome classification"
    ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")) 


#### 17.) Stratified analysis by CONTINENT (this file has been made by allocating each country into a continent manually in excel)

continent_analysis<-read_csv("path/to/continent_analysis.csv")
View(continent_analysis)


continent_totals <- continent_analysis %>%
  filter(Continent %in% c(
    "Asia",
    "Africa",
    "Europe",
    "North America",
    "South America",
    "Oceania"
  )) %>%
  group_by(Continent) %>%
  summarise(
    total_Metabolomic_Absolute_counts = sum(Metabolomic_Absolute_counts, na.rm = TRUE) ) %>%
  mutate( percent_of_total = (total_Metabolomic_Absolute_counts / sum(total_Metabolomic_Absolute_counts)) * 100)

continent_totals
