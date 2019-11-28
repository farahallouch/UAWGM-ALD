# UAWGM-ALD
Using the UAW-GM occupational cohort to explore the association between age at leaving work and alcohol related liver disease (ALD) mortality between 1941 and 2015.

# loading libraries
library(tidyverse)
library(haven)
library(splitstackshape)
library(knitr)
library(survival)
library(survminer)

# reading in original dataset
auto_vs_15 <- read_sas("auto_vs_15.sas7bdat")

# Subsetting dataset
gm_dta <- auto_vs_15 %>% 
  filter(!(NOHIST == 1), 
         !(WH == 0)) %>% # excluding people with < 50% work history (standard)
  select(cod_15, V_ICD, FINRACE, SEX, # selecting relevant variables (ICD codes, ICD 9 or 10, race, sex)
         oddend, STUDYNO, YOB, PLANT, # (subjects with a weird end date for follow-up, study number, year of birth, plant #)
         yod15, YOUT16, yrin16, yrout15_new) %>% # If subject was still at GM for end of exposure assessment, YOUT = 95 (year of death, year of leaving work, year of entering work, year of leaving dataset (either yod or 2016 end of follow-up)
  mutate(yout16.0 = trunc(YOUT16) + 1900, # truncating years of leaving work for P-Y calculation
         yrout15_new.0 = ceiling(yrout15_new) + 1900, # taking ceiling of year of leaving cohort for P-Y calculation
         race = ifelse(FINRACE != 2, 1, 0), # changing unknowns to white = 1 (standard), black = 0
         sex = ifelse(SEX == 2, 0, 1), # changing female to 0, male = 1
         person_years = yrout15_new.0 - yout16.0, # creating P-Y column that starts at the beginning of the year people leave work and ends at the end of the year they die or end of follow up (2016)
         age_out = YOUT16 - YOB, # calculate age at leaving work
         ALD = ifelse(cod_15 == "5710" | # Alcoholic fatty liver ICD9 (2) # ADD K70, K73, K74 with ICD9 EQUIVALENTS (CASE AND DEATON 2015)
                      cod_15 == "5711" | # Acute alcoholic hepatitis ICD9 (0)
                      cod_15 == "5712" | # Alcoholic cirrhosis of liver ICD9 (12)
                        cod_15 == "5713" | # Alcoholic liver damage, unspecified ICD9 (2)
                        cod_15 == "5714" | # Chronic hepatitis ICD9 (1)
                        cod_15 == "5715" | # Cirrhosis of liver without mention of alcohol ICD9 (16)
                        cod_15 == "5716" | # Biliary cirrhosis
                        cod_15 == "5719" | # Unspecified chronic liver disease without mention of alcohol
                        cod_15 == "571x" | # alcohol-related, but missing last nb
                        cod_15 == "K700" | # Alcoholic fatty liver ICD10
                        cod_15 == "K701" | # Alcoholic hepatitis ICD10
                        cod_15 == "K702" | # Alcoholic fibrosis and sclerosis of liver ICD10
                        cod_15 == "K703" | # Alcoholic cirrhosis of liver ICD10
                        cod_15 == "K704" | # Alcoholic hepatic failure ICD10
                        cod_15 == "K709" | # Alcoholic liver disease, unspecified ICD 10
                        cod_15 == "K730" | # Chronic persistent hepatitis, not elsewhere classified
                        cod_15 == "K731" | # Chronic lobular hepatitis, not elsewhere classified
                        cod_15 == "K732" | # Chronic active hepatitis, not elsewhere classified
                        cod_15 == "K738" | # Other chronic hepatitis, not elsewhere classified
                        cod_15 == "K739" | # Chronic hepatitis, unspecified
                        cod_15 == "K740" | # Hepatic fibrosis ICD10
                        cod_15 == "K741" | # Hepatic sclerosis ICD10
                        cod_15 == "K742" | # Hepatic fibrosis with hepatic sclerosis ICD10
                        cod_15 == "K743" | # Primary biliary cirrhosis ICD10
                        cod_15 == "K744" | # Secondary biliary cirrhosis ICD10
                        cod_15 == "K745" | # Biliary cirrhosis, unspecified ICD10
                        cod_15 == "K746", 1, 0)) %>% # Other and unspecified cirrhosis of liver
  select(-FINRACE, # deselect original race and sex coding to make dataset cleaner
         -SEX) %>% 
  mutate(yod15m = ifelse(is.na(yod15), 9999, yod15)) %>% # Removing people whose date of death comes before their date 
  filter(!(yod15m < YOUT16)) %>% 
  select(-(yod15m)) %>% 
  filter(!(person_years <= 0)) # excluding negative person years

# Determining age categories (5) with = number of cases in each 
gm_dta_cox %>% 
  filter(ALD == 1) %>% summarize(twenty = round(quantile(age_out, 0.2)),
                                          fourty = round(quantile(age_out, 0.4)),
                                          sixty = round(quantile(age_out, 0.6)),
                                          eighty = round(quantile(age_out, 0.8)))

gm_dta_cox %>% summarize(min = trunc(min(age_out)),
                     max = ceiling(max(age_out)))

# Adding age out categories to gm_dta
gm_dta <- gm_dta %>% 
  mutate(age_cat_out = as.factor(ifelse(age_out <= 33, "(18, 33]", # min age_cat to 20th percentile 
                                       ifelse(age_out > 33 & age_out <= 42, "(34, 42]", # 20th to 40th percentile
                                       ifelse(age_out > 42 & age_out <= 50, "(42, 50]", # 40th to 60th percentile 
                                       ifelse(age_out > 50 & age_out <= 59, "(50, 59]", "(59, 85]")))))) # 60th to 80th percentile and 80th to max age_out
  
# Relevel age categories
gm_dta$age_cat_out <- relevel(gm_dta$age_cat_out, "(59, 85]") # 59+ as reference

# Cox model data with non-time varying covariates
gm_dta_cox <- gm_dta %>% 
  filter(sex == 1)

# Cox model data with calendar year time-varing covariate
gm_dta_cox_tv <- gm_dta %>% 
  expandRows("person_years", count.is.col = TRUE, drop = FALSE) %>% 
  group_by(STUDYNO) %>%
  mutate(year_obs = row_number(),
         cal_obs = floor(1900 + YOUT16 + (year_obs - 1)),
         ALD = ifelse(ALD == 1 & (1900 + floor(yod15)) == cal_obs, 1, 0),
         cal_obs_cat = cut(cal_obs, breaks = c(1900, seq(1950, 2010, by = 10)))) %>% 
  filter(sex == 1)

# Relevel calendar year categories
gm_dta_cox_tv$cal_obs_cat <- relevel(gm_dta_cox_tv$cal_obs_cat, "(1.97e+03,1.98e+03]")
```
