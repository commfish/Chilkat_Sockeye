# missfill code

# load ----
library(dplyr)
library(tidyr)
source('coho_apportionment/code/missfill_functions.r')

# data ----
read.csv('coho_apportionment/data/chilkat_coho.csv')-> chilkat

# pivot year columns into long format for import into functions
df_long <- chilkat %>%
  pivot_longer(
    cols = -date,    # all year columns now start with 'x'
    names_to = "year",
    values_to = "total_count") %>%
  mutate(year = as.numeric(sub("X", "", year)))

# run functions ----
impute_global(df_long) #imputes all NAs in all years iteratively
impute_local(df_long) #imputes a 10-year rolling imputation (prev & following 5 years)

# read and save output in wide format----
read.csv("coho_apportionment/output/global_imputed.csv") %>%
  distinct(date, year, .keep_all = TRUE) %>%  # remove duplicates
  pivot_wider(
    id_cols = date,
    names_from = year,
    values_from = total_count) %>%
  arrange(date) %>%                           # sort by date
  write.csv("coho_apportionment/output/global_imputed_transform.csv")

read.csv("coho_apportionment/output/local_imputed.csv") %>%
  distinct(date, year, .keep_all = TRUE) %>%  # remove duplicates
  pivot_wider(
    id_cols = date,
    names_from = year,
    values_from = total_count) %>%
  arrange(date) %>%                           # sort by date
  write.csv("coho_apportionment/output/local_imputed_transform.csv")
