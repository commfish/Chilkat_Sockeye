# missfile example code

# load ----
library(dplyr)
library(tidyr)
source('coho_apportionment/example/missfill_functions.r')

# data ----
read.csv('coho_apportionment/example/example.csv')-> situk

# pivot year columns into long format for import into functions
df_long <- situk %>%
  pivot_longer(
    cols = -date,    # all year columns now start with 'x'
    names_to = "year",
    values_to = "total_count") %>%
  mutate(year = as.numeric(sub("X", "", year)))

# run functions ----
impute_global(df_long)
impute_local(df_long)
#impute_local_improved(df_long)

# read and save output in wide format----
read.csv("coho_apportionment/example/global_imputed.csv") %>%
  distinct(date, year, .keep_all = TRUE) %>%  # remove duplicates
  pivot_wider(
    id_cols = date,
    names_from = year,
    values_from = total_count) %>%
  arrange(date) %>%                           # sort by date
  write.csv("coho_apportionment/example/global_imputed_transform.csv")


read.csv("coho_apportionment/example/local_imputed.csv") %>%
  distinct(date, year, .keep_all = TRUE) %>%  # remove duplicates
  pivot_wider(
    id_cols = date,
    names_from = year,
    values_from = total_count) %>%
  arrange(date) %>%                           # sort by date
  write.csv("coho_apportionment/example/local_imputed_transform.csv")
