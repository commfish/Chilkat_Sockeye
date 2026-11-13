# coho apportionment (missfill) for Chilkat sockeye weir
# this file is used to apportion the coho/sockeye counts at the Chilkat Lake weir
# when seine is unavailable for an extended amount of time, it imputes missing daily coho ratios (coho / total weir count) within each season
# script written by Sara Miller (sara.miller@alaska.gov)
# last updated: November 2026

# load ----
library(dplyr)
library(tidyr)
library(readr)
source('coho_apportionment/code/missfill_functions.r')

# settings ----
method <- "global"   # "global" (all years) or "local" (5 years before and after)

# data ----
chilkat <- read.csv('coho_apportionment/data/chilkat_coho_ratio.csv',
                    na.strings = c("", "NA", "-"))

# pivot year columns into long format for import into functions
# (read.csv adds an "X" in front of the year column names, e.g. X2018)
df_long <- chilkat %>% 
  pivot_longer(
    cols = -date,
    names_to = "year",
    values_to = "coho_ratio",
    names_prefix = "X",
    names_transform = list(year = as.integer)) %>%
  arrange(year, date)

# flag days within each year's season (first to last day with data);
# only missing days inside the season are imputed, days outside stay NA
df_long <- df_long %>%
  group_by(year) %>%
  mutate(in_season = date >= min(date[!is.na(coho_ratio)]) &
           date <= max(date[!is.na(coho_ratio)])) %>%
  ungroup()

# report the gaps that will be imputed
df_long %>%
  filter(in_season, is.na(coho_ratio)) %>%
  group_by(year) %>%
  summarise(n_days = n(), first_day = min(date), last_day = max(date)) %>%
  print()

# impute ----
imputed <- switch(method,
                  global = impute_global(df_long, Count_column = "coho_ratio", Season_column = "in_season"),
                  local  = impute_local(df_long,  Count_column = "coho_ratio", Season_column = "in_season"),
                  stop("method must be 'global' or 'local'"))

imputed <- imputed %>% rename(coho_ratio = total_count)

# imputed ratios should be proportions
if (any(imputed$coho_ratio < 0 | imputed$coho_ratio > 1, na.rm = TRUE)) {
  warning("some imputed coho ratios fall outside 0-1")
}

# imputed values
imputed %>% filter(imputed) %>% print(n = Inf)

# save output ----
# long format (includes the imputed TRUE/FALSE flag)
write.csv(imputed, paste0("coho_apportionment/output/", method, "_imputed.csv"),
          row.names = FALSE)

# wide format (one column per year), same layout as the input file
imputed %>%
  pivot_wider(id_cols = date, names_from = year, values_from = coho_ratio) %>%
  arrange(date) %>%
  write.csv(paste0("coho_apportionment/output/", method, "_imputed_transform.csv"),
            row.names = FALSE, na = "")   # blank cells for days with no data, like the input
