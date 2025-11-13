# missfill
library(dplyr)
library(tidyr)
# load ----
source('code/missfill_functions.r')


# data ----
# data inputs are date (mm/dd/yyyy) and weir count
read_csv('data/example.csv')-> situk

df_long <- situk %>%
  pivot_longer(
    cols = -date,               # pivot all columns except 'Date'
    names_to = "year",          # new column for year names
    values_to = "total_count"         # new column for values (0 or blank)
  )

df_long
# run functions ----

# format data
edit(impute_global(df_long))
edit(impute_local(chilkat))
edit(impute_local_improved(chilkat))

