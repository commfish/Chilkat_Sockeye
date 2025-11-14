# these functions were created by Justin Priest and modified by Sara Miller
# assign stat week
statweek <- function(x) {
  as.numeric(format(as.Date(x), "%U")) - as.numeric(format(as.Date(cut(x, "year")), "%U")) + 1
  # function modified from:
  # https://stackoverflow.com/questions/17286983/calculate-statistical-week-starting-1st-january-as-used-in-fisheries-data
  # example usage: dataframe %>% mutate(week = statweek(datecolumn))
}

# unsummarize data
duplicaterows <- function(dataframename, duplicatecolname = "specimen_count", replacenaswithone = FALSE){
  # Use this function if you encounter already summarized data and need to split it into long format
  # For example, you might have ASL data that has a row with Sex=M, Length=687, Specimen Count=3
  # This means that there are 3 fish that are males all with the same length of 687
  # Often though, we will want to run stats on our data which require it to NOT be summarized
  # This function takes summarized data and adds rows based on the count column

  # Using replacenaswithone allows you to decide whether an NA in the count column were really a count of 1
  # Carefully use this because most often an NA is a true zero.

  require(tidyverse)

  # Make an index of which rows will be repeated, and how many times
  .dupcount <- dataframename %>% dplyr::select(duplicatecolname) %>% tibble::deframe()

  # NAs will normally make this fail. We can replace NAs though
  # This replaces NAs with 1. THIS IS A LARGE ASSUMPTION SO BE CAREFUL
  if(sum(is.na(.dupcount) > 0) && replacenaswithone == TRUE){
    .dupcount <- replace_na(.dupcount, 1)
  }

  # Now repeat this for every row to duplicate.
  # A specimen count of 1 will mean the row isn't duplicated; a count of 5, repeats the row 5 times
  dataframename[rep(1:nrow(dataframename), .dupcount), ] %>%
    dplyr::select(-duplicatecolname) # Removes the count row now that it is incorrect!

  # Use like so: duplicaterows(dataframename = newdf, duplicatecolname = "Number.of.Specimens")
  # Thanks to: https://stackoverflow.com/questions/29743691/duplicate-rows-in-a-data-frame-in-r
}

# summarize proportions
count_pct <- function(df) {
  # https://stackoverflow.com/questions/24576515/relative-frequencies-proportions-with-dplyr
  return(
    df %>%
      tally %>%
      mutate(n_pct = 100*n/sum(n))
  )
}

# impute data
# There are four functions here:
# impute_global() which imputes all NAs in all years iteratively
# impute_local() which imputes a 10-year rolling imputation (prev & following 5 years)
# impute_local_improved which is similar to impute_local() but accounts for early years better

impute_global <- function(dfname, Year_column="year", DateName_column="date",
                         outputname = "globalimpute", # Only used if Step 3 turned "on"
                         Count_column = "total_count"){
  ### SUMMARY: Global Impute ###
  # This creates a dataframe that imputes NA values.
  # This algorithm interpolates across rows and columns, following Blick
  # In essence, imputing across rows (years) and columns (streams) allows for an NA in
  #  a year/stream to be informed by typical counts and for that year AND stream
  # This function can be easily modified to auto-create a named dataframe

  # Make sure that all NAs are present (a missing row is NOT same as a row with an NA)


  ### EXAMPLE USAGE ###
  # impute_global(ktn_index, Year_column="year")

  # Step 1: Set up dataframe to impute
  require(dplyr)
  .test <- dfname %>% rename(year = Year_column, date = DateName_column, total_count = Count_column)
  .test <- .test %>% dplyr::select(year, date, total_count)
  .test <- .test %>% mutate(imputed = is.na(total_count))

  # Step 2: Use multiplicative imputation as per Blick, in an iterative procedure
  j=1
  repeat{
    for(i in 1:nrow(.test)){
      .temprow = .test[i,]

      if(.temprow$imputed == TRUE){
        .sumyr = sum((.test %>% filter(year == .temprow$year) )$total_count, na.rm = TRUE)
        .sumrvr = sum((.test %>% filter(date == .temprow$date) )$total_count, na.rm = TRUE)
        .sumall = sum(.test$total_count, na.rm = TRUE)
        .test$total_count[i] = .sumyr * .sumrvr / .sumall
        # this interpolates across rows and columns
      }
    }
    j=j+1
    if(j>100){break} # repeat the above 100 times
  }
  print(.test)

  write.csv(.test, "coho_apportionment/example/global_imputed.csv")
}

impute_local <- function(dfname, Year_column="year", DateName_column="date",
                         Count_column = "total_count"){
  ### SUMMARY: 10-yr Localized Imputation ###
  # This takes a dataframe with NA values and imputes missing data
  # This algorithm uses "local" imputation: only 5 years before and after impute a missing value
  # i.e., only using the preceding 5 years and following 5 years
  # Make sure that all NAs are present (a missing row is NOT same as a row with an NA)

  ### EXAMPLE USAGE ###
  # impute_local(ktn_index, Year_column="year")

  # Step 1: Set up dataframe to impute
  require(dplyr)
  .test <- dfname %>% rename(year = Year_column, date = DateName_column, total_count = Count_column)
  .test <- .test %>% dplyr::select(year, date, total_count)
  .test <- .test %>% mutate(imputed = is.na(total_count))

  # Step 2: Use multiplicative imputation as per Blick, in an iterative procedure

  j=1
  repeat{
    for(i in 1:nrow(.test)){
      .temprow = .test[i,]

      if(.temprow$imputed == TRUE){
        .yr_range = .test %>% filter(between(year, .temprow$year - 5, .temprow$year + 5)) #5 yrs before / after
        .sumyr = sum((.yr_range %>% filter(year == .temprow$year) )$total_count, na.rm = TRUE)
        .sumrvr = sum((.yr_range %>% filter(date == .temprow$date) )$total_count, na.rm = TRUE)
        .sumall = sum(.yr_range$total_count, na.rm = TRUE)
        .test$total_count[i] = .sumyr * .sumrvr / .sumall
        # this is multiplicative imputation as per Blick
      }
    }
    j=j+1
    if(j>100){break} # repeat the above 100 times
  }
  print(.test)

  write.csv(.test, "coho_apportionment/example/local_imputed.csv")

}

impute_local_improved <- function(dfname, Year_column="year", DateName_column="date",
                                  Count_column = "total_count"){
  ### SUMMARY: 10-yr Localized Imputation, improved ###
  # This takes a dataframe with NA values and imputes missing data
  # This algorithm uses "local" imputation: only 5 years before and after impute a missing value
  # i.e., only using the preceding 5 years and following 5 years
  # However this version adds a rule for early years (1987-1996) to use 10 next years (10 yr minimum)
  # Make sure that all NAs are present (a missing row is NOT same as a row with an NA)

  ### EXAMPLE USAGE ###
  # impute_local_improved(ktn_index, Year_column="year")

  # Step 1: Set up dataframe to impute
  require(dplyr)
  .test <- dfname %>% rename(year = Year_column, date = DateName_column, total_count = Count_column)
  .test <- .test %>% dplyr::select(year, date, total_count)
  .test <- .test %>% mutate(imputed = is.na(total_count))

  # Step 2: Use multiplicative imputation as per Blick, in an iterative procedure

  j=1
  repeat{
    for(i in 1:nrow(.test)){
      .temprow = .test[i,]

      if(.test$year < 1997){
        if(.temprow$imputed == TRUE){
          .test_early <- .test %>% filter(year < 1997)
          .sumyr = sum((.test_early %>% filter(year == .temprow$year) )$total_count, na.rm = TRUE)
          .sumrvr = sum((.test_early %>% filter(date == .temprow$date) )$total_count, na.rm = TRUE)
          .sumall = sum(.test_early$total_count, na.rm = TRUE)
          .test$total_count[i] = .sumyr * .sumrvr / .sumall
        } # end early
      }}
    j=j+1
    if(j>50){break} # repeat the above 50 times. Needs to be iterative (imputing depends on other imputed values)
  }# end early

  j=1
  repeat{
    for(i in (nrow(.test %>% filter(year < 1997))+1):nrow(.test)){
      .temprow = .test[i,]
      if(.temprow$imputed == TRUE){
        .yr_range = .test %>% filter(between(year, .temprow$year - 5, .temprow$year + 5)) #5 yrs before / after
        .sumyr = sum((.yr_range %>% filter(year == .temprow$year) )$total_count, na.rm = TRUE)
        .sumrvr = sum((.yr_range %>% filter(date == .temprow$date) )$total_count, na.rm = TRUE)
        .sumall = sum(.yr_range$total_count, na.rm = TRUE)
        .test$total_count[i] = .sumyr * .sumrvr / .sumall
        # this is multiplicative imputation as per Blick
      }
    }
    j=j+1
    if(j>50){break} # repeat the above 50 times. Needs to be iterative (imputing depends on other imputed values)
  } # end late
  print(.test)
  write.csv(.test, "coho_apportionment/example/local_imputed_improved_matrix.csv")
}



