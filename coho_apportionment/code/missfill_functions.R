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
  
  # Using replace nas with one allows you to decide whether an NA in the count column were really a count of 1
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
# There are three functions here:
# impute_global() which imputes all NAs in all years iteratively
# impute_local() which imputes a 10-year rolling imputation (prev & following 5 years)
# impute_local_improved which is similar to impute_local() but accounts for early years better

# global impute #
impute_global <- function(dfname, Year_column = "year", DateName_column = "date",
                          Count_column = "total_count", Season_column = NULL,
                          max_iter = 100, tol = 1e-8) {
  
  # Imputes NAs using multiplicative imputation as per Blick, in an iterative procedure,
  # using all years. Imputing across rows (years) and columns (dates) lets an NA be
  # informed by the typical value for that year AND that date.
  # If Season_column is given (logical column, TRUE = within that year's weir season),
  # only in-season NAs are imputed; out-of-season NAs stay NA and are ignored in the sums.
  # Make sure all date-year combinations are present (a missing row is NOT the same as an NA).
  # Returns a data frame with columns year, date, total_count, imputed.
  
  ### EXAMPLE USAGE ###
  # imputed <- impute_global(df_long, Season_column = "in_season")
  
  require(dplyr)
  
  # Step 1: Set up dataframe to impute
  .test <- dfname %>%
    rename(year = all_of(Year_column), date = all_of(DateName_column),
           total_count = all_of(Count_column))
  
  if (is.null(Season_column)) {
    .test$in_season <- TRUE
  } else {
    .test <- .test %>% rename(in_season = all_of(Season_column))
  }
  
  .test <- .test %>%
    mutate(imputed = is.na(total_count) & in_season) %>%
    select(year, date, total_count, imputed)
  
  .rows <- which(.test$imputed)
  if (length(.rows) == 0) return(.test)
  
  # Step 2: Multiplicative imputation, repeated until values stop changing
  # (iterative because imputed values feed into later imputations)
  for (j in seq_len(max_iter)) {
    .old <- .test$total_count[.rows]
    
    for (i in .rows) {
      .sumyr  <- sum(.test$total_count[.test$year == .test$year[i]], na.rm = TRUE)
      .sumrvr <- sum(.test$total_count[.test$date == .test$date[i]], na.rm = TRUE)
      .sumall <- sum(.test$total_count, na.rm = TRUE)
      .test$total_count[i] <- .sumyr * .sumrvr / .sumall
    }
    
    if (!anyNA(.old) && max(abs(.test$total_count[.rows] - .old)) < tol) {
      message("impute_global converged after ", j, " iterations")
      break
    }
    if (j == max_iter) warning("impute_global did not converge in ", max_iter, " iterations")
  }
  
  return(.test)
}

# 10-yr Localized Imputation #
impute_local <- function(dfname, Year_column = "year", DateName_column = "date",
                         Count_column = "total_count", Season_column = NULL,
                         window = 5, max_iter = 100, tol = 1e-8) {
  
  # Imputes NAs using multiplicative imputation as per Blick, in an iterative procedure,
  # using only years within +/- `window` years of the year being imputed
  # (default: preceding 5 years and following 5 years).
  # If Season_column is given (logical column, TRUE = within that year's weir season),
  # only in-season NAs are imputed; out-of-season NAs stay NA and are ignored in the sums.
  # Make sure all date-year combinations are present (a missing row is NOT the same as an NA).
  # Returns a data frame with columns year, date, total_count, imputed.
  
  ### EXAMPLE USAGE ###
  # imputed <- impute_local(df_long, Season_column = "in_season")
  
  require(dplyr)
  
  # Step 1: Set up dataframe to impute
  .test <- dfname %>%
    rename(year = all_of(Year_column), date = all_of(DateName_column),
           total_count = all_of(Count_column))
  
  if (is.null(Season_column)) {
    .test$in_season <- TRUE
  } else {
    .test <- .test %>% rename(in_season = all_of(Season_column))
  }
  
  .test <- .test %>%
    mutate(imputed = is.na(total_count) & in_season) %>%
    select(year, date, total_count, imputed)
  
  .rows <- which(.test$imputed)
  if (length(.rows) == 0) return(.test)
  
  # Step 2: Multiplicative imputation, repeated until values stop changing
  # (iterative because imputed values feed into later imputations)
  for (j in seq_len(max_iter)) {
    .old <- .test$total_count[.rows]
    
    for (i in .rows) {
      .yr_range <- .test %>%
        filter(between(year, .test$year[i] - window, .test$year[i] + window))
      .sumyr  <- sum(.yr_range$total_count[.yr_range$year == .test$year[i]], na.rm = TRUE)
      .sumrvr <- sum(.yr_range$total_count[.yr_range$date == .test$date[i]], na.rm = TRUE)
      .sumall <- sum(.yr_range$total_count, na.rm = TRUE)
      .test$total_count[i] <- .sumyr * .sumrvr / .sumall
    }
    
    if (!anyNA(.old) && max(abs(.test$total_count[.rows] - .old)) < tol) {
      message("impute_local converged after ", j, " iterations")
      break
    }
    if (j == max_iter) warning("impute_local did not converge in ", max_iter, " iterations")
  }
  
  return(.test)
}

# 10-yr Localized Imputation, improved #
impute_local_improved <- function(dfname, Year_column="year", DateName_column="date",
                                  Count_column = "total_count"){
  
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
      
      if(.temprow$year < 1997){
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
  write.csv(.test, "coho_apportionment/output/local_imputed_improved_matrix.csv")
}

# Depends on dplyr
tickr <- function(
    data, # dataframe
    var, # column of interest
    to # break point definition 
){
  
  VAR <- enquo(var) # makes VAR a dynamic variable
  
  data %>% 
    distinct(!!VAR) %>%
    #    ungroup(!!VAR) %>% 
    mutate(labels = ifelse(!!VAR %in% seq(to * round(min(!!VAR) / to), max(!!VAR), to),
                           !!VAR, "")) %>%
    dplyr::select(breaks = UQ(VAR), labels)
}