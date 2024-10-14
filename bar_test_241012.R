
bar_test_241012 <- function(df, dx_pre, messages=FALSE) {

  require(tidyverse)
  starttime=Sys.time()

  # Verify input
    if(!is.data.frame(df)) stop("First argument must be a dataframe")
    if(NROW(df) == 0) stop("Data file contains no observations. It must contain at least one row")
    if(!is.character(dx_pre)) stop("Second argument must be a character string")
  # Ensure dx_pre is a valid variable name
    if(make.names(dx_pre) != dx_pre) stop("Second argument must be a valid variable name in R")
  # Check if user entered a correct prefix for the diagnosis code variables in the input file
  # Determine how many diagnosis code variables there are in the data
    regex_dx <- paste0("^", dx_pre, "([0-9]+)$")
    dx_colnames <- grep(regex_dx, names(df), value = TRUE)
  # Replace full column name with first capture group and convert to number
    dx_nums <- as.numeric(sub(regex_dx, "\\1", dx_colnames))
    num_dx <- length(dx_nums)
    if(num_dx == 0) stop("No variables with prefix found in data")

  # Make sure df is not a tibble and if it is convert back to regular dataframe
    df <- data.frame(df)

    btab <- i10_map_bar.csv
    btab <- btab[ , c("dx","cell","PsCell")]

  #---------------------------------------------------------------------------------------#
  #  Merge diagnosis code variables with Barell Code reference table to obtain cell name  #
  #  and cell survival probability for each diagnosis code and add them to the data       #
  #---------------------------------------------------------------------------------------#
  for(i in dx_nums){

    if(messages==TRUE){
      message("Determining cell for Diagnosis ", i, " of ", num_dx)
      }

    # Create column name
      dx_name <- paste0(dx_pre, i)

    # Pull just the diagnosis code column of interest
      df_ss <- df[ , dx_name, drop = FALSE]

    # Add row variable for sorting back to original order
      df_ss$n <- 1:NROW(df_ss)

    # Strip out decimal in all codes
      df_ss[ , dx_name] <- sub("\\.", "", df_ss[ , dx_name])

    # Get rid of codes that do not start with a valid character
      i10_valid <- c("S","T")
      df_ss[ , dx_name] <- ifelse(substr(df_ss[,dx_name], 1, 1) %in% c(i10_valid), df_ss[,dx_name], NA)

    # Merge with lookup table for severity
      temp <- merge(df_ss, btab, by.x = dx_name, by.y = "dx", all.x = TRUE, all.y = FALSE, sort = FALSE)

    # Reorder rows after merge
      temp <- temp[order(temp$n), ]

    # Reorder columns and drop dx and n
      temp <- temp[ , c("cell","PsCell")]

    # Rename columns
      names(temp) <- paste0(c("cell_","PsCell_"), i)

    # Add temp columns to dataframe
      df <- .insert_columns(df, dx_name, temp)

    }   #END FOR LOOP (i in dxno)


  #---------------------------------------------------------------------#
  # Add mortality prediction for ICD-10-cm codes from cells #
  #---------------------------------------------------------------------#


    if(messages==TRUE){
      mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
      message("Time elapsed ", mindiff, " minutes")
      message("Calculating mortality predictions")
    }

   coef_df <- select(btab,dx,PsCell)

   # Create hash table
     coef_df <- coef_df[!is.na(coef_df$PsCell), ]
     effect_hash <- coef_df$PsCell
     names(effect_hash) <- coef_df$dx
     calc_mortality_prediction <- function(dx){
       # dx is a character vector of diagnosis codes for one person
       x <- prod(effect_hash[sub("\\.", "", dx)], na.rm = TRUE)
        }
    mat <- as.matrix(df[,grepl(paste0("^", dx_pre), names(df))])
    df$bPSprod <- apply(mat, 1, calc_mortality_prediction)

    # Create hash table
    coef_df <- coef_df[!is.na(coef_df$PsCell), ]
    effect_hash <- coef_df$PsCell
    names(effect_hash) <- coef_df$dx
    calc_mortality_prediction <- function(dx){
      # dx is a character vector of diagnosis codes for one person
      x <- max(effect_hash[sub("\\.", "", dx)], na.rm = TRUE)
      x <- if_else( (x>1|x<0),1,x )
    }
    mat <- as.matrix(df[,grepl(paste0("^", dx_pre), names(df))])
    df$bPSmin <- apply(mat, 1, calc_mortality_prediction)


  # Set rownames
    rownames(df) <- 1:nrow(df)
    if(messages==TRUE){
      mindiff=round(as.double(difftime(Sys.time(),starttime,units="secs"))/60)
      message("Time elapsed ", mindiff, " minutes")
     }

  # Return dataframe
    df

} #END bar_test
