library(tidyverse)
#library(zoo)
library(bfastSpatial)
library(raster)
library(lubridate)
library(lemon)
library(rstatix)
library(plotrix)

# ---------------Data wrangling--------------------------
# Get the location of the file containing the NDVI time series exported from GEE
files <- list.files("GEE_csv/",
           "*.csv",
           full.names = T)

# Read data
df <- read.csv(files, na.strings = "-11000", stringsAsFactors = F)

# Subset the data for each of the 4 scenarios
TF_no <- df %>%
  filter(FrstTyp == "Temperate" & Change == "No")
TF_yes <- df %>%
  filter(FrstTyp == "Temperate" & Change == "Yes")
TDF_no <- df %>%
  filter(FrstTyp == "TropicalDry" & Change == "No")
TDF_yes <- df %>%
  filter(FrstTyp == "TropicalDry" & Change == "Yes")

# Add them as a list and name them
files <- list(TF_no, TF_yes, TDF_no, TDF_yes)
names(files) <- c("TF_no", "TF_yes", "TDF_no", "TDF_yes")

# Loop over the 4 files
resul <- lapply(files, function(j){
  
  # Read data
  df <- j
  
  # See id values
  unique(df$id)
  
  # Remove fid, shape and id columns
  df <- df[,-c(1, (ncol(df)-5):ncol(df))]
  
  # transpose and convert to data fram (from matrix)
  df <- t(df)
  df <- as.data.frame(df)
  
  dates <- row.names(df)
  
  # Add dates as first column
  df <- data.frame(dates = dates,df)
  row.names(df) <- seq(1,nrow(df),1)
  # Transform NDVI values to numeric
  df[,2:ncol(df)] <- as.numeric(unlist(df[,2:ncol(df)]))
  
  # Transform dates in date format
  df2 <- df %>%
    as_tibble() %>%
    mutate_at(vars(dates), function(x) ymd(str_extract(x,"[0-9]{8}"))) %>%
    # Order rows according to dates values
    arrange(dates) %>%
    as.data.frame()


  napercentage <- df2 %>%
    # Transform the characters inside the dates column to Date format
    mutate_at(vars(dates), function(x) ymd(x)) %>%
    # Change data to long format
    pivot_longer(cols = -dates, 
                 names_to = "var",
                 values_to = "value") %>%
    # Extract year, month and day from the date column
    mutate(year = year(dates),
           month = month(dates),
           day = day(dates)) %>%
    # Eliminate dates before 1993, as practically are only NA
    filter(year > 1993) %>%
    # Ungroup (I think this step is unnecessary)
    ungroup() %>%
    # Group by month
    group_by(year, var) %>%
    # Add a column that is set to 0 if value is NA (0 in the data) or not. If not, is set to 1
    mutate(Valor = ifelse(is.na(value), 0, 1)) %>%
    # Count values per month
    count(Valor) %>%
    # Get the number of total obervations
    mutate(total = sum(n)) %>%
    # Get proportion
    mutate(propVal = n / total) %>%
    # Ungroup
    ungroup() %>%
    # Stay with NA only
    filter(Valor == 0) %>%
    # Group by var
    group_by(var) %>%
    # Get proportion of NA
    summarise(propVal_pix = mean(propVal)) %>%
    select(propVal_pix)
  
  # Change column names to Time and NDVIX
  colnames(df2) <- c("Time", paste0("NDVI",seq(1,(ncol(df2)-1))))

  # Apply bfmPixel to extract r2, magnitude, amplitude and stable historic period length
  
  # Create rasters from ndvi time series
  img1 <- lapply(1:nrow(df2), function(i){
    raster(matrix(unlist(df2[i,2:ncol(df2)]),nrow = (ncol(df2)-1)))
  })
  
  # Create a stack
  stackedimg1 <- stack(img1)

  # Get number of cells in raster
  cell_num <- dim(stackedimg1)[1]
  
  # Loop over each cell in raster stack and apply bfmPixel 
  r2_df <- sapply(1:cell_num, function(i){
    resul_TDDF <- bfmPixel(stackedimg1, 
                           dates =  as.Date(df2$Time),
                           # Start Monitoring period
                           start=c(2016, 1), 
                           cell = i,
                           interactive=F)
    
    # Extract r2 and r2 adj values
    r2 <- summary(resul_TDDF$bfm$model)$r.squared
    r2_adj <- summary(resul_TDDF$bfm$model)$adj.r.squared
    
    magnitude <- resul_TDDF$bfm$magnitude
    
    # Copy model information into mod_model
    mod_model <- resul_TDDF$bfm$model
    # plot(resul_TDDF$bfm)
    
    # Set to zero the trend coefficient to "detrend" the Time series
    mod_model$coefficients["trend"] <- 0
    
    # Make predictions over the same data, but "detrended", i.e., trend coef set to zero
    detrend_vals <- predict.lm(mod_model)
    # plot(temp)
    
    # Calculate amplitude as max - min 
    # Remember this is 10 000 based NDVI values
    amplitude <- max(detrend_vals) - min(detrend_vals)
  
    # History length
    history_length = resul_TDDF$bfm$history[2] - resul_TDDF$bfm$history[1]

    # Create a data frame containing extracteds values
    data.frame(r2 = r2, 
               r2_adj = r2_adj,
               amplitude = amplitude,
               history_length = history_length,
               magnitude = magnitude)
  })
  
  # Convert list to data frame
  r2_df <- as.data.frame(r2_df) 
  
  # Calculate mean and SE for r2 and r2 adj
  r2_df %>%
    t() %>%
    as_tibble() %>%
    # Unnest gets values outside of list structure to data.frame style
    unnest(c(r2, r2_adj, amplitude, history_length, magnitude)) %>%
    # Add na percentage
    add_column(napercentage)

})

# Join the list of 4 files into a single data.frame
df_exp <- bind_rows(resul, .id = "file") %>%
  as_tibble() %>%
  # Put file info first
  dplyr::select(file, everything()) 

df4tests <- df_exp %>%
  # Get new column with type of forest
  mutate(forest = str_extract(file, "TDF|TF")) %>%
  # Get new column with yes or no
  mutate(verified = str_extract(file, "yes|no"))

df4tests %>%
  # Rescale magnitude and amplitde to 0-1 NDVI Values instead of 0 - 10 000
  mutate_at(vars(amplitude, magnitude), function(x) x /10000) %>%
  write.csv("DF_logistic_reg.csv", row.names = F)