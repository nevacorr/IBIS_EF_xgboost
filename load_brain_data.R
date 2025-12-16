library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)

reshape_dataframe <- function(df) {
  
  # Drop Combined_ID
  df <- df %>% select(-Combined_ID)
  
  # Identify region columns to pivot
  region_columns <- names(df)[
    !str_detect(names(df), "_VQC|_ExcludeReason|Edited") &
      !names(df) %in% c("DCCID", "Visit")
  ]
  
  # Pivot: make separate columns for each Visit
  df_long <- df %>%
    pivot_longer(
      cols = all_of(region_columns),
      names_to = "region",
      values_to = "value"
    )
  
  # Make a temporary column for Visit with lowercase "v"
  df_long <- df_long %>%
    mutate(Visit_lower = sub("^V", "v", Visit))
  
  df_wide <- df_long %>%
    select(DCCID, Visit_lower, region, value) %>%  # drop VQC/ExcludeReason
    pivot_wider(
      id_cols = DCCID,
      names_from = c(region, Visit_lower),
      values_from = value,
      names_sep="_"
    ) 
  
  return(df_wide)
}

load_infant_subcortical_data <- function(filepath) {
  
  # List all CSV files in the directory that are subcortical (exclude LobeParcel)
  csv_files <- list.files(filepath, pattern = "\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!str_detect(csv_files, "LobeParcel")]
  
  # Load CSVs into a named list of data frames
  dfs <- list()
  for (file in csv_files) {
    fname <- basename(file)
    match <- str_match(fname, "IBIS_v3\\.13_([^_]+)")
    key <- if (!is.na(match[2])) match[2] else fname
    df <- read_csv(file, show_col_types = FALSE)
    
    # Reshape each dataframe
    dfs[[key]] <- reshape_dataframe(df)
  }
  
  # Exclude "ICV" from merging
  dfs_to_merge <- dfs[names(dfs) != "ICV"]
  
  # Merge all dataframes by DCCID
  subcort_merged_df <- reduce(dfs_to_merge, full_join, by = "DCCID")
  
  # Remove columns containing "Edited"
  subcort_merged_df <- subcort_merged_df %>% select(-matches("Edited"))
  # 
  # # Normalize _v12 columns by totTiss_v12
  # if ("totTiss_v12" %in% names(subcort_merged_df)) {
  #   v12_cols <- names(subcort_merged_df)[str_detect(names(subcort_merged_df), "_v12") & 
  #                                          !str_detect(names(subcort_merged_df), "totTiss")]
  #   subcort_merged_df <- subcort_merged_df %>%
  #     mutate(across(all_of(v12_cols), ~ .x / totTiss_v12))
  # }
  # 
  # # Normalize _v24 columns by totTiss_v24
  # if ("totTiss_v24" %in% names(subcort_merged_df)) {
  #   v24_cols <- names(subcort_merged_df)[str_detect(names(subcort_merged_df), "_v24") & 
  #                                          !str_detect(names(subcort_merged_df), "totTiss")]
  #   subcort_merged_df <- subcort_merged_df %>%
  #     mutate(across(all_of(v24_cols), ~ .x / totTiss_v24))
  # }
  
  # Drop totTiss columns 
  subcort_merged_df <- subcort_merged_df %>%
    select(-any_of(c("totTiss_v12", "totTiss_v24")))
  
  # Rename DCCID to CandID
  subcort_merged_df <- subcort_merged_df %>% rename(CandID = DCCID)
  
  return(subcort_merged_df)
}

load_vsa_subcortical_data <- function(filepath, datafilename) {
  
  df <- read_csv(file.path(filepath, datafilename), show_col_types = FALSE)
  
  unwanted_substrings <- c(
    "PassFail", "Exclude_Reason", "score",
    "Visit", "Septum", "Ventricle"
  )
  
  # Remove columns in unwanted substrings, replace . with Nan and convert
  #values to numeric
  df <- df %>%
    select(-any_of(names(df)[str_detect(names(df), paste(unwanted_substrings, collapse = "|"))])) %>%
    mutate(across(where(is.character), ~ na_if(.x, "."))) %>%
    mutate(across(where(is.character), as.numeric))
  
  df <- df %>%
    rename_with(~ ifelse(.x == "CandID", .x, paste0(.x, "_VSA")))
  
  return(df)
}

load_vsa_ct_sa_data <- function(filepath, datafilename, datastr) {
  
  # Read the CSV
  df <- read_csv(file.path(filepath, datafilename), show_col_types = FALSE)
  
  # Define unwanted substrings
  unwanted_substrings <- c("Visit")
  
  # Remove columns containing any of the unwanted substrings, convert '.' to 
  # NA and convert character columns to numeric
  df <- df %>%
    select(-any_of(names(df)[str_detect(names(df), paste(unwanted_substrings, collapse = "|"))])) %>%
    mutate(across(where(is.character), ~ na_if(.x, "."))) %>%
    mutate(across(where(is.character), as.numeric))
  
  # Rename columns except CandID
  df <- df %>%
    rename_with(~ ifelse(.x == "CandID", .x, paste0(.x, "_", datastr, "_VSA")))
  
  return(df)
}

library(dplyr)
library(readr)
library(stringr)

load_and_clean_vsa_dti_data <- function(dir, datafilename) {
  
  # Read CSV
  dti_df <- read_csv(file.path(dir, datafilename), show_col_types = FALSE)
  
  # Drop unwanted columns
  dti_df <- dti_df %>%
    select(-Visit_label, -CandID_Visit, -dMRI_protocol, -FileID) %>%
    select(-matches("Optic|Motor|Fornix|CorticoSpinal|UNC|Reticular|ILF|CT|CF|CC|Temporo|IFOF"))
  
  # Replace '.' and empty strings with NA
  dti_df <- dti_df %>%
    mutate(across(
      where(is.character),
        ~ {
          x <- na_if(.x, ".")
          x <- na_if(x, "")
          as.numeric(x)
        }
      )
    )
  
  # Convert all non-CandID columns to numeric
  dti_df <- dti_df %>%
    mutate(across(-CandID, as.numeric))
  
  return(dti_df)
}

library(dplyr)
library(readr)
library(stringr)
library(fastDummies)

load_and_clean_infant_volume_data_and_all_behavior <- function(filepath, filename) {
  
  df <- read_csv(file.path(filepath, filename), show_col_types = FALSE)
  
  # Encode Sex: Female = 0, Male = 1
  df <- df %>%
    mutate(Sex = recode(Sex, Female = 0, Male = 1))
  
  # One-hot encode Group column (keep NAs)
  df <- df %>%
    fastDummies::dummy_cols(select_columns = "Group",
                            remove_selected_columns = FALSE,
                            ignore_na = FALSE,
                            remove_first_dummy = FALSE)
  
  # Remove the Group_NA column if it exists
  df <- df %>% select(-matches("Group_NA"))
  
  # Find new group columns
  group_cols <- grep("^Group_", names(df), value = TRUE)
  
  # Reorder columns: first 5, then Group dummies, then rest
  group_cols <- grep("^Group_", names(df), value = TRUE)
  other_cols <- setdiff(names(df), group_cols)
  
  new_col_order <- c(other_cols[1:5], group_cols, other_cols[6:length(other_cols)])
  df <- df %>% select(all_of(new_col_order))
  
  return(df)
}

load_and_clean_vsa_volume_data <- function(
    dir, datafilename,
    tissue_dir, tissue_datafilename
) {
  
  df <- read_csv(file.path(dir, datafilename), show_col_types = FALSE)
  
  unwanted_substrings <- c(
    "PassFail", "Cerebellum",
    "ExcludeReason", "score", "Visit"
  )
  
  # Remove columns containing any of the unwanted substrings, convert '.' to 
  # NA and convert character columns to numeric
  df <- df %>%
    select(-any_of(names(df)[str_detect(names(df), paste(unwanted_substrings, collapse = "|"))])) %>%
    mutate(across(where(is.character), ~ na_if(.x, "."))) %>%
    mutate(across(where(is.character), as.numeric))
  
  # Load total intracranial volume and total tissue volume
  tissue_df <- read_csv(file.path(tissue_dir, tissue_datafilename), show_col_types = FALSE) %>%
    select(CandID, ICV_vol, total_Tissue_vol)
  
  # Merge volume and intracranial volume dataframes
  merged_df <- full_join(df, tissue_df, by = "CandID")
  
  # Append string _VSA to all columns that are not CandID
  merged_df <- merged_df %>%
    rename_with(~ ifelse(.x == "CandID", .x, paste0(.x, "_VSA")))
  
  return(merged_df)
}


