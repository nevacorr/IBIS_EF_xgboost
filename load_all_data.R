library(dplyr)
library(purrr)
library(stringr)

load_all_data <- function() {
  
  #############################
  #### Define directories #####
  #############################
  
  working_dir <- getwd()
  
  vol_infant_dir <- "/Users/nevao/R_Projects/IBIS_EF/"
  volume_infant_datafilename <- "final_df_for_xgboost.csv"
  
  #############################
  #### Load infant lobe volume + behavior ####
  #############################
  
  df_infant_dem_lobe <-
    load_and_clean_infant_volume_data_and_all_behavior(
      vol_infant_dir,
      volume_infant_datafilename
    )
  
  #############################
  #### Load school-age lobe volumes ####
  #############################
  
  vol_dir_SA <- "/Users/nevao/Documents/IBIS_EF/source data/Brain_Data/updated imaging_2-27-25/IBISandDS_VSA_Cerebrum_and_LobeParcel_Vols_v01.04_20250221"
  volume_SA_datafilename <- "IBISandDS-VSA_Lobe_Vols_v01.04_20250221.csv"
  
  tot_tiss_dir_SA <- "/Users/nevao/Documents/IBIS_EF/source data/Brain_Data/updated imaging_2-27-25/IBISandDS_VSA_TissueSeg_Vols_v01.04_20250221"
  tot_tiss_SA_datafilename <- "IBISandDS_VSA_TissueSeg_Vols_v01.04_20250221.csv"
  
  df_vsa_lobe <-
    load_and_clean_vsa_volume_data(
      vol_dir_SA,
      volume_SA_datafilename,
      tot_tiss_dir_SA,
      tot_tiss_SA_datafilename
    )
  
  #############################
  #### Load infant subcortical volumes ####
  #############################
  
  subcort_infant_dir <- "/Users/nevao/Documents/IBIS_EF/source data/Brain_Data/IBIS1&2_volumes_v3.13"
  
  df_infant_subcort <-
    load_infant_subcortical_data(subcort_infant_dir) %>%
    as.data.frame() %>%
    mutate(row_id = row_number()) %>%
    select(-row_id)
  
  #############################
  #### Load school-age subcortical volumes ####
  #############################
  
  subcort_vsa_dir <- "/Users/nevao/Documents/IBIS_EF/source data/Brain_Data/updated imaging_2-27-25/IBISandDS_VSA_Subcort_and_LV_Vols_v01.04_20250221"
  subcort_vsa_datafilename <- "IBISandDS_VSA_Subcort_and_LV_Vols_v01.04_20250221.csv"
  
  df_vsa_subcort <-
    load_vsa_subcortical_data(subcort_vsa_dir, subcort_vsa_datafilename)
  
  #############################
  #### Load cortical thickness ####
  #############################
  
  ct_vsa_dir <- "/Users/nevao/Documents/IBIS_EF/source data/Brain_Data/IBISandDS_VSA_SurfaceData_v01.02_20210809"
  ct_vsa_datafilename <- "IBISandDS_VSA_CorticalThickness_DKT_v01.02_20210708.csv"
  
  df_vsa_ct <-
    load_vsa_ct_sa_data(ct_vsa_dir, ct_vsa_datafilename, "CT")
  
  #############################
  #### Load surface area ####
  #############################
  
  sa_vsa_datafilename <- "IBISandDS_VSA_SurfaceArea_DKT_v01.02_20210708.csv"
  
  df_vsa_sa <-
    load_vsa_ct_sa_data(ct_vsa_dir, sa_vsa_datafilename, "SA")
  
  #############################
  #### Load VSA DTI data ####
  #############################
  
  dti_vsa_dir <- paste0(
    "/Users/nevao/Documents/IBIS_EF/source data/Brain_Data/updated imaging_2-27-25/",
    "IBISandDS_VSA_DTI_Siemens_CMRR_v02.02_20250227/Siemens_CMRR"
  )
  
  metric_files <- list(
    FA = "IBISandDS_VSA_DTI_SiemensAndCMRR_FiberAverage_AD_v02.02_20250227.csv",
    AD = "IBISandDS_VSA_DTI_SiemensAndCMRR_FiberAverage_FA_v02.02_20250227.csv",
    MD = "IBISandDS_VSA_DTI_SiemensAndCMRR_FiberAverage_MD_v02.02_20250227.csv",
    RD = "IBISandDS_VSA_DTI_SiemensAndCMRR_FiberAverage_RD_v02.02_20250227.csv"
  )
  
  dti_dfs <- imap(metric_files, function(filename, metric) {
    
    df <- load_and_clean_vsa_dti_data(dti_vsa_dir, filename)
    
    df %>%
      rename_with(
        ~ paste0(metric, "_", .x),
        .cols = -CandID
      )
  })
  
  # Drop duplicate CandID columns except first
  dti_dfs <- c(
    list(dti_dfs[[1]]),
    lapply(dti_dfs[-1], function(df) df %>% select(-CandID))
  )
  
  df_vsa_dti <- bind_cols(dti_dfs)
  
  #############################
  #### Combine all data ####
  #############################

  #--------Two datasets have duplicate CandIDs. For df_infant_dem_lobe, one of the
  #--------ID rows is kept because the rows are identical. For df_vsa_dti
  #--------all of those rows are removed because the same CandID has different data
  
  # Get all duplicated CandID values for df_infant_dem_lobe
  duplicated_ids <- df_infant_dem_lobe$CandID[duplicated(df_infant_dem_lobe$CandID)]
  # Show unique duplicates
  unique_dup <- unique(duplicated_ids[!is.na(duplicated_ids)])
  print(unique_dup)
  
  # Remove exact duplicate rows in df_infant_dem_lobe (4 rows)
  df_infant_dem_lobe <- df_infant_dem_lobe %>%
    distinct()
  # Get all duplicated CandID values
  duplicated_ids <- df_infant_dem_lobe$CandID[duplicated(df_infant_dem_lobe$CandID)]
  # Show unique duplicates
  unique_dup <- unique(duplicated_ids[!is.na(duplicated_ids)])
  print(unique_dup)

  #--------------------------------
  # Get all duplicated CandID values for df_vsa_dti
  duplicated_ids <- df_vsa_dti$CandID[duplicated(df_vsa_dti$CandID)]
  # Show unique duplicates
  unique_dup <- unique(duplicated_ids[!is.na(duplicated_ids)])
  print(unique_dup)

  # Remove 38 subjects with duplicate CandIDs in df_vsa_dti (row values not identical)
  df_vsa_dti <- df_vsa_dti %>%
    filter(!CandID %in% unique_dup)
  # Get all duplicated CandID values
  duplicated_ids <- df_vsa_dti$CandID[duplicated(df_vsa_dti$CandID)]
  # Show unique duplicates
  unique_dup <- unique(duplicated_ids[!is.na(duplicated_ids)])
  print(unique_dup)
  
  #---------------Now combine all dataframes
  
  dfs_list <- list(
    df_infant_dem_lobe,
    df_vsa_lobe,
    df_infant_subcort,
    df_vsa_subcort,
    df_vsa_ct,
    df_vsa_sa,
    df_vsa_dti
  )
  
  dfs_combined <- reduce(
    dfs_list,
    ~ full_join(.x, .y, by = "CandID")
  )
  
  dfs_combined <- dfs_combined %>%
    mutate(CandID = as.integer(CandID))
  
  #############################
  #### Remove rows with no brain data ####
  #############################
  
  non_brain_cols <- c(
    "CandID", "Identifiers", "Combined_ASD_DX", "Risk", "Sex",
    "AB_12_Percent", "AB_24_Percent",
    "BRIEF2_GEC_T_score", "BRIEF2_GEC_raw_score",
    "Flanker_Standard_Age_Corrected", "DCCS_Standard_Age_Corrected",
    "Group", "Group_HR+", "Group_HR-", "Group_LR-"
  )
  
  cols_to_check <- setdiff(colnames(dfs_combined), non_brain_cols)
  
  dfs_all <- dfs_combined %>%
    filter(!if_all(all_of(cols_to_check), is.na))
  
  #############################
  #### Remove rows with no behavior data ####
  #############################
  
  behav_cols <- c(
    "AB_12_Percent", "AB_24_Percent",
    "BRIEF2_GEC_T_score", "BRIEF2_GEC_raw_score",
    "Flanker_Standard_Age_Corrected",
    "DCCS_Standard_Age_Corrected"
  )
  
  df_all_brain_behav <- dfs_all %>%
    filter(!if_all(all_of(behav_cols), is.na))
  
  #############################
  #### Normalize by total tissue ####
  #############################
  
  df_all_brain_behav <- divide_columns_by_tottiss(df_all_brain_behav, df_infant_dem_lobe, "V12")
  df_all_brain_behav <- divide_columns_by_tottiss(df_all_brain_behav, df_infant_dem_lobe, "V24") 
  df_all_brain_behav <- divide_columns_by_tottiss(df_all_brain_behav, df_vsa_lobe, "VSA") 
  df_all_brain_behav <- divide_columns_by_tottiss(df_all_brain_behav, df_infant_subcort, "v12") 
  df_all_brain_behav <- divide_columns_by_tottiss(df_all_brain_behav, df_infant_subcort, "v24") 
  df_all_brain_behav <- divide_columns_by_tottiss(df_all_brain_behav, df_vsa_subcort, "VSA")
  df_all_brain_behav <- divide_columns_by_tottiss(df_all_brain_behav, df_vsa_sa, "VSA")
  
  #############################
  #### Drop tissue / ICV columns ####
  #############################
  
  df_all_brain_behav <- df_all_brain_behav %>%
    select(-matches("Tiss|ICV"))
  
  return(df_all_brain_behav)
}