#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FUNCTIONS TO DOWNLOAD USGS AND SWOT DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### READ ME ####
# Last Modified: 1/23/2026
# Contact: Merritt Harlan (mharlan@usgs.gov)

#This script downloads and processes USGS discharge measurement data and 
#SWOT river slope data. Additionally, Scale Specific Sinuosity data and 
#SWOT River Database (SWORD) data are compiled, and used to estimate 
#discharge resistance parameters using flow law equations. [SWOT: 
#Surface Water Ocean Topography satellite]

#Disclaimer (Provisional): "This software is preliminary or provisional and is 
#subject to revision. It is being provided to meet the need for timely best 
#science. The software has not received final approval by the U.S. Geological 
#Survey (USGS). No warranty, expressed or implied, is made by the USGS or the 
#U.S. Government as to the functionality of the software and related material 
#nor shall the fact of release constitute any such warranty. The software is 
#provided on the condition that neither the USGS nor the U.S. Government shall 
#be held liable for any damages resulting from the authorized or unauthorized 
#use of the software."

# ============================================================
# Load Libraries
# ============================================================
library(RNetCDF)
library(dplyr)
library(tidyr)
library(dataRetrieval)
library(data.table)
library(purrr)
library(httr)
library(jsonlite)
library(future.apply)
library(progressr)
library(lubridate)
library(pracma)
library(sf)
handlers(global = TRUE)

data_folder <- "data/"

# ============================================================
# Helper Functions
# ============================================================

# Convert English-unit ADCP values to SI units
convert_adcp_units <- function(df) {
  df %>%
    mutate(
      Q_cms = discharge_va * 0.0283168466,
      width_m = chan_width * 0.3048,
      depth_m = chan_depth * 0.3048,
      area_m2 = chan_area * 0.09290304,
      velocity_ms = chan_velocity * 0.3048
    )
}

# USGS ADCP retrieval
get_USGS_ADCP <- function(site_id) {
  raw <- readNWISmeas(site_id, tz = "UTC", expanded = TRUE)
  
  raw %>%
    select(
      measurement_dateTime, site_no, gage_height_va, discharge_va,
      chan_discharge, chan_width, chan_area, chan_velocity,
      measured_rating_diff, meas_type, streamflow_method, velocity_method
    ) %>%
    mutate(
      chan_depth = chan_area / chan_width,
      Fr = chan_velocity / sqrt(32.17405 * chan_depth),
      DD = decimal_date(measurement_dateTime)
    ) %>%
    drop_na()
}

# SWOT API retrieval
get_SWOT_reach_data <- function(reach_id = "73270200201", 
                                start_date = "2023-01-01", 
                                end_date = "2025-12-01", 
                                field_length = 24, 
                                fields = "reach_id,time_str,wse,slope,slope_u,slope_r_u,slope2,reach_q,reach_q_b,dark_frac,ice_clim_f,partial_f,n_good_nod,xovr_cal_q,p_length,p_width,p_wse,p_low_slp,xtrk_dist,sword_version,collection_shortname,cycle_id,pass_id,crid"){
  url = paste0("https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries?collection_name=SWOT_L2_HR_RiverSP_D&feature=Reach&feature_id=", 
               reach_id, 
               "&output=csv&start_time=", 
               start_date, "T00:00:00Z", 
               "&end_time=", end_date, "T00:00:00Z", 
               "&fields=", fields) 
  res = GET(url) 
  if(res$status_code == 200){ 
    res_text = fromJSON(rawToChar(res$content)) 
    data = read.csv(text = res_text$results$csv, stringsAsFactors = FALSE) %>% 
      filter(time_str != "no_data") %>% 
      select(1:field_length) %>% 
      mutate(across(where(is.numeric), ~ ifelse(. < -1E+11, NA, .))) 
    if(nrow(data) > 0){
      return(data) 
    } 
  } 
}

# Manning n computation
compute_manning <- function(df) {
  df %>%
    mutate(
      Manning_n = (width_m * (depth_m^1.67) * sqrt(median_slope)) / Q_cms,
      Manning_n_slope_sd = (width_m * (depth_m^1.67) * sqrt(sd_slope)) / Q_cms
    )
}

# ============================================================
# Read Priors + Map USGS Gages to Reach IDs
# ============================================================

# Citation: SWOT Discharge Algorithm Working Group (DAWG). 2024. SWOT Level 4 River Discharge Products Version 1. Ver. 1. PO.DAAC, CA, USA. Dataset accessed 2025-12-03 at https://doi.org/10.5067/SWOT-SOS-V1 
na_priors <- open.nc(file.path(data_folder, "na_sword_v16_SOS_unconstrained_0001_20240726T123358_priors.nc"))

#Accessed from https://www.swordexplorer.com/
SWORD_map <- fread(file.path(data_folder, "SWORD_v17b_to_v16_Translation/v17b/NA_ReachIDs_v17b_vs_v16.csv"))

USGS_group <- grp.inq.nc(na_priors, "USGS")$self

USGS_df <- data.frame(
  Agency_id = var.get.nc(USGS_group, "USGS_id"),
  reach_id  = var.get.nc(USGS_group, "USGS_reach_id"),
  Agency = "USGS",
  stringsAsFactors = FALSE
)

#Query site information (latitude and longitude)
siteNumbers1 <- unique(paste0("USGS-", USGS_df$Agency_id))[1:400]
siteINFO1 <- read_waterdata_monitoring_location(siteNumbers1)
siteNumbers2 <- unique(paste0("USGS-", USGS_df$Agency_id))[400:800]
siteINFO2 <- read_waterdata_monitoring_location(siteNumbers2)
siteNumbers3 <- unique(paste0("USGS-", USGS_df$Agency_id))[800:length(unique(USGS_df$Agency_id))]
siteINFO3 <- read_waterdata_monitoring_location(siteNumbers3)

site_locations = rbind(siteINFO1,siteINFO2, siteINFO3) %>%
  select(monitoring_location_id, monitoring_location_name)

longitude = st_coordinates(site_locations)[, "X"]
latitude = st_coordinates(site_locations)[, "Y"]

gage_locations = cbind(site_locations %>% st_drop_geometry(),longitude, latitude)

# Attach v17 reach IDs
USGS_df <- USGS_df %>%
  rename(v16_reach_id = reach_id) %>%
  mutate(v16_reach_id = as.character(v16_reach_id)) %>%
  left_join(SWORD_map %>% mutate(v16_reach_id = as.character(v16_reach_id)), by = "v16_reach_id") %>%
  filter(!is.na(Agency_id)) %>%
  select(Agency_id, lon, lat, v17_reach_id, v16_reach_id) %>%
  mutate(monitoring_location_id = paste0("USGS-", Agency_id)) %>%
  left_join(gage_locations) %>%
  distinct()

USGS_df %>%
  group_by(Agency_id) %>%
  summarise(count = n()) %>%
  filter(count > 1) %>%
  arrange(desc(count))

#remove duplicate gage/reach pairings from reach splits in SWORD v17b
USGS_df = USGS_df %>% filter(!v17_reach_id %in% c(74100600516, 72554000181, 81130600076))

# ============================================================
# Download USGS ADCP measurements
# ============================================================

site_list <- unique(USGS_df$Agency_id)

with_progress({
  p <- progressor(along = site_list)
  
  res <- future_lapply(
    site_list,
    function(id) {
      p()
      possibly(get_USGS_ADCP, otherwise = NULL)(id)
    }
  )

  res <- res[!sapply(res, is.null)]

})

USGS_gage_measurements <- bind_rows(res) %>%
  convert_adcp_units()

# ============================================================
# Download SWOT Reach Time Series
# ============================================================

plan(multisession, workers = parallel::detectCores() - 2)

reach_list <- unique(USGS_df$v17_reach_id)

SWOT_raw <- future_lapply(paste(reach_list), get_SWOT_reach_data)
SWOT_df  <- bind_rows(SWOT_raw)

plan(sequential)

# ============================================================
# Clean SWOT Data
# ============================================================

SWOT_clean <- SWOT_df %>%
  mutate(v17_reach_id = as.character(reach_id)) %>%
  left_join(SWORD_map %>% mutate(v17_reach_id = as.character(v17_reach_id)), by = "v17_reach_id") %>%
  filter(
    !is.na(wse),
    reach_q < 3,
    slope2 >= 0,
    partial_f == 0,
    ice_clim_f < 2,
    dark_frac < 0.4,
    xovr_cal_q <= 1,
    abs(wse - p_wse) < 10,
    between(abs(xtrk_dist), 10000, 60000),
    reach_q_b <= 4194304
  ) %>%
  group_by(v17_reach_id) %>%
  filter(n() >= 10) %>%
  mutate(outlier_flag = 1:length(slope2) %in% hampel(wse, k = 1, t0 = 3)$ind) %>%
  ungroup() %>%
  filter(outlier_flag == 0)

SWOT_slope <- SWOT_clean %>%
  summarize(
    median_slope = median(slope2),
    sd_slope = sd(slope2),
    .by = v17_reach_id
  )

# ============================================================
# Merge USGS & SWOT and compute hydraulics
# ============================================================

USGS_all <- USGS_gage_measurements %>%
  mutate(Agency_id = site_no) %>%
  left_join(USGS_df %>% mutate(v17_reach_id = as.character(v17_reach_id)), by = "Agency_id") %>%
  left_join(SWOT_slope, by = "v17_reach_id") %>%
  filter(
    !is.na(v17_reach_id),
    !is.na(Q_cms),
    measurement_dateTime >= as.POSIXct("2015-01-01"),
    measured_rating_diff %in% c("Excellent", "Good", "Fair"),
    meas_type != "ICE",
    !is.na(median_slope)
  ) %>%
  compute_manning() %>%
  mutate(
    Q_vtimesarea = velocity_ms * area_m2,
    Q_ratio = Q_vtimesarea / Q_cms,
    slope_check = between(median_slope, 0, 0.01),
    Q_check = between(Q_ratio, 0.90, 1.10),
    Fr_check = between(Fr, 0, 1),
    log10_Manning_n = log10(Manning_n)
  )

USGS_clean <- USGS_all %>%
  filter(
    slope_check,
    Q_check,
    Fr_check,
    between(Manning_n, 0.01, 0.15)
  ) %>%
  group_by(Agency_id) %>%
  filter(n() >= 10) %>%
  mutate(outlier_flag = 1:length(width_m) %in% hampel(width_m, k = 3, t0 = 3)$ind) %>%
  ungroup() %>%
  filter(outlier_flag == 0) %>%
  mutate(SWOT_median_slope = median_slope,
         SWOT_sd_slope = sd_slope,
         gage_id = Agency_id,
         gage_name = monitoring_location_name,
         Q_m3s = Q_cms) %>%
  select(gage_id, gage_name, longitude, latitude, v17_reach_id,
         measurement_dateTime, Q_m3s, width_m, depth_m, area_m2, velocity_ms,
         SWOT_median_slope, SWOT_sd_slope,
         Manning_n, Manning_n_slope_sd, Fr,measured_rating_diff, meas_type, streamflow_method)

# ============================================================
# Aggregated Manning & Regression
# ============================================================

SWORD_data <- read_sf(paste0(data_folder, "na_sword_reaches_v17b.gpkg")) %>%
  mutate(v17_reach_id = as.character(reach_id)) %>%
  filter(v17_reach_id %in% USGS_clean$v17_reach_id) %>%
  select(v17_reach_id, width, width_var, max_width, n_chan_max, slope, wse_var) %>%
  st_drop_geometry()

SWORD_node_data <- read_sf(paste0(data_folder, "na_sword_nodes_v17b.gpkg")) %>%
  mutate(v17_reach_id = as.character(reach_id)) %>%
  filter(v17_reach_id %in% USGS_clean$v17_reach_id) %>%
  st_drop_geometry() %>%
  group_by(v17_reach_id) %>%
  summarize(median_meand_len = median(meand_len),
            median_sinuosity = median(sinuosity))

SWORD_combined <- SWORD_data %>% left_join(SWORD_node_data) %>%
  mutate(SWORD_width_m = width,
         SWORD_width_var_m2 = width_var,
         SWORD_max_width_m = max_width,
         SWORD_n_chan_max = n_chan_max,
         SWORD_slope = slope,
         SWORD_meand_len = median_meand_len,
         SWORD_sinuosity = median_sinuosity,
         SWORD_wse_var_m2 = wse_var) %>%
  select(v17_reach_id, SWORD_width_m, SWORD_width_var_m2, SWORD_max_width_m, 
         SWORD_wse_var_m2, SWORD_slope, SWORD_n_chan_max, SWORD_meand_len, SWORD_sinuosity)

S3_data <- fread(paste0(data_folder, "w10s15-Training_Target_Data.csv")) %>% 
  mutate(v17_reach_id = as.character(reach_id),
         S3_sinuosity = 10^ areaS3_smth,
         S3_emd = emd,
         S3_cosT = cosT,
         S3_peakStride = linStride) %>%
  select(v17_reach_id, S3_sinuosity, S3_emd, S3_cosT, S3_peakStride) 

max_depth <- USGS_clean %>%
  summarize(max_depth = max(depth_m), .by = gage_id)

base_df <- USGS_clean %>%
  left_join(max_depth, by = "gage_id") %>%
  mutate(log10_Ym_Y = log10(max_depth / depth_m),
         log10_Manning_n = log10(Manning_n))

agg_df <- base_df %>%
  summarize(
    max_ADCP_width_m     = max(width_m),
    max_ADCP_depth_m     = max(depth_m),
    max_ADCP_velocity_ms = max(velocity_ms),
    max_ADCP_Q_m3s = max(Q_m3s),
    max_Q_Manning_n = Manning_n[which.max(Q_m3s)],
    ADCP_meas_count  = n(),
    .by = c(gage_id, gage_name, longitude, latitude, v17_reach_id, SWOT_median_slope)
  ) %>%
  filter(ADCP_meas_count >= 20) %>%
  left_join(SWORD_combined) %>%
  left_join(S3_data) %>%
  drop_na()

reg_df <- base_df %>%
  filter(
    is.finite(log10_Ym_Y),
    is.finite(log10_Manning_n)
  ) %>%
  nest(
    .by = c(
      gage_id,
      gage_name,
      v17_reach_id,
      SWOT_median_slope,
      SWOT_sd_slope
    )
  ) %>%
  mutate(
    OK = map_lgl(
      data,
      ~ nrow(.x) >= 2 && sd(.x$log10_Ym_Y) > 0
    )
  ) %>%
  mutate(
    model = map2(
      data, OK,
      function(d, ok) {
        if (!ok) return(NULL)
        lm(log10_Manning_n ~ log10_Ym_Y, data = d)
      }
    )
  ) %>%
  mutate(
    slope     = map_dbl(model, ~ coef(.x)[2], .default = NA_real_),
    intercept = map_dbl(model, ~ coef(.x)[1], .default = NA_real_),
    r2        = map_dbl(model, ~ summary(.x)$r.squared, .default = NA_real_)
  ) %>%
  transmute(
    gage_id,
    v17_reach_id,
    SWOT_median_slope,
    nb_param = 10^intercept,
    x_param  = slope,
    regression_r2 = r2
  )

width_depth_cor <- base_df %>%
  group_by(gage_id, v17_reach_id) %>%
  summarize(
    cor_ADCP_width_depth =
      if (n() >= 2 &&
          sd(width_m, na.rm = TRUE) > 0 &&
          sd(depth_m, na.rm = TRUE) > 0)
        cor(width_m, depth_m, use = "complete.obs", method = "spearman")
    else NA_real_,
    .groups = "drop"
  )

USGS_manning_clean <- agg_df %>%
  left_join(reg_df, by = c("gage_id", "v17_reach_id", "SWOT_median_slope")) %>%
  left_join(width_depth_cor, by = c("gage_id", "v17_reach_id")) %>%
  mutate(nb_ADCP_Qmax = max_Q_Manning_n,
         nb_ADCP_regress = nb_param,
         x_ADCP_regress = x_param) %>%
  mutate(bankfull_depth_ADCPwidth = 0.08 * max_ADCP_width_m ^ (0.39) / (SWOT_median_slope ^ (0.24))) %>%
  mutate(bankfull_depth_SWORDwidth = 0.08 *  SWORD_max_width_m ^ (0.39) / (SWOT_median_slope ^ (0.24))) %>%
  mutate(nb_case1_adcp_max_depth_width = 0.73 * (max_ADCP_depth_m ^ (2/3)) * (SWOT_median_slope ^ (0.19)) / (S3_peakStride ^ (0.32))) %>%
  mutate(nb_case2_adcp_max_width_depth_from_calc = 0.73 * (bankfull_depth_ADCPwidth ^ (2/3)) * (SWOT_median_slope ^ (0.19)) / (S3_peakStride ^ (0.32))) %>%
  mutate(nb_case3_sword_max_width_depth_from_calc = 0.73 * (bankfull_depth_SWORDwidth ^ (2/3)) * (SWOT_median_slope ^ (0.19)) / (S3_peakStride ^ (0.32))) %>%
  mutate(nb_multiple_regress = 10^(-0.7043 + 0.192 * log10(SWOT_median_slope) - 0.0077 * log10(max(0.01, SWORD_wse_var_m2)) - 0.058 * log10(SWORD_n_chan_max) - 0.0264 * log10(S3_peakStride))) %>%
  drop_na() %>%
  mutate(max_ADCP_width_m = round(max_ADCP_width_m, 1),
         max_ADCP_depth_m = round(max_ADCP_depth_m, 1),
         max_ADCP_velocity_ms = round(max_ADCP_velocity_ms, 2),
         max_ADCP_Q_m3s = round(max_ADCP_Q_m3s, 1),
         cor_ADCP_width_depth = round(cor_ADCP_width_depth, 2),
         SWOT_median_slope = round(SWOT_median_slope, 6),
         SWORD_width_m = round(SWORD_width_m, 0),
         SWORD_max_width_m = round(SWORD_max_width_m, 0),
         SWORD_width_var_m2 = round(SWORD_width_var_m2, 0),
         SWORD_wse_var_m2 = round(SWORD_wse_var_m2, 0),
         SWORD_slope = round(SWORD_slope, 3),
         SWORD_meand_len = round(SWORD_meand_len, 0),
         SWORD_sinuosity = round(SWORD_sinuosity, 2),
         S3_peakStride = round(S3_peakStride, 0),
         S3_sinuosity = round(S3_sinuosity, 2),
         S3_emd = round(S3_emd, 2),
         S3_cosT = round(S3_cosT, 2),
         nb_ADCP_Qmax = round(nb_ADCP_Qmax, 3),
         nb_ADCP_regress = round(nb_ADCP_regress, 3),
         x_ADCP_regress = round(x_ADCP_regress, 2),
         regression_r2 = round(regression_r2, 2),
         nb_multiple_regress = round(nb_multiple_regress, 3),
         bankfull_depth_ADCPwidth = round(bankfull_depth_ADCPwidth, 1),
         bankfull_depth_SWORDwidth = round(bankfull_depth_SWORDwidth, 1),
         nb_case1_adcp_max_depth_width = round(nb_case1_adcp_max_depth_width, 3),
         nb_case2_adcp_max_width_depth_from_calc = round(nb_case2_adcp_max_width_depth_from_calc, 3),
         nb_case3_sword_max_width_depth_from_calc = round(nb_case3_sword_max_width_depth_from_calc, 3)) %>%
  select(gage_id, gage_name, longitude, latitude, v17_reach_id, 
         max_ADCP_width_m, max_ADCP_depth_m, max_ADCP_velocity_ms, max_ADCP_Q_m3s, ADCP_meas_count, cor_ADCP_width_depth,
         SWOT_median_slope,
         SWORD_width_m, SWORD_max_width_m, SWORD_width_var_m2, SWORD_wse_var_m2, SWORD_slope, 
         SWORD_n_chan_max, SWORD_meand_len, SWORD_sinuosity,
         S3_peakStride, S3_sinuosity, S3_emd, S3_cosT,
         nb_ADCP_Qmax, nb_ADCP_regress, x_ADCP_regress, regression_r2, nb_multiple_regress, 
         bankfull_depth_ADCPwidth, bankfull_depth_SWORDwidth,
         nb_case1_adcp_max_depth_width, nb_case2_adcp_max_width_depth_from_calc, nb_case3_sword_max_width_depth_from_calc) %>%
  distinct() %>%
  filter(max_ADCP_width_m >= 50) %>%
  filter(!v17_reach_id %in% c("73286000021", "73258000161", "73180300011"))

fwrite(USGS_manning_clean, file.path(data_folder, "Manning_resistance.csv"))

USGS_compiled <- USGS_clean %>% filter(gage_id %in% USGS_manning_clean$gage_id) %>%
  mutate(Q_m3s = round(Q_m3s, 1),
         width_m = round(width_m, 1),
         depth_m = round(depth_m, 1),
         area_m2 = round(area_m2, 1), 
         velocity_ms = round(velocity_ms, 2),
         SWOT_median_slope = round(SWOT_median_slope, 6),
         SWOT_sd_slope = round(SWOT_sd_slope, 6),
         Manning_n = round(Manning_n, 3),
         Manning_n_slope_sd = round(Manning_n_slope_sd, 3),
         Fr = round(Fr, 2)) %>%
  select(gage_id, gage_name, longitude, latitude, v17_reach_id,
         measurement_dateTime, Q_m3s, width_m, depth_m, area_m2,
         velocity_ms, SWOT_median_slope, SWOT_sd_slope, Manning_n,
         Manning_n_slope_sd, Fr, measured_rating_diff, meas_type,
         streamflow_method)

fwrite(USGS_compiled, file.path(data_folder, "ADCP_compiled.csv"))

