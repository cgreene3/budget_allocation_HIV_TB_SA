#Estimating active TB mortalities/TB incidence/HIV prev per 100K to calibrate to
#Uses GBD
#regional population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to KZN_south_africa
indir <- paste0(here(),'/param_files/calculated_param_gen/raw_input_data/GBD')
outdir <- paste0(here(),'/param_files/target_calibration_estimates')

##########read in pop files#################
#data pull date 02/21/2022
calib_years <-c(2005, 2017)

setwd(indir)

#get popestimates
pop_df<-read.csv('pop_estimates_15_59.csv')%>%
  group_by(region_name, year)%>%
  summarise(expected_total_pop = sum(expected_total_pop))

#input TB and HIV mortality rate calibration targets
TBHIV_mort_num_df<-read.csv('disease_mort_num.csv')
SA_TBHIV_mort_num_df<-read.csv('SA_TBHIVmort_num.csv')

#TB incidence rate calibration targets
TB_inc_num_df<-read.csv('tb_incidence_num_df.csv')

#HIV prevalence calibration targets
HIV_prev_num_df<-read.csv('hiv_prev_num.csv')

#TB mort rate calibration estimates 2005 and 2017
region_rate_TBHIV_mort_calibration_df<-TBHIV_mort_num_df%>%
  rename(region_name = location)%>%
  group_by(year, region_name)%>%
  summarise(val_mort_num = sum(val),
            upper_mort_num = sum(upper),
            lower_mort_num = sum(lower))%>%
  left_join(pop_df, by = c('year', 'region_name'))%>%
  mutate(val_rate = (val_mort_num/expected_total_pop)*100000,
         upper_rate = (upper_mort_num/expected_total_pop)*100000,
         lower_rate = (lower_mort_num/expected_total_pop)*100000)%>%
  select(c('region_name', 'year',  
           'val_rate', 'upper_rate', 'lower_rate'))

#TB mort rate calibration estimates 2005 and 2017
allregions_num_TBHIV_mort_calibration<-TBHIV_mort_num_df%>%
  group_by(year)%>%
  summarise(val_mort_num = sum(val),
            upper_mort_num = sum(upper),
            lower_mort_num = sum(lower))

#TB incidence calibration estimates 
region_rate_TB_inc_calibration_df<-TB_inc_num_df%>%
  rename(region_name = location)%>%
  group_by(year, region_name)%>%
  summarise(val_inc_num = sum(val),
            upper_inc_num = sum(upper),
            lower_inc_num = sum(lower))%>%
  left_join(pop_df, by = c('year', 'region_name'))%>%
  mutate(val_rate = (val_inc_num/expected_total_pop)*100000,
         upper_rate = (upper_inc_num/expected_total_pop)*100000,
         lower_rate = (lower_inc_num/expected_total_pop)*100000)%>%
  select(c('region_name', 'year',  
           'val_rate', 'upper_rate', 'lower_rate'))%>%
  mutate(lower_rate = if_else(region_name == 'Western Cape', lower_rate*.9, lower_rate))%>%
  mutate(upper_rate = if_else(region_name == 'Western Cape', upper_rate*1.1, upper_rate))
  

#HIV prev calibration estimates 2005 and 2017
region_rate_HIV_prev_calibration_df<-HIV_prev_num_df%>%
  rename(region_name = location)%>%
  group_by(year, region_name)%>%
  summarise(val_prev_num = sum(val),
            upper_prev_num = sum(upper),
            lower_prev_num = sum(lower))%>%
  left_join(pop_df, by = c('year', 'region_name'))%>%
  mutate(val_rate = (val_prev_num/expected_total_pop)*100000,
         upper_rate = (upper_prev_num/expected_total_pop)*100000,
         lower_rate = (lower_prev_num/expected_total_pop)*100000)%>%
  select(c('region_name', 'year',  
           'val_rate', 'upper_rate', 'lower_rate'))

setwd(outdir)
write.csv(region_rate_TBHIV_mort_calibration_df, 'SA_regional_GBD_TBHIV_mort_rate_calibration_df.csv', row.names = FALSE)
write.csv(allregions_num_TBHIV_mort_calibration, 'SA_all_GBD_TBHIV_mort_num_calibration_df.csv', row.names = FALSE)
write.csv(region_rate_TB_inc_calibration_df, 'SA_regional_GBD_TB_inc_rate_calibration_df.csv', row.names = FALSE)
write.csv(region_rate_HIV_prev_calibration_df, 'SA_regional_GBD_HIV_prev_rate_calibration_df.csv', row.names = FALSE)

