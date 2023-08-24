#last updated aug 23, 2023
#run after program runs (get metrics overtime)

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle'), require, character.only=T)

##GITHUB/LOCAL###
#location where we sent outputs
library(here)
indir_outputs<-paste0(here(), '/results/program_runs/metrics/')
indir_target_calibration_estimates<-paste0(here(), '/param_files/target_calibration_estimates')
indir_input_params<-paste0(here(), '/param_files/input_parameters')
# 
# #location where we send calib analysis results
outdir<-paste0(here(), '/results/calibration_analysis/graphs')

setwd(indir_target_calibration_estimates)
HIV_prev_df<-read.csv('SA_regional_GBD_HIV_prev_rate_calibration_df.csv')
TB_inc_df<-read.csv('SA_regional_GBD_TB_inc_rate_calibration_df.csv')
TBHIV_mort_df<-read.csv('SA_regional_GBD_TBHIV_mort_rate_calibration_df.csv')

setwd(indir_input_params)
accepted_param_sets_ref_df<-read.csv('accepted_param_sets_ref.csv')
region_id_ref_df<-read.csv('region_id_ref.csv')
  

for (i in region_id_ref_df$region_id){
  
  region_name_temp<-region_id_ref_df$region_name[i]
  
  setwd(paste0(indir_outputs, region_name_temp))
  all_files_in_outdir <- list.files(pattern="summarised_eval_metrics_df*")
  outputs_combined_df<-read.csv(all_files_in_outdir[1])%>%
    filter(year < 2018)%>% #before program eval
    select(-c('program_id'))
  
  for(i in 2:length(all_files_in_outdir)){
    temp<-read.csv(all_files_in_outdir[i])%>%
      filter(year < 2018)%>% #before program eval
      select(-c('program_id'))
    outputs_combined_df<-rbind(outputs_combined_df, temp)
  }
  
  outputs_combined_df_reshape<-melt(outputs_combined_df, id = c("year", "general_sim_id", "regional_sim_id"))
  
  calib_TBHIV_mort_rate<-outputs_combined_df_reshape%>%
    filter(grepl('mort', variable) == TRUE)%>%
    group_by(general_sim_id, regional_sim_id, year)%>%
    summarise(value = sum(value))%>%
    mutate(variable = 'TBHIV_mort_per_Y_100k_ppl')
  
  calib_TBinc_rate<-outputs_combined_df_reshape%>%
    filter(grepl('TB_inc', variable) == TRUE)%>%
    group_by(variable, general_sim_id, regional_sim_id, year)%>%
    summarise(value = sum(value))
  
  
  calib_HIVprev_rate<-outputs_combined_df_reshape%>%
    filter(grepl('HIV_prev', variable) == TRUE)%>%
    group_by(variable, general_sim_id, regional_sim_id, year)%>%
    summarise(value = sum(value))
  
  
  
}





outputs_combined_df<-outputs_combined_df%>%
  filter(year <= 2017)



calib_TBHIV_mort_rate<-outputs_combined_df_reshape%>%
  filter(grepl('mort', variable) == TRUE)%>%
  group_by(sim_id, year)%>%
  summarise(value = sum(value))%>%
  left_join(TBHIV_mort_df, by = c('year'))%>%
  mutate(in_confidence_interval = if_else(((value < upper_rate) & (value > lower_rate)),
                                          1, 0))%>%
  filter(year == 2005 | year == 2017)%>%
  mutate(variable = 'TBHIV_mort_rate')

calib_TBinc_rate<-outputs_combined_df_reshape%>%
  filter(grepl('inc', variable) == TRUE)%>%
  group_by(sim_id, year)%>%
  summarise(value = sum(value))%>%
  left_join(TB_inc_df, by = c('year'))%>%
  mutate(in_confidence_interval = if_else(((value < upper_rate) & (value > lower_rate)),
                                          1, 0))%>%
  filter(year == 2005 | year == 2017)%>%
  mutate(variable = "TB_inc_rate")


calib_HIVprev_rate<-outputs_combined_df_reshape%>%
  filter(grepl('hiv_prev', variable) == TRUE)%>%
  group_by(sim_id, year)%>%
  summarise(value = sum(value))%>%
  left_join(HIV_prev_df, by = c('year'))%>%
  mutate(in_confidence_interval = if_else(((value < upper_rate) & (value > lower_rate)),
                                          1, 0))%>%
  filter(year == 2005 | year == 2017)%>%
  mutate(variable = "HIV_prev_rate")


##combine and analyze
results_df<-rbind(calib_TBHIV_mort_rate, calib_TBinc_rate, calib_HIVprev_rate)

results_df_summarised_selection_options<-results_df%>%
  group_by(sim_id)%>%
  summarise(total_in_confidence = sum(in_confidence_interval)
            )%>%
  left_join(sim_calibration_ref_df, by = c('sim_id'))

#if you want to look at the distribution of each parameter for sets that hit 80% of targets
almost_best_calibration_sets_ref_df<-results_df_summarised_selection_options%>%
  filter(total_in_confidence >= 5)

#accepted parameter sets
accepted_calibration_sets_ref_df<-results_df_summarised_selection_options%>%
  filter(total_in_confidence == 6)

#to summarise metrics of accepted points
accepted_calibration_metrics_df<-results_df%>%
  filter(sim_id %in% accepted_calibration_sets_ref_df$sim_id)

#to evaluate metrics that are missed by sets that hit at least 90% of metrics
almost_calibration_metrics_df<-results_df%>%
  filter(sim_id %in% almost_best_calibration_sets_ref_df$sim_id)

setwd(outdir)
write.csv(accepted_calibration_sets_ref_df, 'accepted_calibration_sets_ref_df.csv', row.names = FALSE)
write.csv(almost_best_calibration_sets_ref_df, 'almost_best_calibration_sets_ref_df.csv', row.names = FALSE)
write.csv(almost_calibration_metrics_df, 'almost_calibration_metrics_df.csv', row.names = FALSE)
write.csv(accepted_calibration_metrics_df, 'accepted_calibration_metrics_df.csv', row.names = FALSE)
#write.csv(missing_sims_df, 'missing_sims_df.csv', row.names = FALSE)
#setwd(paste0(here(), '/test/calibration_results/'))
setwd(paste0(indir_input_params, '/general_param_sets_itr'))
#print remaining general parameter sets
general_accepted_calibration_sets_ref<-data.frame(general_parameter_id = unique(accepted_calibration_sets_ref_df$general_sim_id))

write.csv(general_accepted_calibration_sets_ref, 
          paste0('general_accepted_calibration_sets_ref_df',
                 itr_id, '.csv'), 
          row.names = FALSE)