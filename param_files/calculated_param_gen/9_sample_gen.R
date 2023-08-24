#clean workspace
rm(list = ls())
gc()


library(readxl)
library(here)
library(dplyr)
library(ExtDist)
library(lhs)
library(ggplot2)
library(stringr)
library(reshape2)

n_general_samples <- 500
n_regional_samples <- 40

set.seed(as.integer(1)) 

#update parameter file
outdir <- paste0(here(),'/param_files/input_parameters')
indir_params <- paste0(here(),'/param_files/calculated_param_gen/raw_input_data')

setwd(indir_params)
model_params_df<-read_excel('all_regions_SA_model_parameters.xlsx')%>%
  select(c('model_matched_param', 'min', 'max', 'value', 'calibrated', 'regional_specific'))

model_params_df$model_matched_param<-str_replace_all(model_params_df$model_matched_param, c("," = "."))

n_calib_params<-sum(model_params_df$calibrated == 'yes')
n_regional_calib_params<-sum(model_params_df$regional_specific == 'yes') #all regional specific are calibrated


###GENERATE GENERAL PARAMETER SETS####
calib_param_general_df<-model_params_df%>%
  filter(calibrated == 'yes',
         regional_specific == 'no')

calib_param_general_names<-calib_param_general_df$model_matched_param

raw_lhs_df_general<-as.data.frame(randomLHS(n_general_samples, length(calib_param_general_names)))
colnames(raw_lhs_df_general)<-calib_param_general_names

sample_df_general <- data.frame(matrix(nrow = n_general_samples, ncol = 0))

for (mp in calib_param_general_names){
  
  model_params_df_temp<-model_params_df%>%
    filter(model_matched_param == mp)
  
  max_temp<-as.numeric(model_params_df_temp[c('max')])
  
  min_temp<-as.numeric(model_params_df_temp[c('min')])
    
  sample_df_general_temp<-as.data.frame(qunif(unlist(raw_lhs_df_general%>%select(mp)),
                                      min_temp, max_temp))
  colnames(sample_df_general_temp)<-mp
  sample_df_general<-cbind(sample_df_general, sample_df_general_temp)
} 

sample_df_general$general_sim_id <-1:nrow(sample_df_general)

###GENERATE REGIONAL PARAMETER SETS####
calib_param_regional_df<-model_params_df%>%
  filter(calibrated == 'yes',
         regional_specific == 'yes')

calib_param_regional_names<-calib_param_regional_df$model_matched_param

raw_lhs_df_regional<-as.data.frame(randomLHS(n_regional_samples, length(calib_param_regional_names)))
colnames(raw_lhs_df_regional)<-calib_param_regional_names

sample_df_regional <- data.frame(matrix(nrow = n_regional_samples, ncol = 0))

for (mp in calib_param_regional_names){
  
  model_params_df_temp_regional<-model_params_df%>%
    filter(model_matched_param == mp)
  
  max_temp<-as.numeric(model_params_df_temp_regional[c('max')])
  
  min_temp<-as.numeric(model_params_df_temp_regional[c('min')])
  
  sample_df_regional_temp<-as.data.frame(qunif(unlist(raw_lhs_df_regional%>%select(mp)),
                                              min_temp, max_temp))
  colnames(sample_df_regional_temp)<-mp
  sample_df_regional<-cbind(sample_df_regional, sample_df_regional_temp)
} 

sample_df_regional$regional_sim_id <-1:nrow(sample_df_regional)

sample_df<-merge(sample_df_general, 
                 sample_df_regional, 
                 by = character())

#####add in non-calibrated parameters#####

non_calib_param_df<-model_params_df%>%
  filter(calibrated == 'no')

non_calib_param_names<-non_calib_param_df$model_matched_param

for (mp in non_calib_param_names){
  value_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                         c('value')])
  sample_df_temp<-as.data.frame(rep(value_temp, times = nrow(sample_df)))
  colnames(sample_df_temp)<-mp
  sample_df<-cbind(sample_df, sample_df_temp)
  
}

sample_df$sim_id<-seq(1:nrow(sample_df))

##move id to the front
colnames_temp<-colnames(sample_df%>%select(-c('sim_id', 'general_sim_id', 'regional_sim_id')))

sample_df<-sample_df%>%
  select(c('sim_id', 'general_sim_id', 'regional_sim_id',
           colnames_temp))


setwd(outdir)
write.csv(sample_df, 'calibration_sets_df.csv', row.names = FALSE)
