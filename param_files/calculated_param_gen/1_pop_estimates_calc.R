#Estimates populations every year by gender in KZN
#Uses GBD

#for 15 - 59 (all ages in model) used for mortality and hiv incidnece estimate
#15-19 (for calculating aging into model, births perc)

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to resource allocation HIV TB SA
indir_outdir <- paste0(here(), '/param_files/calculated_param_gen/raw_input_data/GBD/')

setwd(indir_outdir)

df<-read.csv('all_cause_mort_num_rate.csv')
regions_list<-unique(df$location)

pop_df_15_59<-data.frame()
pop_df_15_19<-data.frame()

for (i in regions_list){
  #get pop estimates by taking 1/rate * num of deaths (for all cause)
  pop_df_15_59_temp<-df%>%
    filter(location == i)%>%
    select(-c('location', 'cause', 'upper', 'lower', 'measure'))
  
  pop_df_15_59_temp<-dcast(pop_df_15_59_temp, sex+age+year ~ metric, mean, 
                      value.var = 'val') #note base population estimates are the same
  
  ##all ages in model
  pop_df_15_59_temp<-pop_df_15_59_temp%>%
    mutate(Rate = Rate/100000)%>% #is per 100,000
    mutate(pop_val = Number/Rate)%>%
    group_by(year, sex)%>%
    summarise(expected_total_pop = sum(pop_val))%>%
    mutate(region_name = i)
  
  if(nrow(pop_df_15_59) == 0){
    pop_df_15_59<-pop_df_15_59_temp
  } else {
    pop_df_15_59<-rbind(pop_df_15_59, pop_df_15_59_temp)
  }
  
  pop_df_15_19_temp<-df%>%
    filter(location == i)%>%
    select(-c('location', 'cause', 'upper', 'lower', 'measure'))%>%
    filter(age == "15-19 years")
  
  pop_df_15_19_temp<-dcast(pop_df_15_19_temp, sex+age+year ~ metric, mean, 
                value.var = 'val') #note base population estimates are the same
  
  pop_df_15_19_temp<-pop_df_15_19_temp%>%
    mutate(pop_val = (Number/Rate))%>%
    group_by(year, sex)%>%
    summarise(expected_total_pop = sum(pop_val)*100000)%>%
    mutate(region_name = i)
  
  if(nrow(pop_df_15_19) == 0){
    pop_df_15_19<-pop_df_15_19_temp
  } else {
    pop_df_15_19<-rbind(pop_df_15_19, pop_df_15_19_temp)
  }
  
  
}

write.csv(pop_df_15_59, 'pop_estimates_15_59.csv', row.names = FALSE)
write.csv(pop_df_15_19, 'pop_estimates_15_19.csv', row.names = FALSE)

