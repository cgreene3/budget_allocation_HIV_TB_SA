#last updated aug 8, 2023
#ids general parameter sets that are accepted in all regions
#finds max min of parameter sets
#finds all general/region combinations
#ids combined parameter sets that fall within 95% confidence interval of 
#num SA HIV-TB mort

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here'), require, character.only=T)

#read in calib analysis results
indir_calib_analysis<-paste0(here(), '/results/calibration_analysis')
setwd(indir_calib_analysis)
region_id_ref_df<-read.csv('region_id_ref.csv')
region_ids_array<-region_id_ref_df$region_id
region_names_array<-array(region_id_ref_df$region_name)

#read in pop df
setwd(paste0(here(), '/param_files/calculated_param_gen/raw_input_data/GBD/'))
pop_df<-read.csv('pop_estimates_15_59.csv')%>%
  filter(year %in% c(2005, 2017))%>%
  group_by(region_name, year)%>%
  summarise(expected_total_pop = sum(expected_total_pop))%>%
  left_join(region_id_ref_df, by = c('region_name'))

#disease_mort_num<-read.csv('disease_mort_num.csv')%>%
#  filter(year == 2017)%>%
#  summarise(total_disease_mort_max= sum(upper),
#            total_disease_mort_min = sum(lower))

total_disease_mort_max_2017 =  214370.07
total_disease_mort_min_2017 =  126005.94
total_disease_mort_val_2017 = 161751.0834

total_disease_mort_max_2005 =  373674.6207
total_disease_mort_min_2005 =  201329.1864
total_disease_mort_val_2005 = 280659.4793

for(i_temp in region_ids_array){
  if(i_temp == 1){
    setwd(paste0(indir_calib_analysis, '/', region_names_array[i_temp]))
    outputs_combined_df_all<-read.csv('accepted_calibration_sets_ref_df.csv')%>%
      select(c('sim_id', 'general_sim_id', 'regional_sim_id'))
    outputs_combined_df_all$region_id=rep(i_temp, times = nrow(outputs_combined_df_all))
    
    #get metrics
    metrics_combined_df_all<-read.csv('accepted_calibration_metrics_df.csv')%>%
      filter(#variable == 'TBHIV_mort_rate',
             year %in% c(2005, 2017))#%>%
      #select(c('sim_id', 'year', 'variable', 'value'))
    metrics_combined_df_all$region_id=rep(i_temp, times = nrow(outputs_combined_df_all))
  } #else if (i_temp == 9){ #will expand acceptance of TB inc rate.
  #   setwd(paste0(indir_calib_analysis, '/', region_names_array[i_temp]))
  #   outputs_combined_df_temp<-read.csv('almost_best_calibration_sets_ref_df.csv')%>%
  #     select(c('sim_id', 'general_sim_id', 'regional_sim_id'))
  #   outputs_combined_df_temp$region_id=rep(i_temp, times = nrow(outputs_combined_df_temp))
  #   outputs_combined_df_all<-rbind(outputs_combined_df_all, outputs_combined_df_temp)
  #   
  #   #get metrics
  #   metrics_combined_df_temp<-read.csv('almost_calibration_metrics_df.csv')%>%
  #     filter(variable == 'TBHIV_mort_rate',
  #            year %in% c(2005, 2017))%>%
  #     select(c('sim_id', 'year', 'value'))
  #   metrics_combined_df_temp$region_id=rep(i_temp, times = nrow(metrics_combined_df_temp))
  #   metrics_combined_df_all<-rbind(metrics_combined_df_all, metrics_combined_df_temp)
  # } 
  else {
    setwd(paste0(indir_calib_analysis, '/', region_names_array[i_temp]))
    outputs_combined_df_temp<-read.csv('accepted_calibration_sets_ref_df.csv')%>%
      select(c('sim_id', 'general_sim_id', 'regional_sim_id'))
    outputs_combined_df_temp$region_id=rep(i_temp, times = nrow(outputs_combined_df_temp))
    outputs_combined_df_all<-rbind(outputs_combined_df_all, outputs_combined_df_temp)
    
    #get metrics
    metrics_combined_df_temp<-read.csv('accepted_calibration_metrics_df.csv')%>%
      filter(#variable == 'TBHIV_mort_rate',
             year %in% c(2005, 2017))#%>%
      #select(c('sim_id', 'year', 'variable', 'value'))
    metrics_combined_df_temp$region_id=rep(i_temp, times = nrow(metrics_combined_df_temp))
    metrics_combined_df_all<-rbind(metrics_combined_df_all, metrics_combined_df_temp)
  }
}

#number of sets accepted in each region
accepted_sets_count<-outputs_combined_df_all%>%
  group_by(region_id)%>%
  summarise(n_accepted = n())

#get general sim ids accepted at least once in each region
general_sim_ids_accepted_analysis_df<-outputs_combined_df_all%>%
  group_by(region_id, general_sim_id)%>%
  summarise(n_regional_sim_ids_accepted = n())%>%
  group_by(general_sim_id)%>%
  summarise(n_regions_accepted_general_sim_id = n())%>% 
  filter(n_regions_accepted_general_sim_id == length(region_ids_array))

accepted_sets_count_after_general_accepted_filtered<-outputs_combined_df_all%>%
  filter(general_sim_id %in% general_sim_ids_accepted_analysis_df$general_sim_id)%>%
  group_by(region_id)%>%
  summarise(n_accepted = n())

#unique combinations that represent max and min 
#(3 metrics x 2 max min = 6 per region
#x 9 regions x24 general parameter sets should be 2,592 rows)
max_min_param_sets_filtered<-outputs_combined_df_all%>%
  filter(general_sim_id %in% general_sim_ids_accepted_analysis_df$general_sim_id)%>%
  left_join(metrics_combined_df_all, by = c('sim_id', 'region_id'))%>%
  mutate(max_diff = upper_rate - value,
         min_diff = value - lower_rate)%>%
  group_by(general_sim_id, region_id, year, variable)%>%
  mutate(rank_max_diff = rank(-max_diff),
         rank_min_diff = rank(-min_diff))%>%
  filter(rank_max_diff == 1 | rank_min_diff == 1)%>%
  filter(year == 2017)

max_min_param_sets_filtered_distinct<-max_min_param_sets_filtered%>%
  ungroup()%>%
  distinct(general_sim_id, regional_sim_id, region_id)

setwd(here('param_files/input_parameters'))
write.csv(max_min_param_sets_filtered_distinct, 'accepted_param_sets.csv', row.names = FALSE)

number_in_each_region<-max_min_param_sets_filtered_distinct%>%
  group_by(region_id)%>%
  summarise(n_distinct_in_region = n())

#to look up total mort
total_mort_df<-outputs_combined_df_all%>%
  filter(general_sim_id %in% general_sim_ids_accepted_analysis_df$general_sim_id)%>%
  left_join(metrics_combined_df_all, by = c('sim_id', 'region_id'))%>%
  filter(variable == 'TBHIV_mort_rate')%>%
  left_join(pop_df, by = c('region_id', 'region_name', 'year'))%>%
  mutate(total_disease_mort_val = value*(expected_total_pop/100000))

##calculate number of combinations
total_combin_all <-rep(0, times = length(unique(general_sim_ids_accepted_analysis_df$general_sim_id)))

#find all combinations
#2,513,680
for (gen_id_ref in 1:length(general_sim_ids_accepted_analysis_df$general_sim_id)){
  lists_of_regional_id <- vector(mode = "list", length = length(region_ids_array))
  
  gen_id_temp<-general_sim_ids_accepted_analysis_df$general_sim_id[gen_id_ref]
  for(i_temp in region_ids_array){
    #look up sims associated with general sim id and each region
    sims_temp_df<-max_min_param_sets_filtered%>%
      filter(general_sim_id == gen_id_temp,
             region_id == i_temp)
    lists_of_regional_id[[i_temp]]<-unique(sims_temp_df$regional_sim_id)
  }
  total_combin_temp<-1
  for(i in region_ids_array){
      total_combin_temp<-length(lists_of_regional_id[[i]])*total_combin_temp
      }
  total_combin_all[gen_id_ref]<-total_combin_temp
  }

last_time<-Sys.time()

# max_min_combin_df<-outputs_combined_df_all_filtered%>%
#   mutate(diff_max = )
#   group_by(general_sim_id, region_id, year, variable)%>%
#   mutate(max_value = max(value),
#          min_value =min(value))%>%
#   filter(value == min_value | value == max_value)

parameter_set_itr<-1034313#1 

#find all combinations
for (gen_id_temp in general_sim_ids_accepted_analysis_df$general_sim_id[11:length(general_sim_ids_accepted_analysis_df$general_sim_id)]){
  print(gen_id_temp)
  
  lists_of_regional_id <- vector(mode = "list", length = length(region_ids_array))
  
  for(i_temp in region_ids_array){
    #look up sims associated with general sim id and each region
    sims_temp_df<-max_min_param_sets_filtered%>%
      filter(general_sim_id == gen_id_temp,
             region_id == i_temp)
    lists_of_regional_id[[i_temp]]<-unique(sims_temp_df$regional_sim_id)
  }
  #reinitialize param_sets_df
  parameter_sets_df<-data.frame(matrix(ncol = 1 + 1 + length(region_ids_array)+2, nrow = 0))
  colnames(parameter_sets_df)<-c('sim_id', 'general_sim_id',
                                 paste("region_id_", region_ids_array, "_sim_id", sep = ''),
                                 'total_disease_mort_SA_2005',
                                 'total_disease_mort_SA_2017')
  #find all combinations
  lapply(lists_of_regional_id[[1]], function(r1){
    

  #for(r1 in lists_of_regional_id[[1]]){
    total_disease_mort_val_r1 <-total_mort_df%>%
      filter(region_id == 1,
             general_sim_id == gen_id_temp,
             regional_sim_id == r1)
    
    total_disease_mort_val_r1_2005<-total_disease_mort_val_r1%>%
      filter(year == 2005)
    total_disease_mort_val_r1_2005<-total_disease_mort_val_r1_2005$total_disease_mort_val
    
    total_disease_mort_val_r1_2017<-total_disease_mort_val_r1%>%
      filter(year == 2017)
    total_disease_mort_val_r1_2017<-total_disease_mort_val_r1_2017$total_disease_mort_val
    
    
    lapply(lists_of_regional_id[[2]], function(r2){
      
    
    #for(r2 in lists_of_regional_id[[2]]){
      total_disease_mort_val_r2 <-total_mort_df%>%
        filter(region_id == 2,
               general_sim_id == gen_id_temp,
               regional_sim_id == r2)
      
      total_disease_mort_val_r2_2005<-total_disease_mort_val_r2%>%
        filter(year == 2005)
      total_disease_mort_val_r2_2005<-total_disease_mort_val_r2_2005$total_disease_mort_val
      
      total_disease_mort_val_r2_2017<-total_disease_mort_val_r2%>%
        filter(year == 2017)
      total_disease_mort_val_r2_2017<-total_disease_mort_val_r2_2017$total_disease_mort_val
      
      
      lapply(lists_of_regional_id[[3]], function(r3){
        
      
      #for(r3 in lists_of_regional_id[[3]]){
        total_disease_mort_val_r3 <-total_mort_df%>%
          filter(region_id == 3,
                 general_sim_id == gen_id_temp,
                 regional_sim_id == r3)
        
        total_disease_mort_val_r3_2005<-total_disease_mort_val_r3%>%
          filter(year == 2005)
        total_disease_mort_val_r3_2005<-total_disease_mort_val_r3_2005$total_disease_mort_val
        
        total_disease_mort_val_r3_2017<-total_disease_mort_val_r3%>%
          filter(year == 2017)
        total_disease_mort_val_r3_2017<-total_disease_mort_val_r3_2017$total_disease_mort_val
        
        
        
        lapply(lists_of_regional_id[[4]], function(r4){
          
        
        #for(r4 in lists_of_regional_id[[4]]){
          total_disease_mort_val_r4 <-total_mort_df%>%
            filter(region_id == 4,
                   general_sim_id == gen_id_temp,
                   regional_sim_id == r4)
          
          total_disease_mort_val_r4_2005<-total_disease_mort_val_r4%>%
            filter(year == 2005)
          total_disease_mort_val_r4_2005<-total_disease_mort_val_r4_2005$total_disease_mort_val
          
          total_disease_mort_val_r4_2017<-total_disease_mort_val_r4%>%
            filter(year == 2017)
          total_disease_mort_val_r4_2017<-total_disease_mort_val_r4_2017$total_disease_mort_val
          
          lapply(lists_of_regional_id[[5]], function(r5){
            
          
          #for(r5 in lists_of_regional_id[[5]]){
            total_disease_mort_val_r5 <-total_mort_df%>%
              filter(region_id == 5,
                     general_sim_id == gen_id_temp,
                     regional_sim_id == r5)
            
            total_disease_mort_val_r5_2005<-total_disease_mort_val_r5%>%
              filter(year == 2005)
            total_disease_mort_val_r5_2005<-total_disease_mort_val_r5_2005$total_disease_mort_val
            
            total_disease_mort_val_r5_2017<-total_disease_mort_val_r5%>%
              filter(year == 2017)
            total_disease_mort_val_r5_2017<-total_disease_mort_val_r5_2017$total_disease_mort_val
            
            lapply(lists_of_regional_id[[6]], function(r6){
              
            
            #for(r6 in lists_of_regional_id[[6]]){
              total_disease_mort_val_r6 <-total_mort_df%>%
                filter(region_id == 6,
                       general_sim_id == gen_id_temp,
                       regional_sim_id == r6)
              
              total_disease_mort_val_r6_2005<-total_disease_mort_val_r6%>%
                filter(year == 2005)
              total_disease_mort_val_r6_2005<-total_disease_mort_val_r6_2005$total_disease_mort_val
              
              total_disease_mort_val_r6_2017<-total_disease_mort_val_r6%>%
                filter(year == 2017)
              total_disease_mort_val_r6_2017<-total_disease_mort_val_r6_2017$total_disease_mort_val
              
              lapply(lists_of_regional_id[[7]], function(r7){
                
              
              #for(r7 in lists_of_regional_id[[7]]){
                total_disease_mort_val_r7 <-total_mort_df%>%
                  filter(region_id == 7,
                         general_sim_id == gen_id_temp,
                         regional_sim_id == r7)
                
                total_disease_mort_val_r7_2005<-total_disease_mort_val_r7%>%
                  filter(year == 2005)
                total_disease_mort_val_r7_2005<-total_disease_mort_val_r7_2005$total_disease_mort_val
                
                total_disease_mort_val_r7_2017<-total_disease_mort_val_r7%>%
                  filter(year == 2017)
                total_disease_mort_val_r7_2017<-total_disease_mort_val_r7_2017$total_disease_mort_val
                
                
                lapply(lists_of_regional_id[[8]], function(r8){
                  
                
                #for(r8 in lists_of_regional_id[[8]]){
                  total_disease_mort_val_r8 <-total_mort_df%>%
                    filter(region_id == 8,
                           general_sim_id == gen_id_temp,
                           regional_sim_id == r8)
                  
                  total_disease_mort_val_r8_2005<-total_disease_mort_val_r8%>%
                    filter(year == 2005)
                  total_disease_mort_val_r8_2005<-total_disease_mort_val_r8_2005$total_disease_mort_val
                  
                  total_disease_mort_val_r8_2017<-total_disease_mort_val_r8%>%
                    filter(year == 2017)
                  total_disease_mort_val_r8_2017<-total_disease_mort_val_r8_2017$total_disease_mort_val
                  
                  lapply(lists_of_regional_id[[9]], function(r9){
                    
                  
                  #for(r9 in lists_of_regional_id[[9]]){
                    total_disease_mort_val_r9 <-total_mort_df%>%
                      filter(region_id == 9,
                             general_sim_id == gen_id_temp,
                             regional_sim_id == r9)
                    
                    total_disease_mort_val_r9_2005<-total_disease_mort_val_r9%>%
                      filter(year == 2005)
                    total_disease_mort_val_r9_2005<-total_disease_mort_val_r9_2005$total_disease_mort_val
                    
                    total_disease_mort_val_r9_2017<-total_disease_mort_val_r9%>%
                      filter(year == 2017)
                    total_disease_mort_val_r9_2017<-total_disease_mort_val_r9_2017$total_disease_mort_val
                    
                    
                    #calculate total TB and HIV num mort
                    total_disease_mort_SA_2005<-total_disease_mort_val_r1_2005+
                      total_disease_mort_val_r2_2005+
                      total_disease_mort_val_r3_2005+
                      total_disease_mort_val_r4_2005+
                      total_disease_mort_val_r5_2005+
                      total_disease_mort_val_r6_2005+
                      total_disease_mort_val_r7_2005+
                      total_disease_mort_val_r8_2005+
                      total_disease_mort_val_r9_2005
                    
                    
                    total_disease_mort_SA_2017<-total_disease_mort_val_r1_2017+
                      total_disease_mort_val_r2_2017+
                      total_disease_mort_val_r3_2017+
                      total_disease_mort_val_r4_2017+
                      total_disease_mort_val_r5_2017+
                      total_disease_mort_val_r6_2017+
                      total_disease_mort_val_r7_2017+
                      total_disease_mort_val_r8_2017+
                      total_disease_mort_val_r9_2017
                    
                    if(total_disease_mort_SA_2005 >= total_disease_mort_min_2005 &
                       total_disease_mort_SA_2005 <= total_disease_mort_max_2005 &
                       total_disease_mort_SA_2017 >= total_disease_mort_min_2017 &
                       total_disease_mort_SA_2017 <= total_disease_mort_max_2017){
                        sim_id_array_temp_w_region_id<-c(parameter_set_itr,
                                                         gen_id_temp,
                                                         r1,
                                                         r2,
                                                         r3,
                                                         r4,
                                                         r5,
                                                         r6,
                                                         r7,
                                                         r8,
                                                         r9,
                                                         total_disease_mort_SA_2005,
                                                         total_disease_mort_SA_2017)
                        parameter_sets_df<<-rbind(parameter_sets_df, 
                                                 data.frame(t(sim_id_array_temp_w_region_id)))
                        #every 1,000 accepted sims print param set itr
                        if(parameter_set_itr%%10000 == 0){
                          print(parameter_set_itr)
                          print(Sys.time()-last_time)
                          last_time<<-Sys.time()
                        }
                        parameter_set_itr<<-parameter_set_itr+1
                    } #else{
                      #print('not accepted')
                      #}
                    })
                  })
                })
              })
            })
          })
        })
      })
    })
  colnames(parameter_sets_df)<-c('sim_id', 'general_sim_id',
                    paste("region_id_", region_ids_array, "_sim_id", sep = ''),
                    'total_disease_mort_SA_2005',
                    'total_disease_mort_SA_2017')
  setwd(here('param_files/input_parameters/accepted_param_sets'))
  write.csv(parameter_sets_df, file = paste0('param_set_gen_id_',gen_id_temp, '.csv'),
            row.names = FALSE)
}