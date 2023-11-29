#clean workspace
rm(list = ls())
gc()

library(here)
library(dplyr)
library(stringr)

#generates n mort by year if all facility or community based
#generates total n mort if all facility or community based
#generates total cost all facility or community based

setwd(paste0(here(), '/param_files/input_parameters'))
region_id_ref_df<-read.csv('region_id_ref.csv')

setwd(here('/param_files/calculated_param_gen/raw_input_data/GBD'))
pop_df<-read.csv('pop_estimates_15_59.csv')%>%
  group_by(year, region_name)%>%
  summarise(n_pop = sum(expected_total_pop))

setwd(here('/param_files/costs'))
cost_df<-readxl::read_xlsx('cost_df.xlsx')

n_mort_df<-data.frame()
cost_df_all<-data.frame()

###get metrics from TB/HIV model
lapply(region_id_ref_df$region_id, function(region_id_temp){
  
  #get region name
  region_name_temp<-region_id_ref_df$region_name[region_id_temp]
  
  #get all results related to region
  setwd(paste0(here(), '/results/program_runs/metrics/',
               region_name_temp))
  all_files_in_outdir <- list.files(pattern="summarised_eval_metrics_df*")
  
  pop_df_region<-pop_df%>%
    filter(region_name == region_name_temp)
  
  pop_df_region_2017<-pop_df_region%>%
    filter(year == 2017)
  
  ####get n_mort by year#####
  if(region_id_temp == 1){
    
    n_mort_df_temp<-read.csv(all_files_in_outdir[1])%>%
      filter(program_id != 2)%>%
      left_join(pop_df_region, by = c('year'))%>%
      mutate(n_pop = if_else(year <= 2017, n_pop, pop_df_region_2017$n_pop))%>%
      mutate(n_TBHIV_mort_Y_total = ((TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl)*(n_pop/100000)))%>%
      mutate(region_id = region_id_temp)%>%
      select(c('year', 'general_sim_id',	'regional_sim_id',
               'region_id', 'program_id', 
               'n_TBHIV_mort_Y_total'))
    
    n_mort_df<<-n_mort_df_temp
    
    for(file_temp in all_files_in_outdir[2:length(all_files_in_outdir)]){
      
      n_mort_df_temp<-read.csv(file_temp)%>%
        filter(program_id != 2)%>%
        left_join(pop_df_region, by = c('year'))%>%
        mutate(n_pop = if_else(year <= 2017, n_pop, pop_df_region_2017$n_pop))%>%
        mutate(n_TBHIV_mort_Y_total = ((TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl)*(n_pop/100000)))%>%
        mutate(region_id = region_id_temp)%>%
        select(c('year', 'general_sim_id',	'regional_sim_id',
                 'region_id', 'program_id', 
                 'n_TBHIV_mort_Y_total'))
      
      n_mort_df<<-rbind(n_mort_df, n_mort_df_temp)
    }
  }
  
  if(region_id_temp > 1){
    for(file_temp in all_files_in_outdir){
    
      n_mort_df_temp<-read.csv(file_temp)%>%
        filter(program_id != 2)%>%
        left_join(pop_df_region, by = c('year'))%>%
        mutate(n_pop = if_else(year <= 2017, n_pop, pop_df_region_2017$n_pop))%>%
        mutate(n_TBHIV_mort_Y_total = ((TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl)*(n_pop/100000)))%>%
        mutate(region_id = region_id_temp)%>%
        select(c('year', 'general_sim_id',	'regional_sim_id',
                 'region_id', 'program_id', 
                 'n_TBHIV_mort_Y_total'))
      n_mort_df<<-rbind(n_mort_df, n_mort_df_temp)
    }
  }
  
  #get total cost
  if(region_id_temp == 1){
    
    cost_temp_df<-read.csv(all_files_in_outdir[1])%>%
      filter(program_id != 2)%>%
      filter(year >= 2018)%>%
      mutate(n_pop = pop_df_region_2017$n_pop[1])%>%
      left_join(cost_df, by = c('program_id'))%>%
      mutate(ART_init_Y_total = ((ART_init_Y_per_100K_ppl)*(n_pop/100000)),
             TPT_init_Y_total = ((TPT_init_Y_per_100K_ppl)*(n_pop/100000)),
             avg_on_ART_Y_total = ((avg_on_ART_Y_per_100K_ppl)*(n_pop/100000)))%>%
      mutate(ART_cost = avg_on_ART_Y_total*ART_cost_annual,
             TPT_cost = TPT_init_Y_total*TPT_cost_course,
             HIV_test_cost = ART_init_Y_total*HIV_test_per_one_person_found,
             region_id = region_id_temp)%>%
      group_by(general_sim_id,	regional_sim_id,
               region_id, program_id)%>%
      summarise(program_cost = sum(ART_cost)+
                  sum(TPT_cost)+
                  sum(HIV_test_cost))
    
    cost_df_all<<-cost_temp_df
    
    for(file_temp in all_files_in_outdir[2:length(all_files_in_outdir)]){
      cost_temp_df<-read.csv(file_temp)%>%
        filter(program_id != 2)%>%
        filter(year >= 2018)%>%
        mutate(n_pop = pop_df_region_2017$n_pop[1])%>%
        left_join(cost_df, by = c('program_id'))%>%
        mutate(ART_init_Y_total = ((ART_init_Y_per_100K_ppl)*(n_pop/100000)),
               TPT_init_Y_total = ((TPT_init_Y_per_100K_ppl)*(n_pop/100000)),
               avg_on_ART_Y_total = ((avg_on_ART_Y_per_100K_ppl)*(n_pop/100000)))%>%
        mutate(ART_cost = avg_on_ART_Y_total*ART_cost_annual,
               TPT_cost = TPT_init_Y_total*TPT_cost_course,
               HIV_test_cost = ART_init_Y_total*HIV_test_per_one_person_found,
               region_id = region_id_temp)%>%
        group_by(general_sim_id,	regional_sim_id,
                 region_id, program_id)%>%
        summarise(program_cost = sum(ART_cost)+
                    sum(TPT_cost)+
                    sum(HIV_test_cost))
      
      cost_df_all<<-rbind(cost_temp_df, cost_df_all)
    }
    }
    
    if(region_id_temp > 1){
      for(file_temp in all_files_in_outdir){
        cost_temp_df<-read.csv(file_temp)%>%
          filter(program_id != 2)%>%
          filter(year >= 2018)%>%
          mutate(n_pop = pop_df_region_2017$n_pop[1])%>%
          left_join(cost_df, by = c('program_id'))%>%
          mutate(ART_init_Y_total = ((ART_init_Y_per_100K_ppl)*(n_pop/100000)),
                 TPT_init_Y_total = ((TPT_init_Y_per_100K_ppl)*(n_pop/100000)),
                 avg_on_ART_Y_total = ((avg_on_ART_Y_per_100K_ppl)*(n_pop/100000)))%>%
          mutate(ART_cost = avg_on_ART_Y_total*ART_cost_annual,
                 TPT_cost = TPT_init_Y_total*TPT_cost_course,
                 HIV_test_cost = ART_init_Y_total*HIV_test_per_one_person_found,
                 region_id = region_id_temp)%>%
          group_by(general_sim_id,	regional_sim_id,
                   region_id, program_id)%>%
          summarise(program_cost = sum(ART_cost)+
                      sum(TPT_cost)+
                      sum(HIV_test_cost))
        
        cost_df_all<<-rbind(cost_temp_df, cost_df_all)
      }
    }
    
  
})

#id max and min param set in 2017 (n_mort)
setwd(paste0(here(), '/param_files/input_parameters/accepted_param_sets_combin'))
all_files_in_outdir <- list.files(pattern="param_set_gen*")

param_set_all<-read.csv(all_files_in_outdir[1])
  
lapply(all_files_in_outdir[2:length(all_files_in_outdir)], function(file){
  param_set_all<<-rbind(param_set_all, read.csv(file))
})

param_set_min_max_avg<-param_set_all%>%
  group_by()%>%
  mutate(avg_2017 = mean(total_disease_mort_SA_2017),
         diff_2017 = abs(total_disease_mort_SA_2017-avg_2017))%>%
  ungroup()%>%
  filter(total_disease_mort_SA_2017 == max(total_disease_mort_SA_2017)|
           total_disease_mort_SA_2017 == min(total_disease_mort_SA_2017)|
           diff_2017 == min(diff_2017))

rm(param_set_all)


#calculate n mort each year and each program for the max min and avg parameter sets

n_mort_summarised_year_prog_id_all<-data.frame()
cost_summarised_prog_id_all<-data.frame()

lapply(param_set_min_max_avg$sim_id, function(sim_temp){
  
  param_sets_min_max_avg_filtered_sim<-param_set_min_max_avg%>%
    filter(sim_id == sim_temp)
  
  sim_ref_array<-param_sets_min_max_avg_filtered_sim[2:11]
  
  n_mort_df_temp_filtered<-n_mort_df%>%
    filter(general_sim_id == as.integer(sim_ref_array[1]))
  
  cost_df_temp_filtered<-cost_df_all%>%
    filter(general_sim_id == as.integer(sim_ref_array[1]))
  
  n_mort_df_all_for_sim<-data.frame()
  cost_summarised_year_prog_for_sim<-data.frame()
  
  for(region_id_temp in region_id_ref_df$region_id){
    regional_param_set<-as.integer(sim_ref_array[region_id_temp+1])
    n_mort_regional_temp_df<-n_mort_df_temp_filtered%>%
      filter(region_id == region_id_temp)%>%
      filter(regional_sim_id == regional_param_set)
    
    cost_regional_temp_df<-cost_df_temp_filtered%>%
      filter(regional_sim_id == regional_param_set)
      
    if(region_id_temp == 1){
      n_mort_df_all_for_sim<-n_mort_regional_temp_df
      cost_summarised_year_prog_for_sim<-cost_regional_temp_df
    } else{
      n_mort_df_all_for_sim<-rbind(n_mort_df_all_for_sim, n_mort_regional_temp_df)
      cost_summarised_year_prog_for_sim<-rbind(cost_summarised_year_prog_for_sim,
                                              cost_regional_temp_df)
    }
    
  }
  
  n_mort_summarised_year_prog_id<-n_mort_df_all_for_sim%>%
    group_by(year, program_id)%>%
    summarise(n_mort_total = sum(n_TBHIV_mort_Y_total))
  
  cost_summarised_prog_id<-cost_summarised_year_prog_for_sim%>%
    group_by(program_id)%>%
    summarise(cost_total = sum(program_cost))
  
  n_mort_summarised_year_prog_id$sim_id<-rep(sim_temp, 
                                             times = nrow(n_mort_summarised_year_prog_id))
  cost_summarised_prog_id$sim_id<-rep(sim_temp, 
                                      times = nrow(cost_summarised_prog_id))
  
  if(nrow(n_mort_summarised_year_prog_id_all) == 0){
    n_mort_summarised_year_prog_id_all<<-n_mort_summarised_year_prog_id
    cost_summarised_prog_id_all<<-cost_summarised_prog_id
  } else {
    n_mort_summarised_year_prog_id_all<<-rbind(n_mort_summarised_year_prog_id_all,
                                              n_mort_summarised_year_prog_id)
    cost_summarised_prog_id_all<<-rbind(cost_summarised_prog_id_all,
                                        cost_summarised_prog_id)
  }
})


n_mort_summarised_year_prog_id_2017<-n_mort_summarised_year_prog_id_all%>%
  filter(year > 2017)%>%
  group_by(program_id, sim_id)%>%
  summarise(total_mort = sum(n_mort_total))

#cost is right from cost_summarised_prog_id_all

