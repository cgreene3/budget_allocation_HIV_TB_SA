#clean workspace
rm(list = ls())
gc()

library(here)
library(dplyr)
library(stringr)

setwd(paste0(here(), '/param_files/input_parameters'))
region_id_ref_df<-read.csv('region_id_ref.csv')

setwd(here('/param_files/calculated_param_gen/raw_input_data/GBD'))
pop_df<-read.csv('pop_estimates_15_59.csv')%>%
  group_by(year, region_name)%>%
  summarise(n_pop = sum(expected_total_pop))

setwd(here('/param_files/costs'))
cost_df<-readxl::read_xlsx('cost_df.xlsx')

#objective function df
metrics_df_i_p_k<-data.frame()


###get metrics from TB/HIV model
lapply(region_id_ref_df$region_id, function(region_id_temp){
  region_name_temp<-region_id_ref_df$region_name[region_id_temp]
  
  setwd(paste0(here(), '/results/program_runs/metrics/',
               region_name_temp))
  all_files_in_outdir <- list.files(pattern="summarised_eval_metrics_df*")
  
  if(region_id_temp == 1){
    
    pop_df_region<-pop_df%>%
      filter(region_name == region_name_temp)%>%
      filter(year == 2017)
    
    metrics_df_i_p_k_temp<-read.csv(all_files_in_outdir[1])%>%
      filter(year >= 2018)%>%
      mutate(n_pop = pop_df_region$n_pop[1])%>%
      left_join(cost_df, by = c('program_id'))
    
    metrics_df_i_p_k_temp_total_pop<-metrics_df_i_p_k_temp%>%
      mutate(TBHIV_mort_Y_total = ((TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl)*(n_pop/100000)),
             ART_init_Y_total = ((ART_init_Y_per_100K_ppl)*(n_pop/100000)),
             TPT_init_Y_total = ((TPT_init_Y_per_100K_ppl)*(n_pop/100000)),
             avg_on_ART_Y_total = ((avg_on_ART_Y_per_100K_ppl)*(n_pop/100000))
             )%>%
      mutate(ART_cost = avg_on_ART_Y_total*ART_cost_annual,
             TPT_cost = TPT_init_Y_total*TPT_cost_course,
             HIV_test_cost = ART_init_Y_total*HIV_test_per_one_person_found,
             region_id = region_id_temp)%>%
      group_by(general_sim_id,	regional_sim_id,
               region_id, program_id)%>%
      summarise(TBHIV_mort_num_total = sum(TBHIV_mort_Y_total),
                program_cost = sum(ART_cost)+
                  sum(TPT_cost)+
                  sum(HIV_test_cost))
      
    
    metrics_df_i_p_k_temp_rates<-metrics_df_i_p_k_temp%>%
      filter(year == 2027)%>%
      mutate(TBHIV_mort_per_Y_100k_ppl = TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl,
             region_id = region_id_temp)%>%
      select(c('general_sim_id',	'regional_sim_id',
               'region_id', 'program_id', 
               'TB_inc_per_Y_100k_ppl', 'HIV_prev_Y_100k_ppl', 
               'TBHIV_mort_per_Y_100k_ppl'))
    
    metrics_df_i_p_k<<-metrics_df_i_p_k_temp_rates%>%
      left_join(metrics_df_i_p_k_temp_total_pop, by = c('general_sim_id',	'regional_sim_id',
                                                        'region_id', 'program_id'))
      
    for(f in all_files_in_outdir[2:length(all_files_in_outdir)]){
      metrics_df_i_p_k_temp<-read.csv(f)%>%
        filter(year >= 2018)%>%
        mutate(n_pop = pop_df_region$n_pop[1])%>%
        left_join(cost_df, by = c('program_id'))
        
        metrics_df_i_p_k_temp_total_pop<-metrics_df_i_p_k_temp%>%
          mutate(TBHIV_mort_Y_total = ((TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl)*(n_pop/100000)),
                 ART_init_Y_total = ((ART_init_Y_per_100K_ppl)*(n_pop/100000)),
                 TPT_init_Y_total = ((TPT_init_Y_per_100K_ppl)*(n_pop/100000)),
                 avg_on_ART_Y_total = ((avg_on_ART_Y_per_100K_ppl)*(n_pop/100000))
          )%>%
          mutate(ART_cost = avg_on_ART_Y_total*ART_cost_annual,
            TPT_cost = TPT_init_Y_total*TPT_cost_course,
            HIV_test_cost = ART_init_Y_total*HIV_test_per_one_person_found,
            region_id = region_id_temp)%>%
          group_by(general_sim_id,	regional_sim_id,
                   region_id, program_id)%>%
          summarise(TBHIV_mort_num_total = sum(TBHIV_mort_Y_total),
                    program_cost = sum(ART_cost)+
                      sum(TPT_cost)+
                      sum(HIV_test_cost))
        
        
        metrics_df_i_p_k_temp_rates<-metrics_df_i_p_k_temp%>%
          filter(year == 2027)%>%
          mutate(TBHIV_mort_per_Y_100k_ppl = TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl,
                 region_id = region_id_temp)%>%
          select(c('general_sim_id',	'regional_sim_id',
                   'region_id', 'program_id', 
                   'TB_inc_per_Y_100k_ppl', 'HIV_prev_Y_100k_ppl', 
                   'TBHIV_mort_per_Y_100k_ppl'))
        
        metrics_df_i_p_k_temp2<-metrics_df_i_p_k_temp_rates%>%
          left_join(metrics_df_i_p_k_temp_total_pop, by = c('general_sim_id',	'regional_sim_id',
                                                            'region_id', 'program_id'))
        
      metrics_df_i_p_k<<-rbind(metrics_df_i_p_k,
                               metrics_df_i_p_k_temp2)
    }
    } else {
      pop_df_region<-pop_df%>%
        filter(region_name == region_name_temp)%>%
        filter(year == 2017)
      for(f in all_files_in_outdir){
        metrics_df_i_p_k_temp<-read.csv(f)%>%
          filter(year >= 2018)%>%
          mutate(n_pop = pop_df_region$n_pop[1])%>%
          left_join(cost_df, by = c('program_id'))
        
        metrics_df_i_p_k_temp_total_pop<-metrics_df_i_p_k_temp%>%
          mutate(TBHIV_mort_Y_total = ((TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl)*(n_pop/100000)),
                 ART_init_Y_total = ((ART_init_Y_per_100K_ppl)*(n_pop/100000)),
                 TPT_init_Y_total = ((TPT_init_Y_per_100K_ppl)*(n_pop/100000)),
                 avg_on_ART_Y_total = ((avg_on_ART_Y_per_100K_ppl)*(n_pop/100000))
          )%>%
          mutate(ART_cost = avg_on_ART_Y_total*ART_cost_annual,
            TPT_cost = TPT_init_Y_total*TPT_cost_course,
            HIV_test_cost = ART_init_Y_total*HIV_test_per_one_person_found,
            region_id = region_id_temp)%>%
          group_by(general_sim_id,	regional_sim_id,
                   region_id, program_id)%>%
          summarise(TBHIV_mort_num_total = sum(TBHIV_mort_Y_total),
                    program_cost = sum(ART_cost)+
                      sum(TPT_cost)+
                      sum(HIV_test_cost))
        
        
        metrics_df_i_p_k_temp_rates<-metrics_df_i_p_k_temp%>%
          filter(year == 2027)%>%
          mutate(TBHIV_mort_per_Y_100k_ppl = TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl,
                 region_id = region_id_temp)%>%
          select(c('general_sim_id',	'regional_sim_id',
                   'region_id', 'program_id', 
                   'TB_inc_per_Y_100k_ppl', 'HIV_prev_Y_100k_ppl', 
                   'TBHIV_mort_per_Y_100k_ppl'))
        
        metrics_df_i_p_k_temp2<-metrics_df_i_p_k_temp_rates%>%
          left_join(metrics_df_i_p_k_temp_total_pop, by = c('general_sim_id',	'regional_sim_id',
                                                            'region_id', 'program_id'))
        
        metrics_df_i_p_k<<-rbind(metrics_df_i_p_k,
                                 metrics_df_i_p_k_temp2)
    }
  }
})

#create region program id for filter
metrics_df_i_p_k<-metrics_df_i_p_k%>%
  mutate("region_program_id" = paste0(region_id,'_',program_id))

setwd(paste0(here(), '/results/program_runs/metrics_summarised_i_p_k'))
write.csv(metrics_df_i_p_k, 'metrics_summarised.csv', row.names = FALSE)


filtered_metrics_df_i_p_k_melt_p1<-metrics_df_i_p_k%>%
   select(-c(region_program_id))%>%
   filter(program_id == 1)%>%
   melt(., id = c('general_sim_id', 'regional_sim_id', 'region_id', 'program_id'))%>%
   group_by(region_id, variable)%>%
   summarise(mean = mean(value),
             min = min(value),
             max = max(value))%>%
    filter(variable != 'HIV_prev_Y_100k_ppl')%>%
  filter(variable != 'TB_inc_per_Y_100k_ppl')%>%
    mutate(mean_rounded = if_else(variable == 'program_cost',
                           as.character(round(mean/1000000)),
                           if_else(variable == 'TBHIV_mort_num_total',
                                   as.character(round(mean/1000)),
                                   as.character(round(mean)))),
            ll_rounded = if_else(variable == 'program_cost',
                                 as.character(round(min/1000000)),
                                 if_else(variable == 'TBHIV_mort_num_total',
                                         as.character(round(min/1000)),
                                         as.character(round(min)))),
            ul_rounded = if_else(variable == 'program_cost',
                                 as.character(round(max/1000000)),
                                 if_else(variable == 'TBHIV_mort_num_total',
                                         as.character(round(max/1000)),
                                         as.character(round(max)))))%>%
  mutate(confidence_interval_txt= paste0(mean_rounded, ' [', ll_rounded, ', ', ul_rounded, ']'))%>%
  select(c('region_id', 'variable', 'confidence_interval_txt'))%>%
  arrange(variable)


#we removed program  2(program 3 = program 2)
filtered_metrics_df_i_p_k_melt_p2<-metrics_df_i_p_k%>%
  select(-c(region_program_id))%>%
  filter(program_id == 3)%>%
  melt(., id = c('general_sim_id', 'regional_sim_id', 'region_id', 'program_id'))%>%
  group_by(region_id, variable)%>%
  summarise(mean = mean(value),
            min = min(value),
            max = max(value))%>%
  filter(variable != 'HIV_prev_Y_100k_ppl')%>%
  filter(variable != 'TB_inc_per_Y_100k_ppl')%>%
  mutate(mean_rounded = if_else(variable == 'program_cost',
                                as.character(round(mean/1000000)),
                                if_else(variable == 'TBHIV_mort_num_total',
                                        as.character(round(mean/1000)),
                                        as.character(round(mean)))),
         ll_rounded = if_else(variable == 'program_cost',
                              as.character(round(min/1000000)),
                              if_else(variable == 'TBHIV_mort_num_total',
                                      as.character(round(min/1000)),
                                      as.character(round(min)))),
         ul_rounded = if_else(variable == 'program_cost',
                              as.character(round(max/1000000)),
                              if_else(variable == 'TBHIV_mort_num_total',
                                      as.character(round(max/1000)),
                                      as.character(round(max)))))%>%
  mutate(confidence_interval_txt= paste0(mean_rounded, ' [', ll_rounded, ', ', ul_rounded, ']'))%>%
  select(c('region_id', 'variable', 'confidence_interval_txt'))%>%
  arrange(variable)


metrics_df_i_p_k_p1_v_p2<-metrics_df_i_p_k%>%
  select(-c(region_program_id))%>%
  filter(program_id != 2)%>%
  melt(., id = c('general_sim_id', 'regional_sim_id', 'region_id', 'program_id'))%>%
  mutate(value = if_else(program_id == 1, value, -value))%>%
  filter(variable != 'HIV_prev_Y_100k_ppl')%>%
  filter(variable != 'TB_inc_per_Y_100k_ppl')%>%
  group_by(general_sim_id, regional_sim_id, region_id, variable)%>%
  summarise(value = sum(value))%>%
  mutate(value = if_else(variable == 'program_cost', -value, value))%>%
  group_by(region_id, variable)%>%
  summarise(mean = mean(value),
            min = min(value),
            max = max(value))%>%
  mutate(mean_rounded = if_else(variable == 'program_cost',
                                as.character(round(mean/1000000)),
                                if_else(variable == 'TBHIV_mort_num_total',
                                        as.character(round(mean/1000)),
                                        as.character(round(mean)))),
         ll_rounded = if_else(variable == 'program_cost',
                              as.character(round(min/1000000)),
                              if_else(variable == 'TBHIV_mort_num_total',
                                      as.character(round(min/1000)),
                                      as.character(round(min)))),
         ul_rounded = if_else(variable == 'program_cost',
                              as.character(round(max/1000000)),
                              if_else(variable == 'TBHIV_mort_num_total',
                                      as.character(round(max/1000)),
                                      as.character(round(max)))))%>%
  mutate(confidence_interval_txt= paste0(mean_rounded, ' [', ll_rounded, ', ', ul_rounded, ']'))%>%
  select(c('region_id', 'variable', 'confidence_interval_txt'))%>%
  arrange(variable)

conf_intervals_all<-filtered_metrics_df_i_p_k_melt_p1%>%
  left_join(filtered_metrics_df_i_p_k_melt_p2, by = c('region_id', 'variable'))%>%
  left_join(metrics_df_i_p_k_p1_v_p2, by = c('region_id', 'variable'))

colnames(conf_intervals_all)<-c('region_id', 'variable', 'P1', 'P2', 'P1_v_P2')

conf_intervals_all<-conf_intervals_all%>%
  arrange(region_id)%>%
  arrange(variable)

write.csv(conf_intervals_all, 'conf_intervals_all.csv', row.names = FALSE)