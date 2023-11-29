#clean workspace
rm(list = ls())
gc()

library(here)
library(dplyr)
library(stringr)

#generates n mort if all facility or community based

setwd(paste0(here(), '/param_files/input_parameters'))
region_id_ref_df<-read.csv('region_id_ref.csv')

setwd(here('/param_files/calculated_param_gen/raw_input_data/GBD'))
pop_df<-read.csv('pop_estimates_15_59.csv')%>%
  group_by(year, region_name)%>%
  summarise(n_pop = sum(expected_total_pop))

n_mort_df<-data.frame()

###get metrics from TB/HIV model
lapply(region_id_ref_df$region_id, function(region_id_temp){
  region_name_temp<-region_id_ref_df$region_name[region_id_temp]
  
  setwd(paste0(here(), '/results/program_runs/metrics/',
               region_name_temp))
  all_files_in_outdir <- list.files(pattern="summarised_eval_metrics_df*")
  
  pop_df_region<-pop_df%>%
    filter(region_name == region_name_temp)
  
  pop_df_region_2017<-pop_df_region%>%
    filter(year == 2017)
  
  if(region_id_temp == 1){
    
    n_mort_df_temp<-read.csv(all_files_in_outdir[1])%>%
      left_join(pop_df_region, by = c('year'))%>%
      mutate(n_pop = if_else(year <= 2017, n_pop, pop_df_region_2017$n_pop))%>%
      mutate(n_TBHIV_mort_Y_total = ((TB_mort_per_Y_100k_ppl+O_mort_per_Y_100k_ppl)*(n_pop/100000)))%>%
      mutate(region_id = region_id_temp)%>%
      select(c('year', 'general_sim_id',	'regional_sim_id',
               'region_id', 'program_id', 
               'n_TBHIV_mort_Y_total'))
    
    n_mort_df<<-n_mort_df_temp
    
    for(file_temp in all_files_in_outdir[2:length(all_files_in_outdir)]){
      pop_df_region<-pop_df%>%
        filter(region_name == region_name_temp)
      
      n_mort_df_temp<-read.csv(file_temp)%>%
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
})

#id max and min param set in 2017

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

lapply(param_set_min_max_avg$sim_id, function(sim_temp){
  param_set_min_max_avg_temp<-param_set_min_max_avg%>%
    filter(sim_id == sim_temp)
  sim_ref_array<-param_set_min_max_avg_temp[2:11]
  n_mort_df_temp<-n_mort_df%>%
    filter(general_sim_id == as.integer(sim_ref_array[1]))
  
  n_mort_df_all_for_sim<-data.frame()
  
  for(region_id_temp in region_id_ref_df$region_id){
    regional_param_set<-as.integer(sim_ref_array[region_id_temp+1])
    n_mort_regional_temp_df<-n_mort_df_temp%>%
      filter(region_id == region_id_temp)%>%
      filter(regional_sim_id == regional_param_set)
    if(region_id_temp == 1){
      n_mort_df_all_for_sim<-n_mort_regional_temp_df
    } else{
      n_mort_df_all_for_sim<-rbind(n_mort_df_all_for_sim, n_mort_regional_temp_df)
    }
    
  }
  
  n_mort_summarised_year_prog_id<-n_mort_df_all_for_sim%>%
    group_by(year, program_id)%>%
    summarise(n_mort_total = sum(n_TBHIV_mort_Y_total))
  
  n_mort_summarised_year_prog_id$sim_id<-rep(sim_temp, times = nrow(n_mort_summarised_year_prog_id))
  
  if(nrow(n_mort_summarised_year_prog_id_all) == 0){
    n_mort_summarised_year_prog_id_all<<-n_mort_summarised_year_prog_id
  } else {
    n_mort_summarised_year_prog_id_all<<-rbind(n_mort_summarised_year_prog_id_all,
                                              n_mort_summarised_year_prog_id)
  }
})

program_eval_summarised<-n_mort_summarised_year_prog_id_all%>%
  group_by(year, program_id)%>%
  summarise(upper_rate = max(n_mort_total),
            lower_rate = min(n_mort_total),
            val_rate = median(n_mort_total))%>%
  filter(program_id != 2) #we do not consider program 2

program_eval_summarised_connect2<-program_eval_summarised%>%
  filter(year == 2017)

last_yr_calib_est<-program_eval_summarised_connect2$val_rate #so 2017 calib estimate is in black


program_eval_summarised_connect2<-rbind(program_eval_summarised_connect2,
                                        program_eval_summarised_connect2)

program_eval_summarised_connect2$group<-c('If all regions implement standard\nfacility-based ART and TPT (p=1)',
                                          'If all regions implement community-based\nART with TPT (p=2)')

program_eval_summarised$group<-if_else(program_eval_summarised$year <= 2017,
                                       "\nCalibration period\n", 
                                       if_else(program_eval_summarised$program_id == 1,
                                               'If all regions implement standard\nfacility-based ART and TPT (p=1)',
                                               'If all regions implement community-based\nART with TPT (p=2)'))

program_eval_summarised<-rbind(program_eval_summarised,
                               program_eval_summarised_connect2)

sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'extrafont', 'RColorBrewer', 'gridExtra'), require, character.only=T)


colors_for_graphs <- brewer.pal(n = 10, name = "Paired")
colors_for_program_graph_fill<-colors_for_graphs[c(7, 9)]
colors_for_program_graph_line<-colors_for_graphs[c(8, 10)]

library(scales)

program_graph_temp<-ggplot(program_eval_summarised)+
  geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = group))+
  geom_point(aes(x = year, y = val_rate, colour = group, shape = group), size = 2)+
  geom_point(aes(x = 2017, y = last_yr_calib_est), colour = "black", size = 2)+
  theme(text = element_text(size=18, family="Times New Roman"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top", legend.title = element_blank(), legend.text=element_text(size=22),
        plot.title = element_text(hjust = .5, size=20),
        legend.background = element_rect(colour = "lightgrey"),
        legend.box.background = element_rect(colour = "black"))+
  geom_vline(xintercept = 2017, linetype="dashed", 
             color = "darkgrey", size=1.5)+
  annotate("text", x=2012, y=400000, label= "calibration period", size = 7)+
  annotate("text", x=2022, y=400000, label= "intervention period", size = 7)+
  #ggtitle(current_var_eval)+
  scale_x_continuous(name = 'year y', breaks=c(seq(from = 1990, to = 2028, by = 3)))+
  scale_color_manual(values=c("black", colors_for_program_graph_line))+
  scale_fill_manual(values=c('grey75', colors_for_program_graph_fill))+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3),
                     limits = c(0, 501000))+
  ylab('Number of TB/HIV-related deaths')+
  scale_shape_manual(values=c(16, 16, 8, 5))

setwd(paste0(here(), "/results"))

png('n_disease_mort_all_or_none_facil_community.png', width = 1100, height = 400)
print(program_graph_temp)
dev.off()
