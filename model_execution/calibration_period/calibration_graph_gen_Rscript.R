#last updated aug 23, 2023
#run after program runs (get metrics overtime)

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle',
         'extrafont', 
         'RColorBrewer'), require, character.only=T)

#####COLORS FOR PROGRAM GRAPHS#######
colors_for_graphs <- brewer.pal(n = 10, name = "Paired")
colors_for_program_graph_fill<-colors_for_graphs[c(1, 3, 5)]
colors_for_program_graph_line<-colors_for_graphs[c(2, 4, 6)]


#add transparecy
for (i in 1:length(colors_for_program_graph_fill)){
  rbg_temp<-as.vector(col2rgb(colors_for_program_graph_fill[i]))
  mycol <- rgb(rbg_temp[1], rbg_temp[2], rbg_temp[3], 
               max = 255, alpha = 185, 
               names = colors_for_program_graph_fill[i])
  colors_for_program_graph_fill[i]<-mycol
}

#specify graph y axis
HIVTB_mort_rate_max_min_breaks_graph_y_axis<-c(0, 2700, 900)
TBinc_rate_max_min_breaks_graph_y_axis<-c(0, 2700, 900)
HIVprev_max_min_breaks_graph_y_axis<-c(0, 40000, 5000)


##GITHUB/LOCAL###
#location where we sent outputs
library(here)
indir_outputs<-paste0(here(), '/results/program_runs/metrics/')
indir_target_calibration_estimates<-paste0(here(), '/param_files/target_calibration_estimates')
indir_input_params<-paste0(here(), '/param_files/input_parameters')
# 
# #location where we send calib analysis results
outdir<-paste0(here(), '/results/calibration_graphs')

setwd(indir_target_calibration_estimates)
GBD_HIV_prev_df<-read.csv('SA_regional_GBD_HIV_prev_rate_calibration_df.csv')%>%
  mutate(group = "GBD Projections")
GBD_TB_inc_df<-read.csv('SA_regional_GBD_TB_inc_rate_calibration_df.csv')%>%
  mutate(group = "GBD Projections")
GBD_TBHIV_mort_df<-read.csv('SA_regional_GBD_TBHIV_mort_rate_calibration_df.csv')%>%
  mutate(group = "GBD Projections")
GBD_TBHIV_mort_num_df<-read.csv('SA_all_GBD_TBHIV_mort_num_calibration_df.csv')%>%
  mutate(group = "GBD Projections")

setwd(indir_input_params)
#accepted_param_sets_ref_df<-read.csv('accepted_param_sets_ref.csv')
region_id_ref_df<-read.csv('region_id_ref.csv')
  
#read in pop df
setwd(paste0(here(), '/param_files/calculated_param_gen/raw_input_data/GBD/'))
pop_df<-read.csv('pop_estimates_15_59.csv')%>%
  filter(year >= 1990,
         year <= 2017)%>%
  group_by(region_name, year)%>%
  summarise(expected_total_pop = sum(expected_total_pop))%>%
  left_join(region_id_ref_df, by = c('region_name'))


num_TBHIV_mort_df_all<-data.frame()

graph_list<-list()
graph_itr_temp<-0

stats_2017_df<-data.frame()

lapply(region_id_ref_df$region_id, function(region_id){
  
  region_name_temp<-region_id_ref_df$region_name[region_id]
  
  setwd(paste0(indir_outputs, region_name_temp))
  all_files_in_outdir <- list.files(pattern="summarised_eval_metrics_df*")
  outputs_combined_df<-read.csv(all_files_in_outdir[1])%>%
    filter(year == 2017)%>% #before program eval
    select(-c('program_id'))
  
  for(n in 2:length(all_files_in_outdir)){
    temp<-read.csv(all_files_in_outdir[n])%>%
      filter(year == 2017)%>% #before program eval
      select(-c('program_id'))
    outputs_combined_df<-rbind(outputs_combined_df, temp)
  }
  outputs_combined_df$region_id <-rep(region_id, times = nrow(outputs_combined_df))
  outputs_combined_df<-outputs_combined_df%>%
    mutate(TBHIV_mort_per_Y_100k_ppl = TB_mort_per_Y_100k_ppl+ O_mort_per_Y_100k_ppl)%>%
    select(c('region_id', 'TB_inc_per_Y_100k_ppl',
             'HIV_prev_Y_100k_ppl', 'TBHIV_mort_per_Y_100k_ppl'))
  
  if(region_id == 1){
    stats_2017_df<<-outputs_combined_df
  } else {
    stats_2017_df<<-rbind(outputs_combined_df,
                          stats_2017_df)
  }
})

stats_2017_df_2<-stats_2017_df%>%
  group_by(region_id)%>%
  summarise(mean_TB_inc_per_Y_100k_ppl = mean(TB_inc_per_Y_100k_ppl),
            min_TB_inc_per_Y_100k_ppl = min(TB_inc_per_Y_100k_ppl),
            max_TB_inc_per_Y_100k_ppl = max(TB_inc_per_Y_100k_ppl),
            mean_HIV_prev_Y_100k_ppl = mean(HIV_prev_Y_100k_ppl),
            min_HIV_prev_Y_100k_ppl = min(HIV_prev_Y_100k_ppl),
            max_HIV_prev_Y_100k_ppl = max(HIV_prev_Y_100k_ppl),
            mean_TBHIV_mort_per_Y_100k_ppl = mean(TBHIV_mort_per_Y_100k_ppl),
            min_TBHIV_mort_per_Y_100k_ppl = min(TBHIV_mort_per_Y_100k_ppl),
            max_TBHIV_mort_per_Y_100k_ppl = max(TBHIV_mort_per_Y_100k_ppl))

lapply(region_id_ref_df$region_id, function(region_id){
  
  region_name_temp<-region_id_ref_df$region_name[region_id]
  
  setwd(paste0(indir_outputs, region_name_temp))
  all_files_in_outdir <- list.files(pattern="summarised_eval_metrics_df*")
  outputs_combined_df<-read.csv(all_files_in_outdir[1])%>%
    filter(year < 2018)%>% #before program eval
    select(-c('program_id'))
  
  for(n in 2:length(all_files_in_outdir)){
    temp<-read.csv(all_files_in_outdir[n])%>%
      filter(year < 2018)%>% #before program eval
      select(-c('program_id'))
    outputs_combined_df<-rbind(outputs_combined_df, temp)
  }
  
  outputs_combined_df_reshape<-melt(outputs_combined_df, id = c("year", "general_sim_id", "regional_sim_id"))
  
  #get total disease mort
  num_TBHIV_mort_df_temp<-outputs_combined_df_reshape%>%
    filter(grepl('mort', variable) == TRUE)%>%
    group_by(general_sim_id, regional_sim_id, year)%>%
    summarise(value = sum(value))%>% ##sum TB and o mort
    left_join(pop_df%>%
                filter(region_name == region_name_temp), 
              by = c('year'))%>%
    mutate(num_TBHIV_mort = value*(expected_total_pop/100000))%>%
    select(-c('value', 'expected_total_pop', 'region_name'))
  
  if(region_id == 1){
    num_TBHIV_mort_df_all<<-num_TBHIV_mort_df_temp
  } else{
    num_TBHIV_mort_df_all<<-rbind(num_TBHIV_mort_df_all,
                                  num_TBHIV_mort_df_temp)
  }
  
  ####make data frames for calib graphs TB HIV mort#####
  calib_TBHIV_mort_rate<-outputs_combined_df_reshape%>%
    filter(grepl('mort', variable) == TRUE)%>%
    group_by(general_sim_id, regional_sim_id, year)%>%
    summarise(value = sum(value))%>%
    mutate(variable = 'TBHIV_mort_per_Y_100k_ppl')%>%
    group_by(year)%>%
    summarise(val_rate = mean(value),
              upper_rate = max(value),
              lower_rate = min(value))%>%
    mutate(group = "Model Projections")
  
  combined_calib_TBHIV_mort_df<-GBD_TBHIV_mort_df%>%
    filter(year <= 2017,
           region_name == region_name_temp)%>%
    select(-c('region_name'))%>%
    bind_rows(calib_TBHIV_mort_rate)
  
  
  ####generate graph####
  calib_graph_TBHIV_mort_temp<-ggplot(combined_calib_TBHIV_mort_df)+
    geom_ribbon(aes(ymin = lower_rate, 
                    ymax = upper_rate, 
                    x = year, 
                    fill = group))+
    geom_point(aes(x = year, y = val_rate, colour = group, shape = group))+
    scale_fill_manual(values=colors_for_program_graph_fill) +
    scale_colour_manual(values=colors_for_program_graph_line)+
    scale_shape_manual(values=c(0, 19))+
    scale_size_manual(values=c(1,3))+
    geom_segment(aes(x = 2005, y = unlist(GBD_TBHIV_mort_df%>%
                                            filter(year == 2005,
                                                   region_name == region_name_temp)%>%
                                            select(c("lower_rate"))),
                     xend = 2005, yend = unlist(GBD_TBHIV_mort_df%>%
                                                  filter(year == 2005,
                                                         region_name == region_name_temp)%>%
                                                  select(c("upper_rate")))),
                 col = "darkblue", size = 1.2, 
                 arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
    geom_segment(aes(x = 2017, y = unlist(GBD_TBHIV_mort_df%>%
                                            filter(year == 2017,
                                                   region_name == region_name_temp)%>%
                                            select(c("lower_rate"))),
                     xend = 2017, yend = unlist(GBD_TBHIV_mort_df%>%
                                                  filter(year == 2017,
                                                         region_name == region_name_temp)%>%
                                                  select(c("upper_rate")))),
                 col = "darkblue", size = 1.2, 
                 arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
    theme(text = element_text(size=24, family="Times New Roman"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.2, 0.7), 
          legend.title = element_blank(), 
          legend.text=element_text(size=22),
          plot.title = element_text(hjust = .5, size=24),
          #legend.background = element_rect(colour = "lightgrey"),
          #legend.box.background = element_rect(colour = "black")
          )+
    scale_x_continuous(name = 'year y', breaks=c(seq(from = 1990, to = 2017, by = 3)))+
    ylim(HIVTB_mort_rate_max_min_breaks_graph_y_axis[1:2])+
    ylab('TB/HIV disease mortality rate,\n per 100K individuals')+
    annotate("text", x = 2005, 
             y = HIVTB_mort_rate_max_min_breaks_graph_y_axis[2]*.95, 
             label = paste0(region_name_temp, ' (i = ', region_id, ')\n', 'TB and HIV-related mortality rate'),
             size = 10,
             family="Times New Roman")
  
  ####TB incidence#####
  calib_TBinc_rate<-outputs_combined_df_reshape%>%
    filter(grepl('TB_inc', variable) == TRUE)%>%
    group_by(year)%>%
    summarise(val_rate = mean(value),
              upper_rate = max(value),
              lower_rate = min(value))%>%
    mutate(group = "Model Projections")
  
  combined_calib_TB_inc_df<-GBD_TB_inc_df%>%
    filter(year <= 2017,
           region_name == region_name_temp)%>%
    select(-c('region_name'))%>%
    bind_rows(calib_TBinc_rate)
  
  ####generate graph####
  calib_graph_TB_inc_temp<-ggplot(combined_calib_TB_inc_df)+
    geom_ribbon(aes(ymin = lower_rate, 
                    ymax = upper_rate, 
                    x = year, 
                    fill = group))+
    geom_point(aes(x = year, y = val_rate, colour = group, shape = group))+
    scale_fill_manual(values=colors_for_program_graph_fill) +
    scale_colour_manual(values=colors_for_program_graph_line)+
    scale_shape_manual(values=c(0, 19))+
    scale_size_manual(values=c(1,3))+
    geom_segment(aes(x = 2005, y = unlist(GBD_TB_inc_df%>%
                                            filter(year == 2005,
                                                   region_name == region_name_temp)%>%
                                            select(c("lower_rate"))),
                     xend = 2005, yend = unlist(GBD_TB_inc_df%>%
                                                  filter(year == 2005,
                                                         region_name == region_name_temp)%>%
                                                  select(c("upper_rate")))),
                 col = "darkblue", size = 1.2, 
                 arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
    geom_segment(aes(x = 2017, y = unlist(GBD_TB_inc_df%>%
                                            filter(year == 2017,
                                                   region_name == region_name_temp)%>%
                                            select(c("lower_rate"))),
                     xend = 2017, yend = unlist(GBD_TB_inc_df%>%
                                                  filter(year == 2017,
                                                         region_name == region_name_temp)%>%
                                                  select(c("upper_rate")))),
                 col = "darkblue", size = 1.2, 
                 arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
    theme(text = element_text(size=24, family="Times New Roman"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.2, 0.7), legend.title = element_blank(), 
          legend.text=element_text(size=22),
          plot.title = element_text(hjust = .5, size=24, region_name_temp),
          #legend.background = element_rect(colour = "lightgrey"),
          #legend.box.background = element_rect(colour = "black")
          )+
    scale_x_continuous(name = 'year y', breaks=c(seq(from = 1990, to = 2017, by = 3)))+
    ylim(TBinc_rate_max_min_breaks_graph_y_axis[1:2])+
    ylab('TB incidence rate,\n per 100K individuals')+
    annotate("text", x = 2005, 
             y = TBinc_rate_max_min_breaks_graph_y_axis[2]*.95, 
             label = paste0(region_name_temp, ' (i = ', region_id, ')\n TB incidence rate'),
             size = 10,
             family="Times New Roman")
  
  
  ####HIV prevalence rate####
  calib_HIVprev_rate<-outputs_combined_df_reshape%>%
    filter(grepl('HIV_prev', variable) == TRUE)%>%
    group_by(year)%>%
    summarise(val_rate = mean(value),
              upper_rate = max(value),
              lower_rate = min(value))%>%
    mutate(group = "Model Projections")
  
  combined_calib_HIV_prev_df<-GBD_HIV_prev_df%>%
    filter(year <= 2017,
           region_name == region_name_temp)%>%
    select(-c('region_name'))%>%
    bind_rows(calib_HIVprev_rate)
  
  ####generate graph####
  calib_graph_HIV_prev_temp<-ggplot(combined_calib_HIV_prev_df)+
    geom_ribbon(aes(ymin = lower_rate, 
                    ymax = upper_rate, 
                    x = year, 
                    fill = group))+
    geom_point(aes(x = year, y = val_rate, colour = group, shape = group))+
    scale_fill_manual(values=colors_for_program_graph_fill) +
    scale_colour_manual(values=colors_for_program_graph_line)+
    scale_shape_manual(values=c(0, 19))+
    scale_size_manual(values=c(1,3))+
    geom_segment(aes(x = 2005, y = unlist(GBD_HIV_prev_df%>%
                                            filter(year == 2005,
                                                   region_name == region_name_temp)%>%
                                            select(c("lower_rate"))),
                     xend = 2005, yend = unlist(GBD_HIV_prev_df%>%
                                                  filter(year == 2005,
                                                         region_name == region_name_temp)%>%
                                                  select(c("upper_rate")))),
                 col = "darkblue", size = 1.2, 
                 arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
    geom_segment(aes(x = 2017, y = unlist(GBD_HIV_prev_df%>%
                                            filter(year == 2017,
                                                   region_name == region_name_temp)%>%
                                            select(c("lower_rate"))),
                     xend = 2017, yend = unlist(GBD_HIV_prev_df%>%
                                                  filter(year == 2017,
                                                         region_name == region_name_temp)%>%
                                                  select(c("upper_rate")))),
                 col = "darkblue", size = 1.2, 
                 arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
    theme(text = element_text(size=24, family="Times New Roman"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.2, 0.7), 
          legend.title = element_blank(), 
          legend.text=element_text(size=22),
          plot.title = element_text(hjust = .5, size=24)#,
          #legend.background = element_rect(colour = "lightgrey"),
          #legend.box.background = element_rect(colour = "black")
          )+
    scale_x_continuous(name = 'year y', breaks=c(seq(from = 1990, to = 2017, by = 3)))+
    ylim(HIVprev_max_min_breaks_graph_y_axis[1:2])+
    ylab('HIV prevalence rate,\n per 100K individuals')+
    annotate("text", x = 2005, 
             y = HIVprev_max_min_breaks_graph_y_axis[2]*.95, 
             label = paste0(region_name_temp, ' (i = ', region_id, ')\n HIV prevalence rate'),
             size = 10,
             family="Times New Roman")
  
  graph_list[[region_id]]<<-list(calib_graph_TB_inc_temp,
                              calib_graph_HIV_prev_temp,
                              calib_graph_TBHIV_mort_temp)
  
  setwd(paste0(outdir))
  #setwd(paste0(outdir, '/', region_name_temp))
  
  graph_itr_temp<<-graph_itr_temp+1
  file_name2<-paste0(graph_itr_temp,
                     '_TB_inc_rate_calib_graph_region_',
                     region_id,
                     ".png")
  png(file_name2, width = 800, height = 500)
  print(calib_graph_TB_inc_temp)
  dev.off()
  
  graph_itr_temp<<-graph_itr_temp+1
  
  file_name3<-paste0(graph_itr_temp,
                     '_HIV_prev_rate_calib_graph_region_',
                     region_id,
                     ".png")
  png(file_name3, width = 800, height = 500)
  print(calib_graph_HIV_prev_temp)
  dev.off()
  
  graph_itr_temp<<-graph_itr_temp+1
  file_name1<-paste0(graph_itr_temp,
                     '_TBHIV_mort_rate_calib_graph_region_',
                    region_id,
                   ".png")
  png(file_name1, width = 800, height = 500)
  print(calib_graph_TBHIV_mort_temp)
  dev.off()
})


####TOTAL TB HIV num pop####
setwd(paste0(here(), '/param_files/input_parameters/accepted_param_sets_combin'))


all_files_in_outdir <- list.files(pattern="param_set*")
outputs_combined_df<-read.csv(all_files_in_outdir[1])%>%
  select(-c('total_disease_mort_SA_2005',
            'total_disease_mort_SA_2017'))
  
lapply(2:length(all_files_in_outdir), function(n){
  temp<-read.csv(all_files_in_outdir[n])%>%
    select(-c('total_disease_mort_SA_2005',
              'total_disease_mort_SA_2017'))
  outputs_combined_df<<-rbind(outputs_combined_df, temp)
})

outputs_combined_df_reshape<-melt(outputs_combined_df, 
                                   id = c("sim_id",
                                          "general_sim_id"))

outputs_combined_df_reshape<-outputs_combined_df_reshape%>%
  rename("regional_sim_id" = "value")%>%
  mutate(region_id = str_split(variable, "_", simplify = T)[, 3])%>%
  select(-c('variable'))

num_TBHIV_mort_calcs_df_all<-data.frame()

lapply(1990:2017, function(yr){
  print(yr)
  print(Sys.time())
  num_TBHIV_mort_calcs_df<-outputs_combined_df_reshape%>%
    mutate(region_id = as.integer(region_id))%>%
    left_join(num_TBHIV_mort_df_all%>%
                filter(year == yr), by = c("general_sim_id", 
                                           "regional_sim_id",
                                           "region_id"))
  
  num_TBHIV_mort_calcs_df<-num_TBHIV_mort_calcs_df%>%
    group_by(year, sim_id)%>%
    summarise(SA_TBHIV_mort = sum(num_TBHIV_mort))%>%
    group_by(year)%>%
    summarise(val_mort_num = mean(SA_TBHIV_mort),
              upper_mort_num = max(SA_TBHIV_mort),
              lower_mort_num = min(SA_TBHIV_mort))%>%
    mutate(group = "Model Projections")
  
  if(yr == 1990){
    num_TBHIV_mort_calcs_df_all<<-num_TBHIV_mort_calcs_df
  } else {
    num_TBHIV_mort_calcs_df_all<<-rbind(num_TBHIV_mort_calcs_df_all,
                                        num_TBHIV_mort_calcs_df)
  }
  
})
  
num_TBHIV_mort_cals_df_GBD_and_MP<-num_TBHIV_mort_calcs_df_all%>%
  bind_rows(GBD_TBHIV_mort_num_df)%>%
  filter(year <= 2017)

max_min_num_mort<-c(0, 
                    round(max(num_TBHIV_mort_cals_df_GBD_and_MP$upper_mort_num)*1.01,
                          digits = -3))


####generate graph####
library(scales)

SA_num_TBHIV_mort_graph<-ggplot(num_TBHIV_mort_cals_df_GBD_and_MP)+
  geom_ribbon(aes(ymin = lower_mort_num, 
                  ymax = upper_mort_num, 
                  x = year, 
                  fill = group))+
  geom_point(aes(x = year, y = val_mort_num, colour = group, shape = group))+
  scale_fill_manual(values=colors_for_program_graph_fill) +
  scale_colour_manual(values=colors_for_program_graph_line)+
  scale_shape_manual(values=c(0, 19))+
  scale_size_manual(values=c(1,3))+
  geom_segment(aes(x = 2005, y = unlist(GBD_TBHIV_mort_num_df%>%
                                          filter(year == 2005)%>%
                                          select(c("lower_mort_num"))),
                   xend = 2005, yend = unlist(GBD_TBHIV_mort_num_df%>%
                                                filter(year == 2005)%>%
                                                select(c("upper_mort_num")))),
               col = "darkblue", size = 1.2, 
               arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
  geom_segment(aes(x = 2017, y = unlist(GBD_TBHIV_mort_num_df%>%
                                          filter(year == 2017)%>%
                                          select(c("lower_mort_num"))),
                   xend = 2017, yend = unlist(GBD_TBHIV_mort_num_df%>%
                                                filter(year == 2017)%>%
                                                select(c("upper_mort_num")))),
               col = "darkblue", size = 1.2, 
               arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
  theme(text = element_text(size=22, family="Times New Roman"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position=c(0.2, 0.7), legend.title = element_blank(), 
        legend.text=element_text(size=22),
        plot.title = element_text(hjust = .5, size=20),
        #legend.background = element_rect(colour = "lightgrey"),
        #legend.box.background = element_rect(colour = "black")
        )+
  scale_x_continuous(name = 'year y', breaks=c(seq(from = 1990, to = 2017, by = 3)))+
  scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3),
                     limits = c(0, 501000))+
  ylab('Number of TB/HIV-related deaths')+
  annotate("text", x = 1997, 
           y = 480000, 
           label = 'Number of TB and HIV-related\ndeaths in South Africa',
           size = 10,
           family="Times New Roman")

setwd(paste0(outdir, '/', "ALL South Africa"))

png('n_disease_mort.png', width = 800, height = 500)
print(SA_num_TBHIV_mort_graph)
dev.off()
