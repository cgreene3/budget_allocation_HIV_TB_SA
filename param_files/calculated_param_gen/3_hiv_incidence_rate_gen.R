#Estimating HIV incidence rate from GBD
#Uses GBD 2019
#(num HIV incidence)/pop estimate = incidence rate by gender
#SA regional population estimates

#scales up linearly between 1980-1990 
#(when hiv incidence estimates don't exist from GBD)
#GBD from 1990-2017
#scales down according to reduced incidence projections from DO ART HPV-HIV model

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to resource_allocation_HIV_TB for here to work
indir_GBD <- paste0(here(),'/param_files/calculated_param_gen/raw_input_data/GBD')
indir_DRIVE<- paste0(here(), '/param_files/calculated_param_gen/raw_input_data/DRIVE')
outdir <- paste0(here(),'/param_files/input_parameters')
graph_outdir<-paste0(here(),'/param_files/dynamic_param_graphs')

setwd(indir_GBD)

#read in pop estimates
pop_df<-read.csv('pop_estimates_15_59.csv')

#read in incidence rate df
hiv_inc_df<-read.csv('hiv_inc_num.csv')%>%
  group_by(year, sex, region_name)%>%
  summarise(val = sum(val))%>%
  left_join(pop_df, by = c('year', 'sex', 'region_name'))%>%
  mutate(val = (val/expected_total_pop))%>%
  select('year', 'sex', 'val', 'region_name')%>%
  ungroup()

#####function to scale up from 1980 to 1989#####
lapply(unique(hiv_inc_df$region_name), function(i){
  male_1990_inc<-hiv_inc_df%>%
    filter(year == 1990,
           sex == 'Male',
           region_name == i)
  
  female_1990_inc<-hiv_inc_df%>%
    filter(year == 1990,
           sex == 'Female', 
           region_name == i)
  
  male_rate<-male_1990_inc$val/11
  female_rate<-female_1990_inc$val/11
  
  est_yrs<-1980:1989
  
  male_vals<-rep(0, times = length(est_yrs))
  female_vals<-rep(0, times = length(est_yrs))
  
  for (yr in 1:length(est_yrs)){
    male_vals[yr]<-male_rate*yr
    female_vals[yr]<-female_rate*yr
  }
  
  male_1980_1989_df<-data.frame(year = est_yrs,
                                sex = rep('Male', times = length(1980:1989)))
  male_1980_1989_df$val = male_vals
  
  female_1980_1989_df<-data.frame(year = est_yrs,
                                  sex = rep('Female', times = length(1980:1989)))
  female_1980_1989_df$val = female_vals
  
  hiv_rate_1980_1989_df<-rbind(male_1980_1989_df, female_1980_1989_df)
  hiv_rate_1980_1989_df$region_name<-rep(i, times = nrow(hiv_rate_1980_1989_df))
  
  
  hiv_inc_df<<-rbind.data.frame(as.data.frame(hiv_inc_df),
                                as.data.frame(hiv_rate_1980_1989_df))
  
})

hiv_inc_df$program_id <- 1

hiv_inc_df<-hiv_inc_df%>%
  filter(year <= 2017)

#####projections from 2018 to 2028 (evaluation years)####
lapply(unique(hiv_inc_df$region_name), function(i){
  GBD_2017_female_val<-hiv_inc_df%>%
    filter(sex == 'Female',
           year == 2017,
           region_name == i)

  GBD_2017_male_val<-hiv_inc_df%>%
    filter(sex == 'Male',
           year == 2017,
           region_name == i)
  
  setwd(indir_DRIVE)
  intervention_period_change<-read.csv('intervention_hiv_incidence.csv')
  
  incidence_do_art_df_SOC<-intervention_period_change%>%
    filter(year > 2017)%>%
    select(c('year', 'sex', 'Percent_Decrease_SOC'))%>%
    arrange(sex, year)%>%
    group_by(sex)%>%
    mutate(cum_perc_decrease = if_else(year == 2018, 
                                       1+Percent_Decrease_SOC,
                                       cumprod(1+Percent_Decrease_SOC)))%>%
    mutate(val = if_else(sex == "Female", 
                         (cum_perc_decrease*GBD_2017_female_val$val),
                         (cum_perc_decrease*GBD_2017_male_val$val)))%>%
    select(c('year', 'sex', 'val'))%>%
    ungroup()

  incidence_do_art_df_SOC$program_id = 1
  incidence_do_art_df_SOC<-as.data.frame(incidence_do_art_df_SOC)

  incidence_do_art_df_DO_ART<-intervention_period_change%>%
    filter(year > 2017)%>%
    select(c('year', 'sex', 'Percent_Decrease_DO_ART'))%>%
    arrange(sex, year)%>%
    group_by(sex)%>%
    mutate(cum_perc_decrease = if_else(year == 2018, 
                                       1+Percent_Decrease_DO_ART,
                                       cumprod(1+Percent_Decrease_DO_ART)))%>%
    mutate(val = if_else(sex == "Female", 
                         (cum_perc_decrease*GBD_2017_female_val$val),
                         (cum_perc_decrease*GBD_2017_male_val$val)))%>%
    select(c('year', 'sex', 'val'))%>%
    ungroup()

  incidence_do_art_df_DO_ART<-as.data.frame(incidence_do_art_df_DO_ART)
  
  incidence_do_art_df_SOC$region_name = i
  incidence_do_art_df_DO_ART$region_name = i

  incidence_do_art_df_DO_ART_2<-incidence_do_art_df_DO_ART
  incidence_do_art_df_DO_ART_2$program_id = 2

  incidence_do_art_df_DO_ART_3<-incidence_do_art_df_DO_ART
  incidence_do_art_df_DO_ART_3$program_id = 3

  hiv_inc_df<<-rbind(hiv_inc_df, incidence_do_art_df_SOC, incidence_do_art_df_DO_ART_2,
                     incidence_do_art_df_DO_ART_3)

})


hiv_inc_df$max = hiv_inc_df$val*1.25
hiv_inc_df$min = hiv_inc_df$val*.75

hiv_inc_df$program_id<-as.factor(hiv_inc_df$program_id)

setwd(outdir)
write.csv(hiv_inc_df, 'hiv_inc_df.csv', row.names = FALSE)


for(i in 1:length(unique(hiv_inc_df$region_name))){
  
  region_name_temp<-unique(hiv_inc_df$region_name)[i]
  
  #switch order of northern cape and north west
  i<-if_else(i == as.integer(8), as.integer(7), 
             if_else(i == as.integer(7), as.integer(8), as.integer(i)))
  
  #########graphing just val
  hiv_inc_df_graph_male<-hiv_inc_df%>%
    filter(sex == 'Male')%>%
    select(-c('sex'))%>%
    filter(program_id != 3,
           region_name == region_name_temp)%>%
    mutate(program_name = if_else(program_id == 1, "Program 1 (p = 1)", "Program 2 and 3 (p = 2, 3)"))%>%
    filter(year <= 2027)
  
  hiv_inc_plot_male<-ggplot(data=hiv_inc_df_graph_male) +
    geom_line(aes(x=year, y=val, group = factor(program_id)), color="grey")+
    geom_point(aes(x=year, y=val, colour=factor(program_name), 
                   shape = factor(program_name)), size = 1.5)+
    scale_color_manual(values = c("blue", "red"))+
    geom_vline(xintercept = 2017, linetype="dashed", 
               color = "darkgrey", size=1)+
    annotate("text", x=2009.5, y=.045, label= "calibration period", size = 5)+
    annotate("text", x=2023.5, y=.045, label= "evaluation period", size = 5)+
    #labs(title="HIV incidence rate, Males")+
    scale_x_continuous(name = 'y', breaks=seq(from = 1980, to = 2027, by = 4))+
    scale_y_continuous(name = substitute(paste(eta, "_VAL"["1,2,1"]^{"i,p"}, "(y)  for i = ", ii), list(ii=i)),
                         #paste(bquote(eta["1,2,1"]^{"i,p"}), "i = ", i), 
                       limits = c(0, .05), breaks=(seq(0, .05, .01)))+
    theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),
          plot.title = element_text(hjust = 0.5), legend.position = "top",
          legend.key = element_rect(fill = "grey94"),
          legend.background = element_rect(colour = "lightgrey"),
          legend.box.background = element_rect(colour = "black"))+
    scale_shape_manual(values=c(8, 2))
  
  hiv_inc_df_graph_female<-hiv_inc_df%>%
    filter(sex == 'Female')%>%
    select(-c('sex'))%>%
    filter(program_id != 3,
           region_name == region_name_temp)%>%
    mutate(program_name = if_else(program_id == 1, "Program 1 (p = 1)", "Program 2 and 3 (p = 2, 3)"))%>%
    filter(year <= 2027)
  
  hiv_inc_plot_female<-ggplot(data=hiv_inc_df_graph_female) +
    geom_line(aes(x=year, y=val, group = factor(program_id)), color="grey")+
    geom_point(aes(x=year, y=val, colour=factor(program_name), 
                   shape = factor(program_name)), size = 1.3)+
    scale_color_manual(values = c("blue", "red"))+
    geom_vline(xintercept = 2017, linetype="dashed", 
               color = "darkgrey", size=1)+
    annotate("text", x=2009.5, y=.045, label= "calibration period", size = 5)+
    annotate("text", x=2023.5, y=.045, label= "evaluation period", size = 5)+
    #labs(title="HIV incidence rate, Females")+
    scale_x_continuous(name = 'y', breaks=seq(from = 1980, to = 2027, by = 4))+
    scale_y_continuous(name = substitute(paste(eta, "_VAL"["1,2,2"]^{"i,p"}, "(y)  for i = ", ii), list(ii=i)), limits = c(0, .05), breaks=(seq(0, .05, .01)))+
    theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),
          plot.title = element_text(hjust = 0.5), legend.position = "top",
          legend.key = element_rect(fill = "grey94"),
          legend.background = element_rect(colour = "lightgrey"),
          legend.box.background = element_rect(colour = "black"))+
    scale_shape_manual(values=c(8, 2))
  
  setwd(graph_outdir)
  png(paste0("hiv_inc_plot_female", region_name_temp, ".png"), width = 600, height = 440)
  print(hiv_inc_plot_female)
  dev.off()
  
  
  png(paste0("hiv_inc_plot_male", region_name_temp, ".png"), width =600, height = 440)
  print(hiv_inc_plot_male)
  dev.off()
  
}
