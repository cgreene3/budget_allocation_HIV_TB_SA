#Estimating background mortality rate from GBD
#Uses GBD 2019
#(Total Mort - HIV Mort - TB mort)/pop estimate = base mort
#regional SA population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to KZN_south_africa for here to work
indir <- paste0(here(),'/param_files/calculated_param_gen/raw_input_data/GBD')
outdir <- paste0(here(),'/param_files/input_parameters')
graph_outdir<-paste0(here(),'/param_files/dynamic_param_graphs')

setwd(indir)
#read in pop estimates
pop_df<-read.csv('pop_estimates_15_59.csv')
region_id_ref_df<-read.csv('region_id_ref.csv')

base_mort_df_all<-data.frame()

for (i in unique(pop_df$region_name)){
  #read in all cause mort num and summarise
  all_cause_mort_df<-read.csv('all_cause_mort_num_rate.csv')%>%
    filter(location == i)%>%
    filter(metric == 'Number')%>%
    group_by(year, sex)%>%
    summarise(all_cause_val = sum(val))

  disease_mort_df<-read.csv('disease_mort_num.csv')%>%
    filter(location == i)%>%
    group_by(year, sex)%>%
    summarise(disease_val = sum(val))

  base_mort_df<-pop_df%>%
    filter(region_name == i)%>%
    left_join(all_cause_mort_df, by = c('year', 'sex'))%>%
    left_join(disease_mort_df, by = c('year', 'sex'))%>%
    mutate(val = ((all_cause_val-disease_val)/expected_total_pop))%>%
    select('year', 'sex', 'val')

  base_mort_df<-base_mort_df%>%
    filter(year <= 2017)

  base_mort_df_2018_mort_vals_female<-base_mort_df%>%filter(year == 2017, sex == "Female")
  base_mort_df_2018_mort_vals_male<-base_mort_df%>%filter(year == 2017, sex == "Male")


  base_mort_df_2018_2028<-data.frame(year = rep(2018:2028, times = 2),
                                     sex = rep(c("Female", "Male"),
                                               each = length(2018:2028)),
                                     val = rep(c(base_mort_df_2018_mort_vals_female$val,
                                                 base_mort_df_2018_mort_vals_male$val),
                                               each = length(2018:2028)))
  
  base_mort_df<-rbind(base_mort_df, base_mort_df_2018_2028)
  base_mort_df$max = base_mort_df$val*1.25
  base_mort_df$min = base_mort_df$val*.75
  base_mort_df$region_name = rep(i, times = nrow(base_mort_df))
  
  if(nrow(base_mort_df_all) == 0){
    base_mort_df_all<-base_mort_df
  } else{
    base_mort_df_all<-rbind(base_mort_df_all, base_mort_df)
  }
}
  

setwd(outdir)
write.csv(base_mort_df_all, 'base_mort_df.csv', row.names = FALSE)

#make plots for appendix


base_mort_df_graph_male<-base_mort_df_all%>%
  filter(sex == 'Male')%>%
  select(-c('sex'))%>%
  left_join(region_id_ref_df, by = c("region_name"))%>%
  mutate(region_name = paste0(region_name, " (i = ", region_id, ")"))%>%
  arrange(region_id)%>%
  select(-c(region_id))

#without max/min
baseline_mort_plot_male<-ggplot(data=base_mort_df_graph_male%>%filter(year <= 2017), 
       aes(x=year, y=val, group=region_name, colour = region_name)) +
  geom_line(size = 0.5)+
  geom_point(size = 1)+
  #geom_line(color="grey", size = 1)+
  #geom_point(color="blue", size = 2)+
  scale_x_continuous(name = 'y', breaks=seq(from = 1990, to = 2017, by = 4))+
  scale_y_continuous(name = substitute(paste(mu, "_VAL"["1"]^{"i"}, (y))), 
                     limits = c(0, .015), breaks=(seq(0, .015, .005)))+
  theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=3, byrow=TRUE))

base_mort_df_graph_female<-base_mort_df_all%>%
  filter(sex == 'Female')%>%
  select(-c('sex'))%>%
  left_join(region_id_ref_df, by = c("region_name"))%>%
  mutate(region_name = paste0(region_name, " (i = ", region_id, ")"))%>%
  arrange(region_id)%>%
  select(-c(region_id))

#without max/min
baseline_mort_plot_female<-ggplot(data=base_mort_df_graph_female%>%filter(year <= 2017), 
                                aes(x=year, y=val, group=region_name, colour = region_name)) +
  geom_line(size = 0.5)+
  geom_point(size = 1)+
  #geom_line(color="grey", size = 1)+
  #geom_point(color="blue", size = 2)+
  scale_x_continuous(name = 'y', breaks=seq(from = 1990, to = 2017, by = 4))+
  scale_y_continuous(name = substitute(paste(mu, "_VAL"["1"]^{"i"}, (y))), 
                     limits = c(0, .015), breaks=(seq(0, .015, .005)))+
  theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=3, byrow=TRUE))

setwd(graph_outdir)

png("baseline_mort_plot_male.png", width = 560, height = 500)
print(baseline_mort_plot_male)
dev.off()

png("baseline_mort_plot_female.png", width = 560, height = 500)
print(baseline_mort_plot_female)
dev.off()

