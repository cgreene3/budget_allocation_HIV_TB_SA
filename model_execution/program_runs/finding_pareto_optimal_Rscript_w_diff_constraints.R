#clean workspace
rm(list = ls())
gc()

library(here)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(scales)

outdir<-paste0(here(), '/results/pareto_optimal_tables')

#different code to find pareto optimal solutions
###calculate objective functions for each solution set####
n_programs <-2
n_regions<-9

all_solution_sets_df<-expand.grid(rep(list(1:n_programs), n_regions))
colnames(all_solution_sets_df)<-1:n_regions
all_solution_sets_df$solution_id <-1:nrow(all_solution_sets_df)

#combine solution values
setwd(paste0(here(), '/results/program_runs/function_values_by_solution_set'))
all_files_in_dir <- list.files(pattern="function_df*")

function_df_all<-read.csv(all_files_in_dir[1])

lapply(all_files_in_dir[2:length(all_files_in_dir)], function(file){
  function_df_temp<-read.csv(file)
  function_df_all<<-rbind(function_df_all, function_df_temp)
})



###critical val function
critical_val_fun<-function(mean_maintained,
                           mean_candidate,
                           sd_maintained,
                           sd_candidate,
                           conf_level){
  
  obj_z_score<-qnorm(1-((1-conf_level)/2))#two tailed
  
  diff_means<-mean_maintained-mean_candidate
  
  critical_val_temp<-diff_means+
    (obj_z_score*(sqrt(sd_maintained+sd_candidate)/(n_samples^(1/2))))
  
  return(critical_val_temp)
}

pareto_optimal_algorithm_fun<-function(budget,
                                       obj1_conf,
                                       obj2_conf,
                                       cost_conf){
  
  #get cut off for budget
  budget_z_score<-if_else(cost_conf == 0, 0, qnorm(cost_conf))
  
  solution_sets_after_budget<-function_df_all%>%
    filter(variable == 'total_cost_fun')%>%
    mutate(conf_interval_upper = func_mean + (func_sd/sqrt(n_samples))*budget_z_score)%>%
    filter(conf_interval_upper <= budget)
  
  print('n solutions after budget')
  print(nrow(solution_sets_after_budget))
  
  #filter solution sets that do not meet the budget constraint
  function_df_after_budget_filter<-function_df_all%>%
    filter(solution_id %in% solution_sets_after_budget$solution_id)%>%
    group_by(variable)%>%
    mutate(rank_mean = rank(func_mean))
  
  rm(solution_sets_after_budget)
  
  #sort remaining solution sets by rank of objective
  function_df_w_ranks<-dcast(function_df_after_budget_filter%>%
                               select(c('solution_id', 'variable', 'rank_mean')),
                             solution_id ~ variable,
                             fun.aggregate = mean,
                             value.var = 'rank_mean')
    
  function_df_w_ranks<-function_df_w_ranks%>%
    arrange(TBHIV_mort_diff_obj)%>%
    arrange(total_TBHIVmort)%>%
    select(c('solution_id', 'TBHIV_mort_diff_obj', 'total_TBHIVmort'))
  
  maintained_solutions<-c(function_df_w_ranks$solution_id[1])
               
  
  for(candidate_solutions_temp in 2:nrow(function_df_w_ranks)){
    
    candidate_solution_id_temp<-function_df_w_ranks$solution_id[candidate_solutions_temp]
    
    maintained_obj_ranks_df<-function_df_w_ranks%>%
      filter(solution_id %in% maintained_solutions)
    
    
    candidate_solution_df<-function_df_after_budget_filter%>%
      filter(solution_id == candidate_solution_id_temp)
    
    candidate_solution_obj1_mean<-(candidate_solution_df$func_mean[candidate_solution_df$variable == 'total_TBHIVmort'])
    candidate_solution_obj2_mean<-(candidate_solution_df$func_mean[candidate_solution_df$variable == 'TBHIV_mort_diff_obj'])

    candidate_solution_obj1_var<-(candidate_solution_df$func_sd[candidate_solution_df$variable == 'total_TBHIVmort'])^2
    candidate_solution_obj2_var<-(candidate_solution_df$func_sd[candidate_solution_df$variable == 'TBHIV_mort_diff_obj'])^2
    
    #allows to change if based on mean or includes uncertainties
    dominated_or_not_dominated<-rep(0, times = nrow(maintained_obj_ranks_df))
    #0 non-dominated
    #1 dominated with confidence
    #2 dominated without confidence
    
    for(maintained_solutions_temp in 1:nrow(maintained_obj_ranks_df)){
      #first check based on mean
      if(maintained_obj_ranks_df$TBHIV_mort_diff_obj[maintained_solutions_temp] >
         function_df_w_ranks$TBHIV_mort_diff_obj[candidate_solutions_temp]){
        
        dominated_or_not_dominated[maintained_solutions_temp]<-0 #not dominated

      } else {
        
        #check if dominated with confidence
        maintained_solution_id_temp<-maintained_obj_ranks_df$solution_id[maintained_solutions_temp]
        
        maintained_solution_df<-function_df_after_budget_filter%>%
          filter(solution_id == maintained_solution_id_temp)
        
        maintained_solution_obj1_mean<-(maintained_solution_df$func_mean[maintained_solution_df$variable == 'total_TBHIVmort'])
        maintained_solution_obj2_mean<-(maintained_solution_df$func_mean[maintained_solution_df$variable == 'TBHIV_mort_diff_obj'])

        maintained_solution_obj1_var<-(maintained_solution_df$func_sd[maintained_solution_df$variable == 'total_TBHIVmort'])^2
        maintained_solution_obj2_var<-(maintained_solution_df$func_sd[maintained_solution_df$variable == 'TBHIV_mort_diff_obj'])^2

        critical_val_obj1<-critical_val_fun(maintained_solution_obj1_mean,
                                            candidate_solution_obj1_mean,
                                            maintained_solution_obj1_var,
                                            candidate_solution_obj1_var,
                                            obj1_conf)
        critical_val_obj2<-critical_val_fun(maintained_solution_obj2_mean,
                                            candidate_solution_obj2_mean,
                                            maintained_solution_obj2_var,
                                            candidate_solution_obj2_var,
                                            obj2_conf)
        if(critical_val_obj1 <= 0 & 
           critical_val_obj2 <= 0){
          dominated_or_not_dominated[maintained_solutions_temp]<-1 #dominated with confidence
        } else {
          dominated_or_not_dominated[maintained_solutions_temp]<-2 #dominated without confidence
        }
      }
    }
    if(!(1 %in% dominated_or_not_dominated)){
      #if better in one of the other objectives, then non-dominated
      #if(!(2 %in% dominated_or_not_dominated)){
        maintained_solutions<-c(maintained_solutions, function_df_w_ranks$solution_id[candidate_solutions_temp])
      #}
    }
  }
  return(maintained_solutions)
}
  


budget_val<-15000000000
obj1_conf_val_mean<-0
obj2_conf_val_mean<-0
cost_conf_val_mean<-0

obj1_conf_val_w_uncertainty<-0.95
obj2_conf_val_w_uncertainty<-0.95
cost_conf_val_w_uncertainty<-0.95

n_samples<-816

#no budget uncertainty
maintained_solution_all_means<-pareto_optimal_algorithm_fun(budget = budget_val,
                                                            obj1_conf = obj1_conf_val_mean,
                                                            obj2_conf = obj2_conf_val_mean,
                                                            cost_conf = cost_conf_val_mean)
 
maintained_solution_obj_uncertainty_only<-pareto_optimal_algorithm_fun(budget = budget_val,
                                                                obj1_conf = obj1_conf_val_w_uncertainty,
                                                                obj2_conf = obj2_conf_val_w_uncertainty,
                                                                cost_conf = cost_conf_val_mean)

#budget uncertainty
maintained_solution_cost_uncertainty_only<-pareto_optimal_algorithm_fun(budget = budget_val,
                                                                     obj1_conf = obj1_conf_val_mean,
                                                                     obj2_conf = obj2_conf_val_mean,
                                                                     cost_conf = cost_conf_val_w_uncertainty)

maintained_solution_cost_obj_uncertainty<-pareto_optimal_algorithm_fun(budget = budget_val,
                                                                         obj1_conf = obj1_conf_val_w_uncertainty,
                                                                         obj2_conf = obj2_conf_val_w_uncertainty,
                                                                         cost_conf = cost_conf_val_w_uncertainty)
 
print(maintained_solution_all_means)
print(maintained_solution_obj_uncertainty_only)
print(maintained_solution_cost_uncertainty_only)
print(maintained_solution_cost_obj_uncertainty)

budget_z_score<-qnorm(cost_conf_val_w_uncertainty)
obj_z_score<-qnorm(1-(1-obj1_conf_val_w_uncertainty)/2)

####adding in highest and lowest solutions 
all_community_sol<-512 #all community-based
all_facility_sol<-1 #all faciliaty based


###make conf intervals

#not considering cost (not considered in paper)
# pareto_solutions_no_cost_uncertainty_df<-data.frame(solution_id = maintained_solution_obj_uncertainty_only)%>%
#   mutate(obj_uncertainty_only = if_else(solution_id %in% maintained_solution_all_means, 'no', 'yes'))%>%
#   full_join(function_df_all%>%filter(solution_id %in% maintained_solution_obj_uncertainty_only), 
#             by = c('solution_id'))%>%
#   mutate(ll = if_else(variable == 'total_cost_fun', 
#                       0,
#                       func_mean - qnorm(obj1_conf_val_w_uncertainty)*func_sd/sqrt(n_samples)),
#          ul = if_else(variable == 'total_cost_fun', 
#                       func_mean + budget_z_score*(func_sd/sqrt(n_samples)),
#                       func_mean + obj_z_score*(func_sd/sqrt(n_samples))))%>%
#   select(-c('func_sd'))%>%
#   filter(variable %in% c('TBHIV_mort_diff_obj', 'total_TBHIVmort', 'total_cost_fun'))%>%
#   mutate(conf_interval_txt = if_else(variable =='TBHIV_mort_diff_obj',
#                                      paste0(round(func_mean), ' [', round(ll),
#                                             ', ', round(ul), ']'),
#                                      if_else(variable=='total_TBHIVmort', 
#                                              paste0(sprintf('%.3f', func_mean/1000000), 
#                                                     'M [', sprintf('%.3f', ll/1000000),
#                                                     'M, ', sprintf('%.3f', ul/1000000), 'M]'),
#                                              paste0(sprintf('%.3f', func_mean/1000000000), 
#                                                     'B [', sprintf('%.3f', ll/1000000000),
#                                                     'B, ', sprintf('%.3f', ul/1000000000), 'B]'))))
# 
# pareto_solutions_tab_df_no_cost_uncertainty<-pareto_solutions_no_cost_uncertainty_df%>%
#   select(c('solution_id', 'variable', 'conf_interval_txt'))%>%
#   dcast(., solution_id~variable)%>%
#   left_join(all_solution_sets_df, by = c('solution_id'))%>%
#   mutate(obj_uncertainty_only = if_else(solution_id %in% maintained_solution_all_means, 'no', 'yes'))
#   

##considering cost

#used for graph of pareto solutions with all community-based and all facility based
pareto_solutions_graph_with_extremes_df<-data.frame(solution_id = c(maintained_solution_cost_obj_uncertainty,
                                                                       all_community_sol,
                                                                       all_facility_sol))%>%
  mutate(obj_uncertainty_only = if_else(solution_id %in% maintained_solution_cost_uncertainty_only, 'no', 'yes'))%>%
  full_join(function_df_all%>%filter(solution_id %in% c(maintained_solution_cost_obj_uncertainty,
                                                        all_community_sol,
                                                        all_facility_sol)), 
            by = c('solution_id'))%>%
  mutate(ll = if_else(variable == 'total_cost_fun', 
                      0,
                      func_mean - obj_z_score*func_sd/sqrt(n_samples)),
         ul = if_else(variable == 'total_cost_fun', 
                      func_mean + budget_z_score*(func_sd/sqrt(n_samples)),
                      func_mean + obj_z_score*(func_sd/sqrt(n_samples))))%>%
  select(-c('func_sd'))%>%
  filter(variable %in% c('TBHIV_mort_diff_obj', 'total_TBHIVmort', 'total_cost_fun'))
  
#used for pareto table
pareto_solutions_tab_with_extremes_df<-pareto_solutions_graph_with_extremes_df%>%
  mutate(conf_interval_txt = if_else(variable =='TBHIV_mort_diff_obj',
                                     paste0(round(func_mean), ' [', round(ll),
                                            ', ', round(ul), ']'),
                                     if_else(variable=='total_TBHIVmort', 
                                             paste0(sprintf('%.3f', func_mean/1000000), 
                                                    ' [', sprintf('%.3f', ll/1000000),
                                                    ', ', sprintf('%.3f', ul/1000000), ']'),
                                             paste0(sprintf('%.3f', func_mean/1000000000), 
                                                    ' [', sprintf('%.3f', ll/1000000000),
                                                    ', ', sprintf('%.3f', ul/1000000000), ']'))))%>%
  left_join(all_solution_sets_df, by = c('solution_id'))
                                         
                                         
##if table includes with and without cost ucertainty
# pareto_solutions_tab_df_with_cost_uncertainty<-pareto_solutions_tab_df_with_cost_uncertainty%>%
#   select(c('solution_id', 'variable', 'conf_interval_txt'))%>%
#   dcast(., solution_id~variable)%>%
#   left_join(all_solution_sets_df, by = c('solution_id'))%>%
#   mutate(obj_uncertainty_only = if_else(solution_id %in% maintained_solution_cost_uncertainty_only, 'no', 'yes'))

###filter extremes
pareto_solutions_graph_df<-pareto_solutions_graph_with_extremes_df%>%
  filter(!(solution_id %in% c(all_community_sol,
                             all_facility_sol)))


pareto_solutions_tab_df<-pareto_solutions_tab_with_extremes_df%>%
  filter(!(solution_id %in% c(all_community_sol,
                              all_facility_sol)))


setwd(outdir)
write.csv(pareto_solutions_tab_df,
          'pareto_solutions_tab_df.csv',
          row.names = FALSE)
# write.csv(pareto_solutions_tab_df_with_cost_uncertainty, 
#           'pareto_optimal_with_cost_uncertainty_considerations.csv',
#           row.names = FALSE)
# write.csv(pareto_solutions_tab_df_no_cost_uncertainty, 
#           'pareto_optimal_cost_means.csv',
#           row.names = FALSE)


##########GRAPHS#####

#### Pareto graph with extremes #####
pareto_solutions_graph_with_extremes_df2<-pareto_solutions_graph_with_extremes_df%>%
  select(c('solution_id', 'variable', 'func_mean', 'ul', 'll'))

pareto_solutions_graph_with_extremes_df2_yaxis<-pareto_solutions_graph_with_extremes_df2%>%
  filter(variable == 'total_TBHIVmort')%>%
  select(-c('variable'))

colnames(pareto_solutions_graph_with_extremes_df2_yaxis)<-
  c('solution_id', 'func_mean_y', 'ul_y', 'll_y')

pareto_solutions_graph_with_extremes_df2_xaxis<-pareto_solutions_graph_with_extremes_df2%>%
  filter(variable == 'TBHIV_mort_diff_obj')%>%
  select(-c('variable'))

colnames(pareto_solutions_graph_with_extremes_df2_xaxis)<-
  c('solution_id', 'func_mean_x', 'ul_x', 'll_x')

pareto_graph_with_extremes<-pareto_solutions_graph_with_extremes_df2_yaxis%>%
  left_join(pareto_solutions_graph_with_extremes_df2_xaxis, 
            by = c('solution_id'))%>%
  mutate(obj_uncertainty_only = if_else(solution_id %in% maintained_solution_cost_uncertainty_only, 'no', 'yes'))%>%
  arrange(func_mean_y)%>%
  arrange(obj_uncertainty_only)%>%
  mutate(func_mean_y = func_mean_y/1000000,
         ul_y = ul_y/1000000,
         ll_y = ll_y/1000000)

pareto_graph_with_extremes$`Solution ID`<-1:nrow(pareto_graph_with_extremes)
pareto_graph_with_extremes$color_ref<-if_else(pareto_graph_with_extremes$`Solution ID` %in% 1:5, 'blue', 
                                              if_else(pareto_graph_with_extremes$solution_id %in% c(all_community_sol, 
                                                                                                    all_facility_sol), 'red', 
                                                      'black'))


ggplot(pareto_graph_with_extremes, aes(x = func_mean_x, y = func_mean_y, 
                         color = color_ref)
       )+
  geom_errorbarh(aes(xmin = ll_x, xmax = ul_x))+
  geom_errorbar(aes(ymin = ll_y, ymax = ul_y), width = 2)+
  geom_point()+
  #geom_label(aes(label = `Solution ID`))+
  scale_color_manual(values = c('black', 'blue', 'red'))+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Objective 1: TB and HIV-related\nMortality Rate Inequity")+
  ylab("Objective 2: Total TB and HIV-related Deaths\nin Millions")+
  ylim(c(0, 1.5))+
  xlim(c(0,400))
dev.off()  
#stat_ellipse(data = ellipse_graph_df, aes(x = TBHIV_mort_diff_obj,
#                                   y = total_TBHIVmort,
#                                   group = solution_id))


#### Pareto graph without extremes #####

pareto_solutions_graph_df2<-pareto_solutions_graph_df%>%
  select(c('solution_id', 'variable', 'func_mean', 'ul', 'll'))

pareto_solutions_tab_df_with_cost_uncertainty_graph_yaxis<-pareto_solutions_graph_df2%>%
  filter(variable == 'total_TBHIVmort')%>%
  select(-c('variable'))

colnames(pareto_solutions_tab_df_with_cost_uncertainty_graph_yaxis)<-
  c('solution_id', 'func_mean_y', 'ul_y', 'll_y')

pareto_solutions_tab_df_with_cost_uncertainty_graph_xaxis<-pareto_solutions_graph_df2%>%
  filter(variable == 'TBHIV_mort_diff_obj')%>%
  select(-c('variable'))

colnames(pareto_solutions_tab_df_with_cost_uncertainty_graph_xaxis)<-
  c('solution_id', 'func_mean_x', 'ul_x', 'll_x')

pareto_graph<-pareto_solutions_tab_df_with_cost_uncertainty_graph_yaxis%>%
  left_join(pareto_solutions_tab_df_with_cost_uncertainty_graph_xaxis, 
            by = c('solution_id'))%>%
  mutate(obj_uncertainty_only = if_else(solution_id %in% maintained_solution_cost_uncertainty_only, 'no', 'yes'))%>%
  arrange(func_mean_x)%>%
  arrange(obj_uncertainty_only)%>%
  mutate(func_mean_y = func_mean_y/1000000,
         ul_y = ul_y/1000000,
         ll_y = ll_y/1000000)

pareto_graph$`Solution ID`<-1:nrow(pareto_graph)
pareto_graph$color_ref<-if_else(pareto_graph$`Solution ID` %in% 1:5, 'blue', 'black')

ggplot(pareto_graph, aes(x = func_mean_x, y = func_mean_y, 
                         color = color_ref, 
                         label = solution_id)
       )+
  geom_errorbarh(aes(xmin = ll_x, xmax = ul_x))+
  geom_errorbar(aes(ymin = ll_y, ymax = ul_y), width = 2)+
  geom_point()+
  geom_label(aes(label = `Solution ID`))+
  scale_color_manual(values = c('black', 'blue'))+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Objective 1: TB- and HIV-related Mortality Rate Inequity")+
  ylab("Objective 2: Total TB and HIV-related Deaths\nin Millions")#+
  #ylim(c(1.175, 1.25))+
  #xlim(c(220, 310))
dev.off()  
#stat_ellipse(data = ellipse_graph_df, aes(x = TBHIV_mort_diff_obj,
#                                   y = total_TBHIVmort,
#                                   group = solution_id))
#ylim(c(1.15, 1.25))+
#xlim(c(200,300))

###with elipse###
# pareto_solutions_tab_df_with_cost_uncertainty_graph2<-melt(pareto_solutions_tab_df_with_cost_uncertainty_graph2,
#                                                            id = c('solution_id', 'variable'))
# 
# colnames(pareto_solutions_tab_df_with_cost_uncertainty_graph2)<-c('solution_id',
#                                                                   'variable1',
#                                                                   'variable',
#                                                                   'value')




# pareto_solutions_tab_df_with_cost_uncertainty_graph3<-dcast(pareto_solutions_tab_df_with_cost_uncertainty_graph2,
#                                                             solution_id+variable~variable1, mean)
# 
# 
# pareto_solutions_graph_ellipse_equity_obj<-pareto_solutions_tab_df_with_cost_uncertainty_graph3%>%
#   select(c('solution_id', 'variable', 'TBHIV_mort_diff_obj'))%>%
#   left_join(pareto_solutions_tab_df_with_cost_uncertainty_graph3%>%
#               filter(variable == 'func_mean')%>%
#               select(-c('TBHIV_mort_diff_obj', 'variable')), by = c('solution_id'))
# 
# 
# pareto_solutions_graph_ellipse_tot_obj<-pareto_solutions_tab_df_with_cost_uncertainty_graph3%>%
#   select(c('solution_id', 'variable', 'total_TBHIVmort'))%>%
#   left_join(pareto_solutions_tab_df_with_cost_uncertainty_graph3%>%
#               filter(variable == 'func_mean')%>%
#               select(-c('total_TBHIVmort', 'variable')), by = c('solution_id'))%>%
#   select(c(colnames(pareto_solutions_graph_ellipse_equity_obj)))
# 
# ellipse_graph_df<-rbind(pareto_solutions_graph_ellipse_equity_obj,
#                         pareto_solutions_graph_ellipse_tot_obj)

setwd(paste0(here(), '/results/pareto_optimal_graphs'))
png('pareto_graph.png')#, width = 700, height = 500)
ggplot(pareto_graph, aes(x = func_mean_x, y = func_mean_y, 
                         color = color_ref, 
                         label = solution_id))+
  geom_errorbarh(aes(xmin = ll_x, xmax = ul_x))+
  geom_errorbar(aes(ymin = ll_y, ymax = ul_y), width = 2)+
  geom_point()+
  geom_label(aes(label = `Solution ID`))+
  scale_color_manual(values = c('black', 'blue'))+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Objective 1: TB and HIV-related Mortality Rate Inequity")+
  ylab("Objective 2: Total TB and HIV-related Deaths\nin Millions")#+
dev.off()  
#stat_ellipse(data = ellipse_graph_df, aes(x = TBHIV_mort_diff_obj,
  #                                   y = total_TBHIVmort,
  #                                   group = solution_id))
#ylim(c(1.15, 1.25))+
  #xlim(c(200,300))








ellipse_graph_df$solution_id<-as.factor(ellipse_graph_df$solution_id)

point_graph_df<-ellipse_graph_df%>%
  filter(variable == 'func_mean')

ellipse_graph_df<-ellipse_graph_df%>%
  filter(variable != 'func_mean')


ggplot(pareto_solutions_tab_df_with_cost_uncertainty)+
  geom_point(data = point_graph_df,
             aes(x = TBHIV_mort_diff_obj, 
                 y = total_TBHIVmort))+
  geom_errorbar(data = pareto_solutions_tab_df_with_cost_uncertainty%>%
                  filter(variable == 'total_TBHIVmort'),
                aes(ymin = ll,
                    ymax = ul))


# ####OBJ 1 Graph ####
# pareto_solutions_no_cost_uncertainty_df_obj1_graph<-pareto_solutions_no_cost_uncertainty_df%>%
#   filter(variable == 'obj1_fun')%>%
#   select(-c(variable))%>%
#   mutate(labels_xaxis = paste0(#'Solution ID ', 
#     solution_id, ': ', 
#     round(func_mean), ' [', round(ll),
#     ', ', round(ul), ']'))%>%
#   arrange(-func_mean)%>%
#   mutate(ymin_ref = rank(-func_mean)-0.25,
#          ymax_ref = rank(-func_mean)+0.25)%>%
#   mutate(color_y_axis = if_else(obj_uncertainty_only == 'yes',
#                                 '#d95f02', 'black'))
# 
# ####OBJ 2 GRAPH####
# pareto_solutions_no_cost_uncertainty_df_obj2_graph<-pareto_solutions_no_cost_uncertainty_df%>%
#   filter(variable == 'obj2_fun')%>%
#   select(-c(variable))%>%
#   mutate(labels_xaxis = paste0(#'Solution ID ', 
#     solution_id, ': ', 
#     round(func_mean), ' [', round(ll),
#     ', ', round(ul), ']'))%>%
#   arrange(-func_mean)%>%
#   mutate(ymin_ref = rank(-func_mean)-0.25,
#          ymax_ref = rank(-func_mean)+0.25)%>%
#   mutate(color_y_axis = if_else(obj_uncertainty_only == 'yes',
#                                 '#d95f02', 'black'))
# 
# pareto_solutions_w_cost_uncertainty_df_obj2_graph<-pareto_solutions_w_cost_uncertainty_df%>%
#   filter(variable == 'obj2_fun')%>%
#   select(-c(variable))%>%
#   mutate(labels_xaxis = paste0(#'Solution ID ', 
#     solution_id, ': ', 
#     round(func_mean), ' [', round(ll),
#     ', ', round(ul), ']'))%>%
#   arrange(-func_mean)%>%
#   mutate(ymin_ref = rank(-func_mean)-0.25,
#          ymax_ref = rank(-func_mean)+0.25)%>%
#   mutate(color_y_axis = if_else(obj_uncertainty_only == 'yes',
#                                 '#d95f02', 'black'))
# 
# ####OBJ 3 Graph####
# pareto_solutions_no_cost_uncertainty_df_obj3_graph<-pareto_solutions_no_cost_uncertainty_df%>%
#   filter(variable == 'obj3_fun')%>%
#   select(-c(variable))%>%
#   mutate(labels_xaxis = paste0(#'Solution ID ', 
#                                solution_id, ': ', 
#                                sprintf('%.3f', func_mean/1000000), 
#                                'M [', sprintf('%.3f', ll/1000000),
#     'M, ', sprintf('%.3f', ul/1000000), 'M]'))%>%
#   arrange(-func_mean)%>%
#   mutate(ymin_ref = rank(-func_mean)-0.25,
#          ymax_ref = rank(-func_mean)+0.25)%>%
#   mutate(color_y_axis = if_else(obj_uncertainty_only == 'yes',
#                                 '#d95f02', 'black'))
# 
# pareto_solutions_w_cost_uncertainty_df_obj3_graph<-pareto_solutions_w_cost_uncertainty_df%>%
#   filter(variable == 'obj3_fun')%>%
#   select(-c(variable))%>%
#   mutate(labels_xaxis = paste0(#'Solution ID ', 
#     solution_id, ': ', 
#     sprintf('%.3f', func_mean/1000000), 
#     'M [', sprintf('%.3f', ll/1000000),
#     'M, ', sprintf('%.3f', ul/1000000), 'M]'))%>%
#   arrange(-func_mean)%>%
#   mutate(ymin_ref = rank(-func_mean)-0.25,
#          ymax_ref = rank(-func_mean)+0.25)%>%
#   mutate(color_y_axis = if_else(obj_uncertainty_only == 'yes',
#                                 '#d95f02', 'black'))
# 
# ###cost####
# pareto_solutions_no_cost_uncertainty_df_cost_func_graph<-pareto_solutions_no_cost_uncertainty_df%>%
#   filter(variable == 'total_cost_fun')%>%
#   select(-c(variable))%>%
#   mutate(labels_xaxis = paste0(#'Solution ID ', 
#     solution_id, ': ', 
#     sprintf('%.3f', func_mean/1000000000), 
#     'B [', sprintf('%.3f', ll/1000000000),
#     'B, ', sprintf('%.3f', ul/1000000000), 'B]'))%>%
#   arrange(-func_mean)%>%
#   mutate(ymin_ref = rank(-func_mean)-0.25,
#          ymax_ref = rank(-func_mean)+0.25)%>%
#   mutate(color_y_axis = if_else(obj_uncertainty_only == 'yes',
#                                 '#d95f02', 'black'))
# 
# pareto_solutions_w_cost_uncertainty_df_cost_func_graph<-pareto_solutions_w_cost_uncertainty_df%>%
#   filter(variable == 'total_cost_fun')%>%
#   select(-c(variable))%>%
#   mutate(labels_xaxis = paste0(#'Solution ID ', 
#     solution_id, ': ', 
#     sprintf('%.3f', func_mean/1000000000), 
#     'B [', sprintf('%.3f', ll/1000000000),
#     'B, ', sprintf('%.3f', ul/1000000000), 'B]'))%>%
#   arrange(-func_mean)%>%
#   mutate(ymin_ref = rank(-func_mean)-0.25,
#          ymax_ref = rank(-func_mean)+0.25)%>%
#   mutate(color_y_axis = if_else(obj_uncertainty_only == 'yes',
#                                 '#d95f02', 'black'))
# 
# ###write graphs to files
# 
# ##without cost uncertainty###
# setwd(paste0(here(), '/results/pareto_optimal_graphs/without_cost_uncertainty'))
# png('obj1_pareto_graph.png', width = 800*.65, height = 500)
# ggplot(data = pareto_solutions_no_cost_uncertainty_df_obj1_graph)+
#   geom_rect(mapping = aes(ymin = ymin_ref, 
#                           ymax = ymax_ref,
#                           xmin = ll,
#                           xmax = ul,
#                           fill = rank(func_mean)))+
#   geom_point(aes(y = (ymin_ref+ymax_ref)/2, x = func_mean), color = 'red')+
#   theme(text = element_text(size=20, family="Times New Roman"), 
#         #panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         legend.position="none", 
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         plot.title = element_text(hjust = .5, size=20),
#         legend.background = element_rect(colour = "lightgrey"),
#         legend.box.background = element_rect(colour = "black"),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(colour = pareto_solutions_no_cost_uncertainty_df_obj1_graph$color_y_axis))+
#   ylab('')+
#   scale_x_continuous(name = 'TB Incidence Equity')+
#   scale_y_continuous(breaks = c(1:length(pareto_solutions_no_cost_uncertainty_df_obj1_graph$solution_id)), 
#                      labels = pareto_solutions_no_cost_uncertainty_df_obj1_graph$labels_xaxis)
# dev.off()
# 
# png('obj2_pareto_graph.png', width = 800*0.65, height = 500)
# 
# ggplot(data = pareto_solutions_no_cost_uncertainty_df_obj2_graph)+
#   geom_rect(mapping = aes(ymin = ymin_ref, 
#                           ymax = ymax_ref,
#                           xmin = ll,
#                           xmax = ul,
#                           fill = rank(func_mean)))+
#   geom_point(aes(y = (ymin_ref+ymax_ref)/2, x = func_mean), color = 'red')+
#   theme(text = element_text(size=20, family="Times New Roman"), 
#         #panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         legend.position="none", 
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         plot.title = element_text(hjust = .5, size=20),
#         legend.background = element_rect(colour = "lightgrey"),
#         legend.box.background = element_rect(colour = "black"),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(colour = pareto_solutions_no_cost_uncertainty_df_obj2_graph$color_y_axis))+
#   ylab('')+
#   scale_x_continuous(name = 'TB/HIV Mortality Equity')+
#   scale_y_continuous(breaks = c(1:length(pareto_solutions_no_cost_uncertainty_df_obj2_graph$solution_id)), 
#                      labels = pareto_solutions_no_cost_uncertainty_df_obj2_graph$labels_xaxis)
# dev.off()
# 
# png('obj3_pareto_graph.png', 800*0.65, height = 500)
# 
# ggplot(data = pareto_solutions_no_cost_uncertainty_df_obj3_graph)+
#   geom_rect(mapping = aes(ymin = ymin_ref, 
#                           ymax = ymax_ref,
#                           xmin = ll,
#                           xmax = ul,
#                           fill = rank(func_mean)))+
#   geom_point(aes(y = (ymin_ref+ymax_ref)/2, x = func_mean, size = 8), color = 'red')+
#   theme(text = element_text(size=20, family="Times New Roman"), 
#         #panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         legend.position="none", 
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         plot.title = element_text(hjust = .5, size=20),
#         legend.background = element_rect(colour = "lightgrey"),
#         legend.box.background = element_rect(colour = "black"),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(colour = pareto_solutions_no_cost_uncertainty_df_obj3_graph$color_y_axis))+
#   ylab('')+
#   scale_x_continuous(name = 'Total of TB/HIV-related deaths', 
#                      labels = label_number(suffix = " M", scale = 1e-6),
#                      n.breaks = 3,
#                      #limits = c(1280000, 1500000)
#   )+
#   scale_y_continuous(breaks = c(1:length(pareto_solutions_no_cost_uncertainty_df_obj3_graph$solution_id)), 
#                      labels = pareto_solutions_no_cost_uncertainty_df_obj3_graph$labels_xaxis)
# 
# dev.off()
# 
# png('cost_pareto_graph.png', 800*0.65, height = 500)
# ggplot(data = pareto_solutions_no_cost_uncertainty_df_cost_func_graph)+
#   geom_rect(mapping = aes(ymin = ymin_ref, 
#                           ymax = ymax_ref,
#                           xmin = ll,
#                           xmax = ul,
#                           fill = rank(func_mean)))+
#   geom_point(aes(y = (ymin_ref+ymax_ref)/2, x = func_mean, size = 8), color = 'red')+
#   geom_vline(xintercept = budget_val)+
#   theme(text = element_text(size=20, family="Times New Roman"), 
#         #panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         legend.position="none", 
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         plot.title = element_text(hjust = .5, size=20),
#         legend.background = element_rect(colour = "lightgrey"),
#         legend.box.background = element_rect(colour = "black"),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(colour = pareto_solutions_no_cost_uncertainty_df_cost_func_graph$color_y_axis))+
#   ylab('')+
#   scale_x_continuous(name = 'Total Program Cost', 
#                      labels = label_number(suffix = " B", scale = 1e-9)
#   )+
#   scale_y_continuous(breaks = c(1:length(pareto_solutions_no_cost_uncertainty_df_cost_func_graph$solution_id)), 
#                      labels = pareto_solutions_no_cost_uncertainty_df_cost_func_graph$labels_xaxis)
# 
# dev.off()
# 
# ###with cost uncertainty
# setwd(paste0(here(), '/results/pareto_optimal_graphs/with_cost_uncertianty'))
# png('obj1_pareto_graph.png', width = 800*.65, height = 500)
# ggplot(data = pareto_solutions_w_cost_uncertainty_df_obj1_graph)+
#   geom_rect(mapping = aes(ymin = ymin_ref, 
#                           ymax = ymax_ref,
#                           xmin = ll,
#                           xmax = ul,
#                           fill = rank(func_mean)))+
#   geom_point(aes(y = (ymin_ref+ymax_ref)/2, x = func_mean), color = 'red')+
#   theme(text = element_text(size=20, family="Times New Roman"), 
#         #panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         legend.position="none", 
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         plot.title = element_text(hjust = .5, size=20),
#         legend.background = element_rect(colour = "lightgrey"),
#         legend.box.background = element_rect(colour = "black"),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(colour = pareto_solutions_w_cost_uncertainty_df_obj1_graph$color_y_axis))+
#   ylab('')+
#   scale_x_continuous(name = 'TB Incidence Equity')+
#   scale_y_continuous(breaks = c(1:length(pareto_solutions_w_cost_uncertainty_df_obj1_graph$solution_id)), 
#                      labels = pareto_solutions_w_cost_uncertainty_df_obj1_graph$labels_xaxis)
# dev.off()
# 
# png('obj2_pareto_graph.png', width = 800*.65, height = 500)
# 
# ggplot(data = pareto_solutions_w_cost_uncertainty_df_obj2_graph)+
#   geom_rect(mapping = aes(ymin = ymin_ref, 
#                           ymax = ymax_ref,
#                           xmin = ll,
#                           xmax = ul,
#                           fill = rank(func_mean)))+
#   geom_point(aes(y = (ymin_ref+ymax_ref)/2, x = func_mean), color = 'red')+
#   theme(text = element_text(size=20, family="Times New Roman"), 
#         #panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         legend.position="none", 
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         plot.title = element_text(hjust = .5, size=20),
#         legend.background = element_rect(colour = "lightgrey"),
#         legend.box.background = element_rect(colour = "black"),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(colour = pareto_solutions_w_cost_uncertainty_df_obj2_graph$color_y_axis))+
#   ylab('')+
#   scale_x_continuous(name = 'TB/HIV Mortality Equity')+
#   scale_y_continuous(breaks = c(1:length(pareto_solutions_w_cost_uncertainty_df_obj2_graph$solution_id)), 
#                      labels = pareto_solutions_w_cost_uncertainty_df_obj2_graph$labels_xaxis)
# dev.off()
# 
# png('obj3_pareto_graph.png', 
#     width = 800*.65, height = 500)
# 
# ggplot(data = pareto_solutions_w_cost_uncertainty_df_obj3_graph)+
#   geom_rect(mapping = aes(ymin = ymin_ref, 
#                           ymax = ymax_ref,
#                           xmin = ll,
#                           xmax = ul,
#                           fill = rank(func_mean)))+
#   geom_point(aes(y = (ymin_ref+ymax_ref)/2, x = func_mean), color = 'red')+
#   theme(text = element_text(size=20, family="Times New Roman"), 
#         #panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         legend.position="none", 
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         plot.title = element_text(hjust = .5, size=20),
#         legend.background = element_rect(colour = "lightgrey"),
#         legend.box.background = element_rect(colour = "black"),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(colour = pareto_solutions_w_cost_uncertainty_df_obj3_graph$color_y_axis))+
#   ylab('')+
#   scale_x_continuous(name = 'Total of TB/HIV-related deaths', 
#                      labels = label_number(suffix = " M", scale = 1e-6),
#                      n.breaks = 3,
#                      #limits = c(1280000, 1500000)
#   )+
#   
#   scale_y_continuous(breaks = c(1:length(pareto_solutions_w_cost_uncertainty_df_obj3_graph$solution_id)), 
#                      labels = pareto_solutions_w_cost_uncertainty_df_obj3_graph$labels_xaxis)
# 
# dev.off()
# 
# png('cost_pareto_graph.png', width = 800*.65, height = 500)
# ggplot(data = pareto_solutions_w_cost_uncertainty_df_cost_func_graph)+
#   geom_rect(mapping = aes(ymin = ymin_ref, 
#                           ymax = ymax_ref,
#                           xmin = ll,
#                           xmax = ul,
#                           fill = rank(func_mean)))+
#   geom_point(aes(y = (ymin_ref+ymax_ref)/2, x = func_mean), color = 'red')+
#   geom_vline(xintercept = budget_val)+
#   theme(text = element_text(size=20, family="Times New Roman"), 
#         #panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         legend.position="none", 
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         plot.title = element_text(hjust = .5, size=20),
#         legend.background = element_rect(colour = "lightgrey"),
#         legend.box.background = element_rect(colour = "black"),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(colour = pareto_solutions_w_cost_uncertainty_df_cost_func_graph$color_y_axis))+
#   ylab('')+
#   scale_x_continuous(name = 'Total Program Cost', 
#                      labels = label_number(suffix = " B", scale = 1e-9)
#   )+
#   scale_y_continuous(breaks = c(1:length(pareto_solutions_w_cost_uncertainty_df_cost_func_graph$solution_id)), 
#                      labels = pareto_solutions_w_cost_uncertainty_df_cost_func_graph$labels_xaxis)
# 
# dev.off()
# 
# 
# 


pareto_optimal_no_budget_fun<-function(obj1_conf,
                                       obj2_conf){
  
  
  #filter solution sets that do not meet the budget constraint
  function_df_all_temp<-function_df_all%>%
    group_by(variable)%>%
    mutate(rank_mean = rank(func_mean))
  
  #sort remaining solution sets by rank of objective
  function_df_w_ranks<-dcast(function_df_all_temp%>%
                               select(c('solution_id', 'variable', 'rank_mean')),
                             solution_id ~ variable,
                             fun.aggregate = mean,
                             value.var = 'rank_mean')
  
  function_df_w_ranks<-function_df_w_ranks%>%
    arrange(TBHIV_mort_diff_obj)%>%
    arrange(total_TBHIVmort)%>%
    select(c('solution_id', 'TBHIV_mort_diff_obj', 'total_TBHIVmort'))
  
  maintained_solutions<-c(function_df_w_ranks$solution_id[1])
  
  
  for(candidate_solutions_temp in 2:nrow(function_df_w_ranks)){
    
    candidate_solution_id_temp<-function_df_w_ranks$solution_id[candidate_solutions_temp]
    
    maintained_obj_ranks_df<-function_df_w_ranks%>%
      filter(solution_id %in% maintained_solutions)
    
    
    candidate_solution_df<-function_df_all_temp%>%
      filter(solution_id == candidate_solution_id_temp)
    
    candidate_solution_obj1_mean<-(candidate_solution_df$func_mean[candidate_solution_df$variable == 'total_TBHIVmort'])
    candidate_solution_obj2_mean<-(candidate_solution_df$func_mean[candidate_solution_df$variable == 'TBHIV_mort_diff_obj'])
    
    candidate_solution_obj1_var<-(candidate_solution_df$func_sd[candidate_solution_df$variable == 'total_TBHIVmort'])^2
    candidate_solution_obj2_var<-(candidate_solution_df$func_sd[candidate_solution_df$variable == 'TBHIV_mort_diff_obj'])^2
    
    #allows to change if based on mean or includes uncertainties
    dominated_or_not_dominated<-rep(0, times = nrow(maintained_obj_ranks_df))
    #0 non-dominated
    #1 dominated with confidence
    #2 dominated without confidence
    
    for(maintained_solutions_temp in 1:nrow(maintained_obj_ranks_df)){
      #first check based on mean
      if(maintained_obj_ranks_df$TBHIV_mort_diff_obj[maintained_solutions_temp] >
         function_df_w_ranks$TBHIV_mort_diff_obj[candidate_solutions_temp]){
        
        dominated_or_not_dominated[maintained_solutions_temp]<-0 #not dominated
        
      } else {
        
        #check if dominated with confidence
        maintained_solution_id_temp<-maintained_obj_ranks_df$solution_id[maintained_solutions_temp]
        
        maintained_solution_df<-function_df_all_temp%>%
          filter(solution_id == maintained_solution_id_temp)
        
        maintained_solution_obj1_mean<-(maintained_solution_df$func_mean[maintained_solution_df$variable == 'total_TBHIVmort'])
        maintained_solution_obj2_mean<-(maintained_solution_df$func_mean[maintained_solution_df$variable == 'TBHIV_mort_diff_obj'])
        
        maintained_solution_obj1_var<-(maintained_solution_df$func_sd[maintained_solution_df$variable == 'total_TBHIVmort'])^2
        maintained_solution_obj2_var<-(maintained_solution_df$func_sd[maintained_solution_df$variable == 'TBHIV_mort_diff_obj'])^2
        
        critical_val_obj1<-critical_val_fun(maintained_solution_obj1_mean,
                                            candidate_solution_obj1_mean,
                                            maintained_solution_obj1_var,
                                            candidate_solution_obj1_var,
                                            obj1_conf)
        critical_val_obj2<-critical_val_fun(maintained_solution_obj2_mean,
                                            candidate_solution_obj2_mean,
                                            maintained_solution_obj2_var,
                                            candidate_solution_obj2_var,
                                            obj2_conf)
        if(critical_val_obj1 <= 0 & 
           critical_val_obj2 <= 0){
          dominated_or_not_dominated[maintained_solutions_temp]<-1 #dominated with confidence
        } else {
          dominated_or_not_dominated[maintained_solutions_temp]<-2 #dominated without confidence
        }
      }
    }
    if(!(1 %in% dominated_or_not_dominated)){
      #if better in one of the other objectives, then non-dominated
      #if(!(2 %in% dominated_or_not_dominated)){
      maintained_solutions<-c(maintained_solutions, function_df_w_ranks$solution_id[candidate_solutions_temp])
      #}
    }
  }
  return(maintained_solutions)
}

obj1_conf_val_w_uncertainty<-0.95
obj2_conf_val_w_uncertainty<-0.95

obj1_conf_val_mean<-0
obj2_conf_val_mean<-0

n_samples<-816

no_budget_sol<-pareto_optimal_no_budget_fun(obj1_conf_val_w_uncertainty,
                                        obj2_conf_val_w_uncertainty)

no_budget_sol<-pareto_optimal_no_budget_fun(obj1_conf_val_mean,
                                            obj2_conf_val_mean)


#used for graph of pareto solutions with all community-based and all facility based
pareto_solutions_graph_with_extremes_df<-data.frame(solution_id = c(maintained_solution_cost_obj_uncertainty,
                                                                    no_budget_sol))%>%
  mutate(obj_uncertainty_only = if_else(solution_id %in% maintained_solution_cost_uncertainty_only, 'no', 'yes'))%>%
  full_join(function_df_all%>%filter(solution_id %in% c(maintained_solution_cost_obj_uncertainty,
                                                        no_budget_sol)), 
            by = c('solution_id'))%>%
  mutate(ll = if_else(variable == 'total_cost_fun', 
                      0,
                      func_mean - obj_z_score*func_sd/sqrt(n_samples)),
         ul = if_else(variable == 'total_cost_fun', 
                      func_mean + budget_z_score*(func_sd/sqrt(n_samples)),
                      func_mean + obj_z_score*(func_sd/sqrt(n_samples))))%>%
  select(-c('func_sd'))%>%
  filter(variable %in% c('TBHIV_mort_diff_obj', 'total_TBHIVmort', 'total_cost_fun'))

#used for pareto table
pareto_solutions_tab_with_extremes_df<-pareto_solutions_graph_with_extremes_df%>%
  mutate(conf_interval_txt = if_else(variable =='TBHIV_mort_diff_obj',
                                     paste0(round(func_mean), ' [', round(ll),
                                            ', ', round(ul), ']'),
                                     if_else(variable=='total_TBHIVmort', 
                                             paste0(sprintf('%.3f', func_mean/1000000), 
                                                    ' [', sprintf('%.3f', ll/1000000),
                                                    ', ', sprintf('%.3f', ul/1000000), ']'),
                                             paste0(sprintf('%.3f', func_mean/1000000000), 
                                                    ' [', sprintf('%.3f', ll/1000000000),
                                                    ', ', sprintf('%.3f', ul/1000000000), ']'))))%>%
  left_join(all_solution_sets_df, by = c('solution_id'))

##########GRAPHS#####

#### Pareto graph with extremes #####
pareto_solutions_graph_with_extremes_df2<-pareto_solutions_graph_with_extremes_df%>%
  select(c('solution_id', 'variable', 'func_mean', 'ul', 'll'))

pareto_solutions_graph_with_extremes_df2_yaxis<-pareto_solutions_graph_with_extremes_df2%>%
  filter(variable == 'total_TBHIVmort')%>%
  select(-c('variable'))

colnames(pareto_solutions_graph_with_extremes_df2_yaxis)<-
  c('solution_id', 'func_mean_y', 'ul_y', 'll_y')

pareto_solutions_graph_with_extremes_df2_xaxis<-pareto_solutions_graph_with_extremes_df2%>%
  filter(variable == 'TBHIV_mort_diff_obj')%>%
  select(-c('variable'))

colnames(pareto_solutions_graph_with_extremes_df2_xaxis)<-
  c('solution_id', 'func_mean_x', 'ul_x', 'll_x')

pareto_graph_with_extremes<-pareto_solutions_graph_with_extremes_df2_yaxis%>%
  left_join(pareto_solutions_graph_with_extremes_df2_xaxis, 
            by = c('solution_id'))%>%
  mutate(obj_uncertainty_only = if_else(solution_id %in% maintained_solution_cost_uncertainty_only, 'no', 'yes'))%>%
  arrange(func_mean_y)%>%
  arrange(obj_uncertainty_only)%>%
  mutate(func_mean_y = func_mean_y/1000000,
         ul_y = ul_y/1000000,
         ll_y = ll_y/1000000)

pareto_graph_with_extremes$`Solution ID`<-1:nrow(pareto_graph_with_extremes)
pareto_graph_with_extremes$color_ref<-if_else(pareto_graph_with_extremes$`Solution ID` %in% 1:5, 'blue', 
                                              if_else(pareto_graph_with_extremes$solution_id %in% c(no_budget_sol), 'red', 
                                                      'black'))


ggplot(pareto_graph_with_extremes, aes(x = func_mean_x, y = func_mean_y, 
                                       color = color_ref)
)+
  #geom_errorbarh(aes(xmin = ll_x, xmax = ul_x))+
  #geom_errorbar(aes(ymin = ll_y, ymax = ul_y), width = 2)+
  geom_point()+
  #geom_label(aes(label = `Solution ID`))+
  scale_color_manual(values = c('black', 'blue', 'red'))+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Objective 1: TB and HIV-related\nMortality Rate Inequity")+
  ylab("Objective 2: Total TB and HIV-related Deaths\nin Millions")+
  ylim(c(0, 1.5))+
  xlim(c(0,400))
dev.off()  
#stat_ellipse(data = ellipse_graph_df, aes(x = TBHIV_mort_diff_obj,
#                                   y = total_TBHIVmort,
#                                   group = solution_id))
