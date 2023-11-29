#clean workspace
rm(list = ls())
gc()

library(here)
library(dplyr)
library(stringr)
library(reshape2)

#different code to find pareto optimal solutions
###calculate objective functions for each solution set####
n_programs <-3
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

function_df_subset<-function_df_all%>%
  select(c('solution_id', 'variable', 'func_mean', 'func_sd'))

###critical val function
critical_val_fun<-function(mean_maintained,
                           mean_candidate,
                           var_maintained,
                           var_candidate,
                           conf_level){
  
  obj_z_score<-qnorm(1-((1-conf_level)/2))#two tailed
  
  critical_val_temp<-(mean_maintained-mean_candidate)+
    (obj_z_score*(sqrt(var_maintained+var_candidate)))
  
  return(critical_val_temp)
}

pareto_optimal_algorithm_fun<-function(budget,
                                       #obj1_conf,
                                       #obj2_conf,
                                       obj3_conf,
                                       obj4_conf,
                                       cost_conf){
  
  #get cut off for budget
  budget_z_score<-qnorm(cost_conf)
  
  solution_sets_after_budget<-function_df_subset%>%
    filter(variable == 'total_cost')%>%
    mutate(conf_interval_upper = func_mean + (func_sd*budget_z_score))%>%
    filter(conf_interval_upper <= budget)
  
  print('n solutions after budget')
  print(nrow(solution_sets_after_budget))
  
  #filter solution sets that do not meet the budget constraint
  function_df_subset_after_budget_filter<-function_df_subset%>%
    filter(solution_id %in% solution_sets_after_budget$solution_id)%>%
    group_by(variable)%>%
    mutate(rank_mean = rank(func_mean))
  
  rm(solution_sets_after_budget)
  
  #sort remaining solution sets by rank of objective
  function_df_subset_w_ranks<-dcast(function_df_subset_after_budget_filter%>%
                                           select(c('solution_id', 'variable', 'rank_mean')),
                                         solution_id ~ variable,
                                         fun.aggregate = mean,
                                         value.var = 'rank_mean')
    
  function_df_subset_w_ranks<-function_df_subset_w_ranks%>%
    arrange(equity_TB_inc_Y)%>%
    arrange(equity_HIV_prev_Y)%>%
    arrange(equity_TBHIV_mort_Y)%>%
    arrange(total_TBHIV_mort_num)%>%
    select(-c('total_cost'))
  
  maintained_solutions<-c(function_df_subset_w_ranks$solution_id[1])
               
  
  for (candidate_solutions_temp in 2:nrow(function_df_subset_w_ranks)
       ){
    
    candidate_solution_id_temp<-function_df_subset_w_ranks$solution_id[candidate_solutions_temp]
    
    maintained_obj_ranks_df<-function_df_subset_w_ranks%>%
      filter(solution_id %in% maintained_solutions)
    
    
    candidate_solution_df<-function_df_subset_after_budget_filter%>%
      filter(solution_id == candidate_solution_id_temp)
    
    #candidate_solution_obj1_mean<-(candidate_solution_df$func_mean[candidate_solution_df$variable == 'equity_TB_inc_Y'])
    #candidate_solution_obj2_mean<-(candidate_solution_df$func_mean[candidate_solution_df$variable == 'equity_HIV_prev_Y'])
    candidate_solution_obj3_mean<-(candidate_solution_df$func_mean[candidate_solution_df$variable == 'equity_TBHIV_mort_Y'])
    candidate_solution_obj4_mean<-(candidate_solution_df$func_mean[candidate_solution_df$variable == 'total_TBHIV_mort_num'])
    
    #candidate_solution_obj1_var<-(candidate_solution_df$func_sd[candidate_solution_df$variable == 'equity_TB_inc_Y'])^2
    #candidate_solution_obj2_var<-(candidate_solution_df$func_sd[candidate_solution_df$variable == 'equity_HIV_prev_Y'])^2
    candidate_solution_obj3_var<-(candidate_solution_df$func_sd[candidate_solution_df$variable == 'equity_TBHIV_mort_Y'])^2
    candidate_solution_obj4_var<-(candidate_solution_df$func_sd[candidate_solution_df$variable == 'total_TBHIV_mort_num'])^2
    
    
    #allows to change if based on mean or includes uncertainties
    dominated_or_not_dominated<-rep(0, times = nrow(maintained_obj_ranks_df))
    #0 non-dominated
    #1 dominated with confidence
    #2 dominated without confidence
    
    for(maintained_solutions_temp in 1:nrow(maintained_obj_ranks_df)){
      #first check based on mean
      if(#maintained_obj_ranks_df$equity_TB_inc_Y[maintained_solutions_temp] >
         #function_df_subset_w_ranks$equity_TB_inc_Y[candidate_solutions_temp]|
         #maintained_obj_ranks_df$equity_HIV_prev_Y[maintained_solutions_temp] >
         #function_df_subset_w_ranks$equity_HIV_prev_Y[candidate_solutions_temp]|
         maintained_obj_ranks_df$equity_TBHIV_mort_Y[maintained_solutions_temp] >
         function_df_subset_w_ranks$equity_TBHIV_mort_Y[candidate_solutions_temp]){
        
        dominated_or_not_dominated[maintained_solutions_temp]<-0 #not dominated

      } else {
        
        #check if dominated with confidence
        maintained_solution_id_temp<-maintained_obj_ranks_df$solution_id[maintained_solutions_temp]
        
        maintained_solution_df<-function_df_subset_after_budget_filter%>%
          filter(solution_id == maintained_solution_id_temp)
        
        #maintained_solution_obj1_mean<-(maintained_solution_df$func_mean[maintained_solution_df$variable == 'equity_TB_inc_Y'])
        #maintained_solution_obj2_mean<-(maintained_solution_df$func_mean[maintained_solution_df$variable == 'equity_HIV_prev_Y'])
        maintained_solution_obj3_mean<-(maintained_solution_df$func_mean[maintained_solution_df$variable == 'equity_TBHIV_mort_Y'])
        maintained_solution_obj4_mean<-(maintained_solution_df$func_mean[maintained_solution_df$variable == 'total_TBHIV_mort_num'])
        
        #maintained_solution_obj1_var<-(maintained_solution_df$func_sd[maintained_solution_df$variable == 'equity_TB_inc_Y'])^2
        #maintained_solution_obj2_var<-(maintained_solution_df$func_sd[maintained_solution_df$variable == 'equity_HIV_prev_Y'])^2
        maintained_solution_obj3_var<-(maintained_solution_df$func_sd[maintained_solution_df$variable == 'equity_TBHIV_mort_Y'])^2
        maintained_solution_obj4_var<-(maintained_solution_df$func_sd[maintained_solution_df$variable == 'total_TBHIV_mort_num'])^2
        
        #critical_val_obj1<-critical_val_fun(maintained_solution_obj1_mean,
        #                                    candidate_solution_obj1_mean,
        #                                    maintained_solution_obj1_var,
        #                                    candidate_solution_obj1_var,
        #                                    obj1_conf)
        #critical_val_obj2<-critical_val_fun(maintained_solution_obj2_mean,
        #                                    candidate_solution_obj2_mean,
        #                                    maintained_solution_obj2_var,
        #                                    candidate_solution_obj2_var,
        #                                    obj2_conf)
        critical_val_obj3<-critical_val_fun(maintained_solution_obj3_mean,
                                            candidate_solution_obj3_mean,
                                            maintained_solution_obj3_var,
                                            candidate_solution_obj3_var,
                                            obj3_conf)
        critical_val_obj4<-critical_val_fun(maintained_solution_obj4_mean,
                                            candidate_solution_obj4_mean,
                                            maintained_solution_obj4_var,
                                            candidate_solution_obj4_var,
                                            obj4_conf)
        if(#critical_val_obj1 <= 0 & 
           #critical_val_obj2 <= 0 & 
           critical_val_obj3 <= 0 & 
           critical_val_obj4 <= 0){
          dominated_or_not_dominated[maintained_solutions_temp]<-1 #dominated with confidence
        } else {
          dominated_or_not_dominated[maintained_solutions_temp]<-2 #dominated without confidence
        }
      }
    }
    if(!(1 %in% dominated_or_not_dominated)){
      #if better in one of the other objectives, then non-dominated
      #if(!(2 %in% dominated_or_not_dominated)){
        maintained_solutions<-c(maintained_solutions, function_df_subset_w_ranks$solution_id[candidate_solutions_temp])
      #}
    }
  }
  return(maintained_solutions)
}
  


budget_val<-16000000000
obj1_conf_val<-0
obj2_conf_val<-0
obj3_conf_val<-0.1
obj4_conf_val<-0.1
cost_conf_val<-.9



maintained_solution_all<-pareto_optimal_algorithm_fun(budget_val,
                                                      #obj1_conf_val,
                                                      #obj2_conf_val,
                                                      obj3_conf_val,
                                                      obj4_conf_val,
                                                      cost_conf_val)
  