#clean workspace
rm(list = ls())
gc()


sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle'), require, character.only=T)

# # ####HYAK OR GITHUB SPECIFIC CODES TO COMMENT/UNCOMMENT####
#hyak specific code
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {
  stop("only need one parameter", call.=FALSE)
}

#different code to find pareto optimal solutions
###calculate objective functions for each solution set####
n_programs <-2
n_regions<-9

all_solution_sets_df<-expand.grid(rep(list(1:n_programs), n_regions))
colnames(all_solution_sets_df)<-1:n_regions
all_solution_sets_df$solution_id <-1:nrow(all_solution_sets_df)

all_solution_sets_df<-melt(all_solution_sets_df, id = "solution_id")
all_solution_sets_df<-all_solution_sets_df%>%
  rename('region_id' = 'variable')%>%
  rename('program_id' = 'value')%>%
  mutate("region_program_id" = paste0(region_id,'_',program_id))

#32 solution sets per run
#args_temp = 1 #if running local
args_temp<-as.integer(args)
args_start<-1+(args_temp-1)*16
args_end<-args_start+15

outdir<-'/gscratch/icrc/cgreene3/SA_resource_allocation/program_runs/function_values_by_solution_set'

#read in param sets
input_params_dir<-'/gscratch/icrc/cgreene3/SA_resource_allocation/input_parameters'
setwd(input_params_dir)
accepted_param_sets_ref_df<-read.csv('accepted_param_sets_ref.csv')
region_id_ref_df<-read.csv('region_id_ref.csv')

#read in metrics
setwd('/gscratch/icrc/cgreene3/SA_resource_allocation/program_runs/metrics_summarised_i_p_k')
metrics_df_i_p_k<-read.csv('metrics_summarised.csv')


#github
# library(here)
# outdir<-paste0(here(), '/test/obj_func_values')
# input_params_dir<-paste0(here(), '/param_files/input_parameters')
# 
# setwd(input_params_dir)
# accepted_param_sets_ref_df<-read.csv('accepted_param_sets_ref.csv')
# region_id_ref_df<-read.csv('region_id_ref.csv')
# 
# #objective function df
# setwd(paste0(here(), '/results/program_runs/metrics_summarised_i_p_k'))
# metrics_df_i_p_k<-read.csv('metrics_summarised.csv')

function_df<-data.frame()

lapply(unique(all_solution_sets_df$solution_id[args_start:args_end]),
       function(uip){
  
  #get uips
  subset_solution_sets_df<-all_solution_sets_df%>%
    filter(solution_id == uip)
  
  #get metrics
  subset_metrics_solution<-metrics_df_i_p_k%>%
    filter(region_program_id %in% subset_solution_sets_df$region_program_id)
  
  
  #calculate each objective for each sim id
  setwd(paste0(input_params_dir, '/accepted_param_sets_combin'))
  files <- list.files(pattern="param_set_gen*")
  
  solution_vals_all<-data.frame()
  
  lapply(files, function(f){
    
    #get all solution sets combin across all regions accepted
    sim_set_df<-read.csv(f)%>%
      select(-c('total_disease_mort_SA_2005',
                'total_disease_mort_SA_2017'))
    sim_set_df_melt<-melt(sim_set_df,
                          id = c("sim_id",
                                 "general_sim_id"))
    
    sim_set_df_melt<-sim_set_df_melt%>%
      rename("regional_sim_id" = "value")%>%
      mutate(region_id = str_split(variable, "_", simplify = T)[, 3])%>%
      select(-c('variable'))
    
    solution_vals<-sim_set_df_melt%>%
      mutate(region_id = as.integer(region_id))%>%
      left_join(subset_metrics_solution, 
                by = c('general_sim_id', 'regional_sim_id', 'region_id'))%>%
      group_by(sim_id)%>%
      mutate(mean_tb_inc = mean(TB_inc_per_Y_100k_ppl),
             mean_tbhiv_mort = mean(TBHIV_mort_per_Y_100k_ppl))%>%
      mutate(diff_TB_inc = TB_inc_per_Y_100k_ppl-mean_tb_inc,
             diff_TBHIV_mort = TBHIV_mort_per_Y_100k_ppl-mean_tbhiv_mort)%>%
      summarise(obj1_fun = max(diff_TB_inc),
                obj2_fun = max(diff_TBHIV_mort),
                obj3_fun = sum(TBHIV_mort_num_total),
                total_cost_fun = sum(program_cost))
    
    if(f == files[1]){
      solution_vals_all<<-solution_vals
    } else {
      solution_vals_all<<-rbind(solution_vals_all,
                          solution_vals)
    }
  })
  #calculate mean sd and conf intervals
  
  #get zscores
  #two_sided_95<-qnorm(1-(.05/2))
  #two_sided_99<-qnorm(1-(.01/2))
  #two_sided_90<-qnorm(1-(.1/2))
  #one_sided_95<-qnorm(1-.05)
  #one_sided_99<-qnorm(1-.01)
  #one_sided_90<-qnorm(1-.1)

  func_df_temp<-melt(solution_vals_all, id = c("sim_id"))
  func_df_temp<-func_df_temp%>%
    mutate(solution_id = uip)%>%
    group_by(solution_id, variable)%>%
    summarise(func_mean = mean(value),
              func_sd = sd(value))#%>%
   # mutate(conf_interval_lower_95 = if_else(variable == 'total_cost', 0, 
  #                                        (func_mean - (func_sd*two_sided_95))),
   #        conf_interval_upper_95 = if_else(variable == 'total_cost', func_mean + (func_sd*one_sided_95), 
    #                                      (func_mean + (func_sd*two_sided_95))),
    #       conf_interval_lower_90 = if_else(variable == 'total_cost', 0, 
     #                                     (func_mean - (func_sd*two_sided_90))),
      #     conf_interval_upper_90 = if_else(variable == 'total_cost', func_mean + (func_sd*one_sided_90), 
       #                                   (func_mean + (func_sd*two_sided_90))),
      #     conf_interval_lower_99 = if_else(variable == 'total_cost', 0, 
       #                                   (func_mean - (func_sd*two_sided_99))),
        #   conf_interval_upper_99 = if_else(variable == 'total_cost', func_mean + (func_sd*one_sided_99), 
         #                                 (func_mean + (func_sd*two_sided_99))))
    
  if(uip == unique(all_solution_sets_df$solution_id[args_start])){
    function_df<<-func_df_temp
  } else {
    function_df<<-rbind(function_df,
                        func_df_temp)
  }
  print(uip)
  print(Sys.time())
})

setwd(outdir)
write.csv(function_df, 
          paste0('function_df_solution_id_', args_start, '_to_', args_end, '.csv'),
          row.names = FALSE)