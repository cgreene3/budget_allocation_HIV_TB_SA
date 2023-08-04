rm(list = ls())
gc()

library(dplyr)
library(lhs)
library(readxl)
library(deSolve)


#####Set Up and Initialize#####
#read in parameters values and ranges#
model_params_df <- read_excel("~/github/resource_allocation_HIV_TB_SA/param_files/simple_sir_model_example_input_params/simple_sir_model_parameters.xlsx")

#identify parameters that vary over range
n_random_params<-sum(model_params_df$randomness == 'yes')
random_param_names<-model_params_df%>%
  filter(randomness == 'yes')%>%
  select(c('model_matched_param_region_program'))
random_param_names<-as.vector(t(random_param_names))

#####Initialize####
n_regions<-3
n_programs<-2
n_years <- 5
time_steps <- 1/12 #one month

#intialize current iteration at 0
current_iteration<-0

end_algorithm<-"no"
m_samples_per_iteration<-3
p_value_threshold<-0.01
max_samples<-100
max_iteration<-max_samples/m_samples_per_iteration

#identify all initial feasible solutions#
feasible_solutions_df<-expand.grid(rep(list(1:n_programs), n_regions))
colnames(feasible_solutions_df)<-c(paste0("region_",as.character(1:n_regions)))
feasible_solutions_df$solution_id <-1:nrow(feasible_solutions_df)

#initialize maintained solutions at all feasible solutions
maintained_solutions_temp<-feasible_solutions_df$solution_id

#store all samples
all_sample_df<-data.frame(matrix(nrow = 0, ncol = nrow(model_params_df)))
colnames(all_sample_df)<-model_params_df$model_matched_param_region_program

#store all metrics by region, program and sim id
metrics_by_region_program_sim_df<-data.frame(#end_prev = as.numeric(),
                           total_mort = as.numeric(),
                           total_cost = as.numeric(),
                           region_id = as.integer(),
                           program_id = as.integer(),
                           sim_id = as.integer())

#and objective function values for all simulations
by_sim_obj_function_values_df<-data.frame(total_mort_summ = as.numeric(),
                                            #end_prev_mean = as.numeric(),
                                            total_cost = as.numeric(),
                                            sim_id = as.integer(),
                                            solution_id = as.integer())


#####sample gen function#####
sample_gen_fun<-function(){
  current_iteration<<-current_iteration+1
  set.seed(current_iteration)
  raw_lhs_df<-as.data.frame(randomLHS(m_samples_per_iteration, 
                                      n_random_params))
  colnames(raw_lhs_df)<-random_param_names
  
  sample_df<-data.frame(matrix(nrow = m_samples_per_iteration, ncol = 0))
  
  for (mp in model_params_df$model_matched_param_region_program){
    if(model_params_df[model_params_df$model_matched_param_region_program == mp, 
                       c('randomness')] == "yes"){
      
      max_temp<-as.numeric(model_params_df[model_params_df$model_matched_param_region_program == mp, 
                                           c('max')])
      min_temp<-as.numeric(model_params_df[model_params_df$model_matched_param_region_program == mp, 
                                           c('min')])
      
      sample_df_temp<-as.data.frame(qunif(unlist(raw_lhs_df%>%select(c(mp))),
                                          min_temp, max_temp))
      colnames(sample_df_temp)<-mp
      sample_df<-cbind(sample_df, sample_df_temp)
    } else {
      value_temp<-as.numeric(model_params_df[model_params_df$model_matched_param_region_program == mp, 
                                             c('value')])
      sample_df_temp<-as.data.frame(rep(value_temp, times = m_samples_per_iteration))
      colnames(sample_df_temp)<-mp
      sample_df<-cbind(sample_df, sample_df_temp)
    }
  }
  sample_df$iteration_id<-rep(current_iteration, times = nrow(sample_df))
  sample_df$sim_id <-1:nrow(sample_df)+(m_samples_per_iteration*(current_iteration-1))
  rownames(sample_df)<-NULL
  all_sample_df<<-rbind(all_sample_df, sample_df)
}

######SIR MODEL EQUATIONS######
SIR_model_equations_fun <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    dS <- mu*Infected - beta*Susceptible*Infected
    dI <- beta*Susceptible*Infected - gamma*Infected - mu*Infected
    dR <- gamma*Infected
    dMort <- mu*Infected
    return(list(c(dS,dI,dR,dMort)))
  })
}

######GENERATE REGIONAL LEVEL PROGRAM HEALTH OUTCOMES AND COSTS#####
metrics_by_region_program_sim_fun<-function(){
  for (i in 1:n_regions){
    sample_df_temp<-all_sample_df%>%
      filter(iteration_id == current_iteration)
    
    sample_df_region_temp<-sample_df_temp[ , grepl(paste0("region_",i), 
                                                   names(sample_df_temp))|grepl("region_0", 
                                                                                names(sample_df_temp))|grepl("sim_id", names(sample_df_temp))]
    solution_test_temp<-feasible_solutions_df%>%
      filter(solution_id %in% maintained_solutions_temp)
    current_programs<-solution_test_temp[, grepl(paste0("region_",i), 
                                                      names(solution_test_temp))]                                                                                                         
    for (p in unique(current_programs)){
      
      sample_df_region_program_temp<-sample_df_region_temp[, grepl(paste0("program_",p), 
                                                                names(sample_df_region_temp))|grepl("program_0", 
                                                                                             names(sample_df_region_temp))|grepl("sim_id",
                                                                                                                                names(sample_df_region_temp))]
      for (k in sample_df_temp$sim_id){
        beta_temp<-sample_df_region_program_temp%>%
          filter(sim_id == k)%>%
          select(paste0("beta_region_",i, "_program_",p))
        
        beta_temp <- as.numeric(beta_temp)
        
        gamma_temp<-sample_df_region_program_temp%>%
          filter(sim_id == k)%>%
          select(paste0("gamma_region_0_program_0"))
        
        gamma_temp <-as.numeric(gamma_temp)
        
        mu_temp<-sample_df_region_program_temp%>%
          filter(sim_id == k)%>%
          select(paste0("mu_region_", i, "_program_0"))
        
        mu_temp<-as.numeric(mu_temp)
        
        I_init_temp<-sample_df_region_program_temp%>%
          filter(sim_id == k)%>%
          select(paste0("I_init_region_",i, "_program_0"))
        
        I_init_temp <-as.numeric(I_init_temp)
        
        init <- c(Susceptible=1-I_init_temp,Infected=I_init_temp,Recovered=0, Mort = 0)
        parameters <- c(beta=beta_temp,gamma=gamma_temp,mu=mu_temp)
        time <- seq(0,n_years,by=n_years/(2*length(1:n_years)))
        
        out<-ode(y=init,times=time,SIR_model_equations_fun,parms=parameters)
        out_df<-as.data.frame(out)
        metrics_by_region_program_sim_df_temp<-out_df%>%
          summarise(#end_prev = sum(if_else(time == n_years, Infected, 0)),
                    total_mort = sum(Mort))%>%
          mutate(total_cost = if_else(i == 3, if_else(p == 1, runif(1, 0.4, 0.6), runif(1, 2, 2.3)),
                                      if_else(i == 2, if_else(p == 1, runif(1, 0.5, 0.7), runif(1, 1.5, 1.7)),
                                              if_else(p == 1, runif(1, 0.4, 0.6), runif(1, 0.9, 1.1)))))%>%
          mutate(region_id = i,
                 program_id = p,
                 sim_id = k)
        metrics_by_region_program_sim_df<<-rbind(metrics_by_region_program_sim_df, metrics_by_region_program_sim_df_temp)
      }
    }
  }
}


#####CALCULATE MEAN OBJECTIVE FUNCTION VALUES#####
objective_mean_function_value_calcs_fun<-function(){
  
  mean_obj_function_values_df_all_temp<-data.frame(total_mort = as.numeric(),
                                               #end_prev_max = as.numeric(),
                                               total_cost = as.numeric(),
                                               solution_id = as.integer())
  
  maintained_solutions_temp_df<-feasible_solutions_df%>%
    filter(solution_id %in% maintained_solutions_temp)
  
  for (s in 1:nrow(maintained_solutions_temp_df)){
    region_program_in_solution_temp<-c()
    for (i in 1:n_regions){
      #get program implemented in region i under solution s
      program_temp<-maintained_solutions_temp_df[s, i] 
      region_program_in_solution_temp<-c(region_program_in_solution_temp,
                                    paste0("region_",i,"_program_",program_temp))
    }
    #calculate mean objective functions
    #min total mort
    mean_obj_function_values_df_temp<-metrics_by_region_program_sim_df%>%
      mutate(region_program_in_solution = paste0("region_",region_id,"_program_",program_id))%>%
      filter(region_program_in_solution %in% region_program_in_solution_temp)%>%
      group_by(region_id, program_id)%>%
      summarise(exp_total_mort = mean(total_mort),
                #exp_end_prev = mean(end_prev),
                exp_cost = mean(total_cost))%>%
      ungroup()%>%
      summarise(total_mort = sum(exp_total_mort),
                #end_prev_mean = mean(exp_end_prev),
                total_cost = sum(exp_cost))%>%
      mutate(solution_id = maintained_solutions_temp_df$solution_id[s])
    
    mean_obj_function_values_df_all_temp <- rbind(mean_obj_function_values_df_temp, 
                                         mean_obj_function_values_df_all_temp)
  }
  return(mean_obj_function_values_df_all_temp)
}

#####CALCULATE OBJECTIVE FUNCTION VALUES BY SIM ID#####
objective_function_value_by_sim_calcs_fun<-function(){
  
  maintained_solutions_temp_df<-feasible_solutions_df%>%
    filter(solution_id %in% maintained_solutions_temp)
  
  for (s in 1:nrow(maintained_solutions_temp_df)){
    region_program_in_solution_temp<-c()
    for (i in 1:n_regions){
      #get program implemented in region i under solution s
      program_temp<-maintained_solutions_temp_df[s, i] 
      region_program_in_solution_temp<-c(region_program_in_solution_temp,
                                         paste0("region_",i,"_program_",program_temp))
    }
    #calculate objective functions by sim id (only for new sims)
    
    new_samples_temp<-all_sample_df%>%
      filter(iteration_id == current_iteration)
    new_samples_temp<-new_samples_temp$sim_id
    
    by_sim_obj_function_values_df_temp<-metrics_by_region_program_sim_df%>%
      mutate(region_program_in_solution = paste0("region_",region_id,"_program_",program_id))%>%
      filter(region_program_in_solution %in% region_program_in_solution_temp,
             sim_id %in% new_samples_temp)%>%
      group_by(sim_id)%>%
      summarise(total_mort = sum(total_mort),
                #end_prev_mean = mean(end_prev),
                total_cost = sum(total_cost))%>%
      mutate(solution_id = maintained_solutions_temp_df$solution_id[s])
    
    by_sim_obj_function_values_df <<- rbind(by_sim_obj_function_values_df, 
                                            by_sim_obj_function_values_df_temp)
  }
}

#####FIND NON-DOMINATED PARETO OPTIMAL SOLUTIONS#####
rank_algorithm_fun<-function(mean_obj_function_values_df_all_temp,
                                  by_sim_obj_function_values_df){
  
  #sort all the solutions in decreasing order of the first objective function
  obj_function_values_df_current<-mean_obj_function_values_df_all_temp%>%
    mutate(rank_total_mort_min = rank(total_mort),
           #rank_end_prev_mean_min = rank(end_prev_mean),
           rank_total_cost_min = rank(total_cost))%>%
    arrange(rank_total_cost_min)%>%
    #arrange(rank_end_prev_mean_min)%>%
    arrange(rank_total_mort_min)
  
  ###step by step procedure for the algorithm#
  
  #add first element to list
  maintained_solutions<-c(obj_function_values_df_current$solution_id[1])
  end_algorithm_temp<-"yes"
  
  for (cs in 2:nrow(obj_function_values_df_current)){
    maintained_obj_ranks_df<-obj_function_values_df_current%>%
      filter(solution_id %in% maintained_solutions)
    
    dominated_or_not_dominated<-rep(0, times = nrow(maintained_obj_ranks_df))  
    #0 non-dominated
    #1 dominated with confidence
    #2 dominated without confidence
    
    for (ms in 1:nrow(maintained_obj_ranks_df)){
      if(maintained_obj_ranks_df$rank_total_mort_min[ms] >
         obj_function_values_df_current$rank_total_mort_min[cs]|
         #maintained_obj_ranks_df$rank_end_prev_mean_min[ms] >
        #    obj_function_values_df_current$rank_end_prev_mean_min[cs]|
         maintained_obj_ranks_df$rank_total_cost_min[ms] >
              obj_function_values_df_current$rank_total_cost_min[cs]){
        dominated_or_not_dominated[ms]<-0 #not dominated
      } else {
        
        cs_sim_objective_function_values_temp<-by_sim_obj_function_values_df%>%
          filter(solution_id == obj_function_values_df_current$solution_id[cs])
        ms_sim_objective_function_values_temp<-by_sim_obj_function_values_df%>%
          filter(solution_id == maintained_obj_ranks_df$solution_id[ms])
        
        t_test_total_mort<-t.test(ms_sim_objective_function_values_temp$total_mort,
                                  cs_sim_objective_function_values_temp$total_mort,
                                  alternative = "less",
                                  paired = TRUE,
                                  var.equal = FALSE)#,
                                  #method = "bonferroni")
        t_test_total_cost<-t.test(ms_sim_objective_function_values_temp$total_cost,
                                  cs_sim_objective_function_values_temp$total_cost,
                                  alternative = "less",
                                  paired = TRUE,
                                  var.equal = FALSE)#,
                                  #method = "bonferroni")
        if((t_test_total_mort$p.value <= p_value_threshold)&
           (t_test_total_cost$p.value <= p_value_threshold)){
          dominated_or_not_dominated[ms]<-1 #dominated with confidence
        } else {
          dominated_or_not_dominated[ms]<-2 #dominated without confidence
        }
      }
    }
      #maintain, delete, do not stop
      if(sum(dominated_or_not_dominated) == 0){
        maintained_solutions<-c(maintained_solutions, 
                                obj_function_values_df_current$solution_id[cs])
      } else {
        if(!(1 %in% dominated_or_not_dominated)){
          end_algorithm_temp<-"no"
          maintained_solutions<-c(maintained_solutions, 
                                  obj_function_values_df_current$solution_id[cs])
        }
      }
  }
  return(list(maintained_solutions = maintained_solutions,
              end_algorithm = end_algorithm_temp))
}

while (end_algorithm == "no") {
  sample_gen_fun()
  metrics_by_region_program_sim_fun()
  mean_obj_function_values_df_all_temp<-objective_mean_function_value_calcs_fun()
  objective_function_value_by_sim_calcs_fun()
  rank_algorithm_list<-rank_algorithm_fun(mean_obj_function_values_df_all_temp,
                                          by_sim_obj_function_values_df)
  print(rank_algorithm_list)
  maintained_solutions_temp<<-rank_algorithm_list$maintained_solutions
  end_algorithm<<-rank_algorithm_list$end_algorithm  
  if(end_algorithm == "no"){
    end_algorithm<<-if_else(current_iteration+1 <= max_iteration,  "no", "yes")
    if(current_iteration+1 > max_iteration){
      print("reached max iterations")
      print(current_iteration)
    }
  }
}






