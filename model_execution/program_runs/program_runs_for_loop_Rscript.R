#last updated August 22, 2023
#runs warmup from 1940 to end of 2017 (beg of 2018)
#writes csv for model outputs from 1990 to 2027

#clean workspace
rm(list = ls())
gc()

#load packages
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

#args_temp<- as.integer(args[1])
region_id_temp<-as.integer(args[1])

# # #location where input parameters are
indir<-'/gscratch/icrc/cgreene3/SA_resource_allocation/input_parameters'
# # #
# # # #location where want outputs
outdir<-paste0('/gscratch/icrc/cgreene3/SA_resource_allocation/program_runs/metrics')

# # # #if running from github
#library(here)
#if running from local set env to resource_allocation_HIV_TB_SA
#region_id_temp<-1

#location where input parameters are
#indir<-paste0(here(), '/param_files/input_parameters')

#location where want outputs
#outdir<-paste0(here(), '/results/program_runs/metrics')

########Time Horizon and Evaluation intervals (1 month)#####
time_interval <- 1/12

#specify warmup period
start_yr_warmup <- 1940
end_yr_warmup <- 2018 #start of 2018 
TT_warmup<-end_yr_warmup-start_yr_warmup
TT_SET_warmup <- seq(from = 0, to = TT_warmup, by = time_interval) #tau

#specify eval period
start_yr_eval <- 2018
end_yr_eval <- 2028  #start of 2028
TT_eval<-end_yr_eval-start_yr_eval
TT_SET_eval <- seq(from = 0, to = TT_eval, by = time_interval) #tau


#these data frames change depending on combin of calib parameters
setwd(indir)

sim_calibration_ref_df<-read.csv('calibration_sets_df.csv')
accepted_param_sets_ref_df<-read.csv('accepted_param_sets_ref.csv')
region_id_ref_df<-read.csv('region_id_ref.csv')

#make temp dfs
pop_init_df<-data.frame()
mort_rate_df<-data.frame()
birth_perc_df<-data.frame()
ipt_initiation_df<-data.frame()
art_coverage_df<-data.frame()
art_prop_eligible_df<-data.frame()
hiv_incidence_df<-data.frame()

input_data_param_extraction<-function(region_name_temp){
  setwd(indir)
  pop_init_df<- read.csv('pop_init_df_1940.csv')%>%
    filter(region_name == region_name_temp)
  mort_rate_df <- read.csv('base_mort_df.csv')%>%
    filter(region_name == region_name_temp)
  birth_perc_df<-read.csv('birth_perc_df_overtime.csv')%>%
    filter(region_name == region_name_temp)
  ipt_initiation_df <-read.csv('ipt_initiation_df.csv')
  art_coverage_df <-read.csv('art_coverage_df.csv')
  art_prop_eligible_df<-read.csv('ART_prop_eligible_df.csv')
  hiv_incidence_df<-read.csv('hiv_inc_df.csv')%>%
    filter(region_name == region_name_temp)
  
  #clean dataframe column names for consistency
  names(pop_init_df)<-str_replace_all(names(pop_init_df), c(" " = "_" , "-" = "_" ))
  
  #establish compartment and dcompartment ids
  pop_init_df<- pop_init_df%>%
    mutate(compartment_id = paste0("N_", TB_compartment, "_", DR_compartment ,"_", HIV_compartment, "_", G_compartment),
           dcompartment_id = paste0("dN_", TB_compartment, "_", DR_compartment ,"_", HIV_compartment, "_", G_compartment))
  
  return(list(pop_init_df,
              mort_rate_df,
              birth_perc_df,
              ipt_initiation_df,
              art_coverage_df,
              art_prop_eligible_df,
              hiv_incidence_df))
}

################ Define sets ###############################

#3 programs
program_sets<-1:3

#7 TB set description (TB)#
#1:Uninfected
#2:LTBI, infected recently, no TPT
#3: LTBI, infected remotely, no TPT
#4: LTBI, infected recently, post-TPT
#5: LTBI, infected remotely, post-TPT
#6: Active
#7: Recovered/Treated

TB_SET<-1:7

#2 Drug Resistance compartments description (DR)#
#1 : Drug Susceptible
#2 : Multi Drug Resistant

DR_SET<-1:2

#4 HIV compartments description (HIV)#
#1 : HIV Negative
#2 : HIV Positive CD4 > 200 - No ART
#3 : HIV Positive CD4 =<: 200 - No Art
#4 : HIV Positive - ART 

HIV_SET<-1:4

#2 Gender compartments description (G)#
#1: Male
#2: Female

G_SET<-1:2


#define number of compartments
n_compartments<-length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET)

########INITIALIZE POPULATION AND ID LOCATION OF CURRENT POPULATION AMOUNTS######

#N_t_r_h_g - matrix that identifies the location of compartment in 1D array#
N_t_r_h_g_ref <- array(0, dim = c(length(TB_SET), length(DR_SET), length(HIV_SET), length(G_SET)))
count_temp <- 1

for (t in TB_SET){
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        N_t_r_h_g_ref[t,r,h,g] <- count_temp
        count_temp <- count_temp + 1
      }
    }
  }
}

#######Parameter Extraction########

######## Parameters that impact force of infection #######

#beta_g - Number of effective contacts for TB transmission per infectious year#
beta_g <- array(0, dim = length(G_SET))

beta_g_param_extraction<-function(sim_id_current_eval){
  
  beta_g <- array(0, dim = length(G_SET))
  
  for (g in G_SET){
    
    temp_sim_val<-sim_calibration_ref_df%>%
      filter(sim_id == sim_id_current_eval)%>%
      select(paste0("beta_"))
    
    beta_g[g] <- as.numeric(temp_sim_val)
    
  }
  return(beta_g)
}


#phi_h - Relative transmissibility of TB in HIV pops#
phi_h <- array(0, dim = length(HIV_SET))

phi_h_param_extraction<-function(sim_id_current_eval){
  phi_h <- array(0, dim = length(HIV_SET))
  for (h in HIV_SET){
    
    temp_sim_val<-sim_calibration_ref_df%>%
      filter(sim_id == sim_id_current_eval)%>%
      select(paste0("phi_",h,"."))
    
    phi_h[h] <- as.numeric(temp_sim_val)
  }
  return(phi_h)
}

#varepsilon_g - Fraction of new TB infections that are MDR-TB #
varepsilon_g <- array(0, dim = length(G_SET))

varepsilon_g_param_extraction<-function(sim_id_current_eval){
  varepsilon_g <- array(0, dim = length(G_SET))
  
  for (g in G_SET){
    
    temp_sim_val<-sim_calibration_ref_df%>%
      filter(sim_id == sim_id_current_eval)%>%
      select("varepsilon_") #not differentiated by g
    
    varepsilon_g[g] <- as.numeric(temp_sim_val)
  }
  return(varepsilon_g)
}

#increased risk of re-infection for recovered populations#
xi <- 0

xi_param_extraction<-function(sim_id_current_eval){
  
  temp_sim_val<-sim_calibration_ref_df%>%
    filter(sim_id == sim_id_current_eval)%>%
    select("xi_")
  
  xi <-as.numeric(temp_sim_val)
  
  return(xi)
}

#zeta - Indicator that diminishes force of infection due to the partially-protective effect of LTBI infection and acquiring a new TB infection#
zeta<-0

zeta_param_extraction<-function(sim_id_current_eval){
  
  temp_sim_val<-sim_calibration_ref_df%>%
    filter(sim_id == sim_id_current_eval)%>%
    select("zeta_")
  
  zeta <-as.numeric(temp_sim_val)
  return(zeta)
}

#########Parameters that describe TB progression ######

#omega - Rate of moving off of IPT, per year #
temp_sim_val<-sim_calibration_ref_df%>%
  select("omega_")
omega <- as.numeric(unique(temp_sim_val))

#pi_i_t - Base rates of TB progression#
pi_i_t <- array(data = 0, c(length(TB_SET), length(TB_SET)))

pi_i_t_param_extraction<-function(sim_id_current_eval){
  pi_i_t <- array(data = 0, c(length(TB_SET), length(TB_SET)))
  for (t_from in TB_SET){
    for (t_to in TB_SET){
      num_temp = (t_from*10) + t_to
      
      if(paste0("pi_",num_temp,".") %in% colnames(sim_calibration_ref_df)){
        
      
      temp_sim_val<-sim_calibration_ref_df%>%
        filter(sim_id == sim_id_current_eval)%>%
        select(paste0("pi_",num_temp,"."))

        pi_i_t[t_from,t_to] <- as.numeric(temp_sim_val)
      }
    }
    }
  return(pi_i_t)
}

#theta_h - relative risk of TB progression#
theta_h <-array(0, dim = length(HIV_SET))

theta_h_param_extraction<-function(sim_id_current_eval){
  theta_h <-array(0, dim = length(HIV_SET))
  for (h in HIV_SET){
    temp_sim_val<-sim_calibration_ref_df%>%
        filter(sim_id == sim_id_current_eval)%>%
        select(paste0("theta_",h,"."))
      theta_h[h] <- as.numeric(temp_sim_val)
    } 
  return(theta_h)
}

#chi - relative reduced risk of TB progression for those on TPT with LTBI
chi<-0

chi_param_extraction<-function(sim_id_current_eval){
  
  temp_sim_val<-sim_calibration_ref_df%>%
    filter(sim_id == sim_id_current_eval)%>%
    select("chi_")
  
  chi <-as.numeric(temp_sim_val)
  return(chi)
}

#upsilon_h - increased time infectious due to delayed treatment#
upsilon_h <- array(0, dim = length(HIV_SET))

upsilon_h_param_extraction<-function(sim_id_current_eval){
  upsilon_h <- array(0, dim = length(HIV_SET))
  
  for (h in HIV_SET){
    temp_sim_val<-sim_calibration_ref_df%>%
        filter(sim_id == sim_id_current_eval)%>%
        select(paste0("upsilon_",h,"."))
      upsilon_h[h] <- as.numeric(temp_sim_val)
    } 
  return(upsilon_h)
}

#gamma_r -indicator if DR compartment can move onto after IPT#
gamma_r <- c(1,0)

#TPT initiation factor for calibration
kappa_factor<-0

kappa_factor_param_extraction<-function(sim_id_current_eval){

  temp_sim_val<-sim_calibration_ref_df%>%
    filter(sim_id == sim_id_current_eval)%>%
    select(paste0("kappa.factor_"))
    
  kappa_factor <- temp_sim_val[[1]]
  return(kappa_factor)
}

#kappa_t_h_g_p - Rate of IPT initiation, per year#

#policy parameter#
#time varying
kappa_t_h_g <- array(data = 0, c(length(TB_SET), length(HIV_SET), length(G_SET)))

kappa_param_extraction_function<-function(yr){
  
  kappa_t_h_g <- array(data = 0, c(length(TB_SET), length(HIV_SET), length(G_SET)))
  
  for (g in G_SET){
    
    if(yr <= 2017){
      IPT_init_temp<-ipt_initiation_df%>%
        filter(gender_id == g,
               year == yr,
               program_id == 1)
    } else {
      IPT_init_temp<-ipt_initiation_df%>%
        filter(gender_id == g,
               year == 2018,
               program_id == current_program_eval)
    }
    
    
    #initiate from TB compartments 2 and 3, 0 rate otherwise
    kappa_t_h_g[c(2,3), 4, g]<-kappa_factor*(IPT_init_temp$ipt_init_perc)
  }
  return(kappa_t_h_g)
}

#########Parameters that describe HIV progression ######
eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                            length(HIV_SET), 
                            length(G_SET)))

#hiv incidence factor for calibration
eta_factor<-0

eta_factor_param_extraction<-function(sim_id_current_eval){
    
  temp_sim_val<-sim_calibration_ref_df%>%
    filter(sim_id == sim_id_current_eval)%>%
    select("eta.factor_")
    
    eta_factor <- temp_sim_val[[1]]
    return(eta_factor)
}

#art coverage factor for calibration
sigma_factor<-0

sigma_factor_param_extraction<-function(sim_id_current_eval){
  
  temp_sim_val<-sim_calibration_ref_df%>%
    filter(sim_id == sim_id_current_eval)%>%
    select("sigma.factor_")
  
  sigma_factor <- temp_sim_val[[1]]
  return(sigma_factor)
}

#changes 
HIV_transitions_param_func<-function(yr, N_t_r_h_g){
  
  eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                              length(HIV_SET), 
                              length(G_SET)))
  
  #hiv incidence
  #calibrated
  #time varying
  for (g in G_SET){
    hiv_incidence_temp<-hiv_incidence_df%>%
      filter(sex == if_else(g == 1, 'Male', 'Female'),
             year == yr, 
             program_id == current_program_eval)
    
    eta_i_h_g[1,2,g]<-as.numeric(eta_factor*hiv_incidence_temp$val[1])
  }
  
  #CD4 decline
  for (g in G_SET){
    
    temp_sim_val<-sim_calibration_ref_df%>%
      filter(sim_id == sim_id_current_eval)%>%
      select(paste0("eta_23.",g))
    
    eta_i_h_g[2,3,g] <- as.numeric(temp_sim_val)
  
    if(yr >= 2004){
      
      #get current model states (V_h_g in the appendix)
      #time varying based on model states
      hiv_prev_current<-sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4), g]])
      n2_prop_current_yr <- sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 2, g]])/hiv_prev_current
      n3_prop_current_yr <- sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 3, g]])/hiv_prev_current
      n4_prop_current_yr <- sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 4, g]])/hiv_prev_current
      
      #ART coverage (calibrated and time varying) #sigma_g(\tau)
      if(yr <= 2017){
      art_coverage_temp<-art_coverage_df%>%
        filter(sex == if_else(g == 1, 'Male', 'Female'),
               year == yr,
               program_id == 1)
      } else {
        art_coverage_temp<-art_coverage_df%>%
          filter(sex == if_else(g == 1, 'Male', 'Female'),
                 year == 2018,
                 program_id == current_program_eval)
      }
      
      #adjusted with calibration factor
      art_coverage_temp<-sigma_factor*(art_coverage_temp$art_coverage)
      
      #get percent eligible (time varying) 
      if(yr <= 2015){
      prop_eligible_temp<-art_prop_eligible_df%>%
        filter(gender_id == g,
               year == yr)
      prop_eligible<<-prop_eligible_temp$percent_hiv_compartment_2_eligible
      } else {
        prop_eligible <<-1
      }
      
      
      #calculate overall art initation rate (all CD4 counts)
      #equation E in the appendix
      art_initiation_rate_all_numerator<-art_coverage_temp-n4_prop_current_yr
      
      art_initiation_rate_all_denominator<-(prop_eligible*n2_prop_current_yr)+(n3_prop_current_yr)
      
      art_initiation_rate_all<-art_initiation_rate_all_numerator/
        art_initiation_rate_all_denominator
      
      eta_i_h_g[2,4,g]<-art_initiation_rate_all*prop_eligible
      eta_i_h_g[3,4,g]<-art_initiation_rate_all
    }
  }
  return(eta_i_h_g)
} 

#########Parameters for death and aging rates ###########
#alpha_out - Rate of exit from the population#
alpha_out<-sim_calibration_ref_df%>%
  select("alpha.out_")
alpha_out <- as.numeric(unique(alpha_out))

##placeholder for mu_t_h_g
mu_t_h_g<-array(0, dim=c(length(TB_SET), 
                               length(HIV_SET), 
                               length(G_SET)))

#non-disease mort
mu_factor<-0

mu_factor_param_extraction<-function(sim_id_current_eval){
  
  temp_sim_val<-sim_calibration_ref_df%>%
    filter(sim_id == sim_id_current_eval)%>%
    select("mu.factor_")
  
  mu_factor <- temp_sim_val[[1]]
  return(mu_factor)
}

#HIV only increased mort
mu_risk_HIV_other_h<-array(0, dim = length(HIV_SET))

mu_risk_HIV_other_h_param_extraction<-function(sim_id_current_eval){
  mu_risk_HIV_other_h<-array(0, dim = length(HIV_SET))
  for (h in HIV_SET){
      
      temp_sim_val<-sim_calibration_ref_df%>%
        filter(sim_id == sim_id_current_eval)%>%
        select(paste0("risk.other_",h,"."))
      
      mu_risk_HIV_other_h[h] <- as.numeric(temp_sim_val)
      }
  return(mu_risk_HIV_other_h)
}

#TB/HIV increased mort
mu_risk_HIV_TB_h<-array(0, dim = length(HIV_SET))

mu_risk_HIV_TB_h_param_extraction<-function(sim_id_current_eval){
  mu_risk_HIV_TB_h<-array(0, dim = length(HIV_SET))
  for (h in HIV_SET){
    
    temp_sim_val<-sim_calibration_ref_df%>%
      filter(sim_id == sim_id_current_eval)%>%
      select(paste0("risk.TB_",h,"."))
    
    mu_risk_HIV_TB_h[h] <- as.numeric(temp_sim_val)
  }
  return(mu_risk_HIV_TB_h)
}


#mu_t_h_g calcs- mortality rates#
mort_param_func <-function(yr){
  
  mu_t_h_g <- array(0, dim = c(length(TB_SET),
                               length(HIV_SET), 
                               length(G_SET)))
  
  for (g in G_SET){
    
    if(yr < 1990){
      mort_base_temp <- mort_rate_df%>%
        filter(year == 1990,
               sex == if_else(g == 1, 'Male', 'Female'))
      
      mort_base_temp <- as.numeric(mu_factor*mort_base_temp$val[1])
      
    } else {
      mort_base_temp <- mort_rate_df%>%
        filter(year == yr,
               sex == if_else(g == 1, 'Male', 'Female'))
      
      mort_base_temp <- as.numeric(mu_factor*mort_base_temp$val[1])
    } 
    
    for (h in HIV_SET){
      #first update hiv only values
      for (t in c(1:5, 7)){
        if (h > 1){
          mu_t_h_g[t,h,g]<-mu_risk_HIV_other_h[h]*mort_base_temp
        } else {
          #no hiv no tb
          mu_t_h_g[t,1,g]<-mort_base_temp
        }
      }
      #tb and/or hiv values
      mu_t_h_g[6,h,g]<-mu_risk_HIV_TB_h[h]*mort_base_temp
    }
  }
  return(mu_t_h_g)
}

#alpha_in_t_h_g - Proportion of population that enters each compartment each year
alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), 
                                      length(DR_SET), 
                                      length(HIV_SET), 
                                      length(G_SET)))

aging_in_param_func <-function(yr){
  alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), 
                                        length(DR_SET), 
                                        length(HIV_SET), 
                                        length(G_SET)))
  if(yr < 1990){
    birth_perc_df_temp<-birth_perc_df%>%
    filter(year == 1989)
  } else if (yr >= 2018){
    birth_perc_df_temp<-birth_perc_df%>%
      filter(year == 2017)
  } else {
    birth_perc_df_temp<-birth_perc_df%>%
      filter(year == yr)
    }
  
  for (t in TB_SET){
    for (r in DR_SET){
      for (h in HIV_SET){
        for (g in G_SET){
          birth_perc_df_temp2<-birth_perc_df_temp%>%
            filter(TB_compartment == t,
                   DR_compartment == r,
                   HIV_compartment == h,
                   G_compartment == g)
          
          alpha_in_t_r_h_g[t,r,h,g] <- birth_perc_df_temp2$prop_of_pop
        }
      }
    }
  }
  return(alpha_in_t_r_h_g)
}

#############Pre-processing parameter equations, for ease of use in ode solver#######

#total_out_t_r_h - total amount leaving from compartment#
total_out_t_r_h_g <- array(0, dim = length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET))

#equation H in the appendix
total_out_param_func<-function(mu_t_h_g){
  total_out_t_r_h_g <- array(0, dim = length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET))
  count_temp <- 1
  for (t in TB_SET){
    for (r in DR_SET){
      for (h in HIV_SET){
        for (g in G_SET){
          total_out_t_r_h_g[count_temp] <- (mu_t_h_g[t,h,g]*(1-alpha_out))+
            ((1-mu_t_h_g[t,h,g])*alpha_out)+
            (mu_t_h_g[t,h,g]*alpha_out)
          count_temp <- count_temp + 1
        }
      }
    }
  }
  return(total_out_t_r_h_g)
}



#######EQUATIONS THAT DESCRIBE TB AND HIV PROG IN DESOLVE######
tb_hiv_prog_calibration_model <- function(time, N_t_r_h_g, parms){
  
  #initiation delta array, 3 extra for tracking model outputs
  dN_t_r_h_g <- array(0, dim = (n_compartments+5)) 
  
  #####time varying parameters#####
  #Force of Infection Calculations#
  FOI_1_g <- array(0, dim = length(G_SET))
  FOI_2_g <- array(0, dim = length(G_SET))
  
  for (g in G_SET){
    FOI_1_g[g]<-(beta_g[g]*(sum((phi_h)*N_t_r_h_g[N_t_r_h_g_ref[6, 1, HIV_SET,g]])/
                              sum(N_t_r_h_g[1:n_compartments])))
  }
  
  for (g in G_SET){
    FOI_2_g[g] <-(varepsilon_g[g]*FOI_1_g[g])/(1-varepsilon_g[g])
  }
  
  FOI_r <- c(sum(FOI_1_g), sum(FOI_2_g))
  FOI <- sum(FOI_r)
  
  current_yr<-as.integer(time + start_yr_temp)
  
  #ART availability starts in 2004
  if (current_yr >= 2004){
    eta_i_h_g<<-HIV_transitions_param_func(current_yr, N_t_r_h_g)
  }
  
  #parameters that change yearly
  if(current_yr > last_yr){
    
    #after HIV introduced/before ART initiations (that are calculated at each time step)
    if (current_yr > 1980){
      if (current_yr < 2004){
        eta_i_h_g<<-HIV_transitions_param_func(current_yr, N_t_r_h_g)
      }
    }
    
    #deaths
    mu_t_h_g<<-mort_param_func(current_yr)
    total_out_t_r_h_g<<-total_out_param_func(mu_t_h_g)
    
    #entries prop
    alpha_in_t_r_h_g<<-aging_in_param_func(current_yr)
    
    #only start IPT initiations in 2005
    if(current_yr > 2004){
      kappa_t_h_g<<-kappa_param_extraction_function(current_yr)
    }
    last_yr<<-current_yr
    #print(current_yr)
  }
  #print(kappa_t_h_g[2,4,1])
  
  #####Model outputs#####

  ##health metric calculations
  TB_inc_loc<-n_compartments+1
  TB_mort_loc<-n_compartments+2
  o_mort_loc<-n_compartments+3
  tpt_init_loc<-n_compartments+4
  art_init_loc<-n_compartments+5
  hiv_inc_loc<-n_compartments+6
  
  #TB inc metric #from recent no tpt
  dN_t_r_h_g[TB_inc_loc]<-((theta_h[1]*pi_i_t[2,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[2,DR_SET,1,G_SET]]))+
    (theta_h[2]*pi_i_t[2,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[2,DR_SET,2,G_SET]]))+
    (theta_h[3]*pi_i_t[2,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[2,DR_SET,3,G_SET]]))+
    (theta_h[4]*pi_i_t[2,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[2,DR_SET,4,G_SET]]))+
    #remote no TPT
    (theta_h[1]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,DR_SET,1,G_SET]]))+
    (theta_h[2]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,DR_SET,2,G_SET]]))+
    (theta_h[3]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,DR_SET,3,G_SET]]))+
    (theta_h[4]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,DR_SET,4,G_SET]]))+
    #recent, post-TPT
    (chi*theta_h[4]*pi_i_t[2,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1,4,G_SET]]))+
    #remote, post-TPT
    (chi*theta_h[4]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1,4,G_SET]])))
  
  #TB HIV mort
  #active TB
  dN_t_r_h_g[TB_mort_loc]<-((mu_t_h_g[6,1,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6, DR_SET, 1, 1]]))+ #males HIV-
    (mu_t_h_g[6,1,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6, DR_SET, 1, 2]]))+ #females HIV-
    (mu_t_h_g[6,2,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6, DR_SET, 2, 1]]))+ #males HIV 2
    (mu_t_h_g[6,2,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6, DR_SET, 2, 2]]))+ #females HIV 2
    (mu_t_h_g[6,3,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6, DR_SET, 3, 1]]))+ #males HIV 3
    (mu_t_h_g[6,3,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6, DR_SET, 3, 2]]))+ #females HIV 3
    (mu_t_h_g[6,4,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6, DR_SET, 4, 1]]))+ #males HIV 4
    (mu_t_h_g[6,4,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6, DR_SET, 4, 2]])) #females HIV 4
      )
  
  #HIV positive only (can use mortality rates for TB compartment 1 because same)
  #for TB compartment 1, 2....etc. (except for 6)
  dN_t_r_h_g[o_mort_loc]<-(mu_t_h_g[1,2,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 2, 1]]))+ #males HIV 2
    (mu_t_h_g[1,2,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 2, 2]]))+ #females HIV 2
    (mu_t_h_g[1,3,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 3, 1]]))+ #males HIV 3
    (mu_t_h_g[1,3,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 3, 2]]))+ #females HIV 3
    (mu_t_h_g[1,4,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 4, 1]]))+ #males HIV 4
    (mu_t_h_g[1,4,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 4, 2]])) #females HIV 4
  
  #Treatment Rates
  dN_t_r_h_g[tpt_init_loc]<-kappa_t_h_g[1, 4, 1]*sum(N_t_r_h_g[N_t_r_h_g_ref[1, DR_SET, 4, 1]])+ #males from uninfected
    kappa_t_h_g[1, 4, 2]*sum(N_t_r_h_g[N_t_r_h_g_ref[1, DR_SET, 4, 2]])+ #females from uninfected
    kappa_t_h_g[2, 4, 1]*sum(N_t_r_h_g[N_t_r_h_g_ref[2, DR_SET, 4, 1]])+ #males from recent
    kappa_t_h_g[2, 4, 2]*sum(N_t_r_h_g[N_t_r_h_g_ref[2, DR_SET, 4, 2]])+#females from recent
    kappa_t_h_g[3, 4, 1]*sum(N_t_r_h_g[N_t_r_h_g_ref[3, DR_SET, 4, 1]])+ #males from remote
    kappa_t_h_g[3, 4, 2]*sum(N_t_r_h_g[N_t_r_h_g_ref[3, DR_SET, 4, 2]]) #males from remote
  
  #art initiaiton
  dN_t_r_h_g[art_init_loc]<-eta_i_h_g[2, 4, 1]*sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 2, 1]])+ #males from CD4 more
    eta_i_h_g[2, 4, 2]*sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 2, 2]])+ #females from CD4 more
    eta_i_h_g[3, 4, 1]*sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 3, 1]])+ #males from CD4 less
    eta_i_h_g[3, 4, 2]*sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 3, 2]]) #females from CD4 more
  
  #hiv incidence loc
  dN_t_r_h_g[hiv_inc_loc]<-eta_i_h_g[1, 2, 1]*sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 1, 1]])+ #males
    eta_i_h_g[1, 2, 2]*sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 1, 2]]) #females
    
    
  #entries and exits from the population
  B <- sum(total_out_t_r_h_g*N_t_r_h_g[1:n_compartments])
  
  #######TB compartment 1 Equations #########
  #Set DR compartment to 1, since not applicable to drug resistant compartments#
  for (h in HIV_SET){
    for (g in G_SET){
      dN_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]<-(
        sum(alpha_in_t_r_h_g[1,1,h,g]*B) - #entries from births
          (total_out_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) - #exits from aging out and death
          (FOI*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])+ #exits from TB infection 
          (sum(eta_i_h_g[HIV_SET, h, g]*N_t_r_h_g[N_t_r_h_g_ref[1,1,HIV_SET,g]])) - #entries into HIV compartment
          (sum(eta_i_h_g[h,HIV_SET, g])*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) #exit from HIV compartment
      )
    }
  }
  
  #############TB compartment 2 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[2,r,h,g]]<-(
          (alpha_in_t_r_h_g[2,r,h,g]*B) - #entries from births
            (total_out_t_r_h_g[N_t_r_h_g_ref[2,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[2,r,h,g]]) + #exits from aging out and death
            (FOI_r[r]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])+ #infection from compartment 1
            (zeta*FOI_r[r]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,DR_SET,h,g]]))+ #re-infection from compartment 3
            (zeta*FOI_r[r]*N_t_r_h_g[N_t_r_h_g_ref[4,1,h,g]])+ #re-infection from compartment 4
            (zeta*FOI_r[r]*N_t_r_h_g[N_t_r_h_g_ref[5,1,h,g]])+ #re-infection from compartment 5
            (xi*FOI_r[r]*sum(N_t_r_h_g[N_t_r_h_g_ref[7,DR_SET,h,g]])) - #re-infection from compartment 7
            (pi_i_t[2,3]*N_t_r_h_g[N_t_r_h_g_ref[2,r,h,g]]) - #from recent to remote infection
            (gamma_r[r]*kappa_t_h_g[2,h,g]*N_t_r_h_g[N_t_r_h_g_ref[2,r,h,g]]) - #from recent to post-TPT (and adherence)
            (theta_h[h]*pi_i_t[2,6]*N_t_r_h_g[N_t_r_h_g_ref[2,r,h,g]]) + #from recent TB infection to active
            (sum(eta_i_h_g[HIV_SET,h,g]*N_t_r_h_g[N_t_r_h_g_ref[2,r,HIV_SET,g]])) - #entries into HIV compartment
            (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[2,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  #############TB compartment 3 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]<-(
          (alpha_in_t_r_h_g[3,r,h,g]*B) - #entries from births
            (total_out_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #exits from aging out and death
            (pi_i_t[2,3]*N_t_r_h_g[N_t_r_h_g_ref[2,r,h,g]])- #from recent to remote infection
            (zeta*FOI*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]])- #re-infection
            (gamma_r[r]*kappa_t_h_g[3,h,g]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) - #onto TPT (and adherence)
            (theta_h[h]*pi_i_t[3,6]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #from remote LTBI to active
            (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[3,r,HIV_SET,g]])) - #entries into HIV compartment
            (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  #############TB compartment 4 Equations########
  for (g in G_SET){
    for (h in HIV_SET){
      dN_t_r_h_g[N_t_r_h_g_ref[4,1,h,g]] <- (
        (-1*(total_out_t_r_h_g[N_t_r_h_g_ref[4,1,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[4,1,h,g]])) + #exits from aging out and death
          (gamma_r[1]*kappa_t_h_g[2,h,g]*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]) - #entries from recent to post-TPT (DR 1 only)
          (pi_i_t[2,3]*N_t_r_h_g[N_t_r_h_g_ref[4,1,h,g]]) - #exits recent to remote
          (chi*pi_i_t[2,6]*N_t_r_h_g[N_t_r_h_g_ref[4,1,h,g]]) - #exits from TB infection remote, on TPT to active
          (zeta*FOI*N_t_r_h_g[N_t_r_h_g_ref[4,1,h,g]]) + #exits from reinfection
          (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[4,1,HIV_SET,g]])) - #entries into HIV compartment
          (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[4,1,h,g]]) #exit from HIV compartment
      )
    }
  }

  #############TB compartment 5 Equations########
  for (h in HIV_SET){
    for (g in G_SET){
      dN_t_r_h_g[N_t_r_h_g_ref[5,1,h,g]] <-(
        (-1*(total_out_t_r_h_g[N_t_r_h_g_ref[5,1,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[5,1,h,g]])) + #exits from aging out and death
          (gamma_r[1]*kappa_t_h_g[3,h,g]*N_t_r_h_g[N_t_r_h_g_ref[3,1,h,g]]) + #entries from recent to post-TPT (DR 1 only)
          (pi_i_t[2,3]*N_t_r_h_g[N_t_r_h_g_ref[4,1,h,g]]) - #entries recent to remote
          (chi*pi_i_t[3,6]*N_t_r_h_g[N_t_r_h_g_ref[5,1,h,g]]) - #exits from TB infection remote, on TPT to active
          (zeta*FOI*N_t_r_h_g[N_t_r_h_g_ref[5,1,h,g]]) + #exits from reinfection
          (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[5,1,HIV_SET,g]])) - #entries into HIV compartment
          (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[5,1,h,g]]) #exit from HIV compartment
      )
    }
  }
  
  #############TB compartment 6 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]] <- (
          (alpha_in_t_r_h_g[6,r,h,g]*B) - #entries from births
            (total_out_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) + #exits total out
            (theta_h[h]*pi_i_t[2,6]*N_t_r_h_g[N_t_r_h_g_ref[2,r,h,g]]) + #entries from recent LTBI, no TPT to active 
            (theta_h[h]*pi_i_t[3,6]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #entries from remote LTBI, no TPT to active
            (chi*theta_h[h]*pi_i_t[2,6]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) + #entries from recent LTBI, post- TPT to active 
            (chi*theta_h[h]*pi_i_t[3,6]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) + #entries from remote LTBI, post- TPT to active
            (pi_i_t[7,6]*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) - #entries TB relapse rate
            (upsilon_h[h]*pi_i_t[6,7]*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) + #exits from active to recovered
            (sum(eta_i_h_g[HIV_SET, h, g]*N_t_r_h_g[N_t_r_h_g_ref[6,r,HIV_SET,g]])) - #entries into HIV compartment
            (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  #############TB compartment 7 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]] <- (
          (-1*total_out_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) + #total out
            (upsilon_h[h]*pi_i_t[6,7]*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) - #from active to recovered
            (xi*FOI*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) - #exits re-infection from compartment 7
            (pi_i_t[7,6]*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) + #exits relapse rate
            (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[7,r,HIV_SET,g]])) - #entries into HIV compartment
            (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  list(dN_t_r_h_g)
}

########Feed paramters into desolve########
sim_id_current_eval<-0

current_program_eval<-1 
start_yr_temp<-start_yr_warmup
current_yr<-start_yr_warmup
last_yr<-start_yr_warmup-1


#for(region_id_temp in region_id_ref_df$region_id){
  
  region_name_temp <- region_id_ref_df$region_name[region_id_temp]
  init_region_list<-input_data_param_extraction(region_name_temp)
  
  pop_init_df<<-init_region_list[[1]]
  mort_rate_df<<-init_region_list[[2]]
  birth_perc_df<<-init_region_list[[3]]
  ipt_initiation_df<<-init_region_list[[4]]
  art_coverage_df<<-init_region_list[[5]]
  art_prop_eligible_df<<-init_region_list[[6]]
  hiv_incidence_df<<-init_region_list[[7]]
  
  #####N_init - total initial population in each compartment flattened in 1D array for ODE solver #####
  #arrange pop init df, TB-->DR-->HIV-->G
  pop_init_df<-pop_init_df%>%
    arrange(G_compartment)%>%
    arrange(HIV_compartment)%>%
    arrange(DR_compartment)%>%
    arrange(TB_compartment)
  
  #add in calc placeholders
  model_output_names<-c('cum_TB_inc',
                        'cum_TB_mort',
                        'cum_O_mort',
                        'cum_TPT_init',
                        'cum_ART_init',
                        'cum_hiv_inc')
  
  N_init <- pop_init_df$total_pop
  N_init <-c(N_init, rep(0, times = length(model_output_names)))
  
  names(N_init) <- c(pop_init_df$compartment_id, 
                     model_output_names)
  
  region_filtered_accepted_calibration_sets<-accepted_param_sets_ref_df%>%
    filter(region_id == region_id_temp)%>%
    mutate(general_region_id = paste0(general_sim_id, "_", regional_sim_id))
  
  sim_calibration_ref_df_filtered<-sim_calibration_ref_df%>%
    mutate(general_region_id = paste0(general_sim_id, "_", regional_sim_id))%>%
    filter(general_region_id %in% region_filtered_accepted_calibration_sets$general_region_id)
  
  for (sim_temp in sim_calibration_ref_df_filtered$sim_id){
    time_start<-Sys.time()
    #initialize metric data frame
    summarised_eval_metrics_df_all<-data.frame()
    
    current_program_eval <<-1 
    start_yr_temp<<-start_yr_warmup
    current_yr<<-start_yr_warmup
    last_yr<<-start_yr_warmup-1
  
    sim_id_current_eval<-sim_calibration_ref_df$sim_id[sim_temp]
    regional_sim_id<-sim_calibration_ref_df$regional_sim_id[sim_temp]
    general_sim_id<-sim_calibration_ref_df$general_sim_id[sim_temp]
    
    last_yr <<- 1939
    
    ###update calibrated parameters###
    sim_id_current_eval<<-sim_id_current_eval
    beta_g<<-beta_g_param_extraction(sim_id_current_eval)
    phi_h<<-phi_h_param_extraction(sim_id_current_eval)
    varepsilon_g<<-varepsilon_g_param_extraction(sim_id_current_eval)
    xi<<-xi_param_extraction(sim_id_current_eval)
    zeta<<-zeta_param_extraction(sim_id_current_eval)
    pi_i_t<<-pi_i_t_param_extraction(sim_id_current_eval)
    theta_h<<-theta_h_param_extraction(sim_id_current_eval)
    chi<<-chi_param_extraction(sim_id_current_eval)
    upsilon_h<<-upsilon_h_param_extraction(sim_id_current_eval)
    kappa_factor<<-kappa_factor_param_extraction(sim_id_current_eval)
    eta_factor<<-eta_factor_param_extraction(sim_id_current_eval)
    #cannot be above 100% coverage (will be more than 100% coverage at 1.2 under program 3 for women)
    sigma_factor<<-sigma_factor_param_extraction(sim_id_current_eval)
    sigma_factor<<-if_else(sigma_factor >= 1.204, 1.204, sigma_factor)
    mu_factor<<-mu_factor_param_extraction(sim_id_current_eval)
    mu_risk_HIV_other_h<<-mu_risk_HIV_other_h_param_extraction(sim_id_current_eval)
    mu_risk_HIV_TB_h<<-mu_risk_HIV_TB_h_param_extraction(sim_id_current_eval)
    
    ##reset ipt initiations and hiv to 0
    eta_i_h_g<<-array(0, dim=c(length(HIV_SET), 
                               length(HIV_SET), 
                               length(G_SET)))
    kappa_t_h_g<<-array(data = 0, dim=c(length(TB_SET), 
                                        length(HIV_SET),
                                        length(G_SET)))
    
    ####RUN MODEL#####
    out_df<-as.data.frame(ode(times = TT_SET_warmup, y = N_init,
                              func = tb_hiv_prog_calibration_model, method = 'lsoda',
                              parms = NULL))
    
    #get init values for start of each program run
    N_init_start_of_eval<-out_df%>%
      filter(time == max(time))%>%
      select(-c('time'))
    
    ##get metrics
    out_df$tb_prev<-rowSums(out_df[,1+c(N_t_r_h_g_ref[6, DR_SET,HIV_SET,G_SET])])
  
    #calculate hiv prev
    out_df$hiv_prev<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),G_SET])]) 
    out_df$hiv_prev_t_male_per_100K_ppl<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),1])]) 
    out_df$hiv_prev_t_female_per_100K_ppl<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),2])])
    
    #calculate n on ART prev
    out_df$n_art_t_male_per_100K_ppl<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, 4,1])])/out_df$hiv_prev_t_male_per_100K_ppl
    out_df$n_art_t_female_per_100K_ppl<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, 4,2])])/out_df$hiv_prev_t_female_per_100K_ppl
    
    out_df$n_art<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET,4,G_SET])])
    
    #initialize summarised eval metrics
    summarised_eval_metrics_df_all<-out_df%>%
      mutate(year = as.integer(time + start_yr_warmup))%>%
      mutate(TB_inc_per_100k_ppl = lead(cum_TB_inc)-cum_TB_inc,
             TB_mort_per_100k_ppl = lead(cum_TB_mort)-cum_TB_mort,
             O_mort_per_100K_ppl = lead(cum_O_mort)-cum_O_mort,
             TPT_init_per_100K_ppl = lead(cum_TPT_init)-cum_TPT_init,
             ART_init_per_100K_ppl = lead(cum_ART_init)-cum_ART_init,
             HIV_inc_per_100k_ppl = lead(cum_hiv_inc)-cum_hiv_inc)%>%
      group_by(year)%>%
      summarise(TB_inc_per_Y_100k_ppl = sum(TB_inc_per_100k_ppl),
                TB_mort_per_Y_100k_ppl = sum(TB_mort_per_100k_ppl),
                O_mort_per_Y_100k_ppl = sum(O_mort_per_100K_ppl),
                TPT_init_Y_per_100K_ppl = sum(TPT_init_per_100K_ppl),
                ART_init_Y_per_100K_ppl = sum(ART_init_per_100K_ppl),
                HIV_inc_Y_per_100k_ppl = sum(HIV_inc_per_100k_ppl),
                HIV_prev_Y_100k_ppl = mean(hiv_prev),
                ART_coverage_Y_100K_males = mean(n_art_t_male_per_100K_ppl),
                ART_coverage_Y_100K_females = mean(n_art_t_female_per_100K_ppl),
                TB_prev_Y_100K_ppl = mean(tb_prev),
                avg_on_ART_Y_per_100K_ppl = mean(n_art))%>%
      filter(year >= 1990,
             year <= 2017)%>%
      mutate(general_sim_id = general_sim_id,
             regional_sim_id = regional_sim_id,
             program_id = 1)
    
    for (p in program_sets){
      current_program_eval<<-p
      #print(current_program_eval)
      start_yr_temp<<-start_yr_eval #for year calculations
      current_yr <<- start_yr_eval #tracks for dynamic params
      last_yr<<-start_yr_eval-1 #the year before eval period updates
      
      out_df2<-as.data.frame(ode(times = TT_SET_eval, y = unlist(N_init_start_of_eval),
                                 func = tb_hiv_prog_calibration_model, method = 'lsoda',
                                 parms = NULL))
      
      ##get metrics
      #calculate hiv prev
      out_df2$hiv_prev<-rowSums(out_df2[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),G_SET])]) 
      out_df2$hiv_prev_t_male_per_100K_ppl<-rowSums(out_df2[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),1])]) 
      out_df2$hiv_prev_t_female_per_100K_ppl<-rowSums(out_df2[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),2])])
      
      #calculate n on ART prev
      out_df2$n_art_t_male_per_100K_ppl<-rowSums(out_df2[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, 4,1])])/out_df2$hiv_prev_t_male_per_100K_ppl
      out_df2$n_art_t_female_per_100K_ppl<-rowSums(out_df2[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, 4,2])])/out_df2$hiv_prev_t_female_per_100K_ppl
      
      out_df2$hiv_prev<-rowSums(out_df2[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),G_SET])])
      out_df2$art_coverage<-rowSums(out_df2[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET,4,G_SET])])/out_df2$hiv_prev
      out_df2$tb_prev<-rowSums(out_df2[,1+c(N_t_r_h_g_ref[6, DR_SET,HIV_SET,G_SET])])
      
      out_df2$n_art<-rowSums(out_df2[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET,4,G_SET])])
      
      summarised_eval_metrics_df_intervention<-out_df2%>%
        mutate(year = as.integer(time + start_yr_eval))%>%
        mutate(TB_inc_per_100k_ppl = lead(cum_TB_inc)-cum_TB_inc,
               TB_mort_per_100k_ppl = lead(cum_TB_mort)-cum_TB_mort,
               O_mort_per_100K_ppl = lead(cum_O_mort)-cum_O_mort,
               TPT_init_per_100K_ppl = lead(cum_TPT_init)-cum_TPT_init,
               ART_init_per_100K_ppl = lead(cum_ART_init)-cum_ART_init,
               HIV_inc_per_100k_ppl = lead(cum_hiv_inc)-cum_hiv_inc)%>%
        group_by(year)%>%
        summarise(TB_inc_per_Y_100k_ppl = sum(TB_inc_per_100k_ppl),
                  TB_mort_per_Y_100k_ppl = sum(TB_mort_per_100k_ppl),
                  O_mort_per_Y_100k_ppl = sum(O_mort_per_100K_ppl),
                  TPT_init_Y_per_100K_ppl = sum(TPT_init_per_100K_ppl),
                  ART_init_Y_per_100K_ppl = sum(ART_init_per_100K_ppl),
                  HIV_inc_Y_per_100k_ppl = sum(HIV_inc_per_100k_ppl),
                  HIV_prev_Y_100k_ppl = mean(hiv_prev),
                  ART_coverage_Y_100K_males = mean(n_art_t_male_per_100K_ppl),
                  ART_coverage_Y_100K_females = mean(n_art_t_female_per_100K_ppl),
                  TB_prev_Y_100K_ppl = mean(tb_prev),
                  avg_on_ART_Y_per_100K_ppl = mean(n_art))%>%
        filter(year <= 2027)%>%
        mutate(general_sim_id = general_sim_id,
               regional_sim_id = regional_sim_id,
               program_id = p)
      
      
      summarised_eval_metrics_df_all<-rbind(summarised_eval_metrics_df_all,
                                            summarised_eval_metrics_df_intervention)
    }
    
    setwd(paste0(outdir, '/', region_name_temp))
    write.csv(summarised_eval_metrics_df_all, 
              file = paste0('summarised_eval_metrics_df_',
                            region_id_temp, '_',
                            sim_id_current_eval,
                            '.csv'), row.names = FALSE)
    
    completion_time<-Sys.time()-time_start
    print(paste0(region_name_temp, ":", general_sim_id, 
                 "regional sim id", regional_sim_id))
    print(completion_time)
  }
#}
