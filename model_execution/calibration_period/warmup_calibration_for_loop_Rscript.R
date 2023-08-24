#last updated August 7, 2023
#runs warmup from 1940 to end of 2017 (beg of 2018)
#for a specific location
#writes csv for model outputs in calibration years for a specific location
#metrics 

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle'), require, character.only=T)

library(R.utils)

#change to only run remaining general parameter sets
itr_id<-1

#region_name_temp<-"Eastern Cape"
#region_name_temp<-"Free State"
#region_name_temp<-"Gauteng"
#region_name_temp<-"KwaZulu-Natal"
region_name_temp<-"Limpopo"
#region_name_temp<-"Mpumalanga"
#region_name_temp<-"Northern Cape"
#region_name_temp<-"North-West"
#region_name_temp<-"Western Cape"

# # ####HYAK OR GITHUB SPECIFIC CODES TO COMMENT/UNCOMMENT####
#hyak specific code

#set # of args to # of general parameter sets

#args = commandArgs(trailingOnly=TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# } else if (length(args)>1) {
#   stop("only need one parameter", call.=FALSE)
 #}

#args_temp<- as.integer(args[1])

# # #location where input parameters are
#indir<-'/gscratch/icrc/cgreene3/SA_resource_allocation/input_parameters'
# # #
# # # #location where want outputs
#outdir<-paste0('/gscratch/icrc/cgreene3/SA_resource_allocation/calibration_outputs/', 
#               region_name_temp)

# # # #if running from github
library(here)
#if running from local set env to resource_allocation_HIV_TB_SA
args_temp<-259

#setwd(paste0(here(), '/test/calibration_results/'))
#general_accepted_calibration_sets_ref_df<-read.csv(paste0('general_accepted_calibration_sets_ref_df',
#                 itr_id-1, '.csv'))

#location where input parameters are
indir<-paste0(here(), '/param_files/input_parameters')

#location where want outputs
outdir<-paste0(here(), '/test/calibration_outputs/', region_name_temp)

########Time Horizon and Evaluation intervals (1 month)#####
start_yr <- 1940
current_yr <- start_yr #tracks for dynamic params
last_year <- current_yr - 1
end_yr <- 2018
TT<-end_yr-start_yr
time_interval <- 1/12
TT_SET <- seq(from = 0, to = TT, by = time_interval) #tau

#these data frames change depending on combin of calib parameters
setwd(indir)

sim_calibration_ref_df<-read.csv('calibration_sets_df.csv')

setwd(paste0(indir, '/general_param_sets_itr'))
#change depending on first or second calibration run
if(itr_id > 1){
  general_accepted_calibration_sets_ref_df<-read.csv(paste0('general_accepted_calibration_sets_ref_df',
                                                            itr_id-1, '.csv'))
  remaining_general_parameters<-general_accepted_calibration_sets_ref_df$general_parameter_id
} else{
  remaining_general_parameters<-unique(sim_calibration_ref_df$general_sim_id)
}

general_parameter<-remaining_general_parameters[args_temp]

setwd(indir)
sim_calibration_ref_df<-sim_calibration_ref_df%>%
  filter(general_sim_id == general_parameter)
  #filter(sim_id >= sim_id_start,
   #      sim_id <= sim_id_end)
pop_init_df <- read.csv('pop_init_df_1940.csv')%>%
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


################ Define sets ###############################

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
pop_init_df<-pop_init_df%>%
  arrange(G_compartment)%>%
  arrange(HIV_compartment)%>%
  arrange(DR_compartment)%>%
  arrange(TB_compartment)

#initialize population values
N_init <- pop_init_df$total_pop
names(N_init)<-pop_init_df$compartment_id

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
    
    IPT_init_temp<-ipt_initiation_df%>%
      filter(gender_id == g,
             year == yr,
             program_id == 1)
    
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
             year == yr, program_id == 1)
    
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
      art_coverage_temp<-art_coverage_df%>%
        filter(sex == if_else(g == 1, 'Male', 'Female'),
               year == yr,
               program_id == 1)
      
      #adjusted with calibration factor
      art_coverage_temp<-sigma_factor*(art_coverage_temp$art_coverage)
      
      #get percent eligible (time varying) 
      prop_eligible_temp<-art_prop_eligible_df%>%
        filter(gender_id == g,
               year == yr)
      prop_eligible<-prop_eligible_temp$percent_hiv_compartment_2_eligible
      
      
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
  
  current_time<-Sys.time()
  #print(current_time-time_start)
  # if condition with break
  if(difftime(current_time, time_start, units = "secs" )>max_seconds) {
    setwd(outdir)
    write.csv(calib_metrics_df_all_df, file = paste0(region_name_temp, 
                                                     '_calib_metrics_df_general_sim_id',
                                                     args_temp, 'endat', regional_sim_id,
                                                     '.csv'), row.names = FALSE)
    break
  }
  
  #initiation delta array, 3 extra for tracking model outputs
  dN_t_r_h_g <- array(0, dim = (n_compartments+3)) 
  
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
  
  #write year parameter for parameters that change over time
  current_yr <-as.integer(start_yr+time)
  #print(current_yr)
  
  #HIV transitions start in 1980
  if (current_yr >= 2004){
    eta_i_h_g<<-HIV_transitions_param_func(current_yr, N_t_r_h_g)
  }
  
  #parameters that change yearly
  if(current_yr > last_year){
    
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
    last_year<<-current_yr
  }
  
  #####Model outputs#####

  ##health metric calculations
  TB_inc_loc<-n_compartments+1
  TB_mort_loc<-n_compartments+2
  o_mort_loc<-n_compartments+3
  
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
  dN_t_r_h_g[o_mort_loc]<-((mu_t_h_g[1,2,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 2, 1]]))+ #males HIV 2
    (mu_t_h_g[1,2,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 2, 2]]))+ #females HIV 2
    (mu_t_h_g[1,3,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 3, 1]]))+ #males HIV 3
    (mu_t_h_g[1,3,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 3, 2]]))+ #females HIV 3
    (mu_t_h_g[1,4,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 4, 1]]))+ #males HIV 4
    (mu_t_h_g[1,4,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1:5, 7), DR_SET, 4, 2]])) #females HIV 4
  )
  
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
setwd(outdir)
sim_id_current_eval<-0
calib_metrics_df_all_df<-data.frame()
#out_df_all<-data.frame() #to test
n_sims<-1:nrow(sim_calibration_ref_df)
max_seconds<-300

for(sim_temp in n_sims){
  
  time_start<-Sys.time()
  
    sim_id_current_eval<-sim_calibration_ref_df$sim_id[sim_temp]
    regional_sim_id<-sim_calibration_ref_df$regional_sim_id[sim_temp]
    
    print(regional_sim_id)
    last_year <<- 1939
    
    #####N_init - total initial population in each compartment flattened in 1D array for ODE solver #####
    #arrange pop init df, TB-->DR-->HIV-->G
    pop_init_df<-pop_init_df%>%
      arrange(G_compartment)%>%
      arrange(HIV_compartment)%>%
      arrange(DR_compartment)%>%
      arrange(TB_compartment)
    
    N_init <- pop_init_df$total_pop
    
    #add in mort calc placeholders
    model_output_names<-c('cum_TB_inc',
                          'cum_TB_mort',
                          'cum_O_mort')
    
    
    N_init <-c(N_init, rep(0, times = length(model_output_names)))
    
    names(N_init) <- c(pop_init_df$compartment_id, 
                       model_output_names)
    
    
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
    sigma_factor<<-sigma_factor_param_extraction(sim_id_current_eval)
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
    out_df<-as.data.frame(ode(times = TT_SET, y = N_init,
                              func = tb_hiv_prog_calibration_model, method = 'lsoda',
                              parms = NULL))
    
    out_df$hiv_prev<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),G_SET])])
    
    out_df$ART_coverage<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET,4,G_SET])])/out_df$hiv_prev
    
    out_df$tb_prev<-rowSums(out_df[,1+c(N_t_r_h_g_ref[6, DR_SET,HIV_SET,G_SET])])
    
    
    
    out_df<- cbind(year = as.integer(start_yr+out_df$time), sim_id = rep(sim_id_current_eval,
                                                                         times = nrow(out_df)),
                   out_df)
    
    #######summarise calibration metrics######
    subset_out_df<-out_df%>%
      select(c('year', 'time', 
               model_output_names))
    
    model_output_names<-c('cum_TB_inc',
                          'cum_TB_mort',
                          'cum_O_mort')
    
    summarised_metrics<-out_df%>%
      mutate(TB_inc_per_100k_ppl = cum_TB_inc-lag(cum_TB_inc),
             TB_mort_per_100k_ppl = cum_TB_mort-lag(cum_TB_mort),
             O_mort_per_100K_ppl = cum_O_mort - lag(cum_O_mort))%>%
      group_by(year)%>%
      summarise(TB_inc_per_Y_100k_ppl = sum(TB_inc_per_100k_ppl),
                TB_mort_per_Y_100k_ppl = sum(TB_mort_per_100k_ppl),
                O_mort_per_Y_100k_ppl = sum(O_mort_per_100K_ppl),
                hiv_prev_Y_100k_ppl = mean(hiv_prev),
                ART_coverage_Y_100K_ppl = mean(ART_coverage),
                TB_prev_Y_100K_ppl = mean(tb_prev))%>%
      #filter(year >= 1990,
      #       year <= 2017)
      filter(year == 2005 | year == 2017)
    
    
    summarised_calib_metrics_df<-summarised_metrics%>%
      mutate(general_sim_id = args_temp,
             regional_sim_id = regional_sim_id, 
             sim_id = sim_id_current_eval)
    
    if(regional_sim_id == 1){
      calib_metrics_df_all_df<<-summarised_calib_metrics_df
      #out_df_all<<-out_df%>%
      #  filter(year >= 1990,
      #         year <= 2017)
    } else {
      calib_metrics_df_all_df<<-rbind(calib_metrics_df_all_df,summarised_calib_metrics_df)
      #out_df_all<<-rbind(out_df_all, out_df%>%
      #                     filter(year >= 1990,
      #                            year <= 2017))
    }
    completion_time<-Sys.time()-time_start
    print(completion_time)
  #}
}

setwd(outdir)
write.csv(calib_metrics_df_all_df, file = paste0(region_name_temp, 
                                                 '_calib_metrics_df_general_sim_id',
                                                 args_temp,
                                                 '.csv'), row.names = FALSE)
#write.csv(out_df_all, file = 'out_df_all.csv', row.names = FALSE)