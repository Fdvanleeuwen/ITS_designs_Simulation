sim_all <- function(N_person, N_time, N_sim, var, Intercept, Slope, seed, var_time, effect_sizes, treatment_step, treatment_slope, additional_step, additional_step_corr, time_2, model = c("OLS, ML, SEM"), cov = FALSE){
  # The function takes in all arguments needed for the simulation study and outputs the power/bias.
  #set.seed(seed)
  
  All_conditions <- All_conditions(effect_sizes, Slope)
  step_conitions <- All_conditions$step_size_all
  slope_conitions <- All_conditions$slope_size_all

  significant_pre_slope_all <- c()
  significant_step_all <- c()
  significant_slope_all <- c()
  significant_step2_all <- c()
  bias_precent_b1_all <- c()
  bias_precent_b2_all <- c()
  bias_precent_b3_all <- c()
  bias_precent_b4_all <- c()
  Precision_b1_all <- c()
  Precision_b2_all <- c()
  Precision_b3_all <- c()
  Precision_b4_all <- c()
  N_t_all <- c()
  N_pers_all <- c()
  effect_sizes_all <- c()
  
  for (EV in 1:length(step_conitions)){
    significant_pre_slope <- c()
    significant_step <- c()
    significant_slope <- c()
    significant_step__2 <- c()
    bias_precent_b1 <- c()
    bias_precent_b2 <- c()
    bias_precent_b3 <- c()
    bias_precent_b4 <- c()
    Precision_b1 <- c()
    Precision_b2 <- c()
    Precision_b3 <- c()
    Precision_b4 <- c()
    N_t <- c()
    N_pers <- c()
    
    step = step_conitions[EV]*treatment_step
    slope = slope_conitions[EV]*treatment_slope
    additional_step_use = step*0.5*additional_step

    for (l in 1:length(N_time)){
      significant_pre_slope_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      significant_step_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      significant_slope_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      significant_step2_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      bias_precent_b1_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      bias_precent_b2_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      bias_precent_b3_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      bias_precent_b4_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      Precision_b1_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      Precision_b2_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      Precision_b3_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      Precision_b4_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
      N_pers_sub <- rep(NA, length(N_person))
      
      for (p in 1:length(N_person)){
        
        for (s in 1:N_sim){
          forms <- forms_for_sim(N_time[l], var, Intercept, Slope, step, slope, additional_step_use)
          data <- gen_data(forms$form, forms$v_name, forms$var, N_person[p], 
                           forms$Intercept, var_time, additional_step_corr)
          if (model == "ML"){
            power <- run_ML(data$long)
          }
          else if (model == "SEM"){
            power <- run_SEM(data$wide, forms$v_name, cov = cov)
          }
          else{
            power <- run_OLS(data$long, additional_step_corr, time_2)
          }
          
          significant_pre_slope_sub [p,s] <- power$significant_pre_slope
          significant_step_sub[p,s] <- power$significant_step
          significant_slope_sub[p,s] <- power$significant_slope
          significant_step2_sub[p,s] <- power$significant_step2
          
          
          bias_precent_b1_sub[p,s] <- (power$beta_1 - Slope)
          
          if (treatment_step == 0){
            bias_precent_b2_sub[p,s] <- power$beta_2
          }
          else{
            bias_precent_b2_sub[p,s] <- power$beta_2 - step
          }
          
          if (treatment_slope == 0){
            bias_precent_b3_sub[p,s] <- power$beta_3
          }
          else{
            bias_precent_b3_sub[p,s] <- power$beta_3 - slope
          }
          
          if (additional_step_corr == 1){ 
                bias_precent_b4_sub[p,s] <-  power$beta_4 - additional_step_use - step
          }
          else{
            bias_precent_b4_sub[p,s] <-  0
          }
          N_pers_sub[p] <- N_person[p]

        }
        significant_pre_slope2 <- rowMeans(significant_pre_slope_sub)
        significant_step2 <- rowMeans(significant_step_sub)
        significant_slope2 <- rowMeans(significant_slope_sub)
        significant_step_2 <- rowMeans(significant_step2_sub)
        bias_precent_b1_2 <- rowMeans(bias_precent_b1_sub)
        bias_precent_b2_2 <- rowMeans(bias_precent_b2_sub)
        bias_precent_b3_2 <- rowMeans(bias_precent_b3_sub)
        bias_precent_b4_2 <- rowMeans(bias_precent_b4_sub)
        Precision_b1_sub_2 <- (apply(bias_precent_b1_sub, 1, sd, na.rm=TRUE) / sqrt(N_person[p]))
        Precision_b2_sub_2 <- (apply(bias_precent_b2_sub, 1, sd, na.rm=TRUE) / sqrt(N_person[p]))
        Precision_b3_sub_2 <- (apply(bias_precent_b3_sub, 1, sd, na.rm=TRUE) / sqrt(N_person[p]))
        Precision_b4_sub_2 <- (apply(bias_precent_b4_sub, 1, sd, na.rm=TRUE) / sqrt(N_person[p]))
      }
      significant_pre_slope <- append(significant_pre_slope, significant_pre_slope2)
      significant_step <- append(significant_step, significant_step2)
      significant_slope <- append(significant_slope, significant_slope2)
      significant_step__2 <- append(significant_step__2, significant_step_2)
      bias_precent_b1 <- append(bias_precent_b1, bias_precent_b1_2)
      bias_precent_b2 <- append(bias_precent_b2, bias_precent_b2_2)
      bias_precent_b3 <- append(bias_precent_b3, bias_precent_b3_2)
      bias_precent_b4 <- append(bias_precent_b4, bias_precent_b4_2)
      Precision_b1 <- append(Precision_b1, Precision_b1_sub_2)
      Precision_b2 <- append(Precision_b2, Precision_b2_sub_2)
      Precision_b3 <- append(Precision_b3, Precision_b3_sub_2)
      Precision_b4 <- append(Precision_b4, Precision_b4_sub_2)
      N_t <- append(N_t,rep(N_time[l], length(N_person)))
      N_pers <- append(N_pers, N_person)
      
    }
    significant_pre_slope_all <- append(significant_pre_slope_all, significant_pre_slope)
    significant_step_all <- append(significant_step_all, significant_step)
    significant_slope_all <- append(significant_slope_all, significant_slope)
    significant_step2_all <- append(significant_step2_all, significant_step__2)
    bias_precent_b1_all <- append(bias_precent_b1_all, bias_precent_b1)
    bias_precent_b2_all <- append(bias_precent_b2_all, bias_precent_b2)
    bias_precent_b3_all <- append(bias_precent_b3_all, bias_precent_b3)
    bias_precent_b4_all <- append(bias_precent_b4_all, bias_precent_b4)
    Precision_b1_all <- append(Precision_b1_all,Precision_b1)
    Precision_b2_all <- append(Precision_b2_all,Precision_b2)
    Precision_b3_all <- append(Precision_b3_all,Precision_b3)
    Precision_b4_all <- append(Precision_b4_all,Precision_b4)
    N_t_all <- append(N_t_all, N_t)
    N_pers_all <- append(N_pers_all, N_pers)
    effect_sizes_all <- append(effect_sizes_all, rep(effect_sizes[EV],length(N_time)*length(N_person)))
  }

  power <- data_frame(effect_sizes = effect_sizes_all,
                      N_t = N_t_all,
                      N_pers = N_pers_all,
                      pre_slope = significant_pre_slope_all,
                      step = significant_step_all, 
                      slope = significant_slope_all,
                      step2 = significant_step2_all,
                      bias_pre_slope = bias_precent_b1_all,
                      bias_step = bias_precent_b2_all,
                      bias_slope = bias_precent_b3_all,
                      bias_step2 =bias_precent_b4_all,
                      Precision_pre_slope = Precision_b1_all,
                      Precision_step = Precision_b2_all,
                      Precision_slope = Precision_b3_all,
                      Precision_step2 = Precision_b4_all)
  return(power)
}