sim_all <- function(N_person, N_time, N_sim, var, Intercept, Slope, step_size, 
                    slope_size, seed, var_time, model = c("OLS, ML, SEM"), cov = FALSE){
  # The function takes in all arguments needed for the simulation study and outputs the power/bias.
  set.seed(seed)
  significant_pre_slope <- c()
  significant_step <- c()
  significant_slope <- c()
  bias_precent_b1 <- c()
  bias_precent_b2 <- c()
  bias_precent_b3 <- c()
  N_t <- c()
  N_pers <- c()
  
  
  for (l in 1:length(N_time)){
    significant_pre_slope_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
    significant_step_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
    significant_slope_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
    bias_precent_b1_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
    bias_precent_b2_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
    bias_precent_b3_sub <- matrix(NA,nrow=length(N_person),ncol=N_sim)
    N_pers_sub <- rep(NA, length(N_person))
    
    for (p in 1:length(N_person)){
      
      for (s in 1:N_sim){
        forms <- forms_for_sim(N_time[l], var, Intercept, Slope, step_size, slope_size)
        data <- gen_data(forms$form, forms$v_name, forms$var, N_person[p], 
                         forms$Intercept, var_time)
        if (model == "ML"){
          power <- run_ML(data$long)
        }
        else if (model == "SEM"){
          power <- run_SEM(data$wide, forms$v_name, cov = cov)
        }
        else{
          power <- run_OLS(data$long)
        }
        significant_pre_slope_sub [p,s] <- power$significant_pre_slope
        significant_step_sub[p,s] <- power$significant_step
        significant_slope_sub[p,s] <- power$significant_slope
        bias_precent_b1_sub[p,s] <- (power$beta_1 - Slope)/Slope * 100
        bias_precent_b2_sub[p,s] <- (power$beta_2 - step_size + 0.5 * slope_size)/step_size * 100
        bias_precent_b3_sub[p,s] <- (power$beta_3 - slope_size)/slope_size * 100
        N_pers_sub[p] <- N_person[p]
      }
      significant_pre_slope2 <- rowMeans(significant_pre_slope_sub)
      significant_step2 <- rowMeans(significant_step_sub)
      significant_slope2 <- rowMeans(significant_slope_sub)
      bias_precent_b1_2 <- rowMeans(bias_precent_b1_sub)
      bias_precent_b2_2 <- rowMeans(bias_precent_b2_sub)
      bias_precent_b3_2 <- rowMeans(bias_precent_b3_sub)
    }
    significant_pre_slope <- append(significant_pre_slope, significant_pre_slope2)
    significant_step <- append(significant_step, significant_step2)
    significant_slope <- append(significant_slope, significant_slope2)
    bias_precent_b1 <- append(bias_precent_b1, bias_precent_b1_2)
    bias_precent_b2 <- append(bias_precent_b2, bias_precent_b2_2)
    bias_precent_b3 <- append(bias_precent_b3, bias_precent_b3_2)
    N_t <- append(N_t,rep(N_time[l], length(N_person)))
    N_pers <- append(N_pers, N_person)
    
  }
  
  power <- data.frame(N_t = N_t,
                      N_pers = N_pers,
                      pre_slope = significant_pre_slope,
                      step = significant_step, 
                      slope = significant_slope,
                      bias_pre_slope = bias_precent_b1,
                      bias_step = bias_precent_b2,
                      bias_slope = bias_precent_b3)
  return(power)
}