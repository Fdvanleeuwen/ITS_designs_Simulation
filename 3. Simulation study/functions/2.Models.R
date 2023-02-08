# run models 
run_OLS <- function(dd_long){
  # The function runs a OLS segmented regression model and outputs the beta's and t-values 
  #of the intervention dummy and the interaction between this dummy and time. 
  
  # OLS regression
  M1 <- lm(score ~ time*treatment, data = dd_long) 

  # obtain coefficients
  beta_1 <- M1$coefficients[2]
  beta_2 <- M1$coefficients[3]
  beta_3 <- M1$coefficients[4]
  
  # Obtain the t-values
  t_value_pre_slope <- summary(M1)$coefficients[2,3]
  t_value_step_M1 <- summary(M1)$coefficients[3,3]
  t_value_slope_M1 <- summary(M1)$coefficients[4,3]
  
  if (abs(t_value_pre_slope) > 1.96){
    significant_pre_slope <- 1
  }
  else{
    significant_pre_slope <- 0
  }
  
  if (abs(t_value_step_M1) > 1.96){
    significant_step <- 1
  }
  else{
    significant_step <- 0
  }
  
  if (abs(t_value_slope_M1) > 1.96){
    significant_slope <- 1
  }
  else{
    significant_slope <- 0
  }
  
  return(list(significant_pre_slope = significant_pre_slope, significant_step = significant_step,
              significant_slope = significant_slope, beta_1 = beta_1, beta_2 = beta_2,
              beta_3 = beta_3))
}

# run models 
run_ML <- function(dd_long){
  # The function runs a ML segmented regression model and outputs the beta's and t-values 
  #of the intervention dummy and the interaction between this dummy and time. 
  
  # ML regression
  M2 <- lmer(score ~ time*treatment + (1|id), data = dd_long )
  
  # obtain coefficients
  beta_1 <- coef(summary(M2))[,"Estimate"][2]
  beta_2 <- coef(summary(M2))[,"Estimate"][3]
  beta_3 <- coef(summary(M2))[,"Estimate"][4]
  
  # Obtain the t-values
  t_value_pre_slope <- coef(summary(M2))[,"t value"][2]
  t_value_step_M2 <- coef(summary(M2))[,"t value"][3]
  t_value_slope_M2 <- coef(summary(M2))[,"t value"][4]
  
  if (abs(t_value_pre_slope) > 1.96){
    significant_pre_slope <- 1
  }
  else{
    significant_pre_slope <- 0
  }
  
  if (abs(t_value_step_M2) > 1.96){
    significant_step <- 1
  }
  else{
    significant_step <- 0
  }
  
  if (abs(t_value_slope_M2) > 1.96){
    significant_slope <- 1
  }
  else{
    significant_slope <- 0
  }
  
  return(list(significant_pre_slope = significant_pre_slope, significant_step = significant_step,
              significant_slope = significant_slope, beta_1 = beta_1, beta_2 = beta_2,
              beta_3 = beta_3))
}