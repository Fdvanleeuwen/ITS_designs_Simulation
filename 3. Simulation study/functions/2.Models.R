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

# growth model
run_SEM <- function(wide, v_name){
  # Latent growth curve model
  v_name <- X$v_name
  form_i <- c("i =~")
  form_s <- c("s =~")
  
  time = seq(1,length(v_name)+1,1)-(length(v_name)+1)/2
  
  for (i in 2:(length(v_name)+1)){
    
    if (i != length(v_name)+1){
      form_i[i] =  paste("1*",v_name[i-1],"+")
      form_s[i] = paste(time[i-1],"*",v_name[i-1], "+")
      
    }
    else{
      form_i[i] =  paste("1*",v_name[i-1],"\n", sep = "")
      form_s[i] = paste(time[i-1],"*",v_name[i-1],sep = "")
    }
  }
  form_i_c = paste(form_i,collapse=" ")
  form_s_c = paste(form_s,collapse=" ")
  
  final_lgm = paste(c(form_i_c, form_s_c),collapse=" ")
  
  # ML regression
  M3 <- growth(final_lgm, data=test$wide)
  M3_sum <- summary(GCM1_fit)
  
  # obtain coefficients
  #beta_1 <- 
  #beta_2 <- 
  #beta_3 <- 
  
  # Obtain the t-values
  #t_value_pre_slope <- 
  #t_value_step_M2 <- 
  #t_value_slope_M2 <- 
  
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