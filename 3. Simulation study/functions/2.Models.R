# run models 
run_OLS <- function(dd_long, additional_step_corr){
  # The function runs a OLS segmented regression model and outputs the beta's and t-values 
  #of the intervention dummy and the interaction between this dummy and time. 
  
  # OLS regression
  M1 <- lm(score ~ time*treatment, data = dd_long) 

  # obtain coefficients
  if (additional_step_corr == 0){
    beta_1 <- M1$coefficients[2]
    beta_2 <- M1$coefficients[3]
    beta_3 <- M1$coefficients[4] 
    # Obtain the t-values
    t_value_pre_slope <- summary(M1)$coefficients[2,3]
    t_value_step_M1 <- summary(M1)$coefficients[3,3]
    t_value_slope_M1 <- summary(M1)$coefficients[4,3]
  }
  else {
      beta_1 <- M1$coefficients[2]
      beta_2 <- M1$coefficients[3]
      beta_3 <- M1$coefficients[5] 
      # Obtain the t-values
      t_value_pre_slope <- summary(M1)$coefficients[2,3]
      t_value_step_M1 <- summary(M1)$coefficients[3,3]
      t_value_slope_M1 <- summary(M1)$coefficients[5,3]
  }

  
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
run_SEM <- function(wide, v_name, cov = FALSE){
  
  # Latent growth curve model
  form_i <- c("i =~")
  form_s <- c("s =~")
  form_i2 <- c("i2 =~")
  form_s2 <- c("s2 =~")
  
  # create indicators for time and step size
  time = seq(0,length(v_name),1)
  step_indc = c(rep(0, length(v_name)/2), rep(1, length(v_name)/2))
  
  for (g in 2:(length(v_name)+1)){
    
    if (g != length(v_name)+1){
      form_i[g] =  paste("1*",v_name[g-1],"+", sep = "")
      form_s[g] = paste(time[g-1],"*",v_name[g-1], "+",sep = "")
      form_i2[g] = paste(step_indc[g-1],"*",v_name[g-1],"+", sep = "")
    }
    else{
      form_i[g] =  paste("1*",v_name[g-1],"\n", sep = "")
      form_s[g] = paste(time[g-1],"*",v_name[g-1],"\n",sep = "")
      form_i2[g] = paste(step_indc[g-1],"*",v_name[g-1],"\n", sep = "")
    }
    
    if (g <= length(v_name)/2 + 1){
      form_s2[g] = ""
    }
    if (g > length(v_name)/2 +1){
      form_s2[g] = paste(time[g-length(v_name)/2],"*",v_name[g-1], "+",sep = "")
    }
    if (g == length(v_name)+1){
      form_s2[g] = paste(time[g-length(v_name)/2],"*",v_name[g-1],sep = "")
    }
  }
  
  form_i_c = paste(form_i,collapse=" ")
  form_s_c = paste(form_s,collapse=" ")
  form_i2_c = paste(form_i2,collapse=" ")
  form_s2_c = paste(form_s2,collapse=" ")
  form_cov ='\n i ~~ 0*i2
         i ~~ 0*s2
         s ~~ 0*i2
         s ~~ 0*s2'
  
  if (cov == FALSE){
    final_lgm = paste(c(form_i_c, form_s_c, form_i2_c, form_s2_c),collapse=" ")
  }
  
  else{
    final_lgm = paste(c(form_i_c, form_s_c, form_i2_c, form_s2_c, form_cov),collapse=" ")
  }
  
  
  # Piece wise laten growth model regression
  M3 <- growth(final_lgm, data=wide)
  X = parameterEstimates(GCM1_fit)
  
  # obtain coefficients
  
  beta_1 <- X %>% 
    filter(lhs == "s",
           op == "~1") %>% 
    select(est)
  
  beta_1 <- beta_1[1,1]
  
  beta_2 <- X %>% 
    filter(lhs == "i2",
           op == "~1") %>% 
    select(est)
  
  beta_2 <- beta_2[1,1]
  
  beta_3 <- X %>% 
    filter(lhs == "s2",
           op == "~1") %>% 
    select(est)
  
  beta_3 <- beta_3[1,1]
  
  # Obtain the t-values
  z_value_pre_slope <- X %>% 
    filter(lhs == "s",
           op == "~1") %>% 
    select(z)
  
  z_value_pre_slope <- z_value_pre_slope[1,1]
  
  z_value_step_M2 <- X %>% 
    filter(lhs == "s2",
           op == "~1") %>% 
    select(z)
  z_value_step_M2 <- z_value_step_M2[1,1]
  
  z_value_slope_M2 <- X %>% 
    filter(lhs == "i2",
           op == "~1") %>% 
    select(z)
  
  z_value_slope_M2 <- z_value_slope_M2[1,1]
  
  if (abs(z_value_pre_slope) > 2.00){
    significant_pre_slope <- 1
  }
  else{
    significant_pre_slope <- 0
  }
  
  if (abs(z_value_step_M2) > 2.00){
    significant_step <- 1
  }
  else{
    significant_step <- 0
  }
  
  if (abs(z_value_slope_M2) > 2.00){
    significant_slope <- 1
  }
  else{
    significant_slope <- 0
  }
  
  return(list(significant_pre_slope = significant_pre_slope, significant_step = significant_step,
              significant_slope = significant_slope, beta_1 = beta_1, beta_2 = beta_2,
              beta_3 = beta_3))
}