forms_for_sim <- function(N_time, var, Intercept, Slope,step_size, slope_size, additional_step){
  # this is a function to create the formula for generating data
  
  #Initialize variables
  v_name_final <- c()
  form <- c()
  
  # Create variable names
  numbs <- seq(1,N_time,1)
  v_name <- rep("Score", N_time)
  for (k in 1:length(v_name)){
    v_name_final[k] <- paste(v_name[k],numbs[k], sep = "")
  }
  
  # Calculate adjustment to the slope size
  step_adjust <- slope_size*0.5
  
  # Create the formulas
  for (h in 1:(length(numbs)-1)){
    if (h < (N_time/2)){
      form[h] <- paste(v_name_final[h],"+", Slope)
    }
    else if (h == (N_time/2)){
      form[h] <- paste(v_name_final[h],"+", Slope, "+", step_size, "+", step_adjust)
    }
    else{
      form[h] <- paste(v_name_final[h],"+", Slope, "+", slope_size)
    }
  }
  
  # Add an additional step if specified
  if (additional_step != 0){
    time_add_param = ceiling((3/4)*(length(numbs)-1))
    form[time_add_param] = paste(form[time_add_param],"+", additional_step)
  }
  
  var <- var
  
  # Return a list of the formula, variable names, and intercept
  return(list(form = form, var = var, v_name = v_name_final, Intercept = Intercept))
}

gen_data <- function(form, v_name, var_init, N_person, Intercept, var_time, additional_step_corr){
  # This function generates data and changes the format from wide to long.
  
  # Define the data structure
  def <- defData(varname = "Score1", dist = "normal", formula = Intercept,
                 variance = var_init)
  
  # Create variables for each time point using the formula
  for (q in 1:(length(v_name)-1)){
    def <- defData(def, varname = v_name[q+1], dist = "normal", formula = form[q],
                   variance = 0)
  }
  # Generate a data set with N_people
  dd <- genData(N_person, def)
  
  # generate noise 
  M_noise<-matrix(rnorm(N_person*length(v_name), 0, var_time),nrow=N_person)
  
  dd[,2:(length(v_name))] = dd[,2:(length(v_name))]  + M_noise
    
  # change df to long format and add a variable for time and a dummy for treatment
  dd_long <- pivot_longer(data = dd, 
                          cols = c(2:(length(form)+2)))
  
  # Set up the treatment variable and time variable based on whether an additional step was added
  if (additional_step_corr == 0){
    treatment = c(rep(0,(length(form)+1)/2), rep(1,(length(form)+1)/2))
    time_add_param = round((3/4)*(length(form)+1))
  }
  else{
    N <- length(form)+1
    treat_0 <- N/2
    treat_2 <-ceiling(N-treat_0*1.5)
    treat_1 <- N - treat_0 - treat_2
    treatment <- c(rep(0, treat_0), rep(1, treat_1), rep(2, treat_2))
  }
  
  # Add the dummy variables, time and score
  dd_long <- dd_long %>%  
    mutate(treatment = as.factor(rep(treatment, N_person)),
           treatment2 = as.factor(rep(c(rep(0,(length(form)+1)/2), rep(1,(length(form)+1)/2)), N_person)),
           time = rep(seq(1,length(form)+1,1), N_person)-(length(form)+1)/2-0.5,
           score = value)
  
  # return the long and wide dataset
  return(list(long = dd_long, wide = dd))
  
}

All_conditions <- function(effect_sizes, slope){
  
  #Initialize variables
  step_size_all <- c()
  slope_size_all <- c()
  
  # create a lsit for all effect sizes for the step and slope
    for (es in effect_sizes){
      step_size_all = append(step_size_all,es*slope)
      slope_size_all = append(slope_size_all,es*slope)
    }

  return(list(step_size_all = step_size_all, slope_size_all = slope_size_all))
}
