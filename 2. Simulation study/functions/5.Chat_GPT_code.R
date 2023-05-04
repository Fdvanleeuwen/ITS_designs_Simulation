sim_all_GPT <- function(N_person, N_time, N_sim, var, Intercept, Slope, seed, var_time, effect_sizes, treatment_step, treatment_slope, additional_step, additional_step_corr, model = c("OLS, ML, SEM"), cov = FALSE){
  
  All_conditions <- All_conditions(effect_sizes, Slope)
  step_conitions <- All_conditions$step_size_all
  slope_conitions <- All_conditions$slope_size_all
  
  results <- lapply(1:length(step_conitions), function(EV) {
    step = step_conitions[EV] * treatment_step
    slope = slope_conitions[EV] * treatment_slope
    additional_step_use = step * 0.5 * additional_step
    
    lapply(1:length(N_time), function(l) {
      lapply(1:length(N_person), function(p) {
        sapply(1:N_sim, function(s) {
          forms <- forms_for_sim(N_time[l], var, Intercept, Slope, step, slope, additional_step_use)
          data <- gen_data(forms$form, forms$v_name, forms$var, N_person[p], forms$Intercept, var_time, additional_step_corr)
          
          if (model == "ML") {
            power <- run_ML(data$long)
          } else if (model == "SEM") {
            power <- run_SEM(data$wide, forms$v_name, cov = cov)
          } else {
            power <- run_OLS(data$long, additional_step_corr)
          }
          
          list(
            significant_pre_slope = power$significant_pre_slope,
            significant_step = power$significant_step,
            significant_slope = power$significant_slope,
            bias_precent_b1 = (power$beta_1 - Slope),
            bias_precent_b2 = ifelse(treatment_step == 0, power$beta_2, power$beta_2 - step),
            bias_precent_b3 = ifelse(treatment_slope == 0, power$beta_3, power$beta_3 - slope)
          )
        })
      })
    })
  })
  
  power <- data_frame(
    effect_sizes = rep(effect_sizes, each = length(N_time) * length(N_person) * N_sim),
    N_t = rep(N_time, each = length(N_person) * N_sim),
    N_pers = rep(N_person, times = length(N_time) * N_sim),
    pre_slope = unlist(lapply(results, "[[", "significant_pre_slope"), recursive = FALSE),
    step = unlist(lapply(results, "[[", "significant_step"), recursive = FALSE),
    slope = unlist(lapply(results, "[[", "significant_slope"), recursive = FALSE),
    bias_pre_slope = unlist(lapply(results, "[[", "bias_precent_b1"), recursive = FALSE),
    bias_step = unlist(lapply(results, "[[", "bias_precent_b2"), recursive = FALSE),
    bias_slope = unlist(lapply(results, "[[", "bias_precent_b3"), recursive = FALSE),
    Precision_pre_slope = unlist(lapply(results, function(x) sd(x$bias_precent_b1, na.rm = TRUE) / sqrt(length(N_person))), recursive = FALSE),
    Precision_step = unlist(lapply(results, function(x) sd(x$bias_precent_b2, na.rm = TRUE) / sqrt(length(N_person))), recursive = FALSE),
    Precision_slope = unlist(lapply(results, function(x) sd(x$bias_precent_b3, na.rm = TRUE) / sqrt(length(N_person))), recursive = FALSE)
  )
  return(power)
}
  