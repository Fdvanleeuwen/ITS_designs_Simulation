# Check https://www.zachburchill.ml/constant_scales/ for a walkthrough. 

# Given a plot and a scale, get the range of the values used in the scale for that plot
#   Not really intended for users
simple_range_extracter <- function(p, scale) {
  d <- ggplot_build(p)
  d$plot$scales$get_scales(scale$aesthetics)$range$range
}

# Given a bunch of plots and a scale you want to apply to them, this returns that scale,
#   but with limits that encompass the union of all the ranges of values in those plots
# Meant to be for users
get_shared_scale <- function(..., scale) {
  plots <- list(...)
  ranges <- purrr::map(plots, ~simple_range_extracter(., scale))
  single_range <- range(unlist(ranges))
  scale$limits <- single_range
  scale
}

# Given the unquoted variable names of a bunch of plots that have been 
#   saved to the current environment, and a scale that you want to apply
#   to them, this will call `get_shared_scale()` using those plots and scale
#   add that scale to the plots, and then assign these new plots to those 
#   variable names, essentially editing the plots 'in place'
set_scale_union <- function(..., scale) {
  exprs <- rlang::enexprs(...)
  scale <- get_shared_scale(..., scale = scale)
  var_nms <- purrr::map_chr(exprs, rlang::as_name)
  edit_plots_in_place(var_nms, env = parent.frame(),
                      scale = scale)
  # Invisibly return the scale, in case you need it later
  invisible(scale)
}

# Sub-function
edit_plots_in_place <- function(names, env, scale) {
  vars <- rlang::env_has(env = env, nms = names)
  if (!all(vars))
    stop("Environment does not have variables for ",
         paste(names(vars[!vars]), collapse=", "))
  
  purrr:::walk(names, function(nm) {
    og_plot <- rlang::env_get(env, nm = nm)
    message("Changing plot `", nm, "`")
    # Muffles messages about already having scales
    withCallingHandlers(
      assign(x = nm, envir = env,
             value = og_plot + scale),
      message = function(err) {
        if (grepl("already present", err$message))
          invokeRestart("muffleMessage")
      })
  })
}