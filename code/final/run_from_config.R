# ==========================================================================================
# run_from_config.R
#
# Execute model runs from a YAML configuration file.
#
# Usage:
#   source("run_from_config.R")
#   
#   # Run all configurations
#   run_config("runs.yaml")
#   
#   # Run specific runs by name
#   run_config("runs.yaml", runs = c("baseline_env_survival", "sensitivity_sst_only"))
#   
#   # Run by group (defined in YAML)
#   run_config("runs.yaml", group = "baselines")
#   
#   # Dry run (show what would run without executing)
#   run_config("runs.yaml", runs = "baseline_env_survival", dry_run = TRUE)
#   
#   # Continue on error (run all even if some fail)
#   run_config("runs.yaml", on_error = "continue")
#
# ==========================================================================================

library(yaml)
library(tidyverse)
library(here)

# Source the main model-fitting function
source(here("code", "final", "model_fitting_function.R"))

# ==========================================================================================
# Main entry point function
# ==========================================================================================

run_config <- function(
  config_file = here('code','final',"runs.yml"),
  runs = NULL,                              # Names of specific runs to execute
  group = NULL,                             # Group name from config$groups
  dry_run = FALSE,                          # Show what would run without executing
  on_error = "stop",                        # "stop" = fail on error; "continue" = skip & continue
  verbose = TRUE,                           # Print detailed progress
  parallel = FALSE                          # Run in parallel (future: not yet implemented)
) {
  
  # ---- Load and parse configuration ----
  if (verbose) cat("\n", strrep("=", 90), "\n")
  if (verbose) cat("Loading configuration from:", config_file, "\n")
  
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }
  
  config <- yaml::read_yaml(config_file)
  
  if (verbose) cat("Found", length(config$runs), "total configurations\n")
  
  # ---- Determine which runs to execute ----
  run_names <- NULL
  
  if (!is.null(group)) {
    if (!group %in% names(config$groups)) {
      stop("Group '", group, "' not found in config. Available groups: ",
           paste(names(config$groups), collapse = ", "))
    }
    run_names <- config$groups[[group]]
    if (verbose) cat("Running group '", group, "' with", length(run_names), "runs\n")
  } else if (!is.null(runs)) {
    run_names <- runs
    if (verbose) cat("Running", length(run_names), "specified runs\n")
  } else {
    # Run all
    run_names <- purrr::map_chr(config$runs, ~.x$name %||% NA_character_)
    run_names <- run_names[!is.na(run_names)]
    if (verbose) cat("Running all", length(run_names), "configurations\n")
  }
  
  # Validate that requested runs exist
  available_runs <- purrr::map_chr(config$runs, ~.x$name %||% NA_character_)
  available_runs <- available_runs[!is.na(available_runs)]
  
  missing <- setdiff(run_names, available_runs)
  if (length(missing) > 0) {
    stop("The following runs were not found in config: ", paste(missing, collapse = ", "))
  }
  
  if (verbose) cat(strrep("=", 90), "\n\n")
  
  # ---- Prepare run specifications ----
  run_specs <- list()
  for (rname in run_names) {
    run_config_item <- purrr::detect(config$runs, ~(.x$name %||% NA_character_) == rname)
    if (is.null(run_config_item)) {
      stop("Could not find run: ", rname)
    }
    
    # Merge defaults with this run's config
    spec <- config$defaults
    spec <- utils::modifyList(spec, run_config_item)
    spec$name <- rname
    
    run_specs[[rname]] <- spec
  }
  
  # ---- Show what will run ----
  if (verbose) {
    cat("Runs to execute:\n")
    for (i in seq_along(run_specs)) {
      spec <- run_specs[[i]]
      cat(sprintf("  %d. %s", i, spec$name))
      if (!is.null(spec$description)) {
        cat(sprintf(" -- %s", spec$description))
      }
      cat("\n")
    }
    cat("\n")
  }
  
  if (dry_run) {
    cat("\n[DRY RUN MODE - No models will be fitted]\n\n")
    if (verbose) {
      cat("First run configuration (as example):\n")
      first_spec <- run_specs[[1]]
      print(str(first_spec, max.level = 1))
    }
    return(invisible(run_specs))
  }
  
  # ---- Execute runs ----
  cat("Starting model runs...\n")
  if (verbose) cat(strrep("-", 90), "\n\n")
  
  results <- list()
  start_time <- Sys.time()
  
  for (i in seq_along(run_specs)) {
    spec <- run_specs[[i]]
    rname <- spec$name
    
    cat(sprintf("\n[%d/%d] Running: %s", i, length(run_specs), rname))
    if (!is.null(spec$description)) {
      cat(sprintf(" -- %s", spec$description))
    }
    cat("\n")
    
    run_start <- Sys.time()
    
    tryCatch({
      # Convert YAML spec to DoRun() arguments
      args <- config_to_dorun_args(spec)
      
      if (verbose) {
        cat("  Config: env_opt =", args$envOpt, 
            ", Yr1:Yr2 =", args$Yr1, ":", args$Yr2,
            ", Bayes =", args$DoBayes, "\n")
      }
      
      # Call DoRun
      res <- do.call(DoRun, args)
      
      run_time <- difftime(Sys.time(), run_start, units = "mins")
      cat(sprintf("  ✓ Completed in %.1f minutes\n", run_time))
      
      results[[rname]] <- list(
        status = "success",
        time = run_time,
        result = res
      )
      
    }, error = function(e) {
      run_time <- difftime(Sys.time(), run_start, units = "mins")
      cat(sprintf("  ✗ ERROR after %.1f minutes: %s\n", run_time, e$message))
      
      results[[rname]] <<- list(
        status = "error",
        time = run_time,
        error = e$message
      )
      
      if (on_error == "stop") {
        stop("Stopping due to error in run: ", rname)
      } else if (on_error == "continue") {
        warning("Continuing to next run after error in: ", rname)
      }
    })
  }
  
  # ---- Summary ----
  total_time <- difftime(Sys.time(), start_time, units = "mins")
  
  if (verbose) cat(strrep("-", 90), "\n\n")
  
  cat("\n")
  cat(strrep("=", 90), "\n")
  cat("RUN SUMMARY\n")
  cat(strrep("=", 90), "\n")
  
  successes <- sum(purrr::map_lgl(results, ~.x$status == "success"))
  failures <- sum(purrr::map_lgl(results, ~.x$status == "error"))
  
  cat(sprintf("Total: %d runs | Success: %d | Failures: %d | Time: %.1f min\n",
              length(results), successes, failures, total_time))
  
  if (failures > 0) {
    cat("\nFailed runs:\n")
    for (rname in names(results)) {
      if (results[[rname]]$status == "error") {
        cat(sprintf("  • %s: %s\n", rname, results[[rname]]$error))
      }
    }
  }
  
  cat("\nRun log written to: run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt\n", sep = "")
  cat(strrep("=", 90), "\n\n")
  
  # Write log file
  write_run_log(results, run_specs, total_time)
  
  invisible(results)
}

# ==========================================================================================
# Helper: Convert YAML config to DoRun() arguments
# ==========================================================================================
#
# Maps YAML keys (snake_case) to DoRun() parameters (camelCase)
# Handles type conversions (vectors, booleans, NULLs)

config_to_dorun_args <- function(spec) {
  
  # Mapping from YAML field names to DoRun parameter names
  param_map <- list(
    # Run metadata
    code = "Code",
    sens_case = "SensCase",
    
    # Data
    data_file_name = "DataFileName",
    bycatch_file = "ByCatchFile",
    catch_series = "CatchSer",
    
    # Years
    yr1 = "Yr1",
    yr2 = "Yr2",
    yr_s_devs = "YrSDevs",
    
    # Demographics
    nage = "Nage",
    iamat = "IAmat",
    stray_base = "StrayBase",
    sa = "SA",
    sc = "SC",
    time_lag = "TimeLag",
    with_mirror = "WithMirror",
    
    # Model selection
    env_opt = "envOpt",
    env_vars = "envVars",
    env_lag = "envlag",
    spline_k = "splineK",
    
    # Random effects
    rvars = "rvars",
    
    # Density dependence
    dens_dep_opt = "DensDepOpt",
    
    # Parameters
    use_k_prior = "UseKPrior",
    kmax = "Kmax",
    add_cv = "AddCV",
    mix_weights = "MixWeights",
    
    # Numerical
    wght_total = "WghtTotal",
    idirichlet = "Idirichlet",
    max_n = "MaxN",
    seed = "seed",
    sf = "SF",
    
    # Output
    all_plots = "AllPlots",
    full_diag = "FullDiag",
    do_boot = "DoBoot",
    do_bayes = "DoBayes",
    
    # Directory
    subdir = "subdir",
    
    # Advanced
    return_model_obj = "return_model_obj",
    init = "Init"
  )
  
  args <- list()
  
  for (yaml_key in names(param_map)) {
    r_param <- param_map[[yaml_key]]
    
    if (yaml_key %in% names(spec)) {
      val <- spec[[yaml_key]]
      
      # Handle NULL/None values
      if (is.null(val)) {
        args[[r_param]] <- NULL
      } else {
        args[[r_param]] <- val
      }
    }
  }
  
  # Ensure required parameters are present
  if (is.null(args$Code)) {
    stop("Missing required parameter: code")
  }
  if (is.null(args$SensCase)) {
    stop("Missing required parameter: sens_case")
  }
  
  args
}

# ==========================================================================================
# Helper: Write log file
# ==========================================================================================

write_run_log <- function(results, run_specs, total_time) {
  
  log_file <- here('code','run logs',paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  
  lines <- c(
    "Run Configuration Log",
    paste("Generated:", Sys.time()),
    "",
    "Execution Summary",
    strrep("-", 80),
    sprintf("Total runs: %d", length(results)),
    sprintf("Successful: %d", sum(purrr::map_lgl(results, ~.x$status == "success"))),
    sprintf("Failed: %d", sum(purrr::map_lgl(results, ~.x$status == "error"))),
    sprintf("Total time: %.1f minutes", total_time),
    "",
    "Run Details",
    strrep("-", 80)
  )
  
  for (rname in names(results)) {
    res <- results[[rname]]
    spec <- run_specs[[rname]]
    
    lines <- c(lines, "")
    lines <- c(lines, paste0("Run: ", rname))
    if (!is.null(spec$description)) {
      lines <- c(lines, paste0("Description: ", spec$description))
    }
    lines <- c(lines, paste0("Status: ", res$status))
    lines <- c(lines, paste0("Time: ", sprintf("%.1f min", res$time)))
    
    if (res$status == "error") {
      lines <- c(lines, paste0("Error: ", res$error))
    } else {
      lines <- c(lines, 
                 paste0("Config: env_opt=", spec$env_opt, 
                        ", years=", spec$yr1, "-", spec$yr2,
                        ", subdir=", spec$subdir))
    }
  }
  
  writeLines(lines, log_file)
  cat("Log written to:", log_file, "\n")
}

# ==========================================================================================
# Utility: List available runs in a config file
# ==========================================================================================

list_runs <- function(config_file = here('code','final',"runs.yml")){
  config <- yaml::read_yaml(config_file)
  
  cat("\nAvailable runs in", config_file, ":\n")
  cat(strrep("-", 80), "\n")
  
  for (run in config$runs) {
    name <- run$name %||% "(unnamed)"
    desc <- run$description %||% "(no description)"
    cat(sprintf("%-30s %s\n", name, desc))
  }
  
  cat("\n")
  cat("Available groups:\n")
  cat(strrep("-", 80), "\n")
  for (group in names(config$groups)) {
    n_runs <- length(config$groups[[group]])
    cat(sprintf("%-30s (%d runs)\n", group, n_runs))
  }
  cat("\n")
  
  invisible(config)
}

# ==========================================================================================
# Example usage (commented out so script can be sourced safely)
# ==========================================================================================

#' # Example 1: Run all baselines
#' run_config("runs.yaml", group = "baselines")
#' 
#' # Example 2: Run specific runs
#' run_config("runs.yaml", runs = c("baseline_env_survival", "sensitivity_sst_only"))
#' 
#' # Example 3: Dry run to see what would execute
#' run_config("runs.yaml", group = "sensitivity_short", dry_run = TRUE)
#' 
#' # Example 4: Continue on errors
#' run_config("runs.yaml", group = "all", on_error = "continue")
#' 
#' # Example 5: List all available runs
#' list_runs("runs.yaml")
