# ==========================================================================================
# bayes_helpers.R
#
# A "WriteOut() for Bayesian runs": calc_dq() (re-running TMBobj$report() once per
# posterior draw) is by far the slowest step in the Bayesian workflow. WriteBayesOut()
# runs it exactly once per model, keeps only the small subset of report() fields the
# dashboard plots actually consume (dropping large duplicated echoes of input data,
# per-herd intermediate arrays, etc.), grabs a compact set of raw posterior draws needed
# for cross-model comparison plots, and writes the result to a
# "<Code><SensCase>_BayesOut.rds" cache file next to the original _Bayes.rds/_TMB.rds/.rds
# triplet. Subsequent dashboard renders just read that cache (fast) instead of
# recomputing calc_dq() (slow).
#
# Depends on calc_dq()/pull_dq() from plotting_functions_Bayes.R, so that file is
# sourced here directly.
# ==========================================================================================

library(tidyverse)
library(here)
library(rstan)
library(tmbstan)
library(ggalluvial)
library(viridis)

# source(here("code", "final", "plotting_functions_Bayes.R"))

# ------------------------------------------------------------------------------------------
# Project color palette, keyed by model variant so a given model type always gets the same
# color across every comparison plot, regardless of the order runs were discovered in.
# ------------------------------------------------------------------------------------------
MODEL_PALETTE <- c(
  "rS"           = "#E63946",
  "ddOnly"       = "#756BB1",
  "env-survival" = "#238B8B",
  "env-K"        = "#E6781E"
)

# ==========================================================================================
# ---- Discovery, variant detection, and caching -------------------------------------------
# ==========================================================================================

# Same field-presence logic used for the ML dashboard: which report() outputs exist
# tells us which of the four model variants (rS / ddOnly / env-survival / env-K)
# produced a given run.
detect_mtype <- function(TMBout) {
  rn <- names(TMBout$report)
  has_KYr    <- "KYr" %in% rn
  has_envidx <- "env_index" %in% rn
  if (has_KYr && has_envidx)  return("env-K")
  if (has_KYr && !has_envidx) return("ddOnly")
  if (!has_KYr && has_envidx) return("env-survival")
  "rS"
}

# ------------------------------------------------------------------------------------------
# Find every *_Bayes.rds / *_TMB.rds / .rds triplet under base_dir (default
# Diags/final/), skipping any path that matches exclude_pattern (default "tests", so
# a sibling or nested Diags/.../tests/ folder is never picked up).
# ------------------------------------------------------------------------------------------
discover_bayes_runs <- function(base_dir = here("Diags", "final"), exclude_pattern = "tests") {

  bayes_files <- list.files(base_dir, pattern = "_Bayes\\.rds$", recursive = TRUE, full.names = TRUE)
  if (nzchar(exclude_pattern)) {
    bayes_files <- bayes_files[!grepl(exclude_pattern, bayes_files, ignore.case = TRUE)]
  }

  rows <- map(bayes_files, function(bf) {
    d    <- dirname(bf)
    base <- sub("_Bayes\\.rds$", "", basename(bf))
    tmb_path <- file.path(d, paste0(base, "_TMB.rds"))
    out_path <- file.path(d, paste0(base, ".rds"))
    if (!file.exists(tmb_path) || !file.exists(out_path)) {
      message("Skipping ", bf, " -- missing matching _TMB.rds or .rds file")
      return(NULL)
    }
    tibble(
      subdir     = sub(paste0("^", normalizePath(base_dir, mustWork = FALSE), "/?"), "",
                        normalizePath(d, mustWork = FALSE)),
      base_name  = base,
      path_bayes = bf,
      path_tmb   = tmb_path,
      path_out   = out_path
    )
  })

  bind_rows(rows)
}

# Which report() fields to keep, per model variant. FeedK/NbS/NfS/NNS/SurvOutF/MortDiff/
# NfitBreed/NfitFeed/PredMix are common to all four model files; KYr and env_index only
# exist for the variants that use them (see detect_mtype()).
bayes_dq_fields <- function(mtype) {
  flds <- c("NfitFeed", "NfitBreed", "PredMix", "NNS", "NbS", "NfS", "SurvOutF", "MortDiff", "FeedK")
  if (mtype %in% c("ddOnly", "env-K"))       flds <- c(flds, "KYr")
  if (mtype %in% c("env-K", "env-survival")) flds <- c(flds, "env_index")
  unique(flds)
}

# Small scalar posterior parameters worth keeping raw draws of (for cross-model
# comparison plots), on top of "logK"/"logBK" which get checked unconditionally below.
bayes_posterior_pars <- function(mtype) {
  pars <- c("rval")
  if (mtype %in% c("ddOnly", "env-K")) pars <- c(pars, "log_alphaK", "log_betaK")
  if (mtype == "env-K")                pars <- c(pars, "log_Ksigma")
  if (mtype == "env-survival")         pars <- c(pars, "log_sigmaEnv")
  if (mtype == "rS")                   pars <- c(pars, "log_SFsigma")
  pars
}

# ------------------------------------------------------------------------------------------
# WriteBayesOut(): run (or load the cached result of) the expensive derived-quantities
# step for one Bayesian run.
#   - path_bayes/path_tmb/path_out: the three files written by DoRun(..., DoBayes = TRUE)
#   - nsamp: number of posterior draws to use for calc_dq() (matches calc_dq()'s own
#     nsamp argument -- more draws = smoother CIs but slower and a bigger cache file)
#   - force: recompute even if a cache file already exists (e.g. after refitting, or to
#     change nsamp)
# ------------------------------------------------------------------------------------------
WriteBayesOut <- function(path_bayes, path_tmb, path_out, nsamp = 2000, force = FALSE) {

  cache_path <- sub("_Bayes\\.rds$", "_BayesOut.rds", path_bayes)

  if (file.exists(cache_path) && !force) {
    message("Using cached derived quantities: ", cache_path)
    return(read_rds(cache_path))
  }

  message("Computing derived quantities for ", basename(path_bayes),
          " (nsamp = ", nsamp, ") -- this is the slow step, please wait...")

  bayes  <- read_rds(path_bayes)
  TMBobj <- read_rds(path_tmb)
  TMBout <- read_rds(path_out)

  mtype <- detect_mtype(TMBout)

  # --- the expensive step: re-run $report() once per posterior draw -----------------------
  dq <- calc_dq(bayes, TMBobj, nsamp = nsamp)

  # --- keep only the report() fields the dashboard's plotting functions need --------------
  dq_fields <- bayes_dq_fields(mtype)
  dq_arrays <- set_names(dq_fields) |>
    map(\(fld) tryCatch(pull_dq(dq, fld) |> simplify2array(),
                         error = function(e) NULL)) |>
    compact()
  rm(dq); gc()  # the raw per-draw list is what makes this step memory-heavy; drop it once extracted

  # --- small parameter summary table (mean/sd/quantiles/Rhat/n_eff per parameter) ---------
  bayes_param_summary <- summary(bayes)$summary |> as_tibble(rownames = "parameter")

  # --- a handful of raw posterior draws needed for cross-model comparison plots -----------
  wanted_pars <- intersect(bayes_posterior_pars(mtype), bayes_param_summary$parameter)
  vec_pars <- c("logK", "logBK")
  vec_pars <- vec_pars[map_lgl(vec_pars, \(nm) any(grepl(paste0("^", nm, "(\\[|$)"), bayes_param_summary$parameter)))]
  extract_pars <- unique(c(wanted_pars, vec_pars))
  posterior_draws <- if (length(extract_pars) > 0) {
    tryCatch(rstan::extract(bayes, pars = extract_pars), error = function(e) list())
  } else list()

  BayesOut <- list(
    Code       = TMBout$Code,
    SensCase   = TMBout$SensCase,
    mtype      = mtype,
    nsamp      = nsamp,
    input      = TMBout$input,    # small: names, years, catch series, env data, etc.
    report_ML  = TMBout$report,   # ML-fit report -- needed for Qest, ObsMixProp(I)
    sdfixed_ML = TMBout$sdfixed,
    bayes_param_summary = bayes_param_summary,
    posterior_draws      = posterior_draws,
    dq_arrays  = dq_arrays
  )

  write_rds(BayesOut, cache_path)
  message("Cached to ", cache_path)
  BayesOut
}

# Convenience wrapper: discover every run under base_dir and load/compute its
# BayesOut cache in one call. Returns a tibble with a `data` list-column.
load_all_bayes_runs <- function(base_dir = here("Diags", "final"), nsamp = 2000, force = FALSE) {
  idx <- discover_bayes_runs(base_dir)
  if (nrow(idx) == 0) {
    warning("No *_Bayes.rds runs found under ", base_dir)
    return(idx)
  }
  idx$data <- pmap(idx, function(path_bayes, path_tmb, path_out, ...) {
    WriteBayesOut(path_bayes, path_tmb, path_out, nsamp = nsamp, force = force)
  })
  idx |>
    mutate(
      Code     = map_chr(data, "Code"),
      SensCase = map_chr(data, "SensCase"),
      mtype    = map_chr(data, "mtype"),
      label    = paste0(Code, SensCase, " [", mtype, "]")
    )
}

# ------------------------------------------------------------------------------------------
# Evaluate a plotting call defensively so one inapplicable/broken plot doesn't take down
# the whole rendered page.
# ------------------------------------------------------------------------------------------
safe_plot <- function(expr, label = "this plot") {
  tryCatch(
    expr,
    error = function(e) {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0,
                           label = paste0(label, " not available:\n", conditionMessage(e)),
                           size = 3.5) +
        ggplot2::theme_void()
    }
  )
}

# ------------------------------------------------------------------------------------------
# Compact summary row for one cached BayesOut object (Overview table).
# ------------------------------------------------------------------------------------------
summarize_bayes_run <- function(obj) {

  bps <- obj$bayes_param_summary |> filter(parameter != "lp__")

  final_mature <- tryCatch({
    NmatB   <- obj$dq_arrays$NfitBreed          # dims: component(3) x breed x year x sample
    last_yr <- dim(NmatB)[3]
    mature_slice <- NmatB[3, , last_yr, ]       # breed x sample
    totals  <- colSums(mature_slice)
    tibble(median = median(totals),
           low    = unname(quantile(totals, 0.025)),
           upper  = unname(quantile(totals, 0.975)))
  }, error = function(e) tibble(median = NA_real_, low = NA_real_, upper = NA_real_))

  rdraws <- obj$posterior_draws$rval

  tibble(
    Code      = obj$Code,
    SensCase  = obj$SensCase,
    mtype     = obj$mtype,
    nsamp     = obj$nsamp,
    max_Rhat  = suppressWarnings(max(bps$Rhat, na.rm = TRUE)),
    min_n_eff = suppressWarnings(min(bps$n_eff, na.rm = TRUE)),
    r_median  = if (!is.null(rdraws)) median(rdraws) else NA_real_,
    r_low     = if (!is.null(rdraws)) unname(quantile(rdraws, 0.025)) else NA_real_,
    r_upper   = if (!is.null(rdraws)) unname(quantile(rdraws, 0.975)) else NA_real_,
    final_mature_median = final_mature$median,
    final_mature_low    = final_mature$low,
    final_mature_upper  = final_mature$upper
  )
}

# ==========================================================================================
# ---- Adapted single-model plots (read from a cached BayesOut object) ---------------------
# These mirror the logic in plotting_functions_Bayes.R exactly, just swapping
# pull_dq(dqlist, "X") |> simplify2array() for the pre-computed obj$dq_arrays$X, and
# TMBout/bayesobj references for the equivalent pieces cached on the BayesOut object.
# ==========================================================================================

bplot_abundance <- function(obj) {

  yrs <- pluck(obj, "input", "Years")
  numyr <- length(yrs)

  breednames <- obj$input$BreedNames
  feednames  <- obj$input$FeedNames
  fopt <- ifelse("EAL+BER+WGOA" %in% feednames, "F2", "F1")
  bopt <- ifelse("MX_ML" %in% breednames, "B2", "B1")
  hyp <- paste0(bopt, fopt)
  if (fopt == "F1") hyp <- c(hyp, "F1 only")
  if (fopt == "F2") hyp <- c(hyp, "F2 only")
  if (bopt == "B1") hyp <- c(hyp, "B1 only")
  if (bopt == "B2") hyp <- c(hyp, "B2 only")
  hyp <- c(hyp, "All")

  Surveys <- read_csv(here("Diags", "SurveyUse", paste0(obj$Code, obj$SensCase, ".csv")), show_col_types = FALSE)
  Qvec <- obj$report_ML$Qest
  Qvec[1] <- 1
  survd <- Surveys |>
    filter(Hypothesis %in% hyp, Use == "Yes") |>
    rename(year = Year2) |>
    dplyr::select(year, Estimate, CV, Area, Class, Rel, Hypothesis, Class, component = Component) |>
    mutate(sd.estimate = Estimate * CV,
           q = 1 / Qvec[Class],
           rescaled.est = q * Estimate,
           rescaled.upper = q * (Estimate + sd.estimate),
           rescaled.lower = q * (Estimate - sd.estimate)) |>
    mutate(component = ifelse(component == 0, "0+", "1+"))

  NmatB <- obj$dq_arrays$NfitBreed
  NmatTot <- apply(NmatB, c(1, 3, 4), sum)
  quants <- apply(NmatTot, c(1, 2), quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))

  d <- tibble(year = rep(as.integer(yrs), each = 3),
              component = rep(c("0+", "1+", "Mature"), numyr)) |>
    mutate(low      = as.numeric(quants[1, , ]),
           lowmid   = as.numeric(quants[2, , ]),
           median   = as.numeric(quants[3, , ]),
           uppermid = as.numeric(quants[4, , ]),
           upper    = as.numeric(quants[5, , ]))
  surv <- survd |> filter(Area == "Total")
  d <- d |> left_join(surv, by = join_by(year, component))

  ggplot() +
    geom_ribbon(data = d, aes(x = year, ymin = low, ymax = upper, fill = component), alpha = 0.5) +
    geom_ribbon(data = d, aes(x = year, ymin = lowmid, ymax = uppermid, fill = component), alpha = 0.7) +
    geom_line(data = d, aes(year, median, color = component)) +
    geom_pointrange(data = surv, aes(year, rescaled.est, ymax = rescaled.upper, ymin = rescaled.lower),
                     linetype = 2, size = 0.5) +
    scale_fill_manual(values = c("#756BB1", "#238B8B", "#E6781E")) +
    scale_color_manual(values = c("#756BB1", "#238B8B", "#E6781E")) +
    labs(x = "Year", y = "Total Abundance", fill = "Component", color = "Component")
}

bplot_survival <- function(obj) {

  sdat <- obj$dq_arrays$SurvOutF
  SF <- pluck(obj, "input", "SF")
  sdat <- sdat[which(SF == 1), , ]
  quants <- apply(sdat, c(1, 2), quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))

  yrs <- pluck(obj, "input", "Years")
  yrs <- yrs[-length(yrs)]
  numyr <- length(yrs)

  zn <- pluck(obj, "input", "FeedNames") |> as.character()
  zn <- zn[which(SF == 1)]
  numz <- length(zn)

  sdf <- tibble(year = rep(yrs, each = numz)) |>
    mutate(zone = rep(zn, numyr)) |>
    mutate(year = as.integer(year)) |>
    mutate(low      = as.numeric(quants[1, , ]),
           lowmid   = as.numeric(quants[2, , ]),
           median   = as.numeric(quants[3, , ]),
           uppermid = as.numeric(quants[4, , ]),
           upper    = as.numeric(quants[5, , ])) |>
    filter(year > 1994)

  sdf |>
    ggplot() +
    geom_ribbon(aes(year, median, ymin = low, ymax = upper), fill = "#DADAEB") +
    geom_ribbon(aes(year, median, ymin = lowmid, ymax = uppermid), fill = "#756BB1") +
    geom_line(aes(year, median)) +
    geom_hline(yintercept = 0.96, linetype = 2, color = "orange") +
    labs(x = "Year", y = "Survival", title = "Survival") +
    facet_wrap(~zone) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

bplot_mixing <- function(obj) {

  BreedNames <- pluck(obj, "input", "BreedNames") |> as.character()
  FeedNames  <- pluck(obj, "input", "FeedNames") |> as.character()

  d <- obj$dq_arrays$NNS
  d <- d[, , dim(d)[3], ] |> apply(c(1, 2), median)

  df <- crossing(Feed = factor(FeedNames, levels = FeedNames),
                  Breed = factor(BreedNames, levels = BreedNames)) |>
    mutate(abun = as.numeric(d)) |>
    mutate(id = row_number())

  ggplot(df, aes(y = abun, axis1 = Breed, axis2 = Feed)) +
    geom_alluvium(aes(fill = Breed)) +
    geom_stratum(width = 1 / 12, fill = "gray80", color = "black", alpha = 0.5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "black", size = 3) +
    theme_classic() +
    scale_fill_manual(values = viridis_pal(option = "G")(length(BreedNames)), guide = "none") +
    theme(panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_text(size = 16)) +
    labs(y = "Breeding Ground") +
    scale_y_continuous(expand = c(0, 0), sec.axis = sec_axis(~., name = "Feeding Ground")) +
    scale_x_continuous(expand = c(0, 0))
}

bplot_proportions <- function(obj, direction = "B-F") {

  if (!(direction %in% c("B-F", "F-B"))) stop("direction must be B-F or F-B")

  BreedNames <- pluck(obj, "input", "BreedNames") |> as.character()
  FeedNames  <- pluck(obj, "input", "FeedNames") |> as.character()
  d      <- pluck(obj, "report_ML", "ObsMixProp")
  preds  <- obj$dq_arrays$PredMix
  quants <- apply(preds, 1, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))

  dwhich <- pluck(obj, "report_ML", "ObsMixPropI") |>
    as_tibble(.name_repair = "minimal") |>
    set_names(c("direction", "dataset", "breed", "feed"))

  dp <- dwhich |>
    mutate(direction = ifelse(direction == 1, "B-F", "F-B"),
           dataset = ifelse(dataset == 1, "mark-recapture", "genetics"),
           breed = BreedNames[breed],
           feed  = FeedNames[feed]) |>
    unite(labBF, breed, feed, remove = FALSE) |>
    unite(labFB, feed, breed, remove = FALSE)

  dpred <- dp |>
    mutate(low      = as.numeric(quants[1, ]),
           lowmid   = as.numeric(quants[2, ]),
           median   = as.numeric(quants[3, ]),
           uppermid = as.numeric(quants[4, ]),
           upper    = as.numeric(quants[5, ])) |>
    mutate(dataset = "model predictions")

  dobs <- dp |>
    mutate(est = d[, 1], sd = d[, 2]) |>
    mutate(upper = est + 1.96 * sd, low = est - 1.96 * sd)

  wanted <- direction
  d1 <- dobs  |> filter(direction == wanted)
  d2 <- dpred |> filter(direction == wanted)
  xcol <- if (direction == "B-F") "labBF" else "labFB"
  xlab <- if (direction == "B-F") "Breeding to Feeding" else "Feeding to Breeding"

  ggplot() +
    geom_pointrange(data = d1, aes(.data[[xcol]], est, ymin = low, ymax = upper,
                                    shape = dataset, color = dataset, group = dataset),
                     position = position_dodge(width = 0.8)) +
    geom_linerange(data = d2, aes(.data[[xcol]], median, ymin = low, ymax = upper, group = dataset),
                    linewidth = 1, color = "lightblue", position = position_nudge(x = 0.1)) +
    geom_linerange(data = d2, aes(.data[[xcol]], median, ymin = lowmid, ymax = uppermid, group = dataset),
                    linewidth = 1.5, color = "gray50", position = position_nudge(x = 0.1)) +
    geom_point(data = d2, aes(.data[[xcol]], median, color = dataset, shape = dataset, group = dataset),
               position = position_nudge(x = 0.1), size = 3) +
    scale_y_continuous(labels = seq(0, 1, by = 0.2), breaks = seq(0, 1, by = 0.2)) +
    coord_cartesian(ylim = c(0, 1.1)) +
    scale_color_manual(values = c("#E6781E", "#756BB1", "#238B8B")) +
    scale_shape_manual(values = c(16, 1, 17)) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = xlab, y = "Proportion", shape = "", color = "") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9))
}

bplot_fixed_effects <- function(obj) {
  df2 <- obj$bayes_param_summary |>
    filter(parameter != "lp__",
           !grepl("epsEnv", parameter),
           !grepl("Kdev", parameter),
           !grepl("SFdev", parameter)) |>
    mutate(type = case_when(
      grepl("MixPars", parameter)   ~ "Mixing",
      grepl("logBK", parameter)     ~ "Initial Depletion",
      grepl("envParams", parameter) ~ "Environmental",
      grepl("logK", parameter)      ~ "Carrying Capacity",
      TRUE ~ "Other"
    ))
  df2 |>
    ggplot() +
    geom_linerange(aes(parameter, mean, ymin = `2.5%`, ymax = `97.5%`), linewidth = 1, color = "lightblue") +
    geom_pointrange(aes(parameter, mean, ymin = `25%`, ymax = `75%`), linewidth = 1.5, color = "gray50") +
    labs(x = "Parameter", y = "Estimate (95% CI)", title = "Parameter Estimates") +
    facet_wrap(~type, scales = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

# Total-population mortality (natural vs. breeding/feeding catch), ported from
# plot_compare_mort(..., opt = "total", ...) in plotting_functions_Bayes.R.
# type = "raw" | "cumulative" | "rate"
bplot_mortality <- function(obj, type = "raw") {

  yrs <- pluck(obj, "input", "Years")
  yrs <- yrs[-length(yrs)]
  bn <- pluck(obj, "input", "BreedNames") |> as.character()
  fn <- pluck(obj, "input", "FeedNames") |> as.character()

  mdat <- obj$dq_arrays$MortDiff * -1
  mdat[mdat < 0] <- 0

  catchb <- obj$input$CatchB |> t()
  catchbdf <- tibble(catch = as.numeric(obj$input$CatchB), zone = rep(bn, each = length(yrs)),
                      year = rep(yrs, length(bn))) |>
    group_by(zone) |> arrange(year) |> mutate(cume_catch = cumsum(catch)) |> ungroup()
  catchf <- obj$input$CatchF |> t()
  catchfdf <- tibble(catch = as.numeric(obj$input$CatchF), zone = rep(fn, each = length(yrs)),
                      year = rep(yrs, length(fn))) |>
    group_by(zone) |> arrange(year) |> mutate(cume_catch = cumsum(catch)) |> ungroup()
  totcatchbdf <- catchbdf |> group_by(year) |> summarise(totcatchb = sum(catch), totcumeb = sum(cume_catch), .groups = "drop")
  totcatchfdf <- catchfdf |> group_by(year) |> summarise(totcatchf = sum(catch), totcumef = sum(cume_catch), .groups = "drop")

  if (type == "raw") {
    tmdat <- apply(mdat, c(3, 4), sum)
    tquants <- apply(tmdat, 1, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
    tmdf <- tibble(year = yrs) |>
      mutate(low = as.numeric(tquants[1, ]), lowmid = as.numeric(tquants[2, ]),
             median = as.numeric(tquants[3, ]), uppermid = as.numeric(tquants[4, ]),
             upper = as.numeric(tquants[5, ])) |>
      left_join(totcatchbdf, by = "year") |> left_join(totcatchfdf, by = "year")

    p <- tmdf |> filter(year > 1994) |>
      ggplot(aes(x = year)) +
      geom_line(aes(y = totcatchb, color = "Catch: Breeding")) +
      geom_line(aes(y = totcatchf, color = "Catch: Feeding")) +
      geom_ribbon(aes(ymin = low, ymax = upper, fill = "Natural"), alpha = 0.5) +
      geom_ribbon(aes(ymin = lowmid, ymax = uppermid, fill = "Natural"), alpha = 0.7) +
      geom_line(aes(y = median, color = "Natural"), color = "#756BB1") +
      scale_fill_manual(values = c("Catch: Breeding" = "#238B8B", "Catch: Feeding" = "#E6781E", "Natural" = "#756BB1")) +
      scale_color_manual(values = c("Catch: Breeding" = "#238B8B", "Catch: Feeding" = "#E6781E", "Natural" = "#756BB1")) +
      labs(x = "Year", y = "Mortality (ind.)", color = "Source", fill = "")
  }

  if (type == "cumulative") {
    tmdat <- apply(mdat, c(3, 4), sum)
    tcmdat <- apply(tmdat, 2, cumsum)
    tcmquants <- apply(tcmdat, 1, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
    tcmdf <- tibble(year = yrs) |>
      mutate(low = as.numeric(tcmquants[1, ]), lowmid = as.numeric(tcmquants[2, ]),
             median = as.numeric(tcmquants[3, ]), uppermid = as.numeric(tcmquants[4, ]),
             upper = as.numeric(tcmquants[5, ])) |>
      left_join(totcatchbdf, by = "year") |> left_join(totcatchfdf, by = "year")

    p <- tcmdf |> filter(year > 1994) |>
      ggplot(aes(x = year)) +
      geom_line(aes(y = totcumeb, color = "Catch: Breeding")) +
      geom_line(aes(y = totcumef, color = "Catch: Feeding")) +
      geom_ribbon(aes(ymin = low, ymax = upper, fill = "Natural"), alpha = 0.5) +
      geom_ribbon(aes(ymin = lowmid, ymax = uppermid, fill = "Natural"), alpha = 0.7) +
      geom_line(aes(y = median, color = "Natural"), color = "#756BB1") +
      scale_fill_manual(values = c("Catch: Breeding" = "#238B8B", "Catch: Feeding" = "#E6781E", "Natural" = "#756BB1")) +
      scale_color_manual(values = c("Catch: Breeding" = "#238B8B", "Catch: Feeding" = "#E6781E", "Natural" = "#756BB1")) +
      labs(x = "Year", y = "Cumulative Mortality (ind.)", color = "Source", fill = "")
  }

  if (type == "rate") {
    NbS <- obj$dq_arrays$NbS
    NbS <- NbS[, -dim(NbS)[2], ]
    NfS <- obj$dq_arrays$NfS
    NfS <- NfS[, -dim(NfS)[2], ]
    NNS <- obj$dq_arrays$NNS
    NNS <- NNS[, , -dim(NNS)[3], ]
    TotN <- apply(NNS, c(3, 4), sum)

    catchBtot <- colSums(catchb)
    catchrateBtot <- map(1:dim(TotN)[2], \(x) catchBtot / TotN[, x]) |> simplify2array()
    totcatchratebquants <- apply(catchrateBtot, 1, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
    totcatchbrdf <- tibble(year = yrs) |>
      mutate(low = as.numeric(totcatchratebquants[1, ]), lowmid = as.numeric(totcatchratebquants[2, ]),
             median = as.numeric(totcatchratebquants[3, ]), uppermid = as.numeric(totcatchratebquants[4, ]),
             upper = as.numeric(totcatchratebquants[5, ])) |>
      pivot_longer(low:upper, names_to = "quant", values_to = "rate") |> mutate(source = "CatchB")

    catchFtot <- colSums(catchf)
    catchrateFtot <- map(1:dim(TotN)[2], \(x) catchFtot / TotN[, x]) |> simplify2array()
    totcatchratefquants <- apply(catchrateFtot, 1, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
    totcatchfrdf <- tibble(year = yrs) |>
      mutate(low = as.numeric(totcatchratefquants[1, ]), lowmid = as.numeric(totcatchratefquants[2, ]),
             median = as.numeric(totcatchratefquants[3, ]), uppermid = as.numeric(totcatchratefquants[4, ]),
             upper = as.numeric(totcatchratefquants[5, ])) |>
      pivot_longer(low:upper, names_to = "quant", values_to = "rate") |> mutate(source = "CatchF")

    tmdat <- apply(mdat, c(3, 4), sum)
    Totrate <- tmdat / TotN
    Totrquants <- apply(Totrate, 1, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
    trdf <- tibble(year = yrs) |>
      mutate(low = as.numeric(Totrquants[1, ]), lowmid = as.numeric(Totrquants[2, ]),
             median = as.numeric(Totrquants[3, ]), uppermid = as.numeric(Totrquants[4, ]),
             upper = as.numeric(Totrquants[5, ])) |>
      pivot_longer(low:upper, names_to = "quant", values_to = "rate") |> mutate(source = "Natural")

    trdf <- trdf |> bind_rows(totcatchbrdf) |> bind_rows(totcatchfrdf) |>
      pivot_wider(names_from = "quant", values_from = "rate")

    p <- trdf |>
      mutate(source = case_when(source == "Natural" ~ "Natural",
                                 source == "CatchB" ~ "Catch: Breeding",
                                 source == "CatchF" ~ "Catch: Feeding")) |>
      filter(year > 1994) |>
      ggplot(aes(x = year, fill = source)) +
      geom_ribbon(aes(ymin = low, ymax = upper), alpha = 0.5) +
      geom_ribbon(aes(ymin = lowmid, ymax = uppermid), alpha = 0.7) +
      geom_line(aes(y = median, color = source)) +
      scale_fill_manual(values = c("Catch: Breeding" = "#238B8B", "Catch: Feeding" = "#E6781E", "Natural" = "#756BB1")) +
      scale_color_manual(values = c("Catch: Breeding" = "#238B8B", "Catch: Feeding" = "#E6781E", "Natural" = "#756BB1")) +
      labs(x = "Year", y = "Mortality Rate", color = "Source", fill = "Source")
  }

  p
}

# ==========================================================================================
# ---- Cross-model comparison plots (all discovered runs together) -------------------------
# ==========================================================================================

bplot_compare_abundance <- function(objlist, scen.names) {
  yrs <- pluck(objlist[[1]], "input", "Years")
  numyr <- length(yrs)

  dall <- map2(objlist, scen.names, \(obj, y) {
    NmatB <- obj$dq_arrays$NfitBreed
    NmatTot <- apply(NmatB, c(1, 3, 4), sum)
    quants <- apply(NmatTot, c(1, 2), quantile, probs = c(0.025, 0.5, 0.975))
    tibble(year = rep(as.integer(yrs), each = 3),
           component = rep(c("0+", "1+", "Mature"), numyr),
           scenario = y) |>
      mutate(low    = as.numeric(quants[1, , ]),
             median = as.numeric(quants[2, , ]),
             upper  = as.numeric(quants[3, , ]))
  }) |> list_rbind()

  dall |>
    ggplot(aes(year, median, ymin = low, ymax = upper, fill = scenario, color = scenario)) +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~component, nrow = 3, scales = "free_y") +
    scale_fill_manual(values = MODEL_PALETTE) +
    scale_color_manual(values = MODEL_PALETTE) +
    labs(x = "Year", y = "Total Abundance", fill = "Model", color = "Model")
}

bplot_compare_r <- function(objlist, scen.names) {
  dall <- map2(objlist, scen.names, \(obj, y) {
    rdraws <- obj$posterior_draws$rval
    if (is.null(rdraws)) return(NULL)
    tibble(value = rdraws, scenario = y)
  }) |> list_rbind()

  if (nrow(dall) == 0) return(NULL)

  dall |>
    ggplot(aes(value, fill = scenario, color = scenario)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = MODEL_PALETTE) +
    scale_color_manual(values = MODEL_PALETTE) +
    labs(x = "r (rmax)", y = "Density", fill = "Model", color = "Model",
         title = "Posterior estimates of r")
}

bplot_compare_K <- function(objlist, scen.names) {
  dall <- map2(objlist, scen.names, \(obj, y) {
    logK <- obj$posterior_draws$logK
    if (is.null(logK)) return(NULL)
    zn <- pluck(obj, "input", "BreedNames") |> as.character()
    colnames(logK) <- zn
    as_tibble(logK) |>
      pivot_longer(everything(), names_to = "breed", values_to = "value") |>
      mutate(scenario = y)
  }) |> list_rbind()

  if (nrow(dall) == 0) return(NULL)

  dsum <- dall |>
    group_by(breed, scenario) |>
    summarise(median = median(value),
              low   = quantile(value, 0.025),
              upper = quantile(value, 0.975), .groups = "drop")

  dsum |>
    ggplot(aes(breed, median, ymin = low, ymax = upper, color = scenario)) +
    geom_pointrange(position = position_dodge(width = 0.4)) +
    scale_color_manual(values = MODEL_PALETTE) +
    labs(x = "Breeding Ground", y = "ln K", color = "Model", title = "Estimates of K")
}

# Small formatting helper for the Overview / Compare summary tables.
format_bayes_summary_tbl <- function(tbl) {
  tbl |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    rename(
      `Stock Structure` = Code,
      `Sensitivity` = SensCase,
      `Model Type` = mtype,
      `# Posterior Draws` = nsamp,
      `Max R-hat` = max_Rhat,
      `Min ESS` = min_n_eff,
      `r (median)` = r_median,
      `r (2.5%)` = r_low,
      `r (97.5%)` = r_upper,
      `Final Mature Abund. (median)` = final_mature_median,
      `Final Mature Abund. (2.5%)` = final_mature_low,
      `Final Mature Abund. (97.5%)` = final_mature_upper
    )
}
