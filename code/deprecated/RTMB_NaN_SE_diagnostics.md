# Diagnosing NaN Standard Errors in RTMB Models

## Quick Diagnostic Steps

### 1. **Check the Hessian Matrix**
```r
# After fitting
fit <- nlminb(start = parms, objective = obj$fn, gradient = obj$gr)

# Get the Hessian
hess <- obj$he(fit$par)

# Check for issues
eigen_vals <- eigen(hess)$values
print(paste("Min eigenvalue:", min(eigen_vals)))
print(paste("Max eigenvalue:", max(eigen_vals)))
print(paste("Condition number:", max(abs(eigen_vals))/min(abs(eigen_vals))))

# Find problematic parameters
if(any(eigen_vals <= 0)) {
  cat("Hessian is not positive definite!\n")
  cat("Number of non-positive eigenvalues:", sum(eigen_vals <= 0), "\n")
}

# Check for NaN/Inf in Hessian
if(any(is.na(hess)) || any(is.infinite(hess))) {
  cat("Hessian contains NaN or Inf values\n")
  na_locations <- which(is.na(hess), arr.ind = TRUE)
  print(na_locations)
}
```

### 2. **Identify Unidentifiable Parameters (Profile Likelihood)**
```r
# Profile each parameter one at a time
profile_check <- function(obj, fit, par_index, par_name) {
  base_val <- fit$par[par_index]
  test_vals <- seq(base_val - 2, base_val + 2, length = 20)
  
  nll_vals <- sapply(test_vals, function(x) {
    test_par <- fit$par
    test_par[par_index] <- x
    obj$fn(test_par)
  })
  
  plot(test_vals, nll_vals, type = "l", 
       main = paste("Profile:", par_name),
       xlab = par_name, ylab = "Negative Log-Likelihood")
  abline(v = base_val, col = "red", lty = 2)
  
  # Check if flat (unidentifiable)
  range_nll <- max(nll_vals) - min(nll_vals)
  if(range_nll < 0.1) {
    cat(par_name, "appears FLAT (unidentifiable)\n")
    return(FALSE)
  }
  return(TRUE)
}

# Run for all parameters
par_names <- names(fit$par)
for(i in seq_along(fit$par)) {
  profile_check(obj, fit, i, par_names[i])
}
```

### 3. **Check Parameter Correlations**
```r
# Get correlation matrix from Hessian
if(!any(is.na(hess)) && all(eigen(hess)$values > 0)) {
  cov_matrix <- solve(hess)
  cor_matrix <- cov2cor(cov_matrix)
  
  # Find highly correlated parameters (|r| > 0.95)
  high_cor <- which(abs(cor_matrix) > 0.95 & abs(cor_matrix) < 1, arr.ind = TRUE)
  
  if(nrow(high_cor) > 0) {
    cat("\nHighly correlated parameters (|r| > 0.95):\n")
    for(i in 1:nrow(high_cor)) {
      row_idx <- high_cor[i, 1]
      col_idx <- high_cor[i, 2]
      if(row_idx < col_idx) {  # Avoid duplicates
        cat(sprintf("%s <-> %s: r = %.3f\n", 
                    par_names[row_idx], 
                    par_names[col_idx], 
                    cor_matrix[row_idx, col_idx]))
      }
    }
  }
}
```

### 4. **Check for Near-Zero Gradients (Flat Regions)**
```r
# Evaluate gradient at solution
grad <- obj$gr(fit$par)

# Identify parameters with very small gradients
small_grad <- abs(grad) < 1e-6
if(any(small_grad)) {
  cat("\nParameters with near-zero gradients (possibly flat):\n")
  print(par_names[small_grad])
  print(grad[small_grad])
}
```

### 5. **Test with Fixed Parameters**
```r
# Fix suspected problematic parameters and refit
# Example: if logK[1] is problematic

# Create a version with logK[1] fixed
obj_fixed <- MakeADFun(
  data = c(dat, list(logK_1_fixed = fit$par["logK[1]"])),
  parameters = parms_without_logK1,
  map = list(logK = factor(c(NA, 2:length(logK)))),  # NA fixes first element
  random = random_effects,
  DLL = "your_model"
)

fit_fixed <- nlminb(start = parms_without_logK1, 
                    objective = obj_fixed$fn, 
                    gradient = obj_fixed$gr)

# Compare SEs - if they're now fine, logK[1] was the problem
```

## Common Culprits in Your Model

### **1. Random Effects Variance Parameters**
```r
# Check log_SFsigma
# If SFdev values are all very close to zero, sigma is unidentifiable

# Quick check:
cat("SFdev range:", range(fit$par[grep("SFdev", names(fit$par))]), "\n")
cat("log_SFsigma:", fit$par["log_SFsigma"], "\n")
cat("Implied sigma:", exp(fit$par["log_SFsigma"]), "\n")

# If SFdev range is tiny AND sigma is large, there's a problem
```

### **2. Carrying Capacity Parameters (logK)**
```r
# Check if any breeding stocks have no data
# These logK values will be unidentifiable

# Check which stocks have survey data
stocks_with_data <- unique(SurveyI[, 4])  # Assuming column 4 is stock ID
all_stocks <- 1:Nbreed

unobserved <- setdiff(all_stocks, stocks_with_data)
if(length(unobserved) > 0) {
  cat("Stocks with NO direct observations:\n")
  cat("logK for these stocks may be unidentifiable:\n")
  print(paste0("logK[", unobserved, "]"))
}
```

### **3. Mixing Parameters (MixPars)**
```r
# Check if any mixing routes have no data
# Check ObsMixBtoFP and ObsMixFtoBP matrices

# Count observations per mixing route
n_obs_per_route <- sum(ObsMixBtoFP > 0) + sum(ObsMixFtoBP > 0)
cat("Total mixing proportion observations:", n_obs_per_route, "\n")

# If very few, mixing parameters may trade off with K
```

### **4. Initial Depletion (logBK)**
```r
# These are often problematic if:
# - No data in early years
# - Confounded with K

# Check if you have early data
early_years <- min(SurveyI[, 1])  # First year with data
cat("First data year:", early_years, "\n")
cat("Model start year:", Yr1, "\n")
cat("Gap:", early_years - Yr1, "years\n")

# Large gap = logBK poorly identified
```

## Systematic Approach

**Step-by-step debugging:**

```r
# 1. Start simple - fix all random effects
map_no_RE <- list(
  SBdev = factor(rep(NA, length(SBdev))),
  SFdev = factor(rep(NA, length(SFdev))),
  FBdev = factor(rep(NA, length(FBdev)))
)

obj_no_RE <- MakeADFun(data = dat, parameters = parms, 
                        map = map_no_RE, DLL = "model")
fit_no_RE <- nlminb(...)

# Do SEs work now? If yes, problem is in RE structure

# 2. Add REs back one at a time
# First SBdev only
map_SB_only <- list(
  SFdev = factor(rep(NA, length(SFdev))),
  FBdev = factor(rep(NA, length(FBdev)))
)
# ... test

# Then SFdev
# ... test

# This isolates which RE is problematic

# 3. Within problematic RE, test subsets
# E.g., if SFdev is problem, test each feeding ground separately
```

## Practical Fixes

### **Fix 1: Increase Penalties**
```r
# If parameters are weakly identified, strengthen priors/penalties
# In your likelihood:

# Change from:
Penal <- Penal + 0.0001*logK[Ibreed]*logK[Ibreed]

# To stronger penalty:
Penal <- Penal + 0.01*logK[Ibreed]*logK[Ibreed]  # 100x stronger
```

### **Fix 2: Fix Unidentifiable Parameters**
```r
# Use map argument to fix parameters at reasonable values
# Example: fix logBK to start at 50% of K

map_list <- list(
  logBK = factor(rep(NA, Nbreed))  # Fix all logBK
)

# Set starting values
parms$logBK <- rep(log(2), Nbreed)  # Start at 50% K
```

### **Fix 3: Reduce Model Complexity**
```r
# Share parameters across stocks/grounds
# Example: common survival deviations

# Instead of separate SFdev for each feeding ground:
# Use groups
map_list <- list(
  SFdev = factor(c(1, 1, 2, 2, 3, 3))  # Group grounds 1-2, 3-4, 5-6
)
```

### **Fix 4: Add Informative Priors**
```r
# Replace weak penalties with proper priors
# Example for logK:

# Add to likelihood:
LogLikeK <- -sum(dnorm(logK, mean = log(5000), sd = 0.5, log = TRUE))
neglogL <- neglogL + LogLikeK

# This regularizes K around 5000 (adjust to your system)
```

## TMB-Specific Diagnostics

```r
# Check if RTMB is happy
obj$env$parList()  # View parameter structure
obj$report()       # Get reported values

# Check for numerics issues
attributes(obj$fn(fit$par))  # Should show no warnings

# Verify randomeffects structure
if(length(obj$env$random) > 0) {
  cat("Random effects:", names(obj$env$random), "\n")
  cat("Number of RE:", length(obj$env$random), "\n")
}
```

## Red Flags in Your Specific Model

Based on your code, watch for:

1. **log_SFsigma too large** → SFdevs hitting boundaries
2. **logK values > 20** → Unrealistically large populations
3. **MixPars extreme** → Mixing proportions near 0 or 1
4. **logBK all identical** → Model can't distinguish initial states
5. **Convergence at boundary** → Check `fit$message`

## Output to Run

```r
# Comprehensive diagnostic report
diagnostic_report <- function(obj, fit) {
  cat("=" %R% 50, "RTMB DIAGNOSTIC REPORT", "=" %R% 50, "\n\n")
  
  # 1. Convergence
  cat("CONVERGENCE:\n")
  cat("Message:", fit$message, "\n")
  cat("Iterations:", fit$iterations, "\n")
  cat("NLL:", fit$objective, "\n\n")
  
  # 2. Hessian
  cat("HESSIAN DIAGNOSTICS:\n")
  hess <- obj$he(fit$par)
  eigs <- eigen(hess)$values
  cat("Eigenvalue range:", range(eigs), "\n")
  cat("Condition number:", max(abs(eigs))/min(abs(eigs)), "\n")
  cat("Positive definite:", all(eigs > 0), "\n\n")
  
  # 3. Gradient
  cat("GRADIENT CHECK:\n")
  grad <- obj$gr(fit$par)
  cat("Max |gradient|:", max(abs(grad)), "\n")
  cat("Parameters with |grad| > 0.01:\n")
  print(names(fit$par)[abs(grad) > 0.01])
  cat("\n")
  
  # 4. Standard errors
  cat("STANDARD ERRORS:\n")
  sdrep <- sdreport(obj)
  ses <- summary(sdrep)[, "Std. Error"]
  cat("NaN SEs:", sum(is.na(ses)), "out of", length(ses), "\n")
  if(sum(is.na(ses)) > 0) {
    cat("Parameters with NaN SE:\n")
    print(names(ses)[is.na(ses)])
  }
  
  cat("\n", "=" %R% 50, "END REPORT", "=" %R% 50, "\n")
}

diagnostic_report(obj, fit)
```

Run this to get a complete picture of what's going wrong!
