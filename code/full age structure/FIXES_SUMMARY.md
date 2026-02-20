# Code Safety Fixes - Summary

## Issues Fixed to Prevent NA Function Evaluations

### 1. **Division by Zero Protection**

**Line 68 (Mixing matrix normalization):**
- **Problem:** `Total` could be zero when normalizing mixing proportions
- **Fix:** Added `Total <- max(Total, 1e-10)` before division

**Line 89-90 (Survival rate protection):**
- **Problem:** If `SA = 1`, the denominator `(1-SA)` = 0 in equilibrium calculations
- **Fix:** Capped survival: `SA <- min(SA, 0.9999)` and `SC <- min(SC, 0.9999)`

**Line 91 (Equilibrium abundance):**
- **Problem:** `NfecEqn` could be zero, causing division issues in line 194
- **Fix:** Added `NfecEqn <- max(NfecEqn, 1e-10)`

**Lines 465, 487, 520 (Likelihood calculations):**
- **Problem:** Division by zero when:
  - `ndat = 0` when calculating catchability `q`
  - `NbS[Ibreed,YearFeedBreed] = 0` for mixing proportions
  - `NfS[Ifeed,YearFeedBreed] = 0` for mixing proportions
- **Fix:** 
  - `ndat <- max(ndat, 1e-10)`
  - `NbS_safe <- max(NbS[Ibreed,YearFeedBreed], ThreshPop)`
  - `NfS_safe <- max(NfS[Ifeed,YearFeedBreed], ThreshPop)`

**Lines 236, 254 (Density dependence calculations):**
- **Problem:** Division by zero when `Nb[Ibreed,Year]` or `NowNf` are zero, or when K values are zero
- **Fix:** 
  ```r
  Nb[Ibreed,Year] <- max(Nb[Ibreed,Year], ThreshPop)
  BreedK[Ibreed] <- max(BreedK[Ibreed], ThreshPop)
  NowNf <- max(NowNf, ThreshPop)
  FeedK[Ifeed] <- max(FeedK[Ifeed], ThreshPop)
  ```

**Line 286 (Fecundity calculation):**
- **Problem:** Division by zero in `Frac <- NowNf/FeedK_safe`
- **Fix:** `FeedK_safe <- max(FeedK[Ifeed]*MultK[Ifeed,Year], ThreshPop)`

### 2. **Log of Zero/Negative Values**

**Lines 436, 475, 580-582 (Log-normal likelihoods):**
- **Problem:** Taking `log()` of predictions that could be ≤ 0
- **Fix:** Enforced minimum threshold for all abundance predictions:
  ```r
  # In projection loop (lines 221, 231)
  NfitBreed[Icomp,Ibreed,Year] <- max(NfitBreed[Icomp,Ibreed,Year], ThreshPop)
  NfitFeed[Icomp,Ifeed,Year] <- max(NfitFeed[Icomp,Ifeed,Year], ThreshPop)
  
  # In final year (lines 377, 388)
  NfitBreed[Icomp,Ibreed,Nyr+1] <- max(NfitBreed[Icomp,Ibreed,Nyr+1], ThreshPop)
  NfitFeed[Icomp,Ifeed,Nyr+1] <- max(NfitFeed[Icomp,Ifeed,Nyr+1], ThreshPop)
  
  # In likelihood (lines 425, 451)
  Pred <- max(Pred, ThreshPop)
  ```

### 3. **Undefined Variable (SD2)**

**Lines 503, 536:**
- **Problem:** Code uses `SD2` but only defines `SD` 
- **Fix:** Added explicit definition: `SD2 <- SD^2`

### 4. **Exponential Overflow Protection**

**Lines 96-99 (Fecundity calculations):**
- **Problem:** `exp(rval*(Amat+1.0))` could overflow for large `rval` values
- **Fix:** Capped exponential arguments:
  ```r
  exp_term1 <- min(rval*(Amat+1.0), 700)
  exp_term2 <- min(rval*Amat, 700)
  fmax <- 2*(exp(exp_term1)-SA*exp(exp_term2))/(SC*SA^(Amat))
  ```

## Key Safety Principles Applied

1. **Abundance Floor:** All abundance values protected with `ThreshPop = 1.0e-20`
2. **Survival Caps:** Survival rates capped at 0.9999 to prevent division by (1-S) = 0
3. **Division Safety:** All denominators protected with `max(value, threshold)` 
4. **Log Safety:** All logged values guaranteed positive via abundance floors
5. **Exponential Safety:** Large exponential arguments capped at 700 to prevent overflow

## Testing Recommendations

1. **Test with extreme parameters:**
   - Very high survival rates (SA, SC near 1)
   - Very low initial abundances
   - High growth rates (rval)
   - Zero catches/removals

2. **Monitor these outputs for reasonableness:**
   - `fmax` should be < 1
   - All `NfitBreed` and `NfitFeed` should be > 0
   - `Mix` matrix rows should sum to 1
   - `neglogL` should be finite

3. **Add bounds to optimizer:**
   - Consider bounding `logK` to reasonable values
   - Bound `rval` to prevent extreme fecundity
   - Bound survival devs to ±2 or ±3 SD

## Additional Recommendations

Consider adding at the start of the function:
```r
# Validate inputs
if (any(is.na(parms)) || any(is.infinite(parms))) {
  return(1e10)  # Large penalty for invalid parameters
}
```

And before returning neglogL:
```r
# Check for invalid result
if (is.na(neglogL) || is.infinite(neglogL)) {
  cat("WARNING: Invalid likelihood. Check parameter values.\n")
  return(1e10)
}
```
