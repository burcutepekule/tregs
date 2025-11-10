
# Function to calculate Cohen's d
cohens_d = function(x, y) {
  nx = length(x)
  ny = length(y)
  mx = mean(x, na.rm = TRUE)
  my = mean(y, na.rm = TRUE)
  sx = sd(x, na.rm = TRUE)
  sy = sd(y, na.rm = TRUE)
  
  # Pooled standard deviation
  pooled_sd = sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2))
  
  # Cohen's d
  d = (my - mx) / pooled_sd
  
  if (is.na(pooled_sd) || pooled_sd == 0) {
    return(0)  # or NA, depending on how you want to interpret it
  }
  return(d)
}

steady_state_idx = function(x, k = 20, tail_frac = 0.25,
                            tol_abs = 0.05*(150-25),     # 2.5 units
                            tol_sd  = 0.02*(150-25),# 1.25 units
                            tol_slope = 0.005*(150-25))  # 0.125/step
{
  n      = length(x)
  tail_n = ceiling(n*tail_frac)
  x_asym = mean(tail(x, tail_n))
  
  m    = rollapply(x, k, mean, align = "right", fill = NA)
  sd   = rollapply(x, k, sd,   align = "right", fill = NA)
  sl   = rollapply(x, k, function(v) mean(diff(v)), align = "right", fill = NA)
  cand = which(abs(m-x_asym) <= tol_abs & sd <= tol_sd & abs(sl) <= tol_slope)
  
  if (length(cand) == 0) return(NA_integer_)
  cand[1]
}