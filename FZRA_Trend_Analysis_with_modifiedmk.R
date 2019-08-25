# Trend Analysis for ISD Historical Freezing Rain Time Series
# Matt Irish
# May 2019

library(modifiedmk)
library(R.matlab)
library(ggplot2)

# Redefine mmkh with an output of confidence interval:
mmkh_ci = function (x, ci = 0.95) 
{
  x = x
  z = NULL
  z0 = NULL
  pval = NULL
  pval0 = NULL
  S = 0
  Tau = NULL
  essf = NULL
  ci = ci
  if (is.vector(x) == FALSE) {
    stop("Input data must be a vector")
  }
  if (any(is.finite(x) == FALSE)) {
    x <- x[-c(which(is.finite(x) == FALSE))]
    warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }
  n <- length(x)
  if (n < 3) {
    stop("Input vector must contain at least three values")
  }
  V <- rep(NA, n * (n - 1)/2)
  k = 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k = k + 1
      V[k] = (x[j] - x[i])/(j - i)
    }
  }
  slp <- median(V, na.rm = TRUE) #sen's slope
  confidence_interval <- quantile(V,c(.05,.95),na.rm=T) #calculate confidence intervals
  t = 1:length(x)
  xn <- (x[1:n]) - ((slp) * (t)) #detrend the data
  #Calculate S:
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      S = S + sign(x[j] - x[i])
    }
  }
  ro <- acf(rank(xn), lag.max = (n - 1), plot = FALSE)$acf[-1]
  sig <- qnorm((1 + ci)/2)/sqrt(n)
  rof <- rep(NA, length(ro))
  for (i in 1:(length(ro))) {
    if (ro[i] > sig || ro[i] < -sig) {
      rof[i] <- ro[i]
    }
    else {
      rof[i] = 0
    }
  }
  cte <- 2/(n * (n - 1) * (n - 2))
  ess = 0
  for (i in 1:(n - 1)) {
    ess = ess + (n - i) * (n - i - 1) * (n - i - 2) * rof[i]
  }
  essf = 1 + ess * cte
  var.S = n * (n - 1) * (2 * n + 5) * (1/18)
  if (length(unique(x)) < n) {
    aux <- unique(x)
    for (i in 1:length(aux)) {
      tie <- length(which(x == aux[i]))
      if (tie > 1) {
        var.S = var.S - tie * (tie - 1) * (2 * tie + 
                                             5) * (1/18)
      }
    }
  }
  VS = var.S * essf
  if (S == 0) {
    z = 0
    z0 = 0
  }
  if (S > 0) {
    z = (S - 1)/sqrt(VS)
    z0 = (S - 1)/sqrt(var.S)
  }
  else {
    z = (S + 1)/sqrt(VS)
    z0 = (S + 1)/sqrt(var.S)
  }
  pval = 2 * pnorm(-abs(z))
  pval0 = 2 * pnorm(-abs(z0))
  Tau = S/(0.5 * n * (n - 1))
  return(c(`Corrected Zc` = z, `new P-value` = pval, `N/N*` = essf, 
           `Original Z` = z0, `old P.value` = pval0, Tau = Tau, 
           `Sen's slope` = slp, old.variance = var.S, new.variance = VS, confidence_interval= confidence_interval))
}



# Analyze trends
########################################################################################

setwd('/Users/mattirish/Documents/MATLAB/freezing_rain_trends_project')
a_all <- readMat("a_all.mat")
YearFreq_rel <- a_all[["a"]][[11]]
YearFreq_rel_domain_avg <- apply(YearFreq_rel,2, median, na.rm = T)
MonthFreq_yearly_rel_series <-  a_all[["a"]][[12]]
MonthFreq_yearly_rel_series_domain_avg <- apply(MonthFreq_yearly_rel_series,2, median, na.rm = T) #Domain average freezing rain for every month in series
MonthFreq_yearly_rel_series_domain_sum <- apply(MonthFreq_yearly_rel_series,2, sum, na.rm = T) #Domain average freezing rain for every month in series


# Estimate the annual trend for each station:
yearly_trends <-apply(YearFreq_rel, 1, function(x) mmkh_ci(x,c=0.95))
yearly_trends_3lag <-apply(YearFreq_rel, 1, function(x) mmkh3lag(x,c=0.95))
# Estimate the annual trend for the entire region:
yearly_trends_domain_avg <- mmkh(YearFreq_rel_domain_avg,c=0.95)

writeMat("yearly_trends.mat",yearly_trends_data=yearly_trends, yearly_trends_3lag_data=yearly_trends_3lag,fixNames=T)


# Estimate the monthly trend for each station:
pval_monthly_domain_avg <- vector("numeric", 12)
sen_monthly_domain_avg <- vector("numeric", 12)
for (month in 1:5){
  monthly_trends_domain_sum <- mmkh(MonthFreq_yearly_rel_series_domain_avg[seq(month,ncol(MonthFreq_yearly_rel_series), by=12)], c=0.95)
  
  pval_monthly_domain_avg[month] <- monthly_trends_domain_sum[2]
  sen_monthly_domain_avg[month] <- monthly_trends_domain_sum[7]
}


ggplot(yearly_trends[2,], aes(1:length(yearly_trends[2,]), yearly_trends[2,] )) +
  ggline()


## DIY-MODIFIED SEASONAL (HAMED et. al):
# Run on each month separately
# Add Z stats together and var(Z)s
# Compute new overall Z stat for each station
# Compare with standard normal distribution for significance.

Zc <- vector("numeric", 12)
sen <- vector("numeric", 12)
varZc <- 
yearly_trends_seasonal <- vector("numeric", 12)
yearly_trends_seasonal <- vector("numeric", 12)
for (m in 1:97){
for (month in 1:12){
  yearly_trends_seasonal[month] <- mmkh(MonthFreq_yearly_rel_series[m,seq(month,ncol(MonthFreq_yearly_rel_series), by=12)], c=0.95)

}
}