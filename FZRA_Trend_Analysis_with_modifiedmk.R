# Trend Analysis for ISD Historical Freezing Rain Time Series
# Matt Irish
# May 2019

library(modifiedmk)
library(R.matlab)
library(ggplot2)
library(ggridges)

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
  if(all(xn == 1)) xn[1] = 1.001
  ro <- acf(rank(xn), lag.max = (n - 1), plot = FALSE)$acf[-1]
  sig <- qnorm((1 + ci)/2)/sqrt(n)
  rof <- rep(NA, length(ro))
  #warning(print(xn))
  #warning(sprintf("ro is \n %s \n and sig is \n %s \n", ro, sig))
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

# FREQUENCY ############################################################################
a_all <- readMat("a_all.mat")
YearFreq_rel <- a_all[["a"]][[11]]
states <- as.factor(unlist(a_all[["a"]][[15]]))
YearFreq_rel_domain_avg <- apply(YearFreq_rel,2, median, na.rm = T)
MonthFreq_yearly_rel_series <-  a_all[["a"]][[12]]
MonthFreq_yearly_rel_series_domain_avg <- apply(MonthFreq_yearly_rel_series,2, median, na.rm = T) #Domain average freezing rain for every month in series
MonthFreq_yearly_rel_series_domain_sum <- apply(MonthFreq_yearly_rel_series,2, sum, na.rm = T) #Domain average freezing rain for every month in series
MonthFreq_yearly_rel_series_domain_75 <- apply(MonthFreq_yearly_rel_series,2, quantile, na.rm = T) #Domain average freezing rain for every month in series
MonthFreq_yearly_rel_series_domain_75 <- MonthFreq_yearly_rel_series_domain_75[4,]

# Estimate the annual trend for each station:
yearly_trends <-apply(YearFreq_rel, 1, function(x) mmkh_ci(x,c=0.95))
yearly_trends_9614 <-apply(YearFreq_rel[,21:39], 1, function(x) mmkh_ci(x,c=0.95))
# Estimate the annual trend for each state/province:
yearly_trends_state_avg <- data.table()
yearly_trends_state_avg_9614 <- data.table()
for(state in unique(states)){
  YearFreq_rel_state_avg <- apply(YearFreq_rel[states == state,],2, median, na.rm = T)
  yearly_trends_state_avg <- rbind(yearly_trends_state_avg, as.list(c(c(state=state),mmkh_ci(YearFreq_rel_state_avg,c=0.95))))
  yearly_trends_state_avg_9614 <- rbind(yearly_trends_state_avg_9614, as.list(c(c(state=state),mmkh_ci(YearFreq_rel_state_avg[21:39],c=0.95))))
}
# Estimate the annual trend for the entire region:
yearly_trends_domain_avg <- mmkh_ci(YearFreq_rel_domain_avg,c=0.95)
yearly_trends_domain_avg_9614 <- mmkh_ci(YearFreq_rel_domain_avg[21:39],c=0.95)

writeMat("yearly_trends.mat",
         yearly_trends_data=yearly_trends, 
         yearly_trends_9614_data=yearly_trends_9614, 
         yearly_trends_state_avg=yearly_trends_state_avg, 
         yearly_trends_state_avg_9614=yearly_trends_state_avg_9614,         
         yearly_trends_domain_avg=yearly_trends_domain_avg, 
         yearly_trends_domain_avg_9614=yearly_trends_domain_avg_9614,
         fixNames=T)


# Estimate the monthly trend for each station:
months <- c(11,12,1,2,3,4)
pval_monthly_domain_avg <- vector("numeric", length(months))
sen_monthly_domain_avg <- vector("numeric", length(months))
monthly_trends <- data.table(month = months, 
                             sen_median = sen_monthly_domain_avg, p_median = pval_monthly_domain_avg, 
                             sen_quant = sen_monthly_domain_avg, p_quant = pval_monthly_domain_avg)
for (thismonth in months){ #would include October and May but all medians are zero throughout the time series and this throws an error
  monthly_trends_domain_median <- mmkh(MonthFreq_yearly_rel_series_domain_avg[seq(thismonth,ncol(MonthFreq_yearly_rel_series), by=12)], c=0.95)
  monthly_trends_domain_75 <- mmkh(MonthFreq_yearly_rel_series_domain_75[seq(thismonth,ncol(MonthFreq_yearly_rel_series), by=12)], c=0.95)
  
  monthly_trends <- monthly_trends[month == thismonth, sen_median:=monthly_trends_domain_median[7]]
  monthly_trends <- monthly_trends[month == thismonth, p_median:=monthly_trends_domain_median[2]]
  monthly_trends <- monthly_trends[month == thismonth, sen_quant:=monthly_trends_domain_75[7]]
  monthly_trends <- monthly_trends[month == thismonth, p_quant:=monthly_trends_domain_75[2]]
}



ggplot(yearly_trends[2,], aes(1:length(yearly_trends[2,]), yearly_trends[2,] )) +
  ggline()


# DURATION ############################################################################
fzra_durations <- readMat("FZRA_durations.mat")
meandurations <- fzra_durations[["meandurations"]]
mediandurations <- fzra_durations[["mediandurations"]]
mediandurations_domain_avg <- apply(mediandurations,2, median, na.rm = T)
# Estimate the annual DURATION trend for each station:
yearly_trends_duration <-apply(mediandurations, 1, function(x) mmkh_ci(x,c=0.95))
yearly_trends_duration_9614 <-apply(mediandurations[,21:39], 1, function(x) mmkh_ci(x,c=0.95))
# Estimate the annual trend for each state/province:
yearly_trends_duration_state_avg <- data.table()
yearly_trends_duration_state_avg_9614 <- data.table()
for(state in unique(states)){
  YearFreq_rel_duration_state_avg <- apply(mediandurations[states == state,],2, median, na.rm = T)
  yearly_trends_duration_state_avg <- rbind(yearly_trends_duration_state_avg, as.list(c(c(state=state),mmkh_ci(YearFreq_rel_duration_state_avg,c=0.95))))
  yearly_trends_duration_state_avg_9614 <- rbind(yearly_trends_duration_state_avg_9614, as.list(c(c(state=state),mmkh_ci(YearFreq_rel_duration_state_avg[21:39],c=0.95))))
}
# Estimate the annual trend for the entire region:
yearly_trends_duration_domain_avg <- mmkh_ci(mediandurations_domain_avg,c=0.95)
yearly_trends_duration_domain_avg_9614 <- mmkh_ci(mediandurations_domain_avg[21:39],c=0.95)

writeMat("yearly_trends_duration.mat",
         yearly_trends_duration_data=yearly_trends_duration,
         yearly_trends_duration_9614_data=yearly_trends_duration_9614,
         yearly_trends_duration_state_avg=yearly_trends_duration_state_avg, 
         yearly_trends_duration_state_avg_9614=yearly_trends_duration_state_avg_9614,
         yearly_trends_duration_domain_avg=yearly_trends_duration_domain_avg, 
         yearly_trends_duration_domain_avg_9614=yearly_trends_duration_domain_avg_9614,
         fixNames=T)

# INTENSITY ############################################################################
fzra_intensities <- readMat("FZRA_intensities.mat")
pct_islightFZRA <- fzra_intensities[["pct.islightFZRA"]]
pct_islightFZRA_domain_avg <- apply(pct_islightFZRA,2, median, na.rm = T)

# Estimate the annual trend for each station:
yearly_trends_intensity <- apply(pct_islightFZRA, 1, function(x) mmkh_ci(x,c=0.95)) #[,21:39]
yearly_trends_intensity_9614 <- apply(pct_islightFZRA[,21:39], 1, function(x) mmkh_ci(x,c=0.95)) #[,21:39]
#yearly_trends_intensity <- apply(pct_islightFZRA[,21:39], 1, function(x) mmkh(x))
# Estimate the annual trend for each state/province:
yearly_trends_intensity_state_avg <- data.table()
yearly_trends_intensity_state_avg_9614 <- data.table()
for(state in unique(states)){
  YearFreq_rel_intensity_state_avg <- apply(pct_islightFZRA[states == state,],2, median, na.rm = T)
  yearly_trends_intensity_state_avg <- rbind(yearly_trends_intensity_state_avg, as.list(c(c(state=state),mmkh_ci(YearFreq_rel_intensity_state_avg,c=0.95))))
  yearly_trends_intensity_state_avg_9614 <- rbind(yearly_trends_intensity_state_avg_9614, as.list(c(c(state=state),mmkh_ci(YearFreq_rel_intensity_state_avg[21:39],c=0.95))))
}
# Estimate the annual trend for the entire region:
yearly_trends_domain_avg <- mmkh(YearFreq_rel_domain_avg,c=0.95)
# Estimate the annual trend for the entire region:
yearly_trends_intensity_domain_avg <- mmkh_ci(pct_islightFZRA_domain_avg,c=0.95)
yearly_trends_intensity_domain_avg_9614 <- mmkh_ci(pct_islightFZRA_domain_avg[21:39],c=0.95)

writeMat("yearly_trends_intensity.mat",
         yearly_trends_intensity_data=yearly_trends_intensity,
         yearly_trends_intensity_9614_data=yearly_trends_intensity_9614,
         yearly_trends_intensity_state_avg=yearly_trends_intensity_state_avg, 
         yearly_trends_intensity_state_avg_9614=yearly_trends_intensity_state_avg_9614,
         yearly_trends_intensity_domain_avg=yearly_trends_intensity_domain_avg, 
         yearly_trends_intensity_domain_avg_9614=yearly_trends_intensity_domain_avg_9614,
         fixNames=T)
