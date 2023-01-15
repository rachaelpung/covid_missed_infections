source('codes/covid_missInf_load_library.R')
load('output processed/inc_all.RData')

# prob of sero conversion
prob_sero_convert = 0.871

# pdf incubation period
dist_incub = data.table(time_since_inf = 1:21)
dist_incub[, prob_onset_since_inf := dlnorm(time_since_inf, meanlog = 1.63, sdlog = 0.50, log = FALSE)]
dist_incub[, prob_onset_since_inf := prob_onset_since_inf/sum(prob_onset_since_inf)]

# prob dist of being S+ since time of onset
dist_prob_sero_since_onset = read_csv('input/dist_Sero_Positivity.csv')
fit_spline = gam(PERCENTAGE_POSITIVE~s(TIME), data = dist_prob_sero_since_onset)
predict_spline = predict(fit_spline,newdata=data.frame(TIME=0:60))/100
predict_spline[predict_spline>1] = 1
dist_prob_sero_since_onset = data.table(time_since_onset = 0:60, prob_sero_since_onset = as.vector(predict_spline))

# assume antibody level will decline after 2 months (~8 weeks)
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0244126
# https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(21)00354-6/fulltext  

# Not used
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0255208


# interpolate between points
dist_prob_sero_since_onset_interpolate = data.table(time_since_onset = c(60,90,120,150,180,210,240,270,330),
                                                    prob_sero_since_onset = c(1,0.58,0.58,0.58,0.58, 0.56,0.56,0.56,0.19))


fit_spline = smooth.spline(dist_prob_sero_since_onset_interpolate$time_since_onset, 
                           dist_prob_sero_since_onset_interpolate$prob_sero_since_onset,
                           df=5)

predict_spline = predict(fit_spline,61:330)$y
dist_prob_sero_since_onset_interpolate = data.table(time_since_onset = 61:330, prob_sero_since_onset = as.vector(predict_spline))
dist_prob_sero_since_onset = rbind(dist_prob_sero_since_onset,dist_prob_sero_since_onset_interpolate)

# convert time since onset to time since infection
dist_prob_sero_since_onset = dist_prob_sero_since_onset[rep(seq_len(.N), times=21),]
dist_prob_sero_since_onset[, time_inf_to_onset:=rep(1:21, each=.N/21)]
dist_prob_sero_since_onset[, time_since_inf := time_inf_to_onset+time_since_onset]
dist_prob_sero_since_onset[, prob_onset_since_inf := rep(dist_incub[,prob_onset_since_inf], each = .N/21)]
dist_prob_sero_since_onset[, prob_sero_since_inf:=prob_onset_since_inf*prob_sero_since_onset]

# prob dist of being S+ since time of infection if successfully seroconvert
dist_prob_sero_since_inf = dist_prob_sero_since_onset[,.(sum(prob_sero_since_inf)), by=.(time_since_inf)]
dist_prob_sero_since_inf = dist_prob_sero_since_inf[time_since_inf <= 330]
setnames(dist_prob_sero_since_inf, old = 'V1', new = 'prob_sero_since_inf')
dist_prob_sero_since_inf = rbind(data.table(time_since_inf=-14:0, prob_sero_since_inf = 0), dist_prob_sero_since_inf)


# validate against population sero survey
doy_test_start = as.numeric(as.Date("2020-09-07") - as.Date('2019-12-31')) 
doy_test_end   = as.numeric(as.Date("2020-10-31") - as.Date('2019-12-31'))
duration_test = doy_test_end - doy_test_start + 1

# prob of being S+ days of testing
prob_sero_pos = data.table(time_inf = rep(1:doy_test_end, each=duration_test),
                           doy_test = rep(doy_test_start:doy_test_end, times=doy_test_end),
                           prob_doy_test = 1/117,
                           prob_sero_convert = prob_sero_convert)

prob_sero_pos[,time_since_inf:= doy_test-time_inf]
prob_sero_pos[dist_prob_sero_since_inf, prob_sero_since_inf:=prob_sero_since_inf, on=c(time_since_inf='time_since_inf')]
prob_sero_pos[is.na(prob_sero_since_inf),prob_sero_since_inf:=0 ]

prob_sero_pos[,detected:=prob_doy_test*prob_sero_convert*prob_sero_since_inf]


prob_sero_pos = prob_sero_pos[,.(sum(detected)), by=.(time_inf)]$V1


# confidence intervals of observed seroprevalence
binom_CI <- function(events, #events = outcomes
                     trials, #number of individuals, test, etc
                     alpha = 0.05){
  
  n <- trials
  x <- events
  p_hat <- x/n
  
  # Calculate upper and lower limit
  upper_lim <- (p_hat + (qnorm(1-(alpha/2))^2/(2*n)) + qnorm(1-(alpha/2)) * sqrt(((p_hat*(1-p_hat))/n) + (qnorm(1-(alpha/2))^2/(4*n^2))))/(1 + (qnorm(1-(alpha/2))^2/(n)))
  
  lower_lim <- (p_hat + (qnorm(alpha/2)^2/(2*n)) + qnorm(alpha/2) * sqrt(((p_hat*(1-p_hat))/n) + (qnorm(alpha/2)^2/(4*n^2))))/(1 + (qnorm(alpha/2)^2/(n)))
  
  
  # Modification for probabilities close to boundaries
  
  if ((n <= 50 & x %in% c(1, 2)) | (n >= 51 & n <= 100 & x %in% c(1:3))) {
    lower_lim <- 0.5 * qchisq(alpha, 2 * x)/n
  }
  
  if ((n <= 50 & x %in% c(n - 1, n - 2)) | (n >= 51 & n <= 100 & x %in% c(n - (1:3)))) {
    upper_lim <- 1 - 0.5 * qchisq(alpha, 2 * (n - x))/n
  }
  
  out <- c(lower_lim,upper_lim)
  names(out) <- c("lower_CI","upper_CI")
  return(out)
  
}

(2/1578)*100
binom_CI(events = 2, trials = 1578)*100

# 0.1267427
# 
# lower_CI   upper_CI 
# 0.03476431 0.46095281 





# validate against model using linked and unlinked cases
model = inc_all$wild_neg_binom_mod_daily_local_N_unlinked[,1:doy_test_end] +
  inc_all$wild_neg_binom_mod_daily_local_N_linked[,1:doy_test_end] +
  inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,1:doy_test_end] +
  inc_all$wild_neg_binom_mod_daily_local_M_linked[,1:doy_test_end]

model_sero_pos = model %*% prob_sero_pos
quantile(model_sero_pos*100/5300000, probs = c(0.5,0.025,0.975))

# 50%       2.5%      97.5% 
# 0.03087959 0.02106481 0.06106685  
# vs
# 0.1267427 0.03476431 0.46095281

# validate against model using all cases
model = inc_all$wild_neg_binom_total_mod_daily_local_N_unlinked[,1:doy_test_end] +
  inc_all$wild_neg_binom_total_mod_daily_local_N_linked[,1:doy_test_end] +
  inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked[,1:doy_test_end] +
  inc_all$wild_neg_binom_total_mod_daily_local_M_linked[,1:doy_test_end]

model_sero_pos = model %*% prob_sero_pos
quantile(model_sero_pos*100/5300000, probs = c(0.5,0.025,0.975))

# 50%       2.5%      97.5% 
# 0.07778554 0.04822535 0.18350233 
# vs
# 0.1267427 0.03476431 0.46095281

# quantile(model_sero_pos*100/5300000 - (2/1578)*100, probs = c(0.5,0.025,0.975))

























### ASSUME NO WANING IMMUNITY 330 DAYS SINCE INFECTION

# prob of sero conversion
prob_sero_convert = 0.871

# pdf incubation period
dist_incub = data.table(time_since_inf = 1:21)
dist_incub[, prob_onset_since_inf := dlnorm(time_since_inf, meanlog = 1.63, sdlog = 0.50, log = FALSE)]
dist_incub[, prob_onset_since_inf := prob_onset_since_inf/sum(prob_onset_since_inf)]

# prob dist of being S+ since time of onset
dist_prob_sero_since_onset = read_csv('input/dist_Sero_Positivity.csv')
fit_spline = gam(PERCENTAGE_POSITIVE~s(TIME), data = dist_prob_sero_since_onset)
predict_spline = predict(fit_spline,newdata=data.frame(TIME=0:60))/100
predict_spline[predict_spline>1] = 1
dist_prob_sero_since_onset = data.table(time_since_onset = 0:60, prob_sero_since_onset = as.vector(predict_spline))

# assume antibody level will decline after 2 months (~8 weeks)
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0244126
# https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(21)00354-6/fulltext  

# Not used
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0255208


# interpolate between points
dist_prob_sero_since_onset_interpolate = data.table(time_since_onset = c(60,90,120,150,180,210,240,270,330),
                                                    prob_sero_since_onset = c(1,1,1,1,1,1,1,1,1))


fit_spline = smooth.spline(dist_prob_sero_since_onset_interpolate$time_since_onset, 
                           dist_prob_sero_since_onset_interpolate$prob_sero_since_onset,
                           df=5)

predict_spline = predict(fit_spline,61:330)$y
dist_prob_sero_since_onset_interpolate = data.table(time_since_onset = 61:330, prob_sero_since_onset = as.vector(predict_spline))
dist_prob_sero_since_onset = rbind(dist_prob_sero_since_onset,dist_prob_sero_since_onset_interpolate)

# convert time since onset to time since infection
dist_prob_sero_since_onset = dist_prob_sero_since_onset[rep(seq_len(.N), times=21),]
dist_prob_sero_since_onset[, time_inf_to_onset:=rep(1:21, each=.N/21)]
dist_prob_sero_since_onset[, time_since_inf := time_inf_to_onset+time_since_onset]
dist_prob_sero_since_onset[, prob_onset_since_inf := rep(dist_incub[,prob_onset_since_inf], each = .N/21)]
dist_prob_sero_since_onset[, prob_sero_since_inf:=prob_onset_since_inf*prob_sero_since_onset]

# prob dist of being S+ since time of infection if successfully seroconvert
dist_prob_sero_since_inf = dist_prob_sero_since_onset[,.(sum(prob_sero_since_inf)), by=.(time_since_inf)]
dist_prob_sero_since_inf = dist_prob_sero_since_inf[time_since_inf <= 330]
setnames(dist_prob_sero_since_inf, old = 'V1', new = 'prob_sero_since_inf')
dist_prob_sero_since_inf = rbind(data.table(time_since_inf=-14:0, prob_sero_since_inf = 0), dist_prob_sero_since_inf)


# validate against population sero survey
doy_test_start = as.numeric(as.Date("2020-09-07") - as.Date('2019-12-31')) 
doy_test_end   = as.numeric(as.Date("2020-10-31") - as.Date('2019-12-31'))
duration_test = doy_test_end - doy_test_start + 1

# prob of being S+ days of testing
prob_sero_pos = data.table(time_inf = rep(1:doy_test_end, each=duration_test),
                           doy_test = rep(doy_test_start:doy_test_end, times=doy_test_end),
                           prob_doy_test = 1/117,
                           prob_sero_convert = prob_sero_convert)

prob_sero_pos[,time_since_inf:= doy_test-time_inf]
prob_sero_pos[dist_prob_sero_since_inf, prob_sero_since_inf:=prob_sero_since_inf, on=c(time_since_inf='time_since_inf')]
prob_sero_pos[is.na(prob_sero_since_inf),prob_sero_since_inf:=0 ]

prob_sero_pos[,detected:=prob_doy_test*prob_sero_convert*prob_sero_since_inf]


prob_sero_pos = prob_sero_pos[,.(sum(detected)), by=.(time_inf)]$V1

# validate against model using linked and unlinked cases
model = inc_all$wild_neg_binom_mod_daily_local_N_unlinked[,1:doy_test_end] +
  inc_all$wild_neg_binom_mod_daily_local_N_linked[,1:doy_test_end] +
  inc_all$wild_neg_binom_mod_daily_local_M_unlinked[,1:doy_test_end] +
  inc_all$wild_neg_binom_mod_daily_local_M_linked[,1:doy_test_end]

model_sero_pos = model %*% prob_sero_pos
quantile(model_sero_pos*100/5300000, probs = c(0.5,0.025,0.975))

# 50%       2.5%      97.5% 
# 0.04940970 0.03443397 0.09851814  
# vs
# 0.1267427 0.03476431 0.46095281

quantile(model_sero_pos*100/5300000 - (2/1578)*100, probs = c(0.5,0.025,0.975))


# validate against model using all cases
model = inc_all$wild_neg_binom_total_mod_daily_local_N_unlinked[,1:doy_test_end] +
  inc_all$wild_neg_binom_total_mod_daily_local_N_linked[,1:doy_test_end] +
  inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked[,1:doy_test_end] +
  inc_all$wild_neg_binom_total_mod_daily_local_M_linked[,1:doy_test_end]

model_sero_pos = model %*% prob_sero_pos
quantile(model_sero_pos*100/5300000, probs = c(0.5,0.025,0.975))

# 50%      2.5%     97.5% 
# 0.12747591 0.07910323 0.30894649 
# vs
# 0.1267427 0.03476431 0.46095281

quantile(model_sero_pos*100/5300000 - (2/1578)*100, probs = c(0.5,0.025,0.975))


(2/1578)*100
binom_CI(events = 2, trials = 1578)*100

# 0.1267427
# 
# lower_CI   upper_CI 
# 0.03476431 0.46095281 

length(model_sero_pos)
bootstrap_sero_pos = sapply(1:23200, function(x){
  sum(sample(c(1,1, rep(0,1576)), 1578, replace=TRUE))*100/1578
})

quantile(model_sero_pos*100/5300000 - bootstrap_sero_pos, probs = c(0.5,0.025,0.975))
