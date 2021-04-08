library(readr)
library(data.table)
library(fitdistrplus)

set.seed(123)

# load distribution of duration from infection to isolation
distIsolate = data.table(read_csv('input/dist_Isolate.csv')) 
distIsolate = distIsolate[PERIOD == 3,]
x = sample(distIsolate[,DAY], 5000, replace = T, prob = distIsolate[,PROB_ISOLATE])

# load distribution of duration from arrival to isolation
distArrivalIsolate = data.table(read_csv('input/dist_Arrival_Isolate.csv')) 
x = sample(distArrivalIsolate[,DAY_SINCE_ARRIVAL], 5000, replace = T, prob = distArrivalIsolate[,PROB_ISOLATE_DAY_SINCE_ARRIVAL])


plotdist(x, histo = TRUE, demp = TRUE)
descdist(x, discrete=FALSE, boot=500)

fit_nb <- fitdist(x, 'nbinom')
fit_p  <- fitdist(x, 'pois')

summary(fit_nb)
summary(fit_p)

fit_w  <- fitdist(x, 'weibull')
fit_g  <- fitdist(x, 'gamma')
fit_ln <- fitdist(x, 'lnorm')

summary(fit_w)
summary(fit_g)
summary(fit_ln)

par(mfrow=c(2,2))
plot.legend <- c('negativebinomial', 'poisson')
denscomp(list(fit_nb, fit_p), legendtext = plot.legend)
cdfcomp (list(fit_nb, fit_p), legendtext = plot.legend)
qqcomp  (list(fit_nb, fit_p), legendtext = plot.legend)
ppcomp  (list(fit_nb, fit_p), legendtext = plot.legend)

par(mfrow=c(2,2))
plot.legend <- c('Weibull', 'lognormal', 'gamma')
denscomp(list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
cdfcomp (list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_ln), legendtext = plot.legend)

distIsolateFit = data.table(DAY = 0:21,
                            PROB_ISOLATE_WEIBULL = dweibull(0:21, shape = summary(fit_w)$estimate[1], scale = summary(fit_w)$estimate[2]),
                            PROB_ISOLATE_GAMMA = dgamma(0:21, shape = summary(fit_g)$estimate[1], rate = summary(fit_g)$estimate[2]),
                            PROB_ISOLATE_LOG_NORMAL = dlnorm(0:21, meanlog = summary(fit_ln)$estimate[1], sdlog = summary(fit_ln)$estimate[2])
)

distIsolateFit[, `:=` (PROB_ISOLATE_WEIBULL    = PROB_ISOLATE_WEIBULL/sum(PROB_ISOLATE_WEIBULL),
                       PROB_ISOLATE_GAMMA      = PROB_ISOLATE_GAMMA/sum(PROB_ISOLATE_GAMMA),
                       PROB_ISOLATE_LOG_NORMAL = PROB_ISOLATE_LOG_NORMAL/sum(PROB_ISOLATE_LOG_NORMAL))]

distIsolateFitParam = data.table(DIST = c('Weibull', 'Gamma', 'Lognormal'),
                                 ESTIMATE_1 = c(summary(fit_w)$estimate[1], summary(fit_g)$estimate[1], summary(fit_ln)$estimate[1]),
                                 ESTIMATE_2 = c(summary(fit_w)$estimate[2], summary(fit_g)$estimate[2], summary(fit_ln)$estimate[2]))

write.csv(distIsolateFit, 'input/dist_IsolateFit_period3.csv',
          row.names = F)
write.csv(distIsolateFitParam, 'input/dist_IsolateFitParam_period3.csv',
          row.names = F)


distArrivalIsolateFit = data.table(DAY = 0:21,
                                   PROB_ISOLATE_NB = dnbinom(0:21, size = summary(fit_nb)$estimate[1], mu = summary(fit_nb)$estimate[2]),
                                   PROB_ISOLATE_POISSON = dpois(0:21, lambda = summary(fit_p)$estimate[1])
)

distArrivalIsolateFit[, `:=` (PROB_ISOLATE_NB = PROB_ISOLATE_NB/sum(PROB_ISOLATE_NB),
                              PROB_ISOLATE_POISSON = PROB_ISOLATE_POISSON/sum(PROB_ISOLATE_POISSON))]

distArrivalIsolateFitParam = data.table(DIST = c('NegativeBinomial', 'Poisson'),
                                        ESTIMATE_1 = c(summary(fit_nb)$estimate[1], summary(fit_p)$estimate[1]),
                                        ESTIMATE_2 = c(summary(fit_nb)$estimate[2], summary(fit_p)$estimate[2]))


write.csv(distArrivalIsolateFit, 'C:/Users/rachaelpung/Desktop/Secured thumb/Projects/COVID-19/Missed infections/input/dist_Arrival_IsolateFit.csv',
          row.names = F)
write.csv(distArrivalIsolateFitParam, 'C:/Users/rachaelpung/Desktop/Secured thumb/Projects/COVID-19/Missed infections/input/dist_Arrival_IsolateFitParam.csv',
          row.names = F)
