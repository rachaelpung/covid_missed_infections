library(readr)
library(scales)
library(data.table)

# load distribution of generation interval
distGen = data.table(read_csv('input/dist_Gen_Interval.csv')) # gamma distributed, mu 7.5, sigma 3.4

# calculate available portions of the generation interval if arrive after start of infectiousness or isolated before end of infectiousness
distGen[, RIGHT_TRUNCATE_GEN_INTERVAL := sapply(1:.N, function(x){ sum(GEN_INTERVAL[x:.N]) })]  # for all imported cases
distGen[, LEFT_TRUNCATE_GEN_INTERVAL := cumsum(GEN_INTERVAL)] # for all notified and isolated cases 

# load distribution of duration from infection to isolation
distIsolate = data.table(read_csv('input/dist_IsolateFit.csv')) 
distIsolate = distIsolate[PERIOD == 2,]
distIsolate[, PROB_ISOLATE := PROB_ISOLATE_WEIBULL]

# load distribution of duration from infection to travel
distTravel = data.table(read_csv('input/dist_Travel.csv')) 

# load distribution of duration from arrival to isolation
distArrivalIsolate = data.table(read_csv('input/dist_Arrival_Isolate.csv')) 

# truncated R for notified community cases due to early isolation
R_notified_comm = sum(distIsolate[, PROB_ISOLATE]*distGen[, LEFT_TRUNCATE_GEN_INTERVAL])

# truncated R for missed imported cases as they could be infectious before arrival
R_missed_imported = sum(distTravel[, PROB_INFECTION]*distGen[DAY %in% c(1:14), RIGHT_TRUNCATE_GEN_INTERVAL])

# truncated R for notified imported cases as they could be infectious before arrival and due to early isolation
convolution_gen_interval = data.table(DAY_INFECTED_BEFORE_ARRIVAL = rep(1:14, each = 22),
                                      PROB_INFECTION = rep(distTravel[, PROB_INFECTION], each =22),
                                      DAY_ISOLATE_SINCE_ARRIVAL = rep(0:21, times = 14),
                                      PROB_ISOLATE = rep(distArrivalIsolate[,PROB_ISOLATE_DAY_SINCE_ARRIVAL], times = 14))

convolution_gen_interval[, `:=` (START_TRUNCATED_GEN_INTERVAL = DAY_INFECTED_BEFORE_ARRIVAL,
                                 END_TRUNCATED_GEN_INTERVAL  = DAY_INFECTED_BEFORE_ARRIVAL + DAY_ISOLATE_SINCE_ARRIVAL)]

convolution_gen_interval[END_TRUNCATED_GEN_INTERVAL > 21, END_TRUNCATED_GEN_INTERVAL:=21]
convolution_gen_interval[, RIGHT_LEFT_TRUNCATED_GEN_INTERVAL := sapply(1:.N, function(x){ 
  distGen[DAY %in% (START_TRUNCATED_GEN_INTERVAL[x]:END_TRUNCATED_GEN_INTERVAL[x]), sum(GEN_INTERVAL)]   })]

R_notified_imported = sum(convolution_gen_interval[,PROB_INFECTION*PROB_ISOLATE*RIGHT_LEFT_TRUNCATED_GEN_INTERVAL])

# test eigenvalue and eigen vector
# effOffspringParentNotified = 0.8
# effOffspringParentMissed = 0.001
# R = 1.5
# R_notified_comm = R_notified_comm*R
# 
# NGMComm =  matrix(c((1-effOffspringParentMissed)*R,  (1-effOffspringParentNotified)*R_notified_comm,
#                      effOffspringParentMissed*R,      effOffspringParentNotified*R_notified_comm),
#                   nrow = 2, ncol = 2, byrow = T)
# 
# eigen(NGMComm)

# calculate eigenvalues and eigenvector of community transmission for all different parameter combinations
effOffspringParentNotified = seq(0,1,0.05)
effOffspringParentMissed   = seq(0,1,0.05)
dataEigen = data.table(expand.grid(effOffspringParentNotified, effOffspringParentMissed))
setnames(dataEigen, c('effOffspringParentNotified', 'effOffspringParentMissed'))

dataEigen = dataEigen[rep(seq_len(.N), times = 3),]
dataEigen[, R := as.numeric(rep(c(1,1.2,1.5), each = .N/3))]
dataEigen[, R_notified_comm := R_notified_comm*R]

eigen = t(sapply(1:dataEigen[, .N], function(x){
  
  NGMComm =  matrix(c((1-dataEigen[x,effOffspringParentMissed])*dataEigen[x,R],  (1-dataEigen[x,effOffspringParentNotified])*dataEigen[x,R_notified_comm],
                      dataEigen[x,effOffspringParentMissed]*dataEigen[x,R],      dataEigen[x,effOffspringParentNotified]*dataEigen[x,R_notified_comm]),
                    nrow = 2, ncol = 2, byrow = T)
  
  eigenNGMComm = eigen(NGMComm)
  index = which.max(eigenNGMComm$values)
  
  unlink = dataEigen[x,effOffspringParentMissed]*dataEigen[x,R]*eigenNGMComm$vectors[1,index]
  link   = dataEigen[x,effOffspringParentNotified]*dataEigen[x,R_notified_comm]*eigenNGMComm$vectors[2,index]
  
  c(eigenNGMComm$values[index], eigenNGMComm$vectors[,index], unlink, link)
  
}))

dataEigen[, eigenvalue := eigen[,1]]
dataEigen[, eigenvectorratio := eigen[,2]/eigen[,3]]
dataEigen[, unlinklinkratio := eigen[,4]/eigen[,5]]

# generate hex code for eigenvalues
pal <- gradient_n_pal(colours = c("#006837", "#1a9850", "#66bd63", "#a6d96a", "#d9ef8b","#ffffbf",
                                  "#fee08b", "#fdae61", "#fdae61", "#d73027", "#a50026"),
                      values  = c(0.5, 0.6, 0.7, 0.8, 0.9, 1,
                                  1.1, 1.2, 1.3, 1.4, 1.5))
dataEigen[, eigenvalue_hex := pal(eigenvalue)]


# for  plotting of eigenvector ratios with values close to 5
dataInterpolate = data.table(effOffspringParentNotified = c(rep(0.2, 201), rep(0.5, 201), rep(0.8, 201)),
                             effOffspringParentMissed = c(seq(0,0.2,0.001), seq(0,0.2,0.001), seq(0,0.2,0.001)),
                             R = 1.5)

dataInterpolate[, R_notified_comm := R_notified_comm*R]

eigen = t(sapply(1:dataInterpolate[, .N], function(x){
  
  NGMComm =  matrix(c((1-dataInterpolate[x,effOffspringParentMissed])*dataInterpolate[x,R],  (1-dataInterpolate[x,effOffspringParentNotified])*dataInterpolate[x,R_notified_comm],
                      dataInterpolate[x,effOffspringParentMissed]*dataInterpolate[x,R],      dataInterpolate[x,effOffspringParentNotified]*dataInterpolate[x,R_notified_comm]),
                    nrow = 2, ncol = 2, byrow = T)
  
  eigenNGMComm = eigen(NGMComm)
  index = which.max(eigenNGMComm$values)
  
  unlink = dataInterpolate[x,effOffspringParentMissed]*dataInterpolate[x,R]*eigenNGMComm$vectors[1,index]
  link   = dataInterpolate[x,effOffspringParentNotified]*dataInterpolate[x,R_notified_comm]*eigenNGMComm$vectors[2,index]
  
  
  c(eigenNGMComm$values[index], eigenNGMComm$vectors[,index], unlink, link)
  
}))

dataInterpolate[, eigenvalue := eigen[,1]]
dataInterpolate[, eigenvectorratio := eigen[,2]/eigen[,3]]
dataInterpolate[, unlinklinkratio := eigen[,4]/eigen[,5]]

# for plotting line of effective R = 1
dataEff_R_1 = data.table(effOffspringParentNotified = seq(0,1,0.05))
dataEff_R_1 = dataEff_R_1[rep(seq_len(.N), times = 3),]
dataEff_R_1[, R := as.numeric(rep(c(1,1.2,1.5), each = .N/3))]
dataEff_R_1[, truncateGen := R_notified_comm]
dataEff_R_1[, effR := 1]
dataEff_R_1[, effOffspringParentMissed := (((effR/R)-1)*(effOffspringParentNotified*truncateGen - (effR/R)))/((effR/R)-truncateGen) ]


# set up datatable of scenarios with varying number of new infections seeded into population
dataImportScenario = data.table(missed = c(5,80), notified = c(95,20))

dataImportScenario[, `:=` (R = 1.5, effOffspringParentNotified = 0.8)]

dataImportScenario[,`:=` (R_notified_comm = R_notified_comm*R,
                          R_missed_imported = R_missed_imported*R,
                          R_notified_imported = R_notified_imported*R
)]

dataImportScenario = dataImportScenario[rep(seq_len(.N), each = 21), ]
dataImportScenario[, effOffspringParentMissed := rep(seq(0,1,0.05), times = .N/21)]
dataImportScenario[, missed_time_exp := 0]
dataImportScenario[, notified_time_exp := 0]
dataImportScenario[, ratio_missed := rep(c('low', 'high'), each =21)]

for(scenario in 1:dataImportScenario[,.N]){
  
  dataGen = numeric()
  
  missedImport   = dataImportScenario[scenario,missed]
  notifiedImport = dataImportScenario[scenario,notified]
  
  tmp_effOffspringParentMissed   = dataImportScenario[scenario, effOffspringParentMissed]
  tmp_effOffspringParentNotified = dataImportScenario[scenario, effOffspringParentNotified]
  
  tmp_R_missed_imported   = dataImportScenario[scenario, R_missed_imported]
  tmp_R_notified_imported = dataImportScenario[scenario, R_notified_imported]
  tmp_R_notified_comm     = dataImportScenario[scenario, R_notified_comm]
  tmp_R                   = dataImportScenario[scenario, R]
  
  caseImport = matrix(c(missedImport, notifiedImport), nrow = 2, ncol = 1, byrow = T)
  
  dataGen = rbind(dataGen,t(caseImport))
  
  # calculate NGM for imported cases
  NGMImport = matrix(c((1-tmp_effOffspringParentMissed)*tmp_R_missed_imported,  (1-tmp_effOffspringParentNotified)*tmp_R_notified_imported,
                       tmp_effOffspringParentMissed*tmp_R_missed_imported,      tmp_effOffspringParentNotified*tmp_R_notified_imported),
                     nrow = 2, ncol = 2, byrow = T)
  
  caseComm = NGMImport %*% caseImport
  dataGen = rbind(dataGen,t(caseComm))
  
  NGMComm =  matrix(c((1-tmp_effOffspringParentMissed)*tmp_R,  (1-tmp_effOffspringParentNotified)*tmp_R_notified_comm,
                      tmp_effOffspringParentMissed*tmp_R,      tmp_effOffspringParentNotified*tmp_R_notified_comm),
                    nrow = 2, ncol = 2, byrow = T)
  
  for(gen in 1:9){
    
    caseComm = NGMComm %*% caseComm
    dataGen = rbind(dataGen,t(caseComm))
    
  }
  
  dataGen = data.table(dataGen)
  setnames(dataGen, c('missed', 'notified'))
  dataGen[, gen := rep(0:10, times = .N/11)]
  
  dataGen[, missed_ratio := c(NA, missed[2:.N]/missed[1:(.N-1)])]
  dataGen[, notified_ratio := c(NA, notified[2:.N]/notified[1:(.N-1)])]
  
  dataGen[, missed_time_exp := abs(missed_ratio-missed_ratio[.N])]
  dataGen[, notified_time_exp := abs(notified_ratio-notified_ratio[.N])]
  
  weightMissed = (0.1 - dataGen[which(missed_time_exp <= 0.1)[1], missed_time_exp])/(dataGen[which(missed_time_exp <= 0.1)[1]-1, missed_time_exp] - dataGen[which(missed_time_exp <= 0.1)[1], missed_time_exp])
  weightNotified = (0.1 - dataGen[which(notified_time_exp <= 0.1)[1], notified_time_exp])/(dataGen[which(notified_time_exp <= 0.1)[1]-1, notified_time_exp] - dataGen[which(notified_time_exp <= 0.1)[1], notified_time_exp])
  
  if(is.na(weightMissed)){
    dataImportScenario[scenario, missed_time_exp := dataGen[which(missed_time_exp <= 0.1), gen][1]]
  } else if(!is.na(weightMissed)){
    dataImportScenario[scenario, missed_time_exp := dataGen[which(missed_time_exp <= 0.1), gen][1] - weightMissed]
  }
  
  if(is.na(weightNotified)){
    dataImportScenario[scenario, notified_time_exp := dataGen[which(notified_time_exp <= 0.1), gen][1]]
  } else if(!is.na(weightNotified)){
    dataImportScenario[scenario, notified_time_exp := dataGen[which(notified_time_exp <= 0.1), gen][1] - weightNotified]
  }
  
}



dataGen = numeric()

for(scenario in c(3,19,24,40)){
  
  missedImport   = dataImportScenario[scenario,missed]
  notifiedImport = dataImportScenario[scenario,notified]
  
  tmp_effOffspringParentMissed   = dataImportScenario[scenario, effOffspringParentMissed]
  tmp_effOffspringParentNotified = dataImportScenario[scenario, effOffspringParentNotified]
  
  tmp_R_missed_imported   = dataImportScenario[scenario, R_missed_imported]
  tmp_R_notified_imported = dataImportScenario[scenario, R_notified_imported]
  tmp_R_notified_comm     = dataImportScenario[scenario, R_notified_comm]
  tmp_R                   = dataImportScenario[scenario, R]
  
  caseImport = matrix(c(missedImport, notifiedImport), nrow = 2, ncol = 1, byrow = T)
  dataGen = rbind(dataGen,t(caseImport))
  
  # calculate NGM for imported cases
  NGMImport = matrix(c((1-tmp_effOffspringParentMissed)*tmp_R_missed_imported,  (1-tmp_effOffspringParentNotified)*tmp_R_notified_imported,
                       tmp_effOffspringParentMissed*tmp_R_missed_imported,      tmp_effOffspringParentNotified*tmp_R_notified_imported),
                     nrow = 2, ncol = 2, byrow = T)
  
  caseComm = NGMImport %*% caseImport
  dataGen = rbind(dataGen,t(caseComm))
  
  NGMComm =  matrix(c((1-tmp_effOffspringParentMissed)*tmp_R,  (1-tmp_effOffspringParentNotified)*tmp_R_notified_comm,
                      tmp_effOffspringParentMissed*tmp_R,      tmp_effOffspringParentNotified*tmp_R_notified_comm),
                    nrow = 2, ncol = 2, byrow = T)
  
  for(gen in 1:9){
    
    caseComm = NGMComm %*% caseComm
    dataGen = rbind(dataGen,t(caseComm))
    
  }
  
}

dataGen = data.table(dataGen)
setnames(dataGen, c('missed', 'notified'))
dataGen[, gen := rep(0:10, times = .N/11)]
dataGen[, scenario := rep(c(3,19,24,40), each = 11)]


