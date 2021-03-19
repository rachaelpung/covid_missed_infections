genIntDensity <- function(currentParam){
  
  t = 0:21
  
  shape_genInterval = (currentParam$mu_genInterval / currentParam$sigma_genInterval) ^ 2
  rate_genInterval = currentParam$mu_genInterval / (currentParam$sigma_genInterval ^ 2)
  distGen = dgamma(t,shape=shape_genInterval,rate=rate_genInterval)
  
  distGen = distGen/sum(distGen)
  distGen = data.table(DAY = t, 
                       GEN_INTERVAL = distGen)
  
  return(distGen)
  
}

isolateDensity <- function(dataCase){
  
  # extract symptomatic community cases
  dataCaseCommSymp = dataCase[LOCAL_OVERSEAS == 0 & !is.na(DATE_ONSET), .(CASE, DATE_NOTIFICATION, DATE_ONSET, DATE_ISOLATE_FROM_COMM)]
  
  # distribution of incubation period
  t = 1:21
  distIncub = dlnorm(t, meanlog = 1.63, sdlog = 0.50, log = FALSE)
  distIncub = data.table(DAY = t, PROB_ONSET = distIncub/sum(distIncub))
  
  # let 0 be the date of arrival
  distIsolate = data.table()
  
  for(i in 1:dataCaseCommSymp[,.N]){
    
    durationOnsetIso = dataCaseCommSymp[i, DATE_ISOLATE_FROM_COMM] - dataCaseCommSymp[i, DATE_ONSET]
    distIncubCase = copy(distIncub)
    
    if(durationOnsetIso <= -21) next
    if(durationOnsetIso < 0){
      # <0 means onset happens after isolation
      distIncubCase = distIncubCase[-1:durationOnsetIso, ]
      distIncubCase[, PROB_ONSET := PROB_ONSET/sum(PROB_ONSET)]
      
    }
    
    distIncubCase[, DATE_ONSET := dataCaseCommSymp[i, DATE_ONSET]]
    distIncubCase[, DATE_INFECTION := DATE_ONSET - DAY]
    distIncubCase[, DATE_ISOLATE_FROM_COMM := dataCaseCommSymp[i, DATE_ISOLATE_FROM_COMM]]
    distIncubCase[, DURATION_INF_ISO := DATE_ISOLATE_FROM_COMM - DATE_INFECTION]
    
    distIsolate = rbindlist(list(distIsolate, distIncubCase), fill = T, use.names = T)
    
  }
  
  distIsolate = distIsolate[, .(DURATION_INF_ISO, PROB_ONSET)]
  distIsolate = distIsolate[, sum(PROB_ONSET), by = DURATION_INF_ISO]
  setnames(distIsolate, c("DAY", "PROB_ISOLATE"))
  
  distIsolate = distIsolate[DAY <= 21,]
  distIsolate = distIsolate[order(DAY)]
  distIsolate = rbindlist(list(data.table(DAY = 0, PROB_ISOLATE = 0), distIsolate), fill = T, use.names = T)
  
  distIsolate[, PROB_ISOLATE := PROB_ISOLATE/sum(PROB_ISOLATE)]
  distIsolate[, CUM_PROB_ISOLATE := cumsum(distIsolate$PROB_ISOLATE)]
  distIsolate[, PROB_AT_LARGE := 1-CUM_PROB_ISOLATE]
  
  return(distIsolate)
  
}

travelDensity <- function(dataCase){
  
  # extract symptomatic imported case 
  dataCaseImportSymp = dataCase[LOCAL_OVERSEAS == 1 & !is.na(DATE_ONSET), ]
  dataCaseImportSymp[, DAY_ONSET_ARRIVAL := DATE_ONSET - DATE_ARRIVAL]
  
  # distribution of incubation period
  t = 1:60
  distIncub = dlnorm(t, meanlog = 1.63, sdlog = 0.50, log = FALSE)
  distIncub = data.table(DAY = t, PROB_ONSET = distIncub/sum(distIncub))
  
  # let 0 be the date of arrival
  distTimeInf = data.table()
  
  # for each case determine the time of infection and the probability 
  for(i in 1:dataCaseImportSymp[,.N]){  
    
    if(dataCaseImportSymp[i, DAY_ONSET_ARRIVAL] == 0){
      
      distIncubCase = copy(distIncub)
      distIncubCase[, DAY := -1*DAY]
      distTimeInf = rbindlist(list(distTimeInf, distIncubCase), fill = T, use.names = T)
      
      
    } else if(dataCaseImportSymp[i, DAY_ONSET_ARRIVAL] >= 1){
      
      distIncubCase = copy(distIncub)
      distIncubCase = distIncubCase[-(1:dataCaseImportSymp[i, DAY_ONSET_ARRIVAL]),]
      distIncubCase[, DAY := -1:-.N]
      distIncubCase[, PROB_ONSET := PROB_ONSET/sum(PROB_ONSET)]
      
      distTimeInf = rbindlist(list(distTimeInf, distIncubCase), fill = T, use.names = T)
      
    } else if(dataCaseImportSymp[i, DAY_ONSET_ARRIVAL] <= -1){
      
      distIncubCase = copy(distIncub)
      distIncubCase[, DAY := (dataCaseImportSymp[i, DAY_ONSET_ARRIVAL]-1):(dataCaseImportSymp[i, DAY_ONSET_ARRIVAL]-.N)]
      distIncubCase[, PROB_ONSET := PROB_ONSET/sum(PROB_ONSET)]
      
      distTimeInf = rbindlist(list(distTimeInf, distIncubCase), fill = T, use.names = T)
      
    }
  }
  
  distTimeInf = distTimeInf[, sum(PROB_ONSET), by = DAY]
  setnames(distTimeInf, c("DAY", "PROB_INFECTION"))
  
  distTimeInf = distTimeInf[DAY >= -14,]
  distTimeInf[, PROB_INFECTION :=  PROB_INFECTION/sum(PROB_INFECTION)]
  
  return(distTimeInf)
  
}




