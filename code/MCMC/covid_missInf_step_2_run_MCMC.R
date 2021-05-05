library(MASS)
library(readr)
library(data.table)
library(doParallel)

source('codes/v5/MCMC/covid_missInf_functions_simulations.r')

# load distribution of generation interval over time
distGen = data.table(read_csv('input/dist_Gen_Interval.csv')) # gamma distributed, mu 7.5, sigma 3.4
  
# load distribution of duration from infection to isolation
distIsolate = data.table(read_csv('input/dist_Isolate.csv')) 

# load distribution of time of infection for imported cases
distTravel = data.table(read_csv('input/dist_Travel.csv')) 
  
# load observed data
dataCaseImportNotified = data.table(read_csv('input/dataCaseImportNotified.csv')) 
dataIncImport = data.table(read_csv('input/dataIncImport.csv')) 
dataInc = data.table(read_csv('input/dataInc.csv')) 

# load initalParam and simulated spillover infections
initialParam = read_csv("input/initialParam.csv")
dataSpilloverOffspring = data.table(read_csv("output processed/20210223/dataSpilloverOffspringSummary.csv"))
dataSpilloverDateIsolate = data.table(read_csv("output processed/20210223/dataSpilloverDateIsolateSummary.csv"))
  
# start MCMC runs
setSeed = scan('input/setSeed.txt')

## set up clusters
cl = makeCluster(10) #detectCores()
registerDoParallel(cl)

foreach(paramSet = 1:10, .packages = c("data.table", "mgcv", "MASS"), .combine = "c") %dopar%  {
#nrow(initalParam)

# attempt to recover parameters used to simulate the above dataset
set.seed(setSeed[paramSet])
  
# set time periods
iterN = 15000 
timePeriod = 5 # period of interest
nParam = 4 # no of param for the period of interest
  
#  60 (29 Feb),   97 (6 Apr), 170 (18 Jun),  244(31 Aug), 367(1 Jan 2021)
#  18 (18 Jan),   61 (1 Mar),   98 (7 Apr), 171 (19 Jun),  245 (1 Sep)
startTimePeriod = c(18, 61, 98, 171, 245)  
endTimePeriod =  c(60, 97, 170, 244, 367)
  

# load distribution of spillover infections in previous time period
dataSpilloverOffspring = dataSpilloverOffspring[TIMEPERIOD == timePeriod-1]
dataSpilloverOffspring[,TIMEPERIOD := NULL]
  
dataSpilloverDateIsolate = dataSpilloverDateIsolate[TIMEPERIOD == timePeriod-1]
dataSpilloverDateIsolate[,TIMEPERIOD := NULL]

# allocate memory
nameParam = c(paste('mu_R_period_', timePeriod, sep = ''),
              paste('effOffspringParentNotified_period_', timePeriod, sep = ''),
              paste('effOffspringParentMissed_period_', timePeriod, sep = ''),
              paste('ratioImportMissed_period_', timePeriod, sep = ''))

output = list()
output$param = matrix(0, iterN, nParam)
colnames(output$param) = nameParam

output$LL = rep(0,iterN)
output$startTimePeriod = startTimePeriod[timePeriod]  
output$endTimePeriod = endTimePeriod[timePeriod]
output$time = startTimePeriod[timePeriod]:endTimePeriod[timePeriod]
  
output$IncNotifiedUnlinkedObs = matrix(0, 1, length(output$time)) 
output$IncNotifiedLinkedObs   = matrix(0, 1, length(output$time)) 

# store output by day of isolation
output$IncNotifiedUnlinkedSim = matrix(0, iterN, length(output$time)) 
output$IncNotifiedLinkedSim   = matrix(0, iterN, length(output$time)) 

# store output by day of infection
output$IncNotifiedUnlinkedSim_DOI = matrix(0, iterN, length(output$time)) 
output$IncNotifiedLinkedSim_DOI   = matrix(0, iterN, length(output$time))
output$IncMissedUnlinkedSim_DOI   = matrix(0, iterN, length(output$time))
output$IncMissedLinkedSim_DOI     = matrix(0, iterN, length(output$time))

output$SpilloverOffspring   = matrix(0, iterN, 2*21)
output$SpilloverDateIsolate = matrix(0, iterN, 3*21)

# initialise parameters
currentParam = list()
currentParam$param = as.matrix(initialParam[paramSet,((timePeriod-1)*4+1):((timePeriod-1)*4+nParam)])
currentParam$LL = NA
currentParam$time = output$time
currentParam$lenTime = length(output$time)
currentParam$startTimePeriod = output$startTimePeriod
currentParam$endTimePeriod = output$endTimePeriod

# imported cases
dataCaseImportNotified = dataCaseImportNotified[DATE_ARRIVAL %in% currentParam$time]
dataIncImport = dataIncImport[TIME %in% currentParam$time]

distIsolate = distIsolate[PERIOD == timePeriod,]

# simulate daily incidence
currentInc = simCommInc(dataCaseImportNotified, dataIncImport, dataSpilloverOffspring, dataSpilloverDateIsolate,
                        currentParam, distGen, distIsolate, distTravel)

# observed daily incidence
dataInc = dataInc[TIME %in% currentParam$time]

# calculate log likelihood
currentParam$LL = calculateLogLik(dataInc, currentInc$currentIncIsolate, currentParam$lenTime)

# store currentParam and timeSeriesInc into output
output$param[1,] = currentParam$param 
output$LL[1]    = currentParam$LL 

output$IncNotifiedUnlinkedObs[1,] = dataInc$NOTIFIED_COMM_UNLINKED
output$IncNotifiedLinkedObs[1,]   = dataInc$NOTIFIED_COMM_LINKED

output$IncNotifiedUnlinkedSim[1,] = currentInc$currentIncIsolate$NOTIFIED_COMM_UNLINKED
output$IncNotifiedLinkedSim[1,]   = currentInc$currentIncIsolate$NOTIFIED_COMM_LINKED

output$IncNotifiedUnlinkedSim_DOI[1,] = currentInc$currentIncInfection$NOTIFIED_COMM_UNLINKED
output$IncNotifiedLinkedSim_DOI[1,]   = currentInc$currentIncInfection$NOTIFIED_COMM_LINKED
output$IncMissedUnlinkedSim_DOI[1,]   = currentInc$currentIncInfection$MISSED_COMM_UNLINKED
output$IncMissedLinkedSim_DOI[1,]     = currentInc$currentIncInfection$MISSED_COMM_LINKED

output$SpilloverOffspring[1,]     = currentInc$spilloverOffspring$OFFSPRING
output$SpilloverDateIsolate[1,]   = c(currentInc$spilloverDateIsolate$TIME, currentInc$spilloverDateIsolate$NOTIFIED_COMM_UNLINKED, currentInc$spilloverDateIsolate$NOTIFIED_COMM_LINKED)
  

lambda = 0.01
sd = rep(1,nParam); covMatrix = diag(sd)
acceptRate = 0.234 ;  acceptedN = 0

# start MCMC
for(iter in 2:iterN){ #iterN
  
  print(sprintf("Iteration %s of %.f. Acceptance rate %.3f", iter, iterN, acceptRate))

  # MH for current transmission parameter
  oldParam               = currentParam
  oldInc                 = currentInc
  
  currentParam$param[1,1:nParam] = mvrnorm(1, currentParam$param[1,1:nParam], Sigma=lambda*covMatrix)
    
  iterOutput             = mh(oldParam,currentParam, 
                              dataInc, oldInc,
                              dataCaseImportNotified, dataIncImport, dataSpilloverOffspring, dataSpilloverDateIsolate)
    
  currentParam      =  iterOutput[[1]]
  acceptedN         =  acceptedN + iterOutput[[2]]
  currentInc        =  iterOutput[[3]]   

  if(iter >= 100){
    # start adaptation after 100 iterations
    acceptRate = acceptedN/(iter-1)
    lambda = max(0.01,min(1,exp(log(lambda)+(acceptRate-0.234)*0.999^iter))) 
    covMatrix  = cov(output$param[2:iter,1:nParam])

  }
  
  output$param[iter,] = currentParam$param 
  output$LL[iter]    = currentParam$LL 
  
  output$IncNotifiedUnlinkedSim[iter,] = currentInc$currentIncIsolate$NOTIFIED_COMM_UNLINKED
  output$IncNotifiedLinkedSim[iter,]   = currentInc$currentIncIsolate$NOTIFIED_COMM_LINKED
  
  output$IncNotifiedUnlinkedSim_DOI[iter,] = currentInc$currentIncInfection$NOTIFIED_COMM_UNLINKED
  output$IncNotifiedLinkedSim_DOI[iter,]   = currentInc$currentIncInfection$NOTIFIED_COMM_LINKED
  output$IncMissedUnlinkedSim_DOI[iter,]   = currentInc$currentIncInfection$MISSED_COMM_UNLINKED
  output$IncMissedLinkedSim_DOI[iter,]     = currentInc$currentIncInfection$MISSED_COMM_LINKED
  
  output$SpilloverOffspring[iter,]     = currentInc$spilloverOffspring$OFFSPRING
  output$SpilloverDateIsolate[iter,]   = c(currentInc$spilloverDateIsolate$TIME, currentInc$spilloverDateIsolate$NOTIFIED_COMM_UNLINKED, currentInc$spilloverDateIsolate$NOTIFIED_COMM_LINKED)
  
  
  if(iter %% 5000 == 0){
    save(output,file= paste0("output raw/outputTest_", paramSet, "_iter_", iter, "_timeperiod_", timePeriod,  ".rdata"))
  }
  
  
}

#acceptedN/iterN

#save(output,file= paste0("output/outputTest_", paramSet, ".rdata"))

paramSet

}

## stop clusters
stopCluster(cl)
