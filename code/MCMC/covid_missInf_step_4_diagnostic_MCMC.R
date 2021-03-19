library(data.table)


# to delete
load('C:/Users/rachaelpung/Desktop/Secured thumb/Projects/COVID-19/Missed infections/output processed/outputTotalChain_100000_iter_period_5.RData')
outputTotal = copy(outputTotalChain_period_5)

# autocorrelation plot
par(mfrow=c(2,2))

for(param in 1:4){ 
  lag = acf(outputTotal$param[,param], lag.max = 1000, plot = FALSE)$lag
  autocorrelation = acf(outputTotal$param[,param], lag.max = 1000, plot = FALSE)$acf
  plot(lag, autocorrelation, col='#00468BFF', type = 'l')
}

# effective sample size
n = dim(outputTotal$param)[1]
acf = acf(outputTotal$param[,1], lag.max = 10000, plot = FALSE)$acf
acf_sums = acf[1:999] + acf[2:1000]
denominator = 1 + 2*sum(acf[2:(which(acf_sums < 0)[1]+1)])

n_eff = n/denominator

n = length(outputTotal$param[seq(1,100000,10),1])
acf = acf(outputTotal$param[seq(1,100000,10),1], lag.max = 1000, plot = FALSE)$acf
acf_sums = acf[1:999] + acf[2:1000]
denominator = 1 + 2*sum(acf[2:(which(acf_sums < 0)[1]+1)])

n_eff = n/denominator

library(coda)

param = mcmc(outputTotal$param[seq(1,1450000,10),4])
effectiveSize(param)






# Gelman-Rubin convergence diagnostic
compute_Rhat <- function(chains){
  
  chains_split <- rbind(chains[,1: (ncol(chains)/2)],
                        chains[,((ncol(chains)/2)+1):ncol(chains)])
  
  n = dim(chains_split)[2]
  m = dim(chains_split)[1]
  
  Ws = sapply(1:m, function(x){ var(chains_split[x,]) })
  
  W = mean(Ws)
  
  chain_means <- sapply(1:m, function(x){ mean(chains_split[x,]) })
  
  B = var(chain_means)
  
  V = (n-1)/n*W + B
  
  Rhat = sqrt(V/W)
  
  return(Rhat)
}





load('outputTotalChain_1_30000_iter_period_1.RData')
load('outputTotalChain_2_30000_iter_period_1.RData')

chains = rbind(outputTotalChain_1_period_1$param[,1],
             outputTotalChain_2_period_1$param[,1])

total_n = dim(outputTotalChain_1_period_1$param)[1]
terminal_n = seq(10, total_n, 50)

Rhats = sapply(1:length(terminal_n), function(x){
 
  truncated_chains = chains[,1:terminal_n[x]]
  Rhats[i] = compute_Rhat(truncated_chains)
  
})


# acceptance rate
compute_accRate <- function(chain){
  
  n = length(chain)
  acceptance_rate = (uniqueN(chain)-1)/(n-1)
  return(acceptance_rate)
  
}
