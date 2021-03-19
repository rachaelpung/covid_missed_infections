library(grid)
library(readr)
library(data.table)
source('codes/v5/covid_missInf_function_plot.r')

load('output processed/20210223/outputTotalChain_10000_thinned_iter_period_all.RData')
outputTotalChain = data.table(outputTotalChain_period_all$param)

# convert ratio of missed imported cases to average daily missed imported cases
# set time periods
startTimePeriod = c(18, 61, 98, 171, 245)  
endTimePeriod =  c(60, 97, 170, 244, 367)

# load observed data
dataIncImport = data.table(read_csv('input/dataIncImport.csv')) 

outputTotalChain[, avgDailyMissedImport_period_1 := ratioImportMissed_period_1*mean(dataIncImport[TIME %in% c(startTimePeriod[1]:endTimePeriod[1]),SPLINE])]
outputTotalChain[, avgDailyMissedImport_period_2 := ratioImportMissed_period_2*mean(dataIncImport[TIME %in% c(startTimePeriod[2]:endTimePeriod[2]),SPLINE])]
outputTotalChain[, avgDailyMissedImport_period_3 := ratioImportMissed_period_3*mean(dataIncImport[TIME %in% c(startTimePeriod[3]:endTimePeriod[3]),SPLINE])]
outputTotalChain[, avgDailyMissedImport_period_4 := ratioImportMissed_period_4*mean(dataIncImport[TIME %in% c(startTimePeriod[4]:endTimePeriod[4]),SPLINE])]
outputTotalChain[, avgDailyMissedImport_period_5 := ratioImportMissed_period_5*mean(dataIncImport[TIME %in% c(startTimePeriod[5]:endTimePeriod[5]),SPLINE])]

# reorder columns
outputTotalChain = outputTotalChain[,c(1:3,21,5:7,22,9:11,23,13:15,24,17:19,25)]

# plot four parameters
paneller=function(row = 1,column=1, outputTotal)
{
  xlm=c(0.5,5.5)
  if(row == 1 & column == 1) ylm=c(0,2)
  if(row == 1 & column == 2) ylm=c(0,1)
  if(row == 2 & column == 1) ylm=c(0,1)
  if(row == 2 & column == 2) ylm=c(log(0.001),log(300))
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(row == 1 & column == 1){ param = 1; colMedian=CONFIG$colsDark[5]; col50=CONFIG$colsDark1[5]; col95=CONFIG$colsDark3[5]}
  if(row == 1 & column == 2){ param = 2; colMedian=CONFIG$colsDark[3]; col50=CONFIG$colsDark1[3]; col95=CONFIG$colsDark3[3]}
  if(row == 2 & column == 1){ param = 3; colMedian=CONFIG$colsDark[2]; col50=CONFIG$colsDark1[2]; col95=CONFIG$colsDark3[2]}
  if(row == 2 & column == 2){ param = 4; colMedian=CONFIG$colsDark[1]; col50=CONFIG$colsDark1[1]; col95=CONFIG$colsDark3[1]}
  
  # grid lines
  if(row == 1 & column == 1){
    grid.lines(xlm, c(0.5,0.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1,1),     default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1.5,1.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  if((row == 1 & column == 2) | (row == 2 & column == 1)){
    grid.lines(xlm, c(0.25,0.25), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.5,0.5),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.75,0.75), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  if(row == 2 & column == 2){
    grid.lines(xlm, log(c(0.01,0.01)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(0.1,0.1)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(1,1)),       default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(10,10)),     default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(100,100)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
 
  
  
  for(period in 1:5){
    
    median  = median(outputTotal[[param+4*(period-1)]])
    low_50  = quantile(outputTotal[[param+4*(period-1)]], probs = 0.25)
    upp_50  = quantile(outputTotal[[param+4*(period-1)]], probs = 0.75)
    low_95  = quantile(outputTotal[[param+4*(period-1)]], probs = 0.025)
    upp_95  = quantile(outputTotal[[param+4*(period-1)]], probs = 0.975)
    
    if(row==2 & column == 2 & period != 3){
      median  = log(median(outputTotal[[param+4*(period-1)]]))
      low_50  = log(quantile(outputTotal[[param+4*(period-1)]], probs = 0.25))
      upp_50  = log(quantile(outputTotal[[param+4*(period-1)]], probs = 0.75))
      low_95  = log(quantile(outputTotal[[param+4*(period-1)]], probs = 0.025))
      upp_95  = log(quantile(outputTotal[[param+4*(period-1)]], probs = 0.975))
    } else if(row==2 & column == 2 & period == 3){
      median  = NA
      low_50  = NA
      upp_50  = NA
      low_95  = NA
      upp_95  = NA
    }
    
    
    grid.lines(c(period,period), c(low_95,upp_95), default.units = 'native',gp=gpar(col=col95,lwd=2))
    grid.lines(c(period,period), c(low_50,upp_50), default.units = 'native',gp=gpar(col=col50,lwd=2))
    grid.points(period,median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row == 1 & column == 1) grid.yaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5))
  if(row == 1 & column == 2) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,1,by=0.25))
  if(row == 2 & column == 1) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,1,by=0.25))
  if(row == 2 & column == 2) grid.yaxis(at=log(c(0.001, 0.01, 0.1, 1, 10, 100)),label=c(0.001, 0.01, 0.1, 1, 10, 100))
  
  grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nAug 31', 'Sep 1-\nJan 1'),
             gp = gpar(fontsize=8))
  
  if(row == 2) grid.text('Time period',y=unit(-3,'lines'))
  
  if(row == 1 & column == 1) {grid.text('R',x=unit(-3,'lines'),rot=90)
                              grid.text('A',x=unit(-2.5,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  if(row == 1 & column == 2) {grid.text(bquote(~epsilon[link]),x=unit(-3,'lines'),rot=90)  #~epsilon[op]
                              grid.text('B',x=unit(-2.5,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))}    
  if(row == 2 & column == 1) {grid.text(bquote(~epsilon[unlink]),x=unit(-3,'lines'),rot=90) #~epsilon[op*'\'']
                              grid.text('C',x=unit(-2.5,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))}  
  if(row == 2 & column == 2) {grid.text('Avg daily missed importation',x=unit(-3.5,'lines'),rot=90)
                              grid.text('D',x=unit(-2.5,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))}  
  if(row == 2 & column == 2) grid.text('Lockdown of borders',x=unit(6.5,'lines'),y=unit(7,'lines'),rot=90)
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/dot plot.png',height=16,width=16,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2)))

for(c in 1:2){
  for(r in 1:2){
    paneller(r,c,outputTotalChain)
  }
}

popViewport()
popViewport()
dev.off()

rm(paneller)



# load distribution of duration from infection to isolation
distIsolate =  data.table(read_csv('input/dist_IsolateFit.csv')) 

# load distribution of generation interval
distGen = data.table(read_csv('input/dist_Gen_Interval.csv')) # gamma distributed, mu 7.5, sigma 3.4

# calculate available portions of the generation interval if arrive after start of infectiousness or isolated before end of infectiousness
distGen[, RIGHT_TRUNCATE_GEN_INTERVAL := sapply(1:.N, function(x){ sum(GEN_INTERVAL[x:.N]) })]  # for all imported cases
distGen[, LEFT_TRUNCATE_GEN_INTERVAL := cumsum(GEN_INTERVAL)] # for all notified and isolated cases 

# store outputs
outputTotalChain[, `:=` (R_notified_comm_period_1 = 0,R_notified_comm_period_2 = 0,R_notified_comm_period_3 = 0,
                         R_notified_comm_period_4 = 0,R_notified_comm_period_5 = 0,eff_R_period_1 = 0,
                         eff_R_period_2 = 0,eff_R_period_3 = 0,eff_R_period_4 = 0,eff_R_period_5 = 0)]

for(period in 1:5){
  
  # PDF of isolation probabilities since day of infection
  probIsolate = distIsolate[PERIOD == period, PROB_ISOLATE_WEIBULL]
  
  for(row in 1:nrow(outputTotalChain)){
    
    R = outputTotalChain[[(1+(period-1)*4)]][row]
    effOffspringParentNotified = outputTotalChain[[(2+(period-1)*4)]][row] 
    effOffspringParentMissed = outputTotalChain[[(3+(period-1)*4)]][row]
    
    # truncated R for notified community cases due to early isolation
    R_notified_comm = sum(probIsolate*distGen[, LEFT_TRUNCATE_GEN_INTERVAL])*R
    
    outputTotalChain[[(period+20)]][row] = R_notified_comm
    
    NGMComm =  matrix(c((1-effOffspringParentMissed)*R,  (1-effOffspringParentNotified)*R_notified_comm,
                        effOffspringParentMissed*R,      effOffspringParentNotified*R_notified_comm),
                      nrow = 2, ncol = 2, byrow = T)
    
    outputTotalChain[[(period+25)]][row] = max(eigen(NGMComm)$values)
    
  }
}

# plot R
paneller=function(row = 1,column=1, outputTotal)
{
  xlm=c(0.5,5.5)
  ylm=c(0,2)
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(column == 1){ colMedian=CONFIG$colsDark[5]; col50=CONFIG$colsDark1[5]; col95=CONFIG$colsDark3[5]}
  if(column == 2){ colMedian=CONFIG$colsDark[1]; col50=CONFIG$colsDark1[1]; col95=CONFIG$colsDark3[1]}
  if(column == 3){ colMedian=CONFIG$colsDark[4]; col50=CONFIG$colsDark1[4]; col95=CONFIG$colsDark3[4]}
  
  
  # grid lines
  grid.lines(xlm, c(0.5,0.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  grid.lines(xlm, c(1,1),     default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  grid.lines(xlm, c(1.5,1.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  
  for(period in 1:5){
    
    if(column == 1){ param = 1; diff = 4*(period-1)}
    if(column == 2){ param = 20; diff = period}
    if(column == 3){ param = 25; diff = period}
    
    median  = median(outputTotal[[param+diff]])
    low_50  = quantile(outputTotal[[param+diff]], probs = 0.25)
    upp_50  = quantile(outputTotal[[param+diff]], probs = 0.75)
    low_95  = quantile(outputTotal[[param+diff]], probs = 0.025)
    upp_95  = quantile(outputTotal[[param+diff]], probs = 0.975)
    
    grid.lines(c(period,period), c(low_95,upp_95), default.units = 'native',gp=gpar(col=col95,lwd=2))
    grid.lines(c(period,period), c(low_50,upp_50), default.units = 'native',gp=gpar(col=col50,lwd=2))
    grid.points(period,median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.yaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5))
  grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nAug 31', 'Sep 1-\nJan 1'),
             gp = gpar(fontsize=8))
  
  grid.text('Time period',y=unit(-3,'lines'))
  
  if(column == 1) {grid.text('R',x=unit(-3,'lines'),rot=90)
    grid.text('A',x=unit(-3,'lines'),y=unit(11.9,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(column == 2) {grid.text('R of a notified case',x=unit(-3,'lines'),rot=90)  
    grid.text('B',x=unit(-3,'lines'),y=unit(11.9,'lines'),gp=gpar(fontsize=unit(10,'pt')))}    
  if(column == 3) {grid.text('Effective R',x=unit(-3,'lines'),gp=gpar(fontsize=unit(9,'pt')),rot=90) 
    grid.text('C',x=unit(-3,'lines'),y=unit(11.9,'lines'),gp=gpar(fontsize=unit(10,'pt')))}  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/dot R.png',height=8,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3)))

paneller(1,1,outputTotalChain)
paneller(1,2,outputTotalChain)
paneller(1,3,outputTotalChain)

popViewport()
popViewport()
dev.off()