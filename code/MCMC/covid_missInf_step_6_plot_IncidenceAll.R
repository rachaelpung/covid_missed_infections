library(grid)
library(readr)
library(data.table)
source('codes/v5/covid_missInf_function_plot.r')

# load outputs
load('output processed/20210223/outputTotalChain_10000_thinned_iter_period_all.rdata')
outputTotal = copy(outputTotalChain_period_all)

# load observed data
dataIncImport = data.table(read_csv("input/dataIncImport.csv"))
dataInc = data.table(read_csv("input/dataInc.csv"))

# generate plots
paneller=function(row = 1,column=1)
{
  
  xlm=c(0, dataInc[,max(TIME)]+0.5)
  if(row == 1 & column == 1) ylm=c(0,100)  # imported cases
  if(row == 1 & column == 2) ylm=c(0,100)  # linked cases
  if(row == 2 & column == 1) ylm=c(0,50)   # unlinked cases
  if(row == 2 & column == 2) ylm=c(0,2500) # missed cases
  
  innermargins = c(2,3,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  t = dataInc[,TIME]
  if(row == 1 & column == 1){col1 =  CONFIG$colsPastel[1]; col2 = CONFIG$colsDark[1]}
  if(row == 1 & column == 2){col1 =  CONFIG$colsPastel[3]; col2a = CONFIG$colsDark[3]; col2b = CONFIG$colsDark3[3]}
  if(row == 2 & column == 1){col1 =  CONFIG$colsPastel[2]; col2a = CONFIG$colsDark[2]; col2b = CONFIG$colsDark3[2]}
  if(row == 2 & column == 2){col1 =  CONFIG$colsDark[5]; col2a = CONFIG$colsDark2[5]; col2b = CONFIG$colsDark3[5]}
    
  # circuit breaker from 7 Apr to 1 Jun
  grid.polygon(c(98, 98, 153.5,153.5),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#ECECEC',col=NA))
  
  grid.polygon(c(153.5, 153.5, 170.5,170.5),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#F5F5F5',col=NA))
  
  grid.polygon(c(170.5, 170.5, 362,362),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
  
  # plot imported cases
  if(row == 1 & column == 1){
    
    plotImportCommExposureSpline = dataIncImport[,SPLINE]
    plotImportCommExposure = dataIncImport[,NOTIFIED_IMPORT_COMM_EXPOSURE]
    plotImport = dataIncImport[,NOTIFIED_IMPORT_COMM_EXPOSURE] + dataIncImport[,NOTIFIED_IMPORT_ISOLATED]
    
    for(time in 1:length(t)){
      
      if(plotImport[time]>0)
      {
        grid.polygon(t[time]+0.5*c(-1,-1,1,1),
                     c(0,plotImport[time],plotImport[time],0),
                     default.units = 'native',gp=gpar(fill=col1,col=NA)) 
      }
      
      if(plotImportCommExposure[time]>0)
      {
        grid.polygon(t[time]+0.5*c(-1,-1,1,1),
                     c(0,plotImportCommExposure[time],plotImportCommExposure[time],0),
                     default.units = 'native',gp=gpar(fill=col2,col=NA))
      }
      
    }
  }
  
  # plot linked cases
  if(row == 1 & column == 2){
    
    plotIncSim = data.table(TIME = outputTotal$time,
                            MEAN = apply(outputTotal$IncNotifiedLinkedSim,2, mean, na.rm = T),
                            LOWER_QUARTILE = apply(outputTotal$IncNotifiedLinkedSim,2, quantile, probs = 0.025, na.rm = T),
                            UPPER_QUARTILE = apply(outputTotal$IncNotifiedLinkedSim,2, quantile, probs = 0.975, na.rm = T))
    
    plotIncObs = outputTotal$IncNotifiedLinkedObs 
    
    for(time in 1:length(t)){

      if(plotIncObs[time]>0)
      {
        grid.polygon(t[time]+0.5*c(-1,-1,1,1),
                     c(0,plotIncObs[time],plotIncObs[time],0),
                     default.units = 'native',gp=gpar(fill=col1,col=NA)) 
      }
    }

    grid.lines(t,plotIncSim$MEAN,default.units = 'native',gp=gpar(col=col2a)) 
    grid.polygon(c(t,rev(t)),
                 c(plotIncSim$LOWER_QUARTILE,rev(plotIncSim$UPPER_QUARTILE)),
                 default.units = 'native',gp=gpar(col=NA,fill=col2b))
  }
  
  # plot unlinked cases
  if(row == 2 & column == 1){
    
    plotIncSim = data.table(TIME = outputTotal$time,
                            MEAN = apply(outputTotal$IncNotifiedUnlinkedSim,2, mean, na.rm = T),
                            LOWER_QUARTILE = apply(outputTotal$IncNotifiedUnlinkedSim,2, quantile, probs = 0.025, na.rm = T),
                            UPPER_QUARTILE = apply(outputTotal$IncNotifiedUnlinkedSim,2, quantile, probs = 0.975, na.rm = T))
    
    plotIncObs = outputTotal$IncNotifiedUnlinkedObs 
    
    for(time in 1:length(t)){
      
      if(plotIncObs[time]>0)
      {
        grid.polygon(t[time]+0.5*c(-1,-1,1,1),
                     c(0,plotIncObs[time],plotIncObs[time],0),
                     default.units = 'native',gp=gpar(fill=col1,col=NA))
      }
    }
    
    grid.lines(t,plotIncSim$MEAN,default.units = 'native',gp=gpar(col=col2a)) 
    grid.polygon(c(t,rev(t)),
                 c(plotIncSim$LOWER_QUARTILE,rev(plotIncSim$UPPER_QUARTILE)),
                 default.units = 'native',gp=gpar(col=NA,fill=col2b))
  }
  
  # plot missed cases
  if(row == 2 & column == 2){
    
    plotIncSim = data.table(TIME = outputTotal$time,
                            MEAN = apply(outputTotal$IncMissedSim_DOI,2, mean, na.rm = T),
                            LOWER_QUARTILE_95 = apply(outputTotal$IncMissedSim_DOI,2, quantile, probs = 0.025, na.rm = T),
                            UPPER_QUARTILE_95 = apply(outputTotal$IncMissedSim_DOI,2, quantile, probs = 0.975, na.rm = T),
                            LOWER_QUARTILE_50 = apply(outputTotal$IncMissedSim_DOI,2, quantile, probs = 0.25, na.rm = T),
                            UPPER_QUARTILE_50 = apply(outputTotal$IncMissedSim_DOI,2, quantile, probs = 0.75, na.rm = T))
    
    grid.polygon(c(t,rev(t)),
                 c(plotIncSim$LOWER_QUARTILE_50,rev(plotIncSim$UPPER_QUARTILE_50)),
                 default.units = 'native',gp=gpar(col=NA,fill=col2a))
    grid.polygon(c(t,rev(t)),
                 c(plotIncSim$LOWER_QUARTILE_95,rev(plotIncSim$UPPER_QUARTILE_95)),
                 default.units = 'native',gp=gpar(col=NA,fill=col2b))
    grid.lines(t,plotIncSim$MEAN,default.units = 'native',gp=gpar(col=col1))
     
  }
  

  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  if(row == 1) grid.xaxis(at=c(1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367),
                          label=rep('',13))
  if(row == 2) grid.xaxis(at=c(1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367),
                          label = c('2020 Jan 1', 'Feb 1', 'Mar 1', 'Apr 1', 'May 1', 'Jun 1', 'Jul 1', 
                                    'Aug 1', 'Sep 1', 'Oct 1', 'Nov 1', 'Dec 1', '2021 Jan 1'))

  if(row == 1 & column == 1) {grid.yaxis(at=seq(0,100,by=10),label=seq(0,100,by=10))
                              grid.text('A',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  if(row == 1 & column == 2) {grid.yaxis(at=seq(0,100,by=10),label=seq(0,100,by=10))
                              grid.text('B',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  if(row == 2 & column == 1) {grid.yaxis(at=seq(0,50,by=5),label=seq(0,50,by=5))
                              grid.text('C',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))} 
  if(row == 2 & column == 2) {grid.yaxis(at=seq(0,2500,by=250),label=seq(0,2500,by=250))
                              grid.text('D',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}

  if(row == 2) grid.text('Day',y=unit(-3,'lines'))
  if(column == 1 ) grid.text('Incidence',x=unit(-3.5,'lines'),rot=90)

  grid.text('Lockdown',x=unit(12.5,'lines'),y=unit(18.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 1',x=unit(19.5,'lines'),y=unit(18.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 2',x=unit(21.5,'lines'),y=unit(18.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 3',x=unit(44.9,'lines'),y=unit(18.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/incidence.png',height=20,width=40,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2)))

for(c in 1:2){
  for(r in 1:2){
    paneller(r,c)
  }
}

popViewport()
popViewport()
dev.off()

rm(paneller)