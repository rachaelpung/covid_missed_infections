library(grid)
library(readr)
library(data.table)
source('codes/v5/covid_missInf_function_plot.r')

# load distribution of duration from infection to isolation
distIsolate = data.table(read_csv('input/dist_Isolate.csv')) 
distIsolateFit = data.table(read_csv('input/dist_IsolateFit.csv')) 

# load distribution of duration from infection to travel
distTravel = data.table(read_csv('input/dist_Travel.csv')) 

# load distribution of duration from arrival to isolation
distArrivalIsolate = data.table(read_csv('input/dist_Arrival_Isolate.csv')) 
distArrivalIsolateFit = data.table(read_csv('input/dist_Arrival_IsolateFit.csv')) 

# plot R
paneller=function(row = 1,column=1, outputTotal)
{
  xlm=c(0,21)
  ylm=c(0,15)
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  colObs=CONFIG$colsDark3[1]; colFit=CONFIG$colsDark[1]
  pchWeibull=16;pchGamma=4;pchLN=2
  
  if(row==1 & column==1){dataIsolate=distIsolate[PERIOD==1,];dataIsolateFit=distIsolateFit[PERIOD==1,]}
  if(row==1 & column==2){dataIsolate=distIsolate[PERIOD==2,];dataIsolateFit=distIsolateFit[PERIOD==2,]}
  if(row==1 & column==3){dataIsolate=distIsolate[PERIOD==3,];dataIsolateFit=distIsolateFit[PERIOD==3,]}
  if(row==2 & column==1){dataIsolate=distIsolate[PERIOD==4,];dataIsolateFit=distIsolateFit[PERIOD==4,]}
  if(row==2 & column==2){dataIsolate=distIsolate[PERIOD==5,];dataIsolateFit=distIsolateFit[PERIOD==5,]}
  
 
  for(d in dataIsolate[,DAY]){
     grid.polygon(d+0.5*c(-1,-1,1,1),
                  c(0,dataIsolate[DAY==d,PROB_ISOLATE]*100,dataIsolate[DAY==d,PROB_ISOLATE]*100,0),
                  default.units = 'native',gp=gpar(fill=colObs,col='white'))
  }
 
  
    
  grid.lines(dataIsolateFit[,DAY], dataIsolateFit[,PROB_ISOLATE_WEIBULL]*100, default.units = 'native',gp=gpar(col=colFit,lwd=2))
  grid.lines(dataIsolateFit[,DAY], dataIsolateFit[,PROB_ISOLATE_GAMMA]*100, default.units = 'native',gp=gpar(col=colFit,lwd=2))
  grid.lines(dataIsolateFit[,DAY], dataIsolateFit[,PROB_ISOLATE_LOG_NORMAL]*100, default.units = 'native',gp=gpar(col=colFit,lwd=2))
    
  grid.points(dataIsolateFit[,DAY][2:21],(dataIsolateFit[,PROB_ISOLATE_WEIBULL]*100)[2:21],default.units = 'native',pch=pchWeibull,gp=gpar(cex=0.5,col=colFit)) 
  grid.points(dataIsolateFit[,DAY][2:21],(dataIsolateFit[,PROB_ISOLATE_GAMMA]*100)[2:21],default.units = 'native',pch=pchGamma,gp=gpar(cex=0.5,col=colFit))  
  grid.points(dataIsolateFit[,DAY][2:21],(dataIsolateFit[,PROB_ISOLATE_LOG_NORMAL]*100)[2:21],default.units = 'native',pch=pchLN,gp=gpar(cex=0.5,col=colFit))
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.yaxis(at=seq(0,15,by=5),label=seq(0,15,by=5))
  grid.xaxis(at=seq(0,20,by=5),label=seq(0,20,by=5))
  
  if(column==1)grid.text('Probability (%)',x=unit(-3.3,'lines'),rot=90)
  
  if(row==1 & column == 1) grid.text('A',x=unit(-3.3,'lines'),y=unit(13.8,'lines'),gp=gpar(fontsize=unit(10,'pt')))
  if(row==1 & column == 2) grid.text('B',x=unit(-3.3,'lines'),y=unit(13.8,'lines'),gp=gpar(fontsize=unit(10,'pt')))
  if(row==1 & column == 3) grid.text('C',x=unit(-3.3,'lines'),y=unit(13.8,'lines'),gp=gpar(fontsize=unit(10,'pt')))
  if(row==2 & column == 1) grid.text('D',x=unit(-3.3,'lines'),y=unit(13.8,'lines'),gp=gpar(fontsize=unit(10,'pt')))
  if(row==2 & column == 2) grid.text('E',x=unit(-3.3,'lines'),y=unit(13.8,'lines'),gp=gpar(fontsize=unit(10,'pt')))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/distIsolate.png',height=16,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)

paneller(2,1)
paneller(2,2)

popViewport()

grid.text('Days since infection',y=unit(-1,'lines'))

popViewport()
dev.off()

rm(paneller)


# plot distribution of time from infection to arrival
png('figure/distTravel.png',height=8,width=8,units='cm',res=300,pointsize=10)

xlm=c(1,14)
ylm=c(0,15)

pushViewport(plotViewport(c(4,4,4,1),
                          xscale=xlm,yscale=ylm))


colFit=CONFIG$colsDark[1]
grid.lines(distTravel[,DAY]*-1, distTravel[,PROB_INFECTION]*100, default.units = 'native',gp=gpar(col=colFit,lwd=2))
grid.points((distTravel[,DAY]*-1)[2:13],(distTravel[,PROB_INFECTION]*100)[2:13],default.units = 'native',pch=16,gp=gpar(cex=0.5,col=colFit)) 


grid.yaxis(at=seq(0,15,by=5),label=seq(0,15,by=5))
grid.xaxis(at=seq(1,14,by=2),label=seq(1,14,by=2))
grid.text('Probability (%)',x=unit(-3.3,'lines'),rot=90)
grid.text('Days since infection',y=unit(-3.3,'lines'))

grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))


popViewport()
dev.off()


# plot distribution of time from arrival to isolation
png('figure/distArrivalIsolate.png',height=8,width=8,units='cm',res=300,pointsize=10)

xlm=c(-0.5,21.5)
ylm=c(0,20)

pushViewport(plotViewport(c(4,4,4,1),
                          xscale=xlm,yscale=ylm))

colObs=CONFIG$colsDark3[1]; colFit=CONFIG$colsDark[1]
pchNB=16;pchPoisson=4

for(d in distArrivalIsolate[,DAY_SINCE_ARRIVAL]){
  grid.polygon(d+0.5*c(-1,-1,1,1),
               c(0,distArrivalIsolate[DAY_SINCE_ARRIVAL==d,PROB_ISOLATE_DAY_SINCE_ARRIVAL]*100,distArrivalIsolate[DAY_SINCE_ARRIVAL==d,PROB_ISOLATE_DAY_SINCE_ARRIVAL]*100,0),
               default.units = 'native',gp=gpar(fill=colObs,col='white'))
}

grid.lines(distArrivalIsolateFit[,DAY], distArrivalIsolateFit[,PROB_ISOLATE_NB]*100,default.units = 'native',gp=gpar(col=colFit,lwd=2))
grid.points(distArrivalIsolateFit[,DAY],distArrivalIsolateFit[,PROB_ISOLATE_NB]*100,default.units = 'native',pch=pchNB,gp=gpar(cex=0.5,col=colFit)) 

grid.lines(distArrivalIsolateFit[,DAY], distArrivalIsolateFit[,PROB_ISOLATE_POISSON]*100,default.units = 'native',gp=gpar(col=colFit,lwd=2))
grid.points(distArrivalIsolateFit[,DAY],distArrivalIsolateFit[,PROB_ISOLATE_POISSON]*100,default.units = 'native',pch=pchPoisson,gp=gpar(cex=0.5,col=colFit)) 


grid.yaxis(at=seq(0,20,by=5),label=seq(0,20,by=5))
grid.xaxis(at=seq(0,21,by=5),label=seq(0,21,by=5))
grid.text('Probability (%)',x=unit(-3.3,'lines'),rot=90)
grid.text('Days since arrival',y=unit(-3.3,'lines'))

grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))


popViewport()
dev.off()
