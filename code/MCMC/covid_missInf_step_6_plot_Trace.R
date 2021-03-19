library(grid)
library(readr)
library(data.table)

load('output processed/20210223/outputFourChain_period_1.RData')
load('output processed/20210223/outputFourChain_period_2.RData')
load('output processed/20210223/outputFourChain_period_3.RData')
load('output processed/20210223/outputFourChain_period_4.RData')
load('output processed/20210223/outputFourChain_period_5.RData')

outputFourChain_period_1 = data.table(outputFourChain_period_1)
outputFourChain_period_2 = data.table(outputFourChain_period_2)
outputFourChain_period_3 = data.table(outputFourChain_period_3)
outputFourChain_period_4 = data.table(outputFourChain_period_4)
outputFourChain_period_5 = data.table(outputFourChain_period_5)

# convert ratio of missed imported cases to average daily missed imported cases
# set time periods
startTimePeriod = c(18, 61, 98, 171, 245)  
endTimePeriod =  c(60, 97, 170, 244, 367)

# load observed data
dataIncImport = data.table(read_csv('input/dataIncImport.csv')) 

outputFourChain_period_1[, avgDailyMissedImport_period_1 := ratioImportMissed_period_1*mean(dataIncImport[TIME %in% c(startTimePeriod[1]:endTimePeriod[1]),SPLINE])]
outputFourChain_period_2[, avgDailyMissedImport_period_2 := ratioImportMissed_period_2*mean(dataIncImport[TIME %in% c(startTimePeriod[2]:endTimePeriod[2]),SPLINE])]
outputFourChain_period_3[, avgDailyMissedImport_period_3 := ratioImportMissed_period_3*mean(dataIncImport[TIME %in% c(startTimePeriod[3]:endTimePeriod[3]),SPLINE])]
outputFourChain_period_4[, avgDailyMissedImport_period_4 := ratioImportMissed_period_4*mean(dataIncImport[TIME %in% c(startTimePeriod[4]:endTimePeriod[4]),SPLINE])]
outputFourChain_period_5[, avgDailyMissedImport_period_5 := ratioImportMissed_period_5*mean(dataIncImport[TIME %in% c(startTimePeriod[5]:endTimePeriod[5]),SPLINE])]

# reorder columns
outputFourChain_period_1 = outputFourChain_period_1[,c(1:4,6,5)]
outputFourChain_period_2 = outputFourChain_period_2[,c(1:4,6,5)]
outputFourChain_period_3 = outputFourChain_period_3[,c(1:4,6,5)]
outputFourChain_period_4 = outputFourChain_period_4[,c(1:4,6,5)]
outputFourChain_period_5 = outputFourChain_period_5[,c(1:4,6,5)]


paneller=function(row = 1,column=1,inputChain)
{
  iterN = xmax = 30000# inputChain[,.N]/4
  xlm=c(0,xmax)
  if(row == 1) ylm=c(0,5)
  if(row == 2) ylm=c(0,1)
  if(row == 3) ylm=c(0,1)
  if(row == 4) ylm=c(0,20)
  if(row == 5) ylm=c(0,2.5)
  
  if(row == 3 & column == 3) ylm=c(0,0.2)
  if(row == 5 & column == 2) ylm=c(0,300)
  if(row == 5 & column == 4) ylm=c(0,1.5)
  
  innermargins = c(2,3,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  
  iteration = seq_len(iterN)
  col1 = rgb(171/223,1,117/223,alpha=0.7) 
  col2 = '#00468BFF'
  col3 = rgb(1,165/255,0,alpha=0.3)
  col4 = rgb(237/255,0,0,alpha=0.5)
  
  chain1 = inputChain[chain == 1][[row]][1:30000]
  chain2 = inputChain[chain == 2][[row]][1:30000]
  chain3 = inputChain[chain == 3][[row]][1:30000]
  chain4 = inputChain[chain == 4][[row]][1:30000]
  
  grid.lines(iteration, chain1,default.units = 'native',gp=gpar(col=col1))
  grid.lines(iteration, chain2,default.units = 'native',gp=gpar(col=col2))
  grid.lines(iteration, chain3,default.units = 'native',gp=gpar(col=col3))
  grid.lines(iteration, chain4,default.units = 'native',gp=gpar(col=col4))
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.xaxis(at=c(0,15000,30000),label = c(0,15000,30000))
  # if(column %in% c(1,5)) grid.xaxis(at=c(0,15000,30000),label = c(0,15000,30000))
  # if(column %in% c(2:4)) grid.xaxis(at=c(0,150000,300000),label = c(0,150000,300000))
  
  if(row == 1) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  if(row == 2) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,1,by=0.25))
  if(row == 3 & column != 3) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,1,by=0.25))
  if(row == 3 & column == 3) grid.yaxis(at=seq(0,0.2,by=0.05),label=seq(0,0.2,by=0.05))
  if(row == 4) grid.yaxis(at=seq(0,20,by=5),label=seq(0,20,by=5))
  if(row == 5 & column %in% c(1,3,5)) grid.yaxis(at=seq(0,2.5,by=0.5),label=seq(0,2.5,by=0.5))
  if(row == 5 & column == 2) grid.yaxis(at=seq(0,300,by=100),label=seq(0,300,by=100))
  if(row == 5 & column == 4) grid.yaxis(at=seq(0,1.5,by=0.25),label=seq(0,1.5,by=0.25))
  
  if(row == 1 & column == 1) grid.text('Jan 18 - Feb 29', y=unit(8,'lines'))
  if(row == 1 & column == 2) grid.text('Mar 1 - Apr 6', y=unit(8,'lines'))
  if(row == 1 & column == 3) grid.text('Apr 7 - Jun 18', y=unit(8,'lines'))
  if(row == 1 & column == 4) grid.text('Jun 19 - Aug 31', y=unit(8,'lines'))
  if(row == 1 & column == 5) grid.text('Sep 1 - Jan 1', y=unit(8,'lines'))
  
  if(row == 5) grid.text('iteration',y=unit(-3,'lines'))
  
  if(column == 1 & row == 1) grid.text('R',x=unit(-3.5,'lines'),rot=90)
  if(column == 1 & row == 2) grid.text(bquote(~epsilon[link]),x=unit(-3.5,'lines'),rot=90) #~epsilon[op]
  if(column == 1 & row == 3) grid.text(bquote(~epsilon[unlink]),x=unit(-3.5,'lines'),rot=90) #~epsilon[op*'\'']
  if(column == 1 & row == 4) grid.text(bquote(~rho),x=unit(-3.5,'lines'),rot=90)
  if(column == 1 & row == 5) grid.text('Avg daily missed\nimportation',x=unit(-3.5,'lines'),rot=90)
  
  if(column == 3 & row == 4) grid.text('Lockdown of borders', y=unit(4.5,'lines'), gp = gpar(fontsize = 8))
  if(column == 3 & row == 5) grid.text('Lockdown of borders', y=unit(4.5,'lines'), gp = gpar(fontsize = 8))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/trace.png',height=25,width=30,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=5,ncol=5)))

for(c in 1:5){
  
  if(c == 1) inputChain = outputFourChain_period_1
  if(c == 2) inputChain = outputFourChain_period_2
  if(c == 3) inputChain = outputFourChain_period_3
  if(c == 4) inputChain = outputFourChain_period_4
  if(c == 5) inputChain = outputFourChain_period_5
  
  for(r in 1:5){
    paneller(r,c,inputChain)
  }
}

popViewport()
popViewport()
dev.off()

rm(paneller)