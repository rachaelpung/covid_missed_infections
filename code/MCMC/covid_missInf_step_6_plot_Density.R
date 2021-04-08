library(grid)
library(data.table)

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
outputTotalChain = outputTotalChain[,c(1:4,21,5:8,22,9:12,23,13:16,24,17:20,25)]


paneller=function(row = 1,column=1, outputTotal)
{
  if(column == 1){xlm=c(0,2); ylm=c(0,20)} 
  if(column == 2){xlm=c(0,1); ylm=c(0,20)}
  if(column == 2 & row == 4){xlm=c(0,1); ylm=c(0,50)}
  
  if(column == 3){xlm=c(0,1); ylm=c(0,10)}
  if(column == 3 & row == 3){xlm=c(0,0.2); ylm=c(0,300)}
  if(column == 3 & row == 4){xlm=c(0,0.2); ylm=c(0,50)}
  
  if(column == 4){xlm=c(0,20); ylm=c(0,1.5)}
  
  if(column == 5){xlm=c(0,2); ylm=c(0,2)}
  if(column == 5 & row == 2){xlm=c(0,100); ylm=c(0,0.1)}
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  col1 = '#00468BFF'

  
  param = 5*(row-1)+column
  
  chain = density(outputTotal[[param]])
  
  if(row == 3 & column == 4 ) chain = data.frame(x = c(0,0), y = c(0,0))
  if(row == 3 & column == 5 ) chain = data.frame(x = c(0,0), y = c(0,0))
  grid.lines(chain$x, chain$y, default.units = 'native',gp=gpar(col=col1))
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  # y axis
  if(column == 1) grid.yaxis(at=seq(0,20,by=5),label=seq(0,20,by=5))
  if(column == 2 & row != 4) grid.yaxis(at=seq(0,20,by=5),label=seq(0,20,by=5))
  if(column == 2 & row == 4) grid.yaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))
  
  if(column == 3 & row %in% c(1,2,5)) grid.yaxis(at=seq(0,10,by=5),label=seq(0,10,by=5))
  if(column == 3 & row == 3) grid.yaxis(at=seq(0,300,by=100),label=seq(0,300,by=100))
  if(column == 3 & row == 4) grid.yaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))
  
  if(column == 4) grid.yaxis(at=seq(0,1.5,by=0.5),label=seq(0,1.5,by=0.5))
  
  if(column == 5 & row !=2) grid.yaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5))
  if(column == 5 & row ==2) grid.yaxis(at=seq(0,0.1,by=0.025),label=seq(0,0.1,by=0.025))
  
  if(column == 1) grid.xaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5))
  if(column == 2) grid.xaxis(at=seq(0,1,by=0.25),label=seq(0,1,by=0.25))
  if(column == 3 & row %in% c(1,2,5)) grid.xaxis(at=seq(0,1,by=0.25),label=seq(0,1,by=0.25))
  if(column == 3 & row %in% c(3,4)) grid.xaxis(at=seq(0,0.2,by=0.05),label=seq(0,0.2,by=0.05))
  if(column == 4) grid.xaxis(at=seq(0,20,by=5),label=seq(0,20,by=5))
  if(column == 5 & row !=2) grid.xaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5))
  if(column == 5 & row ==2) grid.xaxis(at=seq(0,100,by=25),label=seq(0,100,by=25))
  
  if(row == 5 & column == 1) grid.text('R',y=unit(-3,'lines'))
  if(row == 5 & column == 2) grid.text(bquote(~epsilon[link]),y=unit(-3,'lines'))  #~epsilon[op]
  if(row == 5 & column == 3) grid.text(bquote(~epsilon[unlink]),y=unit(-3,'lines')) #~epsilon[op*'\'']
  if(row == 5 & column == 4) grid.text(bquote(~rho),y=unit(-3,'lines'))
  if(row == 5 & column == 5) grid.text('Avg daily missed\nimportation',y=unit(-3,'lines'))

  if(row == 1 & column == 1) grid.text('Jan 18 - Feb 29', x=unit(-3,'lines'),rot=90)
  if(row == 2 & column == 1) grid.text('Mar 1 - Apr 6', x=unit(-3,'lines'),rot=90)
  if(row == 3 & column == 1) grid.text('Apr 7 - Jun 18', x=unit(-3,'lines'),rot=90)
  if(row == 4 & column == 1) grid.text('Jun 19 - Aug 31', x=unit(-3,'lines'),rot=90)
  if(row == 5 & column == 1) grid.text('Sep 1 - Jan 1', x=unit(-3,'lines'),rot=90)
  
  if(row == 3 & column == 4) grid.text('Lockdown of borders', x=unit(5,'lines'))
  if(row == 3 & column == 5) grid.text('Lockdown of borders', x=unit(5,'lines'))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/density.png',height=25,width=30,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=5,ncol=5)))

for(c in 1:5){
  for(r in 1:5){
    paneller(r,c,outputTotalChain)
  }
}

popViewport()
popViewport()
dev.off()

rm(paneller)



