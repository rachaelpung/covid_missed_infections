library(grid)
library(readr)
library(data.table)

source('codes/v5/covid_missInf_function_plot.r')
source('codes/v5/Eigen/covid_missInf_eigen.r')

# load infections
load('output processed/20210223/outputTotalChain_10000_thinned_iter_period_all.rdata')
outputTotal = copy(outputTotalChain_period_all)


paneller=function(row = 1,column=1)
{
  if(row == 1){ xlm=c(-0.025, 1.025); ylm=c(-0.025, 1.025)}
  if(row == 2){ xlm=c(0, 1); ylm=c(0.5,1.5)} 

  innermargins = c(2,2,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(column == 1) data = dataEigen[R == 1.5,]
  if(column == 2) data = dataEigen[R == 1.2,]
  if(column == 3) data = dataEigen[R == 1.0,]
  
  if(row == 1){
    
    for(i in 1:data[,.N]){

      grid.polygon(data[i,effOffspringParentNotified]+0.025*c(-1,1,1,-1),
                   data[i,effOffspringParentMissed]+0.025*c(-1,-1,1,1),
                   default.units = 'native',gp=gpar(fill=data[i,eigenvalue_hex],col=data[i,eigenvalue_hex]))

      
    }
    
    if(column==2){
      grid.lines(dataEff_R_1[R == 1.2 & effOffspringParentMissed<=1 & effOffspringParentMissed>=0,effOffspringParentNotified],
                 dataEff_R_1[R == 1.2 & effOffspringParentMissed<=1 & effOffspringParentMissed>=0,effOffspringParentMissed],
                 default.units = 'native',gp=gpar(col='#B5B5B5',lty = 'dashed'))
    }
    
    # if(column==1){
    #   grid.lines(dataEff_R_1[R == 1.5 & effOffspringParentMissed<=1 & effOffspringParentMissed>=0,effOffspringParentNotified],
    #              dataEff_R_1[R == 1.5 & effOffspringParentMissed<=1 & effOffspringParentMissed>=0,effOffspringParentMissed],
    #              default.units = 'native',gp=gpar(col='#B5B5B5',lty = 'dashed'))
    # }
    
    
  }
  
  if(row == 2){
    
    eff = c(0.2, 0.5, 0.8)
    col = c( CONFIG$colsDark[1], CONFIG$colsDark[4], CONFIG$colsDark[5])
    
    for(line in 1:3){
      
      x = data[effOffspringParentNotified == eff[line], effOffspringParentMissed]
      y = data[effOffspringParentNotified == eff[line], eigenvalue]
      grid.lines(x,y,default.units = 'native',gp=gpar(col=col[line])) 
      grid.lines(c(0,1),c(1,1),default.units = 'native',gp=gpar(col='#A50026', lty = 'dashed')) 
    }
    
    
  }
  
  popViewport()
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))

  if(row == 1) {
    grid.xaxis(at=seq(0,1,0.25), label=seq(0,1,0.25))
    grid.yaxis(at=seq(0,1,0.25), label=seq(0,1,0.25))
    grid.text(bquote(~epsilon[link]),y=unit(-2.5,'lines'))
  }
  if(row == 2) {
    grid.xaxis(at=seq(0,1,0.25), label=seq(0,1,0.25))
    grid.yaxis(at=seq(0.5,1.5,0.25), label=seq(0.5,1.5,0.25))
    grid.text(bquote(~epsilon[unlink]),y=unit(-2.5,'lines'))
  }
  
  if(row == 1 & column == 1) grid.text(bquote(~epsilon[unlink]),x=unit(-3.5,'lines'),rot=90)
  if(row == 2 & column == 1) grid.text('Effective reproduction number',x=unit(-3.5,'lines'),rot=90,gp=gpar(fontsize=unit(9,'pt')))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  
  
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row == 1 & column == 3){
    
    colourbar(x_bottom_left = 1.2, y_bottom_left = 0.125, y_length = 0.75, x_width = 0.075)
  }
  
  if(row == 2 & column == 3){
    
    col = c( CONFIG$colsDark[1], CONFIG$colsDark[4], CONFIG$colsDark[5])
    grid.text(bquote(~epsilon[link]),x=1.15, y=1.25,default.units = 'native')
    grid.lines(c(1.1,1.15),c(1.125,1.125),default.units = 'native',gp=gpar(col=col[1], lwd = 1.5))
    grid.text('0.2',x=1.25,y=1.125,default.units = 'native',gp=gpar(fontsize=unit(8,'pt')))
    grid.lines(c(1.1,1.15),c(1.0,1.0),default.units = 'native',gp=gpar(col=col[2], lwd = 1.5))
    grid.text('0.5',x=1.25,y=1.0,default.units = 'native',gp=gpar(fontsize=unit(8,'pt')))
    grid.lines(c(1.1,1.15),c(0.875,0.875),default.units = 'native',gp=gpar(col=col[3], lwd = 1.5))
    grid.text('0.8',x=1.25,y=0.875,default.units = 'native',gp=gpar(fontsize=unit(8,'pt')))
    
  }
  
  # caption
  if(row == 1 & column == 1) grid.text('A',x=unit(-3,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(10,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-3,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(10,'pt')))   
  if(row == 1 & column == 3) grid.text('C',x=unit(-3,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(10,'pt')))   
  if(row == 2 & column == 1) grid.text('D',x=unit(-3,'lines'),y=unit(11.7,'lines'),gp=gpar(fontsize=unit(10,'pt'))) 
  if(row == 2 & column == 2) grid.text('E',x=unit(-3,'lines'),y=unit(11.7,'lines'),gp=gpar(fontsize=unit(10,'pt')))  
  if(row == 2 & column == 3) grid.text('F',x=unit(-3,'lines'),y=unit(11.7,'lines'),gp=gpar(fontsize=unit(10,'pt')))   
  
  
  popViewport()
  popViewport()
  
}


colourbar <- function(x_bottom_left = 1.1, y_bottom_left = 0.4, y_length = 0.75, x_width = 0.2){
  
  colour_hex = c("#006837", "#1A9850", "#66BD63", "#A6D96A", "#D9EF8B","#FFFFBF",
                 "#FEE08B", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  cols <- colorRampPalette(colour_hex)
  
  y_start = seq(y_bottom_left, y_bottom_left+y_length,length.out = 126)[-126]
  y_end = seq(y_bottom_left, y_bottom_left+y_length,length.out = 126)[-1]
  
  for(i in 1:125){
    grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_width,x_bottom_left+x_width),
                 c(y_start[i],y_end[i],y_end[i],y_start[i]),
                 default.units = 'native',gp=gpar(fill=cols(125)[i], col=cols(125)[i]))
  }
  
  y_start <-  seq(y_bottom_left, y_bottom_left+y_length,length.out = 5)
  y_start[1] = 0.15
  y_start[5] = 0.85
  label <- c('0.50', '0.75', '1.00', '1.25', '1.50')
  
  for(breaks in 1:5){
    
    if(breaks %in% c(2:4)){
      grid.lines(c(x_bottom_left, x_bottom_left+(x_width/5)),
                 c(y_start[breaks], y_start[breaks]),default.units = 'native',gp=gpar(col='black')) 
       grid.lines(c(x_bottom_left+x_width, x_bottom_left+x_width-(x_width/5)),
                  c(y_start[breaks], y_start[breaks]),default.units = 'native',gp=gpar(col='black')) 
    }
    
    
    grid.text(label[breaks],x=x_bottom_left+x_width+0.05,y=y_start[breaks],gp = gpar(fontsize = 8))
    
  }
  
  grid.text('Effective reproduction number',x=x_bottom_left-0.1, y=y_bottom_left+(y_length/2), 
            gp = gpar(fontsize = 8), rot=90)
  
  grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_width,x_bottom_left+x_width),
               c(y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length,y_bottom_left),
               default.units = 'native',gp=gpar(fill=NA, col='black'))
  
  
}

png('figure/eigenvalue_period2.png',height=15,width=22.5,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,4)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)

paneller(2,1)
paneller(2,2)
paneller(2,3)

popViewport()
popViewport()
dev.off()

rm(paneller)
rm(colourbar)


# plot eigenvector ratios
paneller=function(row=1, column=1)
{
  if(row == 1 & column == 1) {xlm=c(0.1,0.9); ylm=c(0,5)}
  if(row == 1 & column == 2) {xlm=c(0.1,0.9); ylm=c(0,10)}
  if(row == 2 & column == 1) {xlm=c(0.5,5.5); ylm=log(c(0.01,100))}
  if(row == 2 & column == 2) {xlm=c(0.1,0.9); ylm=c(0,10)}
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(c(3,4,1,1),xscale=xlm,yscale=ylm))
  
  if(row == 1 & column == 1){
    for(line in 1:3){
      
      eff = c(0.2, 0.5, 0.8)
      col = c(CONFIG$colsDark[1], CONFIG$colsDark[4], CONFIG$colsDark[5])
      data = dataEigen[R == 1.5 & effOffspringParentMissed >=0.1 & effOffspringParentMissed <=0.9,]
      
      x = data[effOffspringParentNotified == eff[line] & eigenvectorratio <=5 & effOffspringParentMissed >= 0.2, effOffspringParentMissed]
      y = data[effOffspringParentNotified == eff[line] & eigenvectorratio <=5 & effOffspringParentMissed >= 0.2, eigenvectorratio]
      grid.lines(x,y,default.units = 'native',gp=gpar(col=col[line])) 
      
      x = dataInterpolate[effOffspringParentNotified == eff[line] & eigenvectorratio <=5 & effOffspringParentMissed >= 0.1, effOffspringParentMissed]
      y = dataInterpolate[effOffspringParentNotified == eff[line] & eigenvectorratio <=5 & effOffspringParentMissed >= 0.1, eigenvectorratio]
      grid.lines(x,y,default.units = 'native',gp=gpar(col=col[line],lty = 'dashed')) 
      
    }
    
    grid.lines(c(0.1,0.9),c(1,1),default.units = 'native',gp=gpar(col='#A50026', lty = 'dashed')) 
  }
  
  if(row == 1 & column == 2){
    for(line in 1:3){
      
      eff = c(0.2, 0.5, 0.8)
      col = c(CONFIG$colsDark[1], CONFIG$colsDark[4], CONFIG$colsDark[5])
      data = dataEigen[R == 1.5 & effOffspringParentMissed >=0.1 & effOffspringParentMissed <=0.9,]
      
      x = data[effOffspringParentNotified == eff[line] & unlinklinkratio <=10, effOffspringParentMissed]
      y = data[effOffspringParentNotified == eff[line] & unlinklinkratio <=10, unlinklinkratio]
      grid.lines(x,y,default.units = 'native',gp=gpar(col=col[line])) 
      
      x = dataInterpolate[effOffspringParentNotified == eff[line] & unlinklinkratio <=10 & effOffspringParentMissed >=0.1, effOffspringParentMissed]
      y = dataInterpolate[effOffspringParentNotified == eff[line] & unlinklinkratio <=10 & effOffspringParentMissed >=0.1, unlinklinkratio]
      grid.lines(x,y,default.units = 'native',gp=gpar(col=col[line])) 
      
    }
    
  }
  
  if(row == 2 & column == 1){
   
    startTimePeriod = c(18, 61, 98, 171, 245)  
    endTimePeriod =  c(60, 97, 170, 244, 367)
    
    grid.lines(xlm, log(c(0.01,0.01)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(0.1,0.1)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(1,1)),       default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(10,10)),     default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(100,100)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    
    for(period in 1:5){
      
      y_missed = apply(outputTotal$IncMissedSim_DOI[,which(outputTotal$time == startTimePeriod[period]):which(outputTotal$time == endTimePeriod[period])],1,sum)
      y_notified = apply((outputTotal$IncNotifiedUnlinkedSim_DOI[,which(outputTotal$time == startTimePeriod[period]):which(outputTotal$time == endTimePeriod[period])]+
                            outputTotal$IncNotifiedLinkedSim_DOI[,which(outputTotal$time == startTimePeriod[period]):which(outputTotal$time == endTimePeriod[period])])  ,1,sum)
      y_ratio = y_missed/y_notified
                                
      median  = log(median(y_ratio))
      low_50  = log(quantile(y_ratio, probs = 0.25))
      upp_50  = log(quantile(y_ratio, probs = 0.75))
      low_95  = log(quantile(y_ratio, probs = 0.025))
      upp_95  = log(quantile(y_ratio, probs = 0.975))
      
      colMedian=CONFIG$colsDark[4]; col50=CONFIG$colsDark1[4]; col95=CONFIG$colsDark3[4]
      
      grid.lines(c(period-0.2,period-0.2), c(low_95,upp_95), default.units = 'native',gp=gpar(col=col95,lwd=2))
      grid.lines(c(period-0.2,period-0.2), c(low_50,upp_50), default.units = 'native',gp=gpar(col=col50,lwd=2))
      grid.points(period-0.2,median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
      
      y_unlink = apply(outputTotal$IncNotifiedUnlinkedSim_DOI[,which(outputTotal$time == startTimePeriod[period]):which(outputTotal$time == endTimePeriod[period])],1,sum)
      y_link   = apply(outputTotal$IncNotifiedLinkedSim_DOI[,which(outputTotal$time == startTimePeriod[period]):which(outputTotal$time == endTimePeriod[period])],1,sum)
      y_ratio = y_unlink/y_link
      
      median  = log(median(y_ratio))
      low_50  = log(quantile(y_ratio, probs = 0.25))
      upp_50  = log(quantile(y_ratio, probs = 0.75))
      low_95  = log(quantile(y_ratio, probs = 0.025))
      upp_95  = log(quantile(y_ratio, probs = 0.975))
      
      colMedian=CONFIG$colsDark[5]; col50=CONFIG$colsDark1[5]; col95=CONFIG$colsDark3[5]
      
      grid.lines(c(period+0.2,period+0.2), c(low_95,upp_95), default.units = 'native',gp=gpar(col=col95,lwd=2))
      grid.lines(c(period+0.2,period+0.2), c(low_50,upp_50), default.units = 'native',gp=gpar(col=col50,lwd=2))
      grid.points(period+0.2,median,default.units = 'native',pch=17,gp=gpar(cex=0.8,col=colMedian)) 
      
      
      
    }
  
  }
  
  if(row == 2 & column == 2){
    
    data = dataImportScenario[effOffspringParentMissed >= 0.1 & effOffspringParentMissed <=0.9 & ratio_missed == 'low',]
    eff = seq(0.1,0.9,0.05)
    colLine = CONFIG$colsDark[1]
    
    grid.lines(data[,effOffspringParentMissed],data[,notified_time_exp],default.units='native',gp=gpar(col=colLine))
    grid.lines(data[,effOffspringParentMissed],data[,missed_time_exp],default.units='native',gp=gpar(col=colLine,lty='dashed')) 
    
    
    
  }

  
  if(row == 1) grid.xaxis(at=seq(0.1,0.9,0.2),label=seq(0.1,0.9,0.2))
  if(row == 2 & column == 1) grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nAug 31', 'Sep 1-\nJan 1'),
                                        gp = gpar(fontsize=8))
  if(row == 2 & column == 2) grid.xaxis(at=seq(0.1,0.9,0.2), label=seq(0.1,0.9,0.2))
  
  if(row == 1 & column == 1) grid.yaxis(at=seq(0,5,1),label=seq(0,5,1))
  if(row == 1 & column == 2) grid.yaxis(at=seq(0,10,2),label=seq(0,10,2))
  if(row == 2 & column == 1) grid.yaxis(at=log(c(0.01,0.1,1,10,100)),label=c(0.01,0.1,1,10,100))
  if(row == 2 & column == 2) grid.yaxis(at=seq(0,10,by=2),label=seq(0,10,by=2))
  
  if(row == 1) grid.text(bquote(~epsilon[unlink]),y=unit(-3,'lines'))
  if(row == 2 & column == 1) grid.text('Time period',y=unit(-3,'lines'))
  if(row == 2 & column == 2) grid.text(bquote(~epsilon[unlink]),y=unit(-3,'lines'))
  
  if(row == 1 & column == 1) grid.text('Ratio of missed to notified cases',x=unit(-3.5,'lines'),rot=90)
  if(row == 1 & column == 2) grid.text('Ratio of unlinked to linked cases',x=unit(-3.5,'lines'),rot=90)
  if(row == 2 & column == 1) grid.text('Ratio',x=unit(-3.5,'lines'),rot=90)
  if(row == 2 & column == 2) grid.text('Generations to\nexponential growth',x=unit(-3,'lines'),rot=90)
  
  if(row == 1 & column == 1){
    
    col = c( CONFIG$colsDark[1], CONFIG$colsDark[4], CONFIG$colsDark[5])
    grid.text(bquote(~epsilon[link]),x=0.75, y=4.5,default.units = 'native')
    grid.lines(c(0.7,0.75),c(4.15,4.15),default.units = 'native',gp=gpar(col=col[1], lwd = 1.5))
    grid.text('0.2',x=0.83,y=4.15,default.units = 'native',gp=gpar(fontsize=unit(8,'pt')))
    grid.lines(c(0.7,0.75),c(3.8,3.8),default.units = 'native',gp=gpar(col=col[2], lwd = 1.5))
    grid.text('0.5',x=0.83,y=3.8,default.units = 'native',gp=gpar(fontsize=unit(8,'pt')))
    grid.lines(c(0.7,0.75),c(3.45,3.45),default.units = 'native',gp=gpar(col=col[3], lwd = 1.5))
    grid.text('0.8',x=0.83,y=3.45,default.units = 'native',gp=gpar(fontsize=unit(8,'pt')))
    
  }
  
  # caption
  if(row == 1 & column == 1) grid.text('A',x=unit(-3,'lines'),y=unit(11.9,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-3,'lines'),y=unit(11.9,'lines'),gp=gpar(fontsize=unit(12,'pt')))  
  if(row == 2 & column == 1) grid.text('C',x=unit(-3,'lines'),y=unit(11.9,'lines'),gp=gpar(fontsize=unit(12,'pt')))  
  if(row == 2 & column == 2) grid.text('D',x=unit(-3,'lines'),y=unit(11.9,'lines'),gp=gpar(fontsize=unit(12,'pt'))) 
 
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/eigenratio.png',height=16,width=16,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(1,0,0,0)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2)))


paneller(1,1)
paneller(1,2)
paneller(2,1)
paneller(2,2)


popViewport()
popViewport()
dev.off()

rm(paneller)

# plot number of cases by generations
paneller=function(row = 1,column=1)
{
  if(row == 1 & column == 1) {xlm=c(0, 10); ylm=c(0,log(1700))}
  if(row == 1 & column == 2) {xlm=c(0, 10); ylm=c(0,log(1700))}
  if(row == 2){xlm=c(0.1,0.9); ylm=c(0,10)}
  
  innermargins = c(3,4,1,1)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row == 1){
    
    data = dataGen
    gen = 0:10
    col = c(CONFIG$colsDark2[3], CONFIG$colsDark1[3], CONFIG$colsDark[3],
            CONFIG$colsDark2[2], CONFIG$colsDark1[2], CONFIG$colsDark[2],
            CONFIG$colsDark2[1], CONFIG$colsDark1[1], CONFIG$colsDark[1])
    
    if(column == 1){
      col = c(CONFIG$colsDark[1], CONFIG$colsDark2[1])
      
      s=3
      grid.lines(gen,log(data[scenario==s,notified]),default.units='native',gp=gpar(col=col[1]))
      grid.lines(gen,log(data[scenario==s,missed]),default.units='native',gp=gpar(col=col[1],lty='dashed'))
      grid.points(gen[2:10],log(data[scenario==s,notified])[2:10],default.units = 'native',pch=16,gp=gpar(cex=0.5,col=col[1])) 
      grid.points(gen[2:10],log(data[scenario==s,missed])[2:10],default.units = 'native',pch=16,gp=gpar(cex=0.5,col=col[1])) 
      
      s=19
      grid.lines(gen,log(data[scenario==s,notified]),default.units='native',gp=gpar(col=col[2]))
      grid.lines(gen,log(data[scenario==s,missed]),default.units='native',gp=gpar(col=col[2],lty='dashed'))
      grid.points(gen[2:10],log(data[scenario==s,notified])[2:10],default.units = 'native',pch=4,gp=gpar(cex=0.5,col=col[2])) 
      grid.points(gen[2:10],log(data[scenario==s,missed])[2:10],default.units = 'native',pch=4,gp=gpar(cex=0.5,col=col[2])) 
      
    }
   
    
    
    if(column == 2){ 
      col = c(CONFIG$colsDark[3], CONFIG$colsDark2[3])
      
      s=24
      grid.lines(gen,log(data[scenario==s,notified]),default.units='native',gp=gpar(col=col[1]))
      grid.lines(gen,log(data[scenario==s,missed]),default.units='native',gp=gpar(col=col[1],lty='dashed'))
      grid.points(gen[2:10],log(data[scenario==s,notified])[2:10],default.units = 'native',pch=16,gp=gpar(cex=0.5,col=col[1])) 
      grid.points(gen[2:10],log(data[scenario==s,missed])[2:10],default.units = 'native',pch=16,gp=gpar(cex=0.5,col=col[1])) 
      
      s=40
      grid.lines(gen,log(data[scenario==s,notified]),default.units='native',gp=gpar(col=col[2]))
      grid.lines(gen,log(data[scenario==s,missed]),default.units='native',gp=gpar(col=col[2],lty='dashed'))
      grid.points(gen[2:10],log(data[scenario==s,notified])[2:10],default.units = 'native',pch=4,gp=gpar(cex=0.5,col=col[2])) 
      grid.points(gen[2:10],log(data[scenario==s,missed])[2:10],default.units = 'native',pch=4,gp=gpar(cex=0.5,col=col[2])) 
      
    }
    
    
  }
  
  if(row == 2){
    
    data = dataImportScenario[effOffspringParentMissed >= 0.1 & effOffspringParentMissed <=0.9]
    eff = seq(0.1,0.9,0.05)
    colLine = c(CONFIG$colsDark[1], CONFIG$colsDark[3])
    
    if(column == 1){ 
      
      data = data[ratio_missed == 'low',]
      grid.lines(data[,effOffspringParentMissed],data[,notified_time_exp],default.units='native',gp=gpar(col=colLine[1]))
      grid.lines(data[,effOffspringParentMissed],data[,missed_time_exp],default.units='native',gp=gpar(col=colLine[1],lty='dashed')) 
      
    }
    if(column == 2){ 
      
      data = data[ratio_missed == 'high',]
      grid.lines(data[,effOffspringParentMissed],data[,notified_time_exp],default.units='native',gp=gpar(col=colLine[2]))
      grid.lines(data[,effOffspringParentMissed],data[,missed_time_exp],default.units='native',gp=gpar(col=colLine[2],lty='dashed')) 
      
    }
    
  }
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  if(row == 1) grid.xaxis(at=seq(0,10,2), label=seq(0,10,2))
  if(row == 2) grid.xaxis(at=seq(0.1,0.9,0.2), label=seq(0.1,0.9,0.2))
  

  if(row == 1 & column == 1) grid.yaxis(at=log(c(1,10,100,1000)),label=c(1,10,100,1000))
  if(row == 1 & column == 2) grid.yaxis(at=log(c(1,10,100,1000)),label=c(1,10,100,1000))
  if(row == 2) grid.yaxis(at=seq(0,10,by=2),label=seq(0,10,by=2))
  
  if(row == 1) grid.text('Generation',y=unit(-3,'lines'))
  if(row == 2) grid.text(bquote(~epsilon[unlink]),y=unit(-3,'lines'))
  
  if(column == 1 & row == 1) grid.text('No. of cases',x=unit(-3.5,'lines'),rot=90)
  if(column == 1 & row == 2) grid.text('Generations to\nexponential growth',x=unit(-3,'lines'),rot=90)
  
  
  # caption
  if(row == 1 & column == 1) grid.text('A',x=unit(-3,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-3,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))  
  if(row == 2 & column == 1) grid.text('C',x=unit(-3,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))  
  if(row == 2 & column == 2) grid.text('D',x=unit(-3,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt'))) 
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}

png('figure/missedImports.png',height=16,width=16,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,0,0)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2)))

paneller(1,1)
paneller(1,2)

paneller(2,1)
paneller(2,2)

popViewport()
popViewport()
dev.off()

rm(paneller)
