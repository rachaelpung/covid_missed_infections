source('codes/covid_missInf_load_library.R')
source('codes/covid_missInf_plot_colour.R')
source('codes/covid_missInf_functions.R')

load('output processed/theta_all.RData')

# plot R estimates
paneller=function(row = 1,column=1)
{
  if(column %in% c(1,2)) xlm=c(0.5,5.5)
  if(column == 3) xlm=c(0.5,4.5)
  if(column == 1) ylm=c(0,2)
  if(column == 2) ylm=c(0,2)
  if(column == 3) ylm=c(0,5)
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  param = 2; colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]
  if(column == 1){variant_fit = 'wild_neg_binom_period_'; period=1:5}
  if(column == 2){variant_fit = 'wild_neg_binom_total_period_'; period=1:5}
  if(column == 3){variant_fit = 'delta_neg_binom_total_period_'; period=1:4}
  
  # grid lines
  if(row == 1 & column %in% c(1,2)){
    grid.lines(xlm, c(0.5,0.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1,1), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1.5,1.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2,2), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  

  
  if(row == 1 & column == 3){
    grid.lines(xlm, c(1,1), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2,2), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(3,3), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(4,4), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(5,5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  
  for(period_num in period){
    
    extract = paste(variant_fit, period_num, sep ='')
    
    data = data.table(median = median(theta_all[[extract]][[param]]),
                      upper95 =  quantile(theta_all[[extract]][[param]], probs = 0.975),
                      lower95 =  quantile(theta_all[[extract]][[param]], probs = 0.025),
                      upper50 =  quantile(theta_all[[extract]][[param]], probs = 0.75),
                      lower50 =  quantile(theta_all[[extract]][[param]], probs = 0.25))
    
    
    
    grid.lines(c(period_num,period_num), c(data$lower95,data$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period_num,period_num), c(data$lower50,data$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period_num,data$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(column %in% c(1,2)) grid.yaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5))
  if(column == 3) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  
  if(column %in% c(1,2)) grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
                                    gp = gpar(fontsize=8))
  if(column == 3) grid.xaxis(at=seq(1,4,by=1),label=c('Apr 1-\nMay 12', 'May 13-\nJun 30', 'Jul 1-\nJul 17', 'Jul 18-\nAug 18'),
                             gp = gpar(fontsize=8))
  
  
  if(row == 1 & column == 1) {grid.text(bquote('R = R'[m]),x=unit(-3,'lines'),rot=90)
    grid.text('A',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  
  if(row == 1 & column == 2) {grid.text('B',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  if(row == 1 & column == 3) {grid.text('C',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/param_R_supp.png',height=8,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)

grid.text('Time period',y=unit(-1,'lines'))

popViewport()
popViewport()
dev.off()

rm(paneller)

# plot R estimates delta vaccine
paneller=function(row = 1,column=1)
{
  xlm=c(0.5,4.5)
  ylm=c(0,6)
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  param = 2; colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]
  variant_fit = 'delta_neg_binom_total_vaccine_period_'; period=1:4
  
  # grid lines
  grid.lines(xlm, c(1,1), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  grid.lines(xlm, c(2,2), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  grid.lines(xlm, c(3,3), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  grid.lines(xlm, c(4,4), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  grid.lines(xlm, c(5,5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  grid.lines(xlm, c(6,6), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  
  
  
  for(period_num in period){
    
    extract = paste(variant_fit, period_num, sep ='')
    
    data = data.table(median = median(theta_all[[extract]][[param]]),
                      upper95 =  quantile(theta_all[[extract]][[param]], probs = 0.975),
                      lower95 =  quantile(theta_all[[extract]][[param]], probs = 0.025),
                      upper50 =  quantile(theta_all[[extract]][[param]], probs = 0.75),
                      lower50 =  quantile(theta_all[[extract]][[param]], probs = 0.25))
    
    grid.lines(c(period_num,period_num), c(data$lower95,data$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period_num,period_num), c(data$lower50,data$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period_num,data$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.yaxis(at=seq(0,6,by=1),label=seq(0,6,by=1))
  
  grid.xaxis(at=seq(1,4,by=1),label=c('Apr 1-\nMay 12', 'May 13-\nJun 30', 'Jul 1-\nJul 17', 'Jul 18-\nAug 18'),
                             gp = gpar(fontsize=8))
  
  
  grid.text(bquote('R = R'[m]),x=unit(-3,'lines'),rot=90)
  # grid.text('A',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))

  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/param_R_vaccine_supp.png',height=8,width=8,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=1)))

paneller(1,1)


grid.text('Time period',y=unit(-1,'lines'))

popViewport()
popViewport()
dev.off()

rm(paneller)
