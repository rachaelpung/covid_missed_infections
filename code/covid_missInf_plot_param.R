source('codes/covid_missInf_load_library.R')
source('codes/covid_missInf_plot_colour.R')

load('input/obs_data.RData')
load('input/time_period.RData')
load('output processed/theta_all.RData')


# plot param estimates 1X3 & R
paneller=function(row = 1,column=1)
{
  if(row!=3) xlm=c(0.5,5.5)
  if(row==3) xlm=c(0.5,4.5)
  
  if(row == 1 & column == 1) ylm=c(log(0.001),log(100))
  if(row == 1 & column == 2) ylm=c(0,1)
  if(row == 1 & column == 3) ylm=c(0,1)
  
  if(row == 2) ylm=c(0,2)
  if(row == 3) ylm=c(0,5)
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(row == 1 & column == 1){param = 1; colMedian=CONFIG$cols[1]; colCI50=CONFIG$colsLight2[1]; colCI95=CONFIG$colsLight3[1]}
  if(row == 1 & column == 2){param = 3; colMedian=CONFIG$cols[3]; colCI50=CONFIG$colsLight2[3]; colCI95=CONFIG$colsLight3[3]}
  if(row == 1 & column == 3){param = 4; colMedian=CONFIG$cols[2]; colCI50=CONFIG$colsLight2[2]; colCI95=CONFIG$colsLight3[2]}
  
  if(row != 1 & column == 1){param = 2; colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]}
  if(row != 1 & column == 2){param = 8; colMedian=CONFIG$cols[6]; colCI50=CONFIG$colsLight2[6]; colCI95=CONFIG$colsLight3[6]}
  if(row != 1 & column == 3){param = 9; colMedian=CONFIG$cols[4]; colCI50=CONFIG$colsLight2[4]; colCI95=CONFIG$colsLight3[4]}
  
  if(row != 3){variant_fit = 'wild_neg_binom_period_'; period=1:5}
  if(row == 3){variant_fit = 'delta_neg_binom_total_period_'; period=1:4}
  
  
  # grid lines
  if((row == 1 & column == 2) | (row == 1 & column == 3)){
    grid.lines(xlm, c(0.25,0.25), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.5,0.5),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.75,0.75), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  if(row == 1 & column == 1){
    grid.lines(xlm, log(c(0.01,0.01)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(0.1,0.1)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(1,1)),       default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(10,10)),     default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(100,100)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  if(row == 2){
    grid.lines(xlm, c(0.5,0.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1,1), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1.5,1.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2,2), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    # grid.lines(xlm, c(2.5,2.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  if(row == 3){
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
    
    if(row ==1 & column == 1 & period_num != 3){
      # convert factor of missed imported cases to average missed imported cases
      
      mean_import = mean(obs_data$daily_arrival_import_N[time_period[period == period_num & variant == 'W']$doy_start:time_period[period == period_num & variant == 'W']$doy_end])
      data = data * mean_import
      data = log(data)
    } else if(row==1 & column == 1 & period_num == 3){
      data[!is.na(data)] =NA
    }
    
   
    
    grid.lines(c(period_num,period_num), c(data$lower95,data$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period_num,period_num), c(data$lower50,data$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period_num,data$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row == 1 & column == 1) grid.yaxis(at=log(c(0.001, 0.01, 0.1, 1, 10, 100)),label=c(0.001, 0.01, 0.1, 1, 10, 100))
  if(row == 1 & column == 2) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 1 & column == 3) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 2) grid.yaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5))
  if(row == 3) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  
  
  if(row != 3) grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 
                                                   'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
                          gp = gpar(fontsize=8))
  if(row == 3) grid.xaxis(at=seq(1,4,by=1),label=c('Apr 1-\nMay 12', 'May 13-\nJun 30', 'Jul 1-\nJul 17', 
                                                   'Jul 18-\nAug 18'),
                          gp = gpar(fontsize=8))
  
  
  if(row == 1 & column == 1) {grid.text('Avg daily missed importation',x=unit(-3.5,'lines'),rot=90)}  
  if(row == 1 & column == 1) {grid.text('Lockdown of borders',x=unit(6.5,'lines'),y=unit(7,'lines'),rot=90)}
  if(row == 1 & column == 2) {grid.text(bquote(~epsilon[ct] ~ ' (%)'),x=unit(-3,'lines'),rot=90)}  #~epsilon[link]  
  if(row == 1 & column == 3) {grid.text(bquote(~epsilon[cf] ~ ' (%)'),x=unit(-3,'lines'),rot=90)} #~epsilon[unlink]  
  
  if(row != 1 & column == 1) {grid.text(bquote('R = R'[m]),x=unit(-3,'lines'),rot=90)}  #'R'[missed]
  if(row != 1 & column == 2) {grid.text(bquote('R'[n]),x=unit(-3,'lines'),rot=90)}    #'R'[notified]
  if(row != 1 & column == 3) {grid.text(bquote('R'[eff]),x=unit(-3,'lines'),gp=gpar(fontsize=unit(9,'pt')),rot=90)}  
  
  if(row == 1 & column ==1){grid.text('A',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 1 & column ==2){grid.text('B',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 1 & column ==3){grid.text('C',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 2 & column ==1){grid.text('D',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 2 & column ==2){grid.text('E',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 2 & column ==3){grid.text('F',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 3 & column ==1){grid.text('G',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 3 & column ==2){grid.text('H',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 3 & column ==3){grid.text('I',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/param_all_fig_2.png',height=24,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=3)))

for(c in 1:3){
  for(r in 1:3){
    paneller(r,c)
  }
}

grid.text('Time period',y=unit(-1,'lines'))

popViewport()
popViewport()
dev.off()

rm(paneller)


# plot param estimates 2x2
paneller=function(row = 1,column=1)
{
  xlm=c(0.5,6.5)
  if(row == 1 & column == 1) ylm=c(log(0.001),log(100))
  if(row == 1 & column == 2) ylm=c(0,1)
  if(row == 2 & column == 1) ylm=c(0,1)
  if(row == 2 & column == 2) ylm=c(0,2)
  
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(row == 1 & column == 1){param = 1; colMedian=CONFIG$cols[1]; colCI50=CONFIG$colsLight2[1]; colCI95=CONFIG$colsLight3[1]}
  if(row == 1 & column == 2){param = 3; colMedian=CONFIG$cols[3]; colCI50=CONFIG$colsLight2[3]; colCI95=CONFIG$colsLight3[3]}
  if(row == 2 & column == 1){param = 4; colMedian=CONFIG$cols[2]; colCI50=CONFIG$colsLight2[2]; colCI95=CONFIG$colsLight3[2]}
  if(row == 2 & column == 2){param = 2; colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]}
 
  variant_fit = 'wild_neg_binom_period_'; period=1:6

  # grid lines
  if((row == 1 & column == 2) | (row == 2 & column == 1)){
    grid.lines(xlm, c(0.25,0.25), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.5,0.5),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.75,0.75), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  if(row == 1 & column == 1){
    grid.lines(xlm, log(c(0.01,0.01)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(0.1,0.1)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(1,1)),       default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(10,10)),     default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(100,100)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  if(row == 2 & column == 2){
    grid.lines(xlm, c(0.5,0.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1,1), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1.5,1.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2,2), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    # grid.lines(xlm, c(2.5,2.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  
  for(period_num in period){
    
    extract = paste(variant_fit, period_num, sep ='')
    
    data = data.table(median = median(theta_all[[extract]][[param]]),
                      upper95 =  quantile(theta_all[[extract]][[param]], probs = 0.975),
                      lower95 =  quantile(theta_all[[extract]][[param]], probs = 0.025),
                      upper50 =  quantile(theta_all[[extract]][[param]], probs = 0.75),
                      lower50 =  quantile(theta_all[[extract]][[param]], probs = 0.25))
    
    if(row ==1 & column == 1 & period_num != 3){
      # convert factor of missed imported cases to average missed imported cases
      
      mean_import = mean(obs_data$daily_arrival_import_N[time_period[period == period_num & variant == 'W']$doy_start:time_period[period == period_num & variant == 'W']$doy_end])
      data = data * mean_import
      data = log(data)
    } else if(row==1 & column == 1 & period_num == 3){
      data[!is.na(data)] =NA
    }
    
    
    
    grid.lines(c(period_num,period_num), c(data$lower95,data$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period_num,period_num), c(data$lower50,data$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period_num,data$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row == 1 & column == 1) grid.yaxis(at=log(c(0.001, 0.01, 0.1, 1, 10, 100)),label=c(0.001, 0.01, 0.1, 1, 10, 100))
  if(row == 1 & column == 2) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 2 & column == 1) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 2 & column == 2) grid.yaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5))
  
  
  grid.xaxis(at=seq(1,6,by=1),label=c('2020 Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 
                                                   'Jun 19-\nJul 12', 'Jul 13-\nDec 31', '2021 Jan 1-\nMar 31'),
                          gp = gpar(fontsize=5))
  
  
  if(row == 1 & column == 1) {grid.text('Avg daily missed importation',x=unit(-3.5,'lines'),rot=90)}  
  if(row == 1 & column == 1) {grid.text('Lockdown of borders',x=unit(5.5,'lines'),y=unit(7,'lines'),rot=90)}
  if(row == 1 & column == 2) {grid.text(bquote(~epsilon[ct] ~ ' (%)'),x=unit(-3,'lines'),rot=90)}  #~epsilon[link]  
  if(row == 2 & column == 1) {grid.text(bquote(~epsilon[cf] ~ ' (%)'),x=unit(-3,'lines'),rot=90)} #~epsilon[unlink]  
  
  if(row == 2 & column == 2) {grid.text(bquote('R = R'[m]),x=unit(-3,'lines'),rot=90)}  #'R'[missed]
  
  if(row == 1 & column ==1){grid.text('A',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 1 & column ==2){grid.text('B',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 2 & column ==1){grid.text('C',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 2 & column ==2){grid.text('D',x=unit(-3,'lines'),y=unit(15,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/param_wild_extend_supp.png',height=16,width=16,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2)))

for(c in 1:2){
  for(r in 1:2){
    paneller(r,c)
  }
}

grid.text('Time period',y=unit(-1,'lines'))

popViewport()
popViewport()
dev.off()
