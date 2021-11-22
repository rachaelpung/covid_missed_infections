source('codes/v6/covid_missInf_load_library.R')
source('codes/v6/covid_missInf_plot_colour.R')

# load observed data
load('input/obs.data.RData')

# load time period
time.period = read_excel('data/param.xlsx', sheet = 'time.period.param')
time.period = data.table(time.period)
setnames(time.period, c('period', 'date.start', 'date.end'))

time.period[, variant:=gsub('(.*)\\d', '\\1',period)]
time.period[, period:=gsub(variant,'',period), by=seq_len(nrow(time.period))]
time.period[, period:=as.numeric(period)]

time.period[, doy.start := as.numeric(as.Date(date.start)-as.Date('2019-12-31'))]
time.period[, doy.end := as.numeric(as.Date(date.end)-as.Date('2019-12-31'))]
time.period[, duration := doy.end-doy.start+1]


# # load posterior
# period=c(1,2,3,4,5)
# 
# folder = c('wild neg binom', 'wild neg binom total', 'delta neg binom total')
# 
# theta.all = list()
# theta.period = list()
# 
# for(folder.num in 1:length(folder)){
#   
#   print(folder[folder.num])
#   
#   for(period.num in period){
#     
#     list.file = dir(paste('output raw/20210906 test runs/', folder[folder.num], '/', sep =''), pattern = paste('_period_', period.num, sep=''))
#     list.file
#     
#     for(f in 1:length(list.file)){
#       
#       load(paste('output raw/20210906 test runs/', folder[folder.num], '/', list.file[f], sep=''))
#       
#       # burn-in first 5000 and thinned every 50
#       retain = seq(5001,60000,50)
#       theta.period[[f]] = as.data.table(store$theta[retain,])
#       theta.period[[f]]$chain = ceiling(f/4)
#       theta.period[[f]]$period = period.num
#       
#       print(f)
#     }
#     
#     # theta.all[[period.num+(folder.num-1)*5]] = rbindlist(theta.period, use.names = T)
#     name.list = gsub(' ','.',paste(folder[folder.num], '.period.', period.num, sep=''))
#     theta.all[[name.list]] = rbindlist(theta.period, use.names = T)
#     
#     
#     print(period.num)
#   }
# }  
# 
# rm(theta.period)
# 
# save(theta.all, file='output processed/theta.all.RData')

load('output processed/theta.all.RData')

# plot param estimates 2X2
paneller=function(row = 1,column=1)
{
  xlm=c(0.5,5.5)
  if(row == 1 & column == 1) ylm=c(0,2.5)
  if(row == 1 & column == 2) ylm=c(0,1)
  if(row == 2 & column == 1) ylm=c(0,1)
  if(row == 2 & column == 2) ylm=c(log(0.001),log(100))
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(row == 1 & column == 1){param = 2; colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]}
  if(row == 1 & column == 2){param = 3; colMedian=CONFIG$cols[3]; colCI50=CONFIG$colsLight2[3]; colCI95=CONFIG$colsLight3[3]}
  if(row == 2 & column == 1){param = 4; colMedian=CONFIG$cols[2]; colCI50=CONFIG$colsLight2[2]; colCI95=CONFIG$colsLight3[2]}
  if(row == 2 & column == 2){param = 1; colMedian=CONFIG$cols[1]; colCI50=CONFIG$colsLight2[1]; colCI95=CONFIG$colsLight3[1]}
  
  # grid lines
  if(row == 1 & column == 1){
    grid.lines(xlm, c(0.5,0.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1,1), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1.5,1.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2,2), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2.5,2.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    
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
  
  
  
  for(period.num in 1:5){
    
    extract = paste('wild.neg.binom.period.', period.num, sep ='')
    
    data = data.table(median = median(theta.all[[extract]][[param]]),
                      upper95 =  quantile(theta.all[[extract]][[param]], probs = 0.975),
                      lower95 =  quantile(theta.all[[extract]][[param]], probs = 0.025),
                      upper50 =  quantile(theta.all[[extract]][[param]], probs = 0.75),
                      lower50 =  quantile(theta.all[[extract]][[param]], probs = 0.25))
    
    if(row==2 & column == 2 & period.num != 3){
      # convert factor of missed imported cases to average missed imported cases
      
      mean.import = mean(obs.data$daily.arrival.import.N[time.period[period == period.num & variant == 'W']$doy.start:time.period[period == period.num & variant == 'W']$doy.end])
      data = data * mean.import
      data = log(data)
    } else if(row==2 & column == 2 & period.num == 3){
      data[!is.na(data)] =NA
    }
    
    grid.lines(c(period.num,period.num), c(data$lower95,data$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period.num,period.num), c(data$lower50,data$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period.num,data$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row == 1 & column == 1) grid.yaxis(at=seq(0,2.5,by=0.5),label=seq(0,2.5,by=0.5))
  if(row == 1 & column == 2) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 2 & column == 1) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 2 & column == 2) grid.yaxis(at=log(c(0.001, 0.01, 0.1, 1, 10, 100)),label=c(0.001, 0.01, 0.1, 1, 10, 100))
  
  grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
             gp = gpar(fontsize=8))
  
  if(row == 2) grid.text('Time period',y=unit(-3,'lines'))
  
  if(row == 1 & column == 1) {grid.text('R',x=unit(-3,'lines'),rot=90)
    grid.text('A',x=unit(-2.5,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  if(row == 1 & column == 2) {grid.text(bquote(~epsilon[link] ~ ' (%)'),x=unit(-3,'lines'),rot=90)  #~epsilon[op]
    grid.text('B',x=unit(-2.5,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))}    
  if(row == 2 & column == 1) {grid.text(bquote(~epsilon[unlink] ~ ' (%)'),x=unit(-3,'lines'),rot=90) #~epsilon[op*'\'']
    grid.text('C',x=unit(-2.5,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))}  
  if(row == 2 & column == 2) {grid.text('Avg daily missed importation',x=unit(-3.5,'lines'),rot=90)
    grid.text('D',x=unit(-2.5,'lines'),y=unit(11.6,'lines'),gp=gpar(fontsize=unit(12,'pt')))}  
  if(row == 2 & column == 2) grid.text('Lockdown of borders',x=unit(6.5,'lines'),y=unit(7,'lines'),rot=90)
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/param_wild_link_unlink.png',height=16,width=16,units='cm',res=300,pointsize=10)
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




# plot param estimates 1X3
paneller=function(row = 1,column=1)
{
  xlm=c(0.5,5.5)
  if(row == 1 & column == 1) ylm=c(log(0.001),log(100))
  if(row == 1 & column == 2) ylm=c(0,1)
  if(row == 1 & column == 3) ylm=c(0,1)

  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(row == 1 & column == 1){param = 1; colMedian=CONFIG$cols[1]; colCI50=CONFIG$colsLight2[1]; colCI95=CONFIG$colsLight3[1]}
  if(row == 1 & column == 2){param = 3; colMedian=CONFIG$cols[3]; colCI50=CONFIG$colsLight2[3]; colCI95=CONFIG$colsLight3[3]}
  if(row == 1 & column == 3){param = 4; colMedian=CONFIG$cols[2]; colCI50=CONFIG$colsLight2[2]; colCI95=CONFIG$colsLight3[2]}
  
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
  
  
  
  for(period.num in 1:5){
    
    extract = paste('wild.neg.binom.period.', period.num, sep ='')
    
    data = data.table(median = median(theta.all[[extract]][[param]]),
                      upper95 =  quantile(theta.all[[extract]][[param]], probs = 0.975),
                      lower95 =  quantile(theta.all[[extract]][[param]], probs = 0.025),
                      upper50 =  quantile(theta.all[[extract]][[param]], probs = 0.75),
                      lower50 =  quantile(theta.all[[extract]][[param]], probs = 0.25))
    
    if(row==1 & column == 1 & period.num != 3){
      # convert factor of missed imported cases to average missed imported cases
      
      mean.import = mean(obs.data$daily.arrival.import.N[time.period[period == period.num & variant == 'W']$doy.start:time.period[period == period.num & variant == 'W']$doy.end])
      data = data * mean.import
      data = log(data)
    } else if(row==1 & column == 1 & period.num == 3){
      data[!is.na(data)] =NA
    }
    
    grid.lines(c(period.num,period.num), c(data$lower95,data$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period.num,period.num), c(data$lower50,data$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period.num,data$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row == 1 & column == 1) grid.yaxis(at=log(c(0.001, 0.01, 0.1, 1, 10, 100)),label=c(0.001, 0.01, 0.1, 1, 10, 100))
  if(row == 1 & column == 2) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 1 & column == 3) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  
  grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
             gp = gpar(fontsize=8))
  
  if(row == 1 & column == 1) {grid.text('Avg daily missed importation',x=unit(-3.5,'lines'),rot=90)
    grid.text('A',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))}  
  if(row == 1 & column == 1) grid.text('Lockdown of borders',x=unit(6.5,'lines'),y=unit(7,'lines'),rot=90)
  if(row == 1 & column == 2) {grid.text(bquote(~epsilon[link] ~ ' (%)'),x=unit(-3,'lines'),rot=90)  #~epsilon[op]
    grid.text('B',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))}    
  if(row == 1 & column == 3) {grid.text(bquote(~epsilon[unlink] ~ ' (%)'),x=unit(-3,'lines'),rot=90) #~epsilon[op*'\'']
    grid.text('C',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))}  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/param_wild_link_unlink_3_param.png',height=8,width=24,units='cm',res=300,pointsize=10)
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




