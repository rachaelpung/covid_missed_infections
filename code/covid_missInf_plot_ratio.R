source('codes/v6/covid_missInf_load_library.R')
source('codes/v6/covid_missInf_plot_colour.R')

load('output processed/inc.all.RData')

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

# plot eigenvector ratios
paneller=function(row=1, column=1)
{

  xlm=c(0.5,5.5)
  # if(column == 1){ylm=log(c(0.01,10))}
  # if(column == 2){ylm=log(c(0.01,1000))}
  
  if(column == 1){ylm=c(0,1)}
  if(column == 2){ylm=c(0,1)}
 
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(c(3,4,1,1),xscale=xlm,yscale=ylm))
  
  if(row == 1 & column ==1){variant.fit = 'wild.neg.binom.mod';        extract.start=time.period[variant == 'W']$doy.start; extract.end=time.period[variant == 'W']$doy.end}
  if(row == 1 & column ==2){variant.fit = 'delta.neg.binom.total.mod'; extract.start=time.period[variant == 'D']$doy.start - 456; extract.end=time.period[variant == 'D']$doy.end - 456}
  
  # grid lines
  # grid.lines(xlm, log(c(0.01,0.01)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  # grid.lines(xlm, log(c(0.1,0.1)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  # grid.lines(xlm, log(c(1,1)),       default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  # grid.lines(xlm, log(c(10,10)),     default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  # grid.lines(xlm, log(c(100,100)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  # grid.lines(xlm, log(c(1000,1000)),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  
  grid.lines(xlm, c(0.25,0.25), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  grid.lines(xlm, c(0.5,0.5),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  grid.lines(xlm, c(0.75,0.75), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    
  
  for(period.num in 1:5){
    
    extract.fit = grep(variant.fit, names(inc.all))
    data = inc.all[extract.fit]

    for(i in 1:length(data)){data[[i]] = data[[i]][, extract.start[period.num]:extract.end[period.num]]}

    missed = apply(data[[3]]+ data[[4]],1,sum)
    notified = apply(data[[1]]+ data[[2]],1,sum)
    # ratio = missed/notified
    ratio = missed/(notified+missed)

    ratio = data.table(median = median(ratio),
                       upper95 = quantile(ratio, probs = 0.975),
                       lower95 = quantile(ratio, probs = 0.025),
                       upper50 = quantile(ratio, probs = 0.75),
                       lower50 = quantile(ratio, probs = 0.25))

    # ratio = log(ratio)

    colMedian=CONFIG$cols[4]; colCI50=CONFIG$colsLight2[4]; colCI95=CONFIG$colsLight3[4]

    grid.lines(c(period.num-0.2,period.num-0.2), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period.num-0.2,period.num-0.2), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period.num-0.2,ratio$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian))

    
    unlinked = apply(data[[1]],1,sum)
    linked   = apply(data[[2]],1,sum)
    # ratio = unlinked/linked
    ratio = unlinked/(linked+unlinked)

    ratio = data.table(median = median(ratio),
                       upper95 = quantile(ratio, probs = 0.975),
                       lower95 = quantile(ratio, probs = 0.025),
                       upper50 = quantile(ratio, probs = 0.75),
                       lower50 = quantile(ratio, probs = 0.25))

    # ratio = log(ratio)

    colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]

    grid.lines(c(period.num+0.2,period.num+0.2), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period.num+0.2,period.num+0.2), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period.num+0.2,ratio$median,default.units = 'native',pch=17,gp=gpar(cex=0.8,col=colMedian))

   
  }
  
  
  if(column == 1) grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 
                                                     'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
                             gp = gpar(fontsize=8))
  if(column == 2) grid.xaxis(at=seq(1,5,by=1),label=c('Apr 1-\nMay 15', 'May 16-\nJun 13', 'Jun 14-\nJul 5', 
                                                     'Jul 6-\nJul 11', 'Jul 12-\nAug 18'),
                             gp = gpar(fontsize=8))
  
  # if(column == 1) grid.yaxis(at=log(c(0.01,0.1,1,10)),label=c(0.01,0.1,1,10))
  # if(column == 2) grid.yaxis(at=log(c(0.01,0.1,1,10,100,1000)),label=c(0.01,0.1,1,10,100,1000))
  
  if(column == 1) grid.yaxis(at=seq(0,1,0.25),label=seq(0,100,25))
  if(column == 2) grid.yaxis(at=seq(0,1,0.25),label=seq(0,100,25))

  if(column == 1) grid.text('Proportion (%)',x=unit(-3.5,'lines'),rot=90)


  
  # caption
  if(row == 1 & column == 1) grid.text('A',x=unit(-3,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))
  if(row == 1 & column == 2) grid.text('B',x=unit(-3,'lines'),y=unit(11.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))  
 
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/proportion.png',height=8,width=16,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(1,0,0,0)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))

paneller(1,1)
paneller(1,2)

grid.text('Time period',y=unit(-0.5,'lines'))

popViewport()
popViewport()
dev.off()

rm(paneller)