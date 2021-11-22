source('codes/v6/covid_missInf_load_library.R')
source('codes/v6/covid_missInf_plot_colour.R')

# load posterior
period=c(1,2,3,4,5)
period.variant='W'

theta.period.all = list()
theta.period = list()

for(period.num in period){
  
  list.file = dir('output raw/20210906 test runs/wild neg binom total/', pattern = paste(period.variant, '_period_', period.num, sep=''))
  list.file
  
  for(f in 1:length(list.file)){
    
    load(paste('output raw/20210906 test runs/wild neg binom total/', list.file[f], sep=''))
    
    # burn-in first 5000 
    retain = seq(5001,60000,1)
    theta.period[[f]] = as.data.table(store$theta[retain,])
    theta.period[[f]]$chain = ceiling(f/4)
    theta.period[[f]]$period = period.num
    
    print(f)
  }
  
  theta.period.all[[period.num]] = rbindlist(theta.period, use.names = T)
  
  print(period.num)
}

rm(theta.period)

# load observed data
load('input/obs.data.RData')

# load time period
param.time.period = read_excel('data/param.xlsx', sheet = 'time.period.param')
param.time.period = data.table(param.time.period)
setnames(param.time.period, c('period', 'date.start', 'date.end'))

param.time.period = param.time.period[grep(period.variant, period)]
param.time.period[, period:=gsub(period.variant,'',period)]
param.time.period[, period:=as.numeric(period)]

param.time.period[, doy.start := as.numeric(as.Date(date.start)-as.Date('2019-12-31'))]
param.time.period[, doy.end := as.numeric(as.Date(date.end)-as.Date('2019-12-31'))]
param.time.period[, duration := doy.end-doy.start+1]


# plot
paneller=function(row = 1,column=1)
{
  niter = xmax = 275000
  xlm=c(0,xmax)
  if(row == 1) ylm=c(0,5)
  if(row == 2) ylm=c(0,5)
  if(row == 3) ylm=c(0,1)
  if(row == 4) ylm=c(0,1)
  if(row == 5) ylm=c(0,5)
  
  if(row == 1 & column == 2) ylm=c(0,75)
  if(row == 1 & column == 4) ylm=c(0,25)
  if(row == 1 & column == 5) ylm=c(0,50)
  
  
  innermargins = c(2,3,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  
  iteration = seq_len(niter)
  chain.colour = c(CONFIG$colsLight2[4],CONFIG$colsLight2[3],CONFIG$colsLight2[8],CONFIG$colsLight2[9])
  
  for(chain.num in 1:4){
    data = theta.period.all[[column]][chain==chain.num][[row]]
    
    # convert to average daily missed imports
    if(row == 1){
      data = data * mean(obs.data$daily.arrival.import.N[param.time.period$doy.start[column]:param.time.period$doy.end[column]])
    }
    
    grid.lines(iteration, data, default.units = 'native',gp=gpar(col=chain.colour[chain.num]))
  }
  

  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.xaxis(at=c(0,100000,200000),label = c(0,100000,200000))
 
  
  if(row == 1 & column %in% c(1,3)) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  if(row == 1 & column == 2) grid.yaxis(at=seq(0,75,by=25),label=seq(0,75,by=25))
  if(row == 1 & column == 4) grid.yaxis(at=seq(0,25,by=5),label=seq(0,25,by=5))
  if(row == 1 & column == 5) grid.yaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))

  if(row == 2) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))

  
  if(row == 3) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 4) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 5) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  
  if(row == 1 & column == 1) grid.text('Jan 18 - Feb 29', y=unit(8,'lines'))
  if(row == 1 & column == 2) grid.text('Mar 1 - Apr 6', y=unit(8,'lines'))
  if(row == 1 & column == 3) grid.text('Apr 7 - Jun 18', y=unit(8,'lines'))
  if(row == 1 & column == 4) grid.text('Jun 19 - Jul 12', y=unit(8,'lines'))
  if(row == 1 & column == 5) grid.text('Jul 13 - Dec 31', y=unit(8,'lines'))
  
  if(column == 1 & row == 1) grid.text('Avg daily missed\nimportation',x=unit(-3.5,'lines'),rot=90)
  if(column == 1 & row == 2) grid.text('R',x=unit(-3.5,'lines'),rot=90)
  if(column == 1 & row == 3) grid.text(bquote(~epsilon[link]~' (%)'),x=unit(-3.5,'lines'),rot=90) #~epsilon[op]
  if(column == 1 & row == 4) grid.text(bquote(~epsilon[unlink]~' (%)'),x=unit(-3.5,'lines'),rot=90) #~epsilon[op*'\'']
  if(column == 1 & row == 5) grid.text('k\'',x=unit(-3.5,'lines'),rot=90)
  
  if(column == 3 & row == 1) grid.text('Lockdown of borders', y=unit(4.5,'lines'), gp = gpar(fontsize = 8))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/trace_wild_total.png',height=25,width=30,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=5,ncol=5)))

for(r in 1:5){
  for(c in 1:5){
    paneller(r,c)
    
    print(c(r,c))
  }
}


grid.text('iteration',y=unit(-3,'lines'))


popViewport()
popViewport()
dev.off()

rm(paneller)







################## delta


# load posterior
period=c(1,2,3,4,5)
period.variant='D'

theta.period.all = list()
theta.period = list()

for(period.num in period){
  
  list.file = dir('output raw/20210906 test runs/delta neg binom total/', pattern = paste(period.variant, '_period_', period.num, sep=''))
  list.file
  
  for(f in 1:length(list.file)){
    
    load(paste('output raw/20210906 test runs/delta neg binom total/', list.file[f], sep=''))
    
    # burn-in first 5000 
    retain = seq(5001,60000,1)
    theta.period[[f]] = as.data.table(store$theta[retain,])
    theta.period[[f]]$chain = ceiling(f/5)
    theta.period[[f]]$period = period.num
    
    print(f)
  }
  
  theta.period.all[[period.num]] = rbindlist(theta.period, use.names = T)
  
  print(period.num)
}

rm(theta.period)

# load observed data
load('input/obs.data.RData')

# load time period
param.time.period = read_excel('data/param.xlsx', sheet = 'time.period.param')
param.time.period = data.table(param.time.period)
setnames(param.time.period, c('period', 'date.start', 'date.end'))

param.time.period = param.time.period[grep(period.variant, period)]
param.time.period[, period:=gsub(period.variant,'',period)]
param.time.period[, period:=as.numeric(period)]

param.time.period[, doy.start := as.numeric(as.Date(date.start)-as.Date('2019-12-31'))]
param.time.period[, doy.end := as.numeric(as.Date(date.end)-as.Date('2019-12-31'))]
param.time.period[, duration := doy.end-doy.start+1]


# plot
paneller=function(row = 1,column=1)
{
  niter = xmax = 275000
  xlm=c(0,xmax)
  if(row == 1) ylm=c(0,50)
  if(row == 2) ylm=c(0,5)
  if(row == 3) ylm=c(0,1)
  if(row == 4) ylm=c(0,1)
  if(row == 5) ylm=c(0,5)
  
  
  if(row == 2 & column == 4) ylm=c(0,30)
  
  innermargins = c(2,3,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  
  iteration = seq_len(niter)
  chain.colour = c(CONFIG$colsLight2[4],CONFIG$colsLight2[3],CONFIG$colsLight2[8],CONFIG$colsLight2[9])
  
  for(chain.num in 1:4){
    data = theta.period.all[[column]][chain==chain.num][[row]]
    
    # convert to average daily missed imports
    if(row == 1){
      data = data * mean(obs.data$daily.arrival.import.N[param.time.period$doy.start[column]:param.time.period$doy.end[column]])
    }
    
    grid.lines(iteration, data, default.units = 'native',gp=gpar(col=chain.colour[chain.num]))
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.xaxis(at=c(0,100000,200000),label = c(0,100000,200000))
  
  
  if(row == 1 & column == 1) grid.yaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))
  if(row == 1 & column != 1) grid.yaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))

  
  if(row == 2 & column != 4) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  if(row == 2 & column == 4) grid.yaxis(at=seq(0,30,by=5),label=seq(0,30,by=5))
  
  if(row == 3) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 4) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 5) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  
  if(row == 1 & column == 1) grid.text('Apr 1 - May 15', y=unit(8,'lines'))
  if(row == 1 & column == 2) grid.text('May 16 - Jun 13', y=unit(8,'lines'))
  if(row == 1 & column == 3) grid.text('Jun 14 - Jul 5', y=unit(8,'lines'))
  if(row == 1 & column == 4) grid.text('Jul 6 - Jul 11', y=unit(8,'lines'))
  if(row == 1 & column == 5) grid.text('Jul 13 - Aug 18', y=unit(8,'lines'))
  
  if(column == 1 & row == 1) grid.text('Avg daily missed\nimportation',x=unit(-3.5,'lines'),rot=90)
  if(column == 1 & row == 2) grid.text('R',x=unit(-3.5,'lines'),rot=90)
  if(column == 1 & row == 3) grid.text(bquote(~epsilon[link]~' (%)'),x=unit(-3.5,'lines'),rot=90) #~epsilon[op]
  if(column == 1 & row == 4) grid.text(bquote(~epsilon[unlink]~' (%)'),x=unit(-3.5,'lines'),rot=90) #~epsilon[op*'\'']
  if(column == 1 & row == 5) grid.text('k\'',x=unit(-3.5,'lines'),rot=90)
  
  # if(column == 3 & row == 1) grid.text('Lockdown of borders', y=unit(4.5,'lines'), gp = gpar(fontsize = 8))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/trace_delta_total.png',height=25,width=30,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=5,ncol=5)))

for(r in 1:5){
  for(c in 1:5){
    paneller(r,c)
    
    print(c(r,c))
  }
}

grid.text('iteration',y=unit(-3,'lines'))


popViewport()
popViewport()
dev.off()

rm(paneller)
