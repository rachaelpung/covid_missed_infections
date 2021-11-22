source('codes/v6/covid_missInf_load_library.R')
source('codes/v6/covid_missInf_plot_colour.R')
source('codes/v6/covid_missInf_functions.R')

load('output processed/theta.all.RData')

# # convert R missed into R notified, R effective, R notified import, R missed import
# 
# # load observed data
# load('input/obs.data.RData')
# 
# # load time period
# time.period = read_excel('data/param.xlsx', sheet = 'time.period.param')
# time.period = data.table(time.period)
# setnames(time.period, c('period', 'date.start', 'date.end'))
# 
# time.period[, variant:=gsub('(.*)\\d', '\\1',period)]
# time.period[, period:=gsub(variant,'',period), by=seq_len(nrow(time.period))]
# time.period[, period:=as.numeric(period)]
# 
# time.period[, doy.start := as.numeric(as.Date(date.start)-as.Date('2019-12-31'))]
# time.period[, doy.end := as.numeric(as.Date(date.end)-as.Date('2019-12-31'))]
# time.period[, duration := doy.end-doy.start+1]
# 
# param.W = list(time.period.matrix = time.period[variant =='W'])
# param.D = list(time.period.matrix = time.period[variant =='D'])
# 
# # fixed disease transmission parameters
# # generation interval
# param.W$gen.mean.log = param.D$gen.mean.log = log(5.4)
# param.W$gen.sd.log = param.D$gen.sd.log = 0.4
# 
# # incubation period
# param.W$incub.mean.log = param.D$incub.mean.log = 1.63
# param.W$incub.sd.log = param.D$incub.sd.log = 0.50
# 
# # load the distribution of duration from infection to isolation
# generation.matrix.local.N.wild = matrixGenerationLocal(param.W, obs.data, notified = T)
# generation.matrix.local.N.delta = matrixGenerationLocal(param.D, obs.data, notified = T)
# 
# # load the distribution of duration from infection to travel to isolation
# generation.matrix.imported.M.wild = matrixGenerationImport(param.W, obs.data, notified = F)
# generation.matrix.imported.M.delta = matrixGenerationImport(param.D, obs.data, notified = F)
# generation.matrix.imported.N.wild = matrixGenerationImport(param.W, obs.data, notified = T)
# generation.matrix.imported.N.delta = matrixGenerationImport(param.D, obs.data, notified = T)
# 
# 
# # reformat data to extract the distribution of duration from infection to isolation in respective time periods
# wild.period.start = param.W$time.period.matrix$doy.start-17
# wild.period.end = (param.W$time.period.matrix$doy.start-17 +13)
# delta.period.start = param.D$time.period.matrix$doy.start-456
# delta.period.end = (param.D$time.period.matrix$doy.start-456 + 13)
# 
# period = 1:5
# folder = c('wild neg binom', 'wild neg binom total', 'delta neg binom total')
# prop.generation.at.large.local.N = list()
# prop.generation.at.large.imported.N = list()
# prop.generation.at.large.imported.M = list()
# 
# for(folder.num in 1:length(folder)){
#   print(folder[folder.num])
#   
#   for(period.num in period){
#     name.list = gsub(' ','.',paste(folder[folder.num], '.period.', period.num, sep=''))
#     
#     if(grepl('wild',folder[folder.num])){
#       
#       prop.generation.at.large.local.N[[name.list]] = sum(generation.matrix.local.N.wild[wild.period.start[period.num],wild.period.start[period.num]:wild.period.end[period.num]])
#       prop.generation.at.large.imported.N[[name.list]] = sum(generation.matrix.imported.N.wild[wild.period.start[period.num],wild.period.start[period.num]:wild.period.end[period.num]])
#       prop.generation.at.large.imported.M[[name.list]] = sum(generation.matrix.imported.M.wild[wild.period.start[period.num],wild.period.start[period.num]:wild.period.end[period.num]])
#       
#     } else if(grepl('delta',folder[folder.num])){
#       
#       prop.generation.at.large.local.N[[name.list]] = sum(generation.matrix.local.N.delta[delta.period.start[period.num],delta.period.start[period.num]:delta.period.end[period.num]])
#       prop.generation.at.large.imported.N[[name.list]] = sum(generation.matrix.imported.N.delta[delta.period.start[period.num],delta.period.start[period.num]:delta.period.end[period.num]])
#       prop.generation.at.large.imported.M[[name.list]] = sum(generation.matrix.imported.M.delta[delta.period.start[period.num],delta.period.start[period.num]:delta.period.end[period.num]])
#       
#     }
#   }
# }
# 
# # compute the R of notified cases and effective R
# eigenNGM <- function(theta){
# 
#   sapply(1:nrow(theta), function(x){
# 
#     NGMComm = matrix(c((1-theta[[4]][x])*theta[[2]][x],  (1-theta[[3]][x])*theta[[8]][x],
#                           theta[[4]][x]*theta[[2]][x],      theta[[3]][x]*theta[[8]][x]),
#                      nrow = 2, ncol = 2, byrow = T)
# 
#     # missed offspring missed infector    missed offspring notified infector
#     # notified offspring missed infector    notified offspring notified infector
# 
#     max(eigen(NGMComm)$values)
# 
#   })
# 
# }
# 
# for(i in 1:length(theta.all)){
# 
#   col.add.name =  paste(c('R.notified.period.', 'R.eff.period.', 'R.notified.imported.period.', 'R.missed.imported.period.'), unique(theta.all[[i]]$period), sep='')
#   
#   theta.all[[i]][[col.add.name[1]]] = prop.generation.at.large[[i]]*theta.all[[i]][[2]]
#   theta.all[[i]][[col.add.name[2]]] = eigenNGM(theta.all[[i]])
#   
#   theta.all[[i]][[col.add.name[3]]] = prop.generation.at.large.imported.N[[i]]*theta.all[[i]][[2]]
#   theta.all[[i]][[col.add.name[4]]] = prop.generation.at.large.imported.M[[i]]*theta.all[[i]][[2]]
# 
#   print(i)
# }
# 
# save(theta.all, file='output processed/theta.all.RData')


# plot effective R
paneller=function(row = 1,column=1, outputTotal)
{
  xlm=c(0.5,5.5)
  if(row == 1) ylm=c(0,2.5)
  if(row == 2) ylm=log(c(0.01,30))
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(column == 1){param = 2; colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]}
  if(column == 2){param = 8; colMedian=CONFIG$cols[1]; colCI50=CONFIG$colsLight2[1]; colCI95=CONFIG$colsLight3[1]}
  if(column == 3){param = 9; colMedian=CONFIG$cols[4]; colCI50=CONFIG$colsLight2[4]; colCI95=CONFIG$colsLight3[4]}
  
  
  if(row == 1){variant.fit = 'wild.neg.binom.period.'}
  if(row == 2){variant.fit = 'delta.neg.binom.total.period.'}
  
  # grid lines
  if(row == 1){
    grid.lines(xlm, c(0.5,0.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1,1), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1.5,1.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2,2), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2.5,2.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  if(row == 2){
    grid.lines(xlm, log(c(0.1,0.1)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(1,1)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, log(c(10,10)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  
  for(period.num in 1:5){
    
    extract = paste(variant.fit, period.num, sep ='')
    
    data = data.table(median = median(theta.all[[extract]][[param]]),
                      upper95 =  quantile(theta.all[[extract]][[param]], probs = 0.975),
                      lower95 =  quantile(theta.all[[extract]][[param]], probs = 0.025),
                      upper50 =  quantile(theta.all[[extract]][[param]], probs = 0.75),
                      lower50 =  quantile(theta.all[[extract]][[param]], probs = 0.25))
    
    if(row == 2){ data = log(data)}
    
    grid.lines(c(period.num,period.num), c(data$lower95,data$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period.num,period.num), c(data$lower50,data$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period.num,data$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(row == 1) grid.yaxis(at=seq(0,2.5,by=0.5),label=seq(0,2.5,by=0.5))
  if(row == 2) grid.yaxis(at=log(c(0.01,0.1,1,10)),label=c(0.01,0.1,1,10))
  
  if(row == 1) grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 
                                                   'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
                          gp = gpar(fontsize=8))
  if(row == 2) grid.xaxis(at=seq(1,5,by=1),label=c('Apr 1-\nMay 15', 'May 16-\nJun 13', 'Jun 14-\nJul 5', 
                                                   'Jul 6-\nJul 11', 'Jul 12-\nAug 18'),
                          gp = gpar(fontsize=8))
  

  if(column == 1) {grid.text(bquote('R = R'[missed]),x=unit(-3,'lines'),rot=90)}
  if(column == 2) {grid.text(bquote('R'[notified]),x=unit(-3,'lines'),rot=90)}    
  if(column == 3) {grid.text(bquote('R'[eff]),x=unit(-3,'lines'),gp=gpar(fontsize=unit(9,'pt')),rot=90)}  
  
  if(row == 1 & column ==1){grid.text('A',x=unit(-3,'lines'),y=unit(13.5,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 1 & column ==2){grid.text('B',x=unit(-3,'lines'),y=unit(13.5,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 1 & column ==3){grid.text('C',x=unit(-3,'lines'),y=unit(13.5,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 2 & column ==1){grid.text('D',x=unit(-3,'lines'),y=unit(13.5,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 2 & column ==2){grid.text('E',x=unit(-3,'lines'),y=unit(13.5,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  if(row == 2 & column ==3){grid.text('F',x=unit(-3,'lines'),y=unit(13.5,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}



png('figure/param_R_eff.png',height=16,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=3)))

paneller(1,1)
paneller(1,2)
paneller(1,3)

paneller(2,1)
paneller(2,2)
paneller(2,3)

grid.text('Time period',y=unit(-1,'lines'))

popViewport()
popViewport()
dev.off()




# plot R estimates
paneller=function(row = 1,column=1)
{
  xlm=c(0.5,5.5)
  if(column %in% c(1,2)) ylm=c(0,2.5)
  # if(column == 3) ylm=log(c(0.1,30))
  
  innermargins = c(2,2,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  param = 2; colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]
  if(column == 1){variant.fit = 'wild.neg.binom.period.'}
  if(column == 2){variant.fit = 'wild.neg.binom.total.period.'}
  # if(column == 3){variant.fit = 'delta.neg.binom.total.period.'}
  
  # grid lines
  if(row == 1 & column %in% c(1,2)){
    grid.lines(xlm, c(0.5,0.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1,1), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(1.5,1.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2,2), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(2.5,2.5), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  # if(row == 1 & column == 3){
  #   grid.lines(xlm, log(c(1,1)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  #   grid.lines(xlm, log(c(10,10)), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  # }
  
  
  for(period.num in 1:5){
    
    extract = paste(variant.fit, period.num, sep ='')
    
    data = data.table(median = median(theta.all[[extract]][[param]]),
                      upper95 =  quantile(theta.all[[extract]][[param]], probs = 0.975),
                      lower95 =  quantile(theta.all[[extract]][[param]], probs = 0.025),
                      upper50 =  quantile(theta.all[[extract]][[param]], probs = 0.75),
                      lower50 =  quantile(theta.all[[extract]][[param]], probs = 0.25))
    
    # if(column == 3){ data = log(data)}
    
    grid.lines(c(period.num,period.num), c(data$lower95,data$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
    grid.lines(c(period.num,period.num), c(data$lower50,data$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
    grid.points(period.num,data$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian)) 
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(column %in% c(1,2)) grid.yaxis(at=seq(0,2.5,by=0.5),label=seq(0,2.5,by=0.5))
  if(column == 3) grid.yaxis(at=log(c(0.1,1,10)),label=c(0.1,1,10))
  
  if(column %in% c(1,2)) grid.xaxis(at=seq(1,5,by=1),label=c('Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
                                    gp = gpar(fontsize=8))
  # if(column == 3) grid.xaxis(at=seq(1,5,by=1),label=c('Apr 1-\nMay 15', 'May 16-\nJun 13', 'Jun 14-\nJul 5', 'Jul 6-\nJul 11', 'Jul 12-\nAug 18'),
  #                            gp = gpar(fontsize=8))
  
  
  if(row == 1 & column == 1) {grid.text(bquote('R = R'[missed]),x=unit(-3,'lines'),rot=90)
    grid.text('A',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  
  if(row == 1 & column == 2) {grid.text('B',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  # if(row == 1 & column == 3) {grid.text('C',x=unit(-2.5,'lines'),y=unit(10,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}


png('figure/param_R.png',height=8,width=16,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))

paneller(1,1)
paneller(1,2)
# paneller(1,3)

grid.text('Time period',y=unit(-1,'lines'))

popViewport()
popViewport()
dev.off()

rm(paneller)
