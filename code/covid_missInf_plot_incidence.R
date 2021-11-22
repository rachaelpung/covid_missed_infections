source('codes/v6/covid_missInf_load_library.R')
source('codes/v6/covid_missInf_plot_colour.R')

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
# extract.wild = seq_len(sum(time.period[variant =='W']$duration))
# extract.delta = seq_len(sum(time.period[variant =='D']$duration))
# 
# 
# # load processed/thinned & burned posterior
# period=c(1,2,3,4,5)
# 
# folder = c('wild neg binom', 'wild neg binom total', 'delta neg binom total')
# 
# inc.all = list() # incidence of all wild and delta cases
# inc.variant = list() # temporary variable, variant specific incidence
# inc.variant.period = list() # temporary variable, variant and period specific incidence
# 
# for(folder.num in 1:length(folder)){
#   
#   print(folder[folder.num])
#   
#   for(period.num in period){
#     
#     
#     list.file = dir(paste('output processed/', folder[folder.num], '/', sep =''), pattern = as.character(period.num))
#     list.file
#     
#     load(paste('output processed/', folder[folder.num], '/', list.file, sep=''))
#     
#     # remove theta and likelihood estimates
#     store.period[[1]] = store.period[[2]] = NULL
#     
#     
#     if(period.num == 1){
#       inc.variant = store.period
#     } else{
#       
#       inc.variant.period = store.period
#       
#       for(i in 1:length(inc.variant.period)){
#         
#         
#         inc.variant.period[[i]][,1:13] = sweep(inc.variant.period[[i]][,1:13], 2, inc.variant.period.median[[i]], FUN='-')
#         extract = (length(inc.variant[[i]])-12) : length(inc.variant[[i]])
#         
#         
#         inc.variant[[i]][,extract] = inc.variant[[i]][,extract,with=FALSE] + inc.variant.period[[i]][,1:13]
#         inc.variant[[i]] = cbind(inc.variant[[i]], inc.variant.period[[i]][,14:length(inc.variant.period[[i]])])
#         
#         if (any(inc.variant[[i]] < 0)){
#           print("Negative values detected")
#           break
#         }
#       }
#     }
#     
#     
#     
#     
#     # extract median estimates to substract from the stored outputs of next time period
#     inc.variant.period.median = store.period
#     
#     for(i in 1:length(inc.variant.period.median)){
#       inc.variant.period.median[[i]] = apply(inc.variant.period.median[[i]], 2, median)
#       extract = (length(inc.variant.period.median[[i]])-12) : length(inc.variant.period.median[[i]])
#       inc.variant.period.median[[i]] = inc.variant.period.median[[i]][extract]
#     }
#     
#     
#     
#     print(period.num)
#   }
#   
#   name.list = paste(gsub(' ','.',folder[folder.num]), names(inc.variant), sep=".") 
#   
#   names(inc.variant) = name.list
#   
#   inc.all = c(inc.all, inc.variant)
#   
# }
# 
# rm(inc.variant, inc.variant.period, inc.variant.period.median, period, period.num)
# 
# 
# 
# # add zeros to wild type incidence
# zero.mat = matrix(0, nrow = dim(inc.all[[1]])[1], ncol = time.period[variant =='W']$doy.start[1]-1)
# 
# for(i in 1:length(inc.all)){
#   
#   inc.all[[i]] = as.matrix(inc.all[[i]])
#   
#   if(grepl('wild',names(inc.all)[i])){
#     
#     inc.all[[i]] = inc.all[[i]][,extract.wild]
#     inc.all[[i]] = cbind(zero.mat,inc.all[[i]])
#     
#   } else if(grepl('delta',names(inc.all)[i])){
#     
#     inc.all[[i]] = inc.all[[i]][,extract.delta]
#     
#   }
#   
#   colnames(inc.all[[i]]) = NULL
# }
# 
# rm(i, zero.mat)
# 
# save(inc.all, file='output processed/inc.all.RData')

load('output processed/inc.all.RData')


# generate incidence plot for wild type fitted using linked and unlinked cases
paneller=function(row = 1,column=1)
{
  
  xlm=c(1, 367)
  if(row == 1 & column == 1) ylm=c(0,100)  # imported cases
  if(row == 1 & column == 2) ylm=c(0,100)  # linked cases
  if(row == 2 & column == 1) ylm=c(0,50)   # unlinked cases
  if(row == 2 & column == 2) ylm=c(0,500) # missed cases
  
  innermargins = c(2,3,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  time = 1:366
  if(row == 1 & column == 1){colObs.all =  CONFIG$colsLight2[1]; colObs.no.shn =  CONFIG$cols[1]}
  if(row == 1 & column == 2){colObs = CONFIG$colsLight2[3]; colMedian =  CONFIG$cols[3]; colCI = CONFIG$colsLight3[3]}
  if(row == 2 & column == 1){colObs = CONFIG$colsLight2[2]; colMedian =  CONFIG$cols[2]; colCI = CONFIG$colsLight3[2]}
  if(row == 2 & column == 2){colMedian =  CONFIG$cols[5]; colCI95 = CONFIG$colsLight3[5]; colCI50 = CONFIG$colsLight2[5]}
  
  # circuit breaker from 7 Apr to 1 Jun
  grid.polygon(c(98, 98, 153.5,153.5),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#ECECEC',col=NA))
  
  grid.polygon(c(153.5, 153.5, 170.5,170.5),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#F5F5F5',col=NA))
  
  grid.polygon(c(170.5, 170.5, 362,362),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
  
  # plot imported cases
  if(row == 1 & column == 1){
    
    data.obs.all = obs.data$daily.arrival.import.N
    data.obs.no.shn = obs.data$daily.arrival.import.N.no.shn
    
    for(t in 1:length(time)){
      
      if(data.obs.all[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data.obs.all[t],data.obs.all[t],0),
                     default.units = 'native',gp=gpar(fill=colObs.all,col=NA)) 
      }
      
      if(data.obs.no.shn[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data.obs.no.shn[t],data.obs.no.shn[t],0),
                     default.units = 'native',gp=gpar(fill=colObs.no.shn,col=NA))
      }
      
    }
  }
  
  # plot linked cases
  if(row == 1 & column == 2){
    
    data.model = data.table(median = apply(inc.all$wild.neg.binom.mod.daily.local.N.linked.by.doy.isolate, 2, median),
                            upper = apply(inc.all$wild.neg.binom.mod.daily.local.N.linked.by.doy.isolate,2, quantile, 0.975),
                            lower = apply(inc.all$wild.neg.binom.mod.daily.local.N.linked.by.doy.isolate,2, quantile, 0.025))
    
    data.obs = obs.data$daily.local.N.linked
    
    for(t in 1:length(time)){
      
      if(data.obs[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data.obs[t],data.obs[t],0),
                     default.units = 'native',gp=gpar(fill=colObs,col=NA)) 
      }
    }
    
    grid.lines(time,data.model$median,default.units = 'native',gp=gpar(col=colMedian)) 
    grid.polygon(c(time,rev(time)),
                 c(data.model$lower,rev(data.model$upper)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI))
  }
  
  # plot unlinked cases
  if(row == 2 & column == 1){
    
    data.model = data.table(median = apply(inc.all$wild.neg.binom.mod.daily.local.N.unlinked.by.doy.isolate, 2, median),
                            upper =  apply(inc.all$wild.neg.binom.mod.daily.local.N.unlinked.by.doy.isolate,2, quantile, 0.975),
                            lower = apply(inc.all$wild.neg.binom.mod.daily.local.N.unlinked.by.doy.isolate,2, quantile, 0.025))
    
    data.obs = obs.data$daily.local.N.unlinked
    
    for(t in 1:length(time)){
      
      if(data.obs[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data.obs[t],data.obs[t],0),
                     default.units = 'native',gp=gpar(fill=colObs,col=NA)) 
      }
    }
    
    grid.lines(time,data.model$median,default.units = 'native',gp=gpar(col=colMedian)) 
    grid.polygon(c(time,rev(time)),
                 c(data.model$lower,rev(data.model$upper)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI))
  }
  
  # plot missed cases
  if(row == 2 & column == 2){
    
    data.model = inc.all$wild.neg.binom.mod.daily.local.M.unlinked + inc.all$wild.neg.binom.mod.daily.local.M.linked
    
    data.model = data.table(median = apply(data.model, 2, median),
                            upper95 =  apply(data.model,2, quantile, 0.975),
                            lower95 = apply(data.model,2, quantile, 0.025),
                            upper50 =  apply(data.model,2, quantile, 0.75),
                            lower50 = apply(data.model,2, quantile, 0.25))
    
    grid.polygon(c(time,rev(time)),
                 c(data.model$lower50,rev(data.model$upper50)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI50))
    grid.polygon(c(time,rev(time)),
                 c(data.model$lower95,rev(data.model$upper95)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI95))
    grid.lines(time,data.model$median,default.units = 'native',gp=gpar(col=colMedian))
    
  }
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  if(row == 1) grid.xaxis(at=c(1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367),
                          label=rep('',13))
  if(row == 2) grid.xaxis(at=c(1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367),
                          label = c('2020 Jan 1', 'Feb 1', 'Mar 1', 'Apr 1', 'May 1', 'Jun 1', 'Jul 1', 
                                    'Aug 1', 'Sep 1', 'Oct 1', 'Nov 1', 'Dec 1', '2021 Jan 1'))
  
  if(row == 1 & column == 1) {grid.yaxis(at=seq(0,100,by=10),label=seq(0,100,by=10))
    grid.text('A',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  if(row == 1 & column == 2) {grid.yaxis(at=seq(0,100,by=10),label=seq(0,100,by=10))
    grid.text('B',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  if(row == 2 & column == 1) {grid.yaxis(at=seq(0,50,by=5),label=seq(0,50,by=5))
    grid.text('C',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))} 
  if(row == 2 & column == 2) {grid.yaxis(at=seq(0,500,by=100),label=seq(0,500,by=100))
    grid.text('D',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  
  if(row == 2) grid.text('Day',y=unit(-3,'lines'))
  if(column == 1 ) grid.text('Incidence',x=unit(-3.5,'lines'),rot=90)
  
  grid.text('Partial lockdown',x=unit(12.5,'lines'),y=unit(17,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 1',x=unit(19.5,'lines'),y=unit(18.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 2',x=unit(21.5,'lines'),y=unit(18.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 3',x=unit(44.9,'lines'),y=unit(18.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/incidence_wild_link_unlink.png',height=20,width=40,units='cm',res=300,pointsize=10)
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


# generate incidence plot for wild type fitted using total cases
paneller=function(row = 1,column=1)
{
  
  xlm=c(1, 367)
  if(row == 1 & column == 1) ylm=c(0,150)  # total notified cases
  if(row == 1 & column == 2) ylm=c(0,500) # missed cases
  
  innermargins = c(2,3,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  time = 1:366
  if(row == 1 & column == 1){colObs = CONFIG$colsLight2[4]; colMedian =  CONFIG$cols[4]; colCI = CONFIG$colsLight3[4]}
  if(row == 1 & column == 2){colMedian =  CONFIG$cols[5]; colCI95 = CONFIG$colsLight3[5]; colCI50 = CONFIG$colsLight2[5]}
  
  # circuit breaker from 7 Apr to 1 Jun
  grid.polygon(c(98, 98, 153.5,153.5),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#ECECEC',col=NA))
  
  grid.polygon(c(153.5, 153.5, 170.5,170.5),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#F5F5F5',col=NA))
  
  grid.polygon(c(170.5, 170.5, 362,362),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
  
  
  # plot all cases
  if(row == 1 & column == 1){
    
    data.model = inc.all$wild.neg.binom.total.mod.daily.local.N.linked.by.doy.isolate + 
                 inc.all$wild.neg.binom.total.mod.daily.local.N.unlinked.by.doy.isolate
    
    data.model = data.table(median = apply(data.model, 2, median),
                            upper = apply(data.model, 2, quantile, 0.975),
                            lower = apply(data.model, 2, quantile, 0.025))
    
    data.obs = obs.data$daily.local.N
    
    for(t in 1:length(time)){
      
      if(data.obs[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data.obs[t],data.obs[t],0),
                     default.units = 'native',gp=gpar(fill=colObs,col=NA)) 
      }
    }
    
    grid.lines(time,data.model$median,default.units = 'native',gp=gpar(col=colMedian)) 
    grid.polygon(c(time,rev(time)),
                 c(data.model$lower,rev(data.model$upper)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI))
  }
  
  
  # plot missed cases
  if(row == 1 & column == 2){
    
    data.model = inc.all$wild.neg.binom.total.mod.daily.local.M.unlinked + 
                 inc.all$wild.neg.binom.total.mod.daily.local.M.linked
    
    data.model = data.table(median = apply(data.model, 2, median),
                            upper95 =  apply(data.model,2, quantile, 0.975),
                            lower95 = apply(data.model,2, quantile, 0.025),
                            upper50 =  apply(data.model,2, quantile, 0.75),
                            lower50 = apply(data.model,2, quantile, 0.25))
    
    grid.polygon(c(time,rev(time)),
                 c(data.model$lower50,rev(data.model$upper50)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI50))
    # grid.polygon(c(time,rev(time)),
    #              c(data.model$lower95,rev(data.model$upper95)),
    #              default.units = 'native',gp=gpar(col=NA,fill=colCI95))
    grid.lines(time,data.model$median,default.units = 'native',gp=gpar(col=colMedian))
    
  }
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  

  if(row == 1) grid.xaxis(at=c(1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367),
                          label = c('2020 Jan 1', 'Feb 1', 'Mar 1', 'Apr 1', 'May 1', 'Jun 1', 'Jul 1', 
                                    'Aug 1', 'Sep 1', 'Oct 1', 'Nov 1', 'Dec 1', '2021 Jan 1'))
  
  if(row == 1 & column == 1) {grid.yaxis(at=seq(0,150,by=50),label=seq(0,150,by=50))
    grid.text('A',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  if(row == 1 & column == 2) {grid.yaxis(at=seq(0,500,by=100),label=seq(0,500,by=100))
    grid.text('B',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  
  if(row == 1) grid.text('Day',y=unit(-3,'lines'))
  if(column == 1) grid.text('Incidence',x=unit(-3.5,'lines'),rot=90)
  
  grid.text('Partial lockdown',x=unit(12.5,'lines'),y=unit(15,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 1',x=unit(19.5,'lines'),y=unit(16.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 2',x=unit(21.5,'lines'),y=unit(16.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 3',x=unit(44.9,'lines'),y=unit(16.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/incidence_wild_total.png',height=10,width=40,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))

paneller(1,1)
paneller(1,2)

popViewport()
popViewport()
dev.off()

rm(paneller)





# generate plots
paneller=function(row = 1,column=1)
{
  
  xlm=c(457, 596)
  if(row == 1 & column == 1) ylm=c(0,250)  # total notified cases
  if(row == 1 & column == 2) ylm=c(0,3500) # missed cases
  
  innermargins = c(2,3,1,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  time = 457:596
  if(row == 1 & column == 1){colObs = CONFIG$colsLight2[4]; colMedian =  CONFIG$cols[4]; colCI = CONFIG$colsLight3[4]}
  if(row == 1 & column == 2){colMedian =  CONFIG$cols[5]; colCI95 = CONFIG$colsLight3[5]; colCI50 = CONFIG$colsLight2[5]}
  
  # phase 2 heighten alert 16 May to 13 Jun
  grid.polygon(c(502, 502, 530, 530),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
  
  # phase 3 heighten alert 14 Jun to 21 July
  # grid.polygon(c(530.5, 530.5, 568.5, 568.5),
  #              c(0,ylm[2],ylm[2],0),
  #              default.units = 'native',gp=gpar(fill=NA,col=NA))
  
  # phase 2 heighten alert 22 Jul to 18 Aug
  grid.polygon(c(569, 569, 596, 596),
               c(0,ylm[2],ylm[2],0),
               default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
  
  
  # plot all cases
  if(row == 1 & column == 1){
    
    data.model = inc.all$delta.neg.binom.total.mod.daily.local.N.linked.by.doy.isolate + 
                 inc.all$delta.neg.binom.total.mod.daily.local.N.unlinked.by.doy.isolate
    
    data.model = data.table(median = apply(data.model, 2, median),
                            upper = apply(data.model, 2, quantile, 0.975),
                            lower = apply(data.model, 2, quantile, 0.025))
    
    data.obs = obs.data$daily.local.N[time]
    
    for(t in 1:length(time)){
      
      if(data.obs[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data.obs[t],data.obs[t],0),
                     default.units = 'native',gp=gpar(fill=colObs,col=NA)) 
      }
    }
    
    grid.lines(time,data.model$median,default.units = 'native',gp=gpar(col=colMedian)) 
    grid.polygon(c(time,rev(time)),
                 c(data.model$lower,rev(data.model$upper)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI))
  }
  
  
  # plot missed cases
  if(row == 1 & column == 2){
    
    data.model = inc.all$delta.neg.binom.total.mod.daily.local.M.unlinked + 
                 inc.all$delta.neg.binom.total.mod.daily.local.M.linked
    
    data.model = data.table(median = apply(data.model, 2, median),
                            upper95 =  apply(data.model,2, quantile, 0.975),
                            lower95 = apply(data.model,2, quantile, 0.025),
                            upper50 =  apply(data.model,2, quantile, 0.75),
                            lower50 = apply(data.model,2, quantile, 0.25))
    
    grid.polygon(c(time,rev(time)),
                 c(data.model$lower50,rev(data.model$upper50)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI50))
    # grid.polygon(c(time,rev(time)),
    #              c(data.model$lower95,rev(data.model$upper95)),
    #              default.units = 'native',gp=gpar(col=NA,fill=colCI95))
    grid.lines(time,data.model$median,default.units = 'native',gp=gpar(col=colMedian))
    
  }
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  if(row == 1) grid.xaxis(at=c(457, 487, 518, 548, 579),
                          label = c('2021 Apr 1', 'May 1', 'Jun 1', 'Jul 1', 'Aug 1'))
  
  if(row == 1 & column == 1) {grid.yaxis(at=seq(0,250,by=50),label=seq(0,250,by=50))
    grid.text('A',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  if(row == 1 & column == 2) {grid.yaxis(at=seq(0,3500,by=500),label=seq(0,3500,by=500))
    grid.text('B',x=unit(-3,'lines'),y=unit(15.5,'lines'),gp=gpar(fontsize=unit(12,'pt')))}
  
  if(row == 1) grid.text('Day',y=unit(-3,'lines'))
  if(column == 1) grid.text('Incidence',x=unit(-3.5,'lines'),rot=90)
  
  grid.text('Phase 3',x=unit(1,'lines'),y=unit(16.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 2 heighten alert',x=unit(15.3,'lines'),y=unit(13.8,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 3 heighten alert',x=unit(24.5,'lines'),y=unit(13.8,'lines'), gp = gpar(fontsize = 9),rot=90)
  grid.text('Phase 2 heighten alert',x=unit(37.3,'lines'),y=unit(13.8,'lines'), gp = gpar(fontsize = 9),rot=90)
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/incidence_delta_total.png',height=10,width=40,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))

paneller(1,1)
paneller(1,2)

popViewport()
popViewport()
dev.off()

rm(paneller)
