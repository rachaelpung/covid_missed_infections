source('codes/covid_missInf_load_library.R')
source('codes/covid_missInf_plot_colour.R')

load('input/obs_data.RData')
load('input/time_period.RData')

# informative prior
load('output processed/20220521/inc_all.RData')
load('output processed/20220521/theta_all.RData')
inc_inform = copy(inc_all) 
theta_inform = copy(theta_all) 
rm(inc_all, theta_all)

# uninformative prior
load('output processed/20220818/inc_all.RData')
load('output processed/20220818/theta_all.RData')

inc_uninform = copy(inc_all) 
theta_uninform = copy(theta_all) 
rm(inc_all, theta_all)

names(inc_uninform) = names(inc_inform)[13:18]
names(theta_uninform) = names(theta_inform)[11:14]

# wild neg binom link and unlink incidence
paneller=function(row = 1,column=1)
{
  
  if(!(row==3 & column==1)) xlm=c(1, 367)
  if(row==3 & column==1) xlm=c(0.5,5.5)
  
  if(row==1 & column==1) ylm=c(0,100)  # linked cases
  if(row==1 & column==2) ylm=c(0,50)   # unlinked cases
  if(row==2 & column==1) ylm=c(0,300)  # missed cases 50CI
  if(row==2 & column==2) ylm=c(0,2000)  # missed cases 95CI
  if(row==3 & column ==1) ylm=c(0,1)    # prop of linked to all cases vs prop of missed to all infections
  
  innermargins = c(2,3,1.5,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  time = 1:366
  if(row==3 & column==1){variant_fit = 'wild_neg_binom_mod'; extract_start=time_period[variant == 'W']$doy_start; extract_end=time_period[variant == 'W']$doy_end}
  
  # set colour
  if(row==1 & column==1){colMedian =  c('#1C4E1B',CONFIG$cols[6]); colCI = c(CONFIG$colsLight1[3], CONFIG$colsLight1[6])}
  if(row==1 & column==2){colMedian =  c('#A73203',CONFIG$cols[6]); colCI = c(CONFIG$colsLight1[2], CONFIG$colsLight1[6])}
  if(row==2){colMedian =  c('#3F2845',CONFIG$cols[6]); colCI95 = c(CONFIG$colsLight3[5],CONFIG$colsLight3[6]); colCI50 = c(CONFIG$colsLight2[5],CONFIG$colsLight2[6])}
  

  # shaded polygon for 2020
  if(row!=3){
    # circuit breaker from 7 Apr to 1 Jun, c(98, 98, 153.5,153.5)
    # phase 1, c(153.5, 153.5, 170.5,170.5)
    # phase 2, c(170.5, 170.5, 362.5,362.5)
    # phase 3, c(362.5,362.5,367,367)
    
    # based on time periods and not policy dates
    grid.polygon(c(98, 98, 153.5,153.5),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#ECECEC',col=NA))
    
    grid.polygon(c(153.5, 153.5, 170.5,170.5),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#F5F5F5',col=NA))
    
    grid.polygon(c(170.5, 170.5, 194.5,194.5),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
    
    grid.polygon(c(194.5,194.5,367,367),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
  }
  
  # grid lines
  if(row==3){
    grid.lines(xlm, c(0.25,0.25), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.5,0.5),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.75,0.75), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  

  # plot linked cases
  if(row==1 & column==1){
    
    # informative and uninformative prior
    for(p in c(2,1)){
      if(p==1) inc_all=inc_inform
      if(p==2) inc_all=inc_uninform
      
      data_model = data.table(median = apply(inc_all$wild_neg_binom_mod_daily_local_N_linked_by_doy_isolate, 2, median),
                              upper = apply(inc_all$wild_neg_binom_mod_daily_local_N_linked_by_doy_isolate,2, quantile, 0.975),
                              lower = apply(inc_all$wild_neg_binom_mod_daily_local_N_linked_by_doy_isolate,2, quantile, 0.025))
      
      grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian[p])) 
      
      grid.polygon(c(time,rev(time)),
                   c(data_model$lower,rev(data_model$upper)),
                   default.units = 'native',gp=gpar(col=NA,fill=colCI[p]))
      
    }
    
  }
  
  # plot unlinked cases
  if(row==1 & column==2){
    
    for(p in c(2,1)){
      if(p==1) inc_all=inc_inform
      if(p==2) inc_all=inc_uninform
      
      data_model = data.table(median = apply(inc_all$wild_neg_binom_mod_daily_local_N_unlinked_by_doy_isolate, 2, median),
                              upper =  apply(inc_all$wild_neg_binom_mod_daily_local_N_unlinked_by_doy_isolate,2, quantile, 0.975),
                              lower = apply(inc_all$wild_neg_binom_mod_daily_local_N_unlinked_by_doy_isolate,2, quantile, 0.025))
      
      
      
      grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian[p])) 
      grid.polygon(c(time,rev(time)),
                   c(data_model$lower,rev(data_model$upper)),
                   default.units = 'native',gp=gpar(col=NA,fill=colCI[p]))
    }
   
  }
  
  # plot missed cases
  if(row==2){
    
    for(p in c(2,1)){
      if(p==1) inc_all=inc_inform
      if(p==2) inc_all=inc_uninform
      data_model = inc_all$wild_neg_binom_mod_daily_local_M_unlinked + inc_all$wild_neg_binom_mod_daily_local_M_linked
      
      data_model = data.table(median = apply(data_model, 2, median),
                              upper95 =  apply(data_model,2, quantile, 0.975),
                              lower95 = apply(data_model,2, quantile, 0.025),
                              upper50 =  apply(data_model,2, quantile, 0.75),
                              lower50 = apply(data_model,2, quantile, 0.25))
      if(column==1){
        grid.polygon(c(time,rev(time)),
                     c(data_model$lower50,rev(data_model$upper50)),
                     default.units = 'native',gp=gpar(col=NA,fill=colCI50[p]))
      }
    
      if(column==2){
        grid.polygon(c(time,rev(time)),
                     c(data_model$lower95,rev(data_model$upper95)),
                     default.units = 'native',gp=gpar(col=NA,fill=colCI95[p]))
      }

      grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian[p]))
    
    }
  }
  
  # plot proportions
  if(row==3){
    
    grid.text('Wild-type',x=unit(0.025,'npc'),unit(0.95,'npc'),just = 'left')
    if(variant_fit =='wild_neg_binom_mod') period=1:5
    
    for(period_num in period){
      
      # uninformed prior
      inc_all=inc_uninform
      extract_fit = grep(variant_fit, names(inc_all))
      data = inc_all[extract_fit]
      
      for(i in 1:length(data)){data[[i]] = data[[i]][, extract_start[period_num]:extract_end[period_num]]}
      
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
      
      colMedian=CONFIG$cols[6]; colCI50=CONFIG$colsLight2[6]; colCI95=CONFIG$colsLight3[6]
      
      grid.lines(c(period_num+0.125,period_num+0.125), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      grid.lines(c(period_num+0.125,period_num+0.125), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      grid.points(period_num+0.125,ratio$median,default.units = 'native',pch=15,gp=gpar(cex=0.8,col=colMedian))
      
      
      
      # informed prior
      inc_all=inc_inform
      extract_fit = grep(variant_fit, names(inc_all))
      data = inc_all[extract_fit]
      
      for(i in 1:length(data)){data[[i]] = data[[i]][, extract_start[period_num]:extract_end[period_num]]}
      
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
      
      colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]
      
      grid.lines(c(period_num,period_num)-0.125, c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      grid.lines(c(period_num,period_num)-0.125, c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      grid.points(period_num-0.125,ratio$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian))
      
      
      # # observed 
      # unlinked = apply(data[[1]],1,sum)
      # linked   = apply(data[[2]],1,sum)
      # # ratio = unlinked/linked
      # ratio = unlinked/(linked+unlinked)
      # 
      # ratio = data.table(median = median(ratio),
      #                    upper95 = quantile(ratio, probs = 0.975),
      #                    lower95 = quantile(ratio, probs = 0.025),
      #                    upper50 = quantile(ratio, probs = 0.75),
      #                    lower50 = quantile(ratio, probs = 0.25))
      # 
      # # ratio = log(ratio)
      # 
      # colMedian=CONFIG$cols[4]; colCI50=CONFIG$colsLight2[4]; colCI95=CONFIG$colsLight3[5]
      # 
      # grid.lines(c(period_num-0.25,period_num-0.25), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      # grid.lines(c(period_num-0.25,period_num-0.25), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      # grid.points(period_num-0.25,ratio$median,default.units = 'native',pch=17,gp=gpar(cex=0.8,col=colMedian))
      # 
      
    }
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  # plot legend
  space.x=0.15
  if(row %in% c(1,2)){
    
    if(row==2 & column==1) colCI=colCI95
    if(row==2 & column==2) colCI=colCI50
    
    grid.lines(c(0.75,0.775)-space.x, c(0.9575,0.9575), default.units = 'npc',gp=gpar(col=colMedian[1],lwd=2))
    grid.polygon(c(0.75,0.775,0.775,0.75)-space.x, c(0.885,0.885,0.91,0.91), default.units = 'npc',gp=gpar(col=NA,fill=colCI[1]))
    grid.text('Median (inform)',x=unit(0.785-space.x,'npc'),unit(0.9575,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(!(row==2 & column==1))grid.text('95% CI (inform)',x=unit(0.785-space.x,'npc'),unit(0.8975,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(row==2 & column==1)grid.text('50% CI (inform)',x=unit(0.785-space.x,'npc'),unit(0.8975,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    
    grid.lines(c(0.75,0.775)-space.x, c(0.8375,0.8375), default.units = 'npc',gp=gpar(col=colMedian[2],lwd=2))
    grid.polygon(c(0.75,0.775,0.775,0.75)-space.x, c(0.765,0.765,0.79,0.79), default.units = 'npc',gp=gpar(col=NA,fill=colCI[2]))
    grid.text('Median (non-inform)',x=unit(0.785-space.x,'npc'),unit(0.8375,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(!(row==2 & column==1))grid.text('95% CI (non-inform)',x=unit(0.785-space.x,'npc'),unit(0.7775,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(row==2 & column==1)grid.text('50% CI (non-inform)',x=unit(0.785-space.x,'npc'),unit(0.7775,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    if(row==1 & column==1){
      x_bottom_left = 0.60; y_bottom_left = 0.6375; y_length = 0.03; x_length = 0.3
      colour_hex = c('white', '#FFFFFF', '#FAFAFA', '#F5F5F5', '#ECECEC')
      cols <- colorRampPalette(colour_hex)
      val=c(11,8,5,2,0)
      
      x_start = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-126]
      x_end = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-1]
      
      for(i in 1:125){
        grid.polygon(c(x_start[i],x_end[i],x_end[i],x_start[i]),
                     c(y_bottom_left,y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length),
                     default.units = 'npc',gp=gpar(fill=cols(125)[i], col=cols(125)[i]))
      }
      
      x_start <-  seq(x_bottom_left, x_bottom_left+x_length,length.out=5)
      label <- c('>10', '8', '5', '2', '0')
      
      for(breaks in 1:length(x_start)){
        
        if(breaks %in% c(2:(length(x_start)-1))){
          grid.lines(c(x_start[breaks], x_start[breaks]),
                     c(y_bottom_left, y_bottom_left+(y_length/5)),default.units = 'npc',gp=gpar(col='black'))
          grid.lines(c(x_start[breaks], x_start[breaks]),
                     c(y_bottom_left+y_length, y_bottom_left+y_length-(y_length/5)),default.units = 'npc',gp=gpar(col='black'))
        }
        
        grid.text('Visitor restrictions (pax/day)',x=x_bottom_left,y=0.7275,gp = gpar(fontsize = 6), default.units = 'npc',just = 'left')
        
        grid.text(label[breaks],x=x_start[breaks],y=0.69,gp = gpar(fontsize = 6), default.units = 'npc')
        
      }
      
      grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_length,x_bottom_left+x_length),
                   c(y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length,y_bottom_left),
                   default.units = 'npc',gp=gpar(fill=NA, col='black'))
    }
    
    
  }
  
  if(row==3){
    
    colMedian.1=CONFIG$cols[4]; colCI50.1=CONFIG$colsLight2[4]; colCI95.1=CONFIG$colsLight3[4]
    colMedian.2=CONFIG$cols[5]; colCI50.2=CONFIG$colsLight2[5]; colCI95.2=CONFIG$colsLight3[5]
    colMedian.3=CONFIG$cols[6]; colCI50.3=CONFIG$colsLight2[6]; colCI95.3=CONFIG$colsLight3[6]
    
    space.x=0.02
    # grid.lines(c(0.6,0.6)+space.x, c(0.17,0.23), default.units = 'npc',gp=gpar(col=colCI95.1,lwd=2))
    # grid.lines(c(0.6,0.6)+space.x, c(0.185,0.215), default.units = 'npc',gp=gpar(col=colCI50.1,lwd=2))
    # grid.points(0.6+space.x,0.2,default.units = 'npc',pch=15,gp=gpar(cex=0.5,col=colMedian.1))
    
    grid.lines(c(0.6,0.6)+space.x, c(0.095,0.155), default.units = 'npc',gp=gpar(col=colCI95.2,lwd=2))
    grid.lines(c(0.6,0.6)+space.x, c(0.11,0.14), default.units = 'npc',gp=gpar(col=colCI50.2,lwd=2))
    grid.points(0.6+space.x,0.125,default.units = 'npc',pch=16,gp=gpar(cex=0.5,col=colMedian.2))
    
    grid.lines(c(0.6,0.6)+space.x, c(0.02,0.08), default.units = 'npc',gp=gpar(col=colCI95.3,lwd=2))
    grid.lines(c(0.6,0.6)+space.x, c(0.035,0.065), default.units = 'npc',gp=gpar(col=colCI50.3,lwd=2))
    grid.points(0.6+space.x,0.05,default.units = 'npc',pch=17,gp=gpar(cex=0.5,col=colMedian.3))
    
    # grid.text('Unlinked / All cases',x=unit(0.625+space.x,'npc'),unit(0.2,'npc'),just = 'left',gp=gpar(fontsize=unit(6,'pt')))
    
    grid.text('Missed / All infections ',x=unit(0.635+space.x,'npc'),unit(0.14,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('(Inform)',x=unit(0.635+space.x,'npc'),unit(0.11,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    grid.text('Missed / All infections',x=unit(0.635+space.x,'npc'),unit(0.065,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('(Non-inform)',x=unit(0.635+space.x,'npc'),unit(0.035,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    
    
  }
  
  
  # secondary axis  
  if(row!=3){
    # grid.xaxis(at=c(60.5, 97.5, 170.5, 194.5),label=FALSE,main=FALSE)
    
    period=c(60.5, 97.5, 170.5, 194.5)
    for(p in 1:length(period)){
      grid.lines(c((period[p]-xlm[1])/(xlm[2]-xlm[1]),(period[p]-xlm[1])/(xlm[2]-xlm[1])), c(0.975,1),default.units = 'npc')
    }
    
    if(row==1){
      grid.text('Wuhan travellers',x=unit(10,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('Global',x=unit(70,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('outbreak',x=unit(85,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      
      grid.text('Partial lockdown',x=unit(107,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('Gradual reopen',x=unit(180,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('Further reopen',x=unit(204,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('Time periods for wild-type',y=unit(1,'npc')+unit(0.5,'lines'))
    }
  }
  
  if(row==3){
    grid.text('Outbreak metrics',y=unit(1,'npc')+unit(0.5,'lines'))
  }
  
  
  # main axis ticks
  if(row!=3) grid.xaxis(at=c(1, 92, 183, 275, 367), label = c('2020\nJan 1', 'Apr 1', 'Jul 1', 'Oct 1', '2021\nJan 1'))
  
  if(row==3) grid.xaxis(at=seq(1,5,by=1),label=c('2020 Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
                        gp = gpar(fontsize=8))
  
  # axis label
  if(row!=3) grid.text('Incidence',x=unit(-3.5,'lines'),rot=90)
  if(row==3) grid.text('Proportion (%)',x=unit(-3.5,'lines'),rot=90)
  
  # plot labels
  if(row==1 & column==1) {grid.yaxis(at=seq(0,100,by=10),label=seq(0,100,by=10))
    grid.text('A',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row==1 & column==2) {grid.yaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))
    grid.text('B',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row==2 & column==1) {grid.yaxis(at=seq(0,300,by=50),label=seq(0,300,by=50))
    grid.text('C',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))} 
  if(row==2 & column==2) {grid.yaxis(at=seq(0,2000,by=500),label=seq(0,2000,by=500))
    grid.text('D',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row==3 & column==1) {grid.yaxis(at=seq(0,1,0.25),label=seq(0,100,25))
    grid.text('E',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}

png('figure/incidence_wild_noninform_prior_supp.png',height=24,width=20,units='cm',res=300,pointsize=10)
# pdf('figure/incidence_wild_and_delta.pdf',height=9.6,width=8)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=2)))

paneller(1,1)
paneller(1,2)
paneller(2,1)
paneller(2,2)
paneller(3,1)

grid.text('Time',y=unit(-1,'lines'),gp = gpar(fontsize = 13))


popViewport()
popViewport()
dev.off()

rm(paneller)

# wild neg binom link and unlink density plots
paneller=function(row = 1,column=1)
{
  if(row==1) xlm=c(0,5)
  if(row==2) xlm=c(0,2)
  if(row==3) xlm=c(0,1)
  if(row==4) xlm=c(0,1)
  if(row==5) xlm=c(0,5)
  
  if(row == 1 & column == 2) xlm=c(0,200)
  if(row == 1 & column == 4) xlm=c(0,50)
  if(row == 1 & column == 5) xlm=c(0,50)
  
  ylm=c(0,5)
  if(row==1 & column==2) ylm=c(0,0.15)
  if(row==1 & column==4) ylm=c(0,2)
  if(row==2 & column==3) ylm=c(0,30)
  if(row==2 & column==5) ylm=c(0,10)
  if(row==3 & column==1) ylm=c(0,20)
  if(row==3 & column==2) ylm=c(0,10)
  if(row==3 & column==3) ylm=c(0,15)
  if(row==3 & column==5) ylm=c(0,10)
  if(row==4 & column==3) ylm=c(0,30)
  if(row==4 & column==5) ylm=c(0,10)
  if(row==5 & column==2) ylm=c(0,10)
  if(row==5 & column==3) ylm=c(0,10)

  
  innermargins = c(2,2,1,1)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(row==1){param = 1; colLines=c(CONFIG$cols[1],CONFIG$cols[6])}
  if(row==2){param = 2; colLines=c(CONFIG$cols[5],CONFIG$cols[6])}
  if(row==3){param = 3; colLines=c(CONFIG$cols[3],CONFIG$cols[6])}
  if(row==4){param = 4; colLines=c(CONFIG$cols[2],CONFIG$cols[6])}
  if(row==5){param = 5; colLines=c(CONFIG$cols[4],CONFIG$cols[6])}
  
  variant_fit = 'wild_neg_binom_period_'
  period_num=column
  
  extract = paste(variant_fit, period_num, sep ='')
  
  for(p in 1:2){
    
    if(p==1) theta_all = theta_inform
    if(p==2) theta_all = theta_uninform
    
    data = theta_all[[extract]][[param]]
    
    if(row==1 & period_num!=3){
      # convert factor of missed imported cases to average missed imported cases

      mean_import = mean(obs_data$daily_arrival_import_N[time_period[period == period_num & variant == 'W']$doy_start:time_period[period == period_num & variant == 'W']$doy_end])
      data = data * mean_import

    } else if(row==1 & period_num==3){
      data[!is.na(data)]=NA
    }
    
    if(!(row==1 & period_num==3)){
      data = density(data)
      grid.lines(data$x, data$y, default.units = 'native',gp=gpar(col=colLines[p],lwd=2))
    }
    
  }
  
  
  if(column==1){
    space.x=0.7
    grid.lines(c(0.75,0.775)-space.x, c(0.95,0.95), default.units = 'npc',gp=gpar(col=colLines[1],lwd=2))
    grid.text('Inform',x=unit(0.785-space.x,'npc'),unit(0.95,'npc'),just = 'left',gp=gpar(fontsize=unit(10,'pt')))

    grid.lines(c(0.75,0.775)-space.x, c(0.85,0.85), default.units = 'npc',gp=gpar(col=colLines[2],lwd=2))
    grid.text('Non-inform',x=unit(0.785-space.x,'npc'),unit(0.85,'npc'),just = 'left',gp=gpar(fontsize=unit(10,'pt')))

  }
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))

  
  if(row==1 & column%in%c(1,3,5)) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1))
  if(row==1 & column==4) grid.yaxis(at=seq(0,2,0.5),label = seq(0,2,0.5))
  if(row==1 & column==2) grid.yaxis(at=seq(0,0.15,0.05),label = seq(0,0.15,0.05))
  
  if(row==2 & column %in%c(1,2,4)) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1),gp = gpar(fontsize = 12))
  if(row==2 & column==3) grid.yaxis(at=seq(0,30,10),label = seq(0,30,10),gp = gpar(fontsize = 12))
  if(row==2 & column==5) grid.yaxis(at=seq(0,10,2),label = seq(0,10,2),gp = gpar(fontsize = 12))
  
  if(row==3 & column==4) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1))
  if(row==3 & column%in%c(2,5)) grid.yaxis(at=seq(0,10,2),label = seq(0,10,2))
  if(row==3 & column==3) grid.yaxis(at=seq(0,15,5),label = seq(0,15,5))
  if(row==3 & column==1) grid.yaxis(at=seq(0,20,5),label = seq(0,20,5))

  if(row==4 & column%in%c(1,2,4)) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1))
  if(row==4 & column==3) grid.yaxis(at=seq(0,30,10),label = seq(0,30,10))
  if(row==4 & column==5) grid.yaxis(at=seq(0,10,2),label = seq(0,10,2))

  if(row==5 & column%in%c(1,4,5)) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1))
  if(row==5 & column%in%c(2,3)) grid.yaxis(at=seq(0,10,2),label = seq(0,10,2))

  
  if(row == 1 & column %in% c(1,3)) grid.xaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  if(row == 1 & column == 2) grid.xaxis(at=seq(0,200,by=50),label=seq(0,200,by=50))
  if(row == 1 & column == 4) grid.xaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))
  if(row == 1 & column == 5) grid.xaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))
  
  if(row == 2) grid.xaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5))
  if(row == 3) grid.xaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 4) grid.xaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
  if(row == 5) grid.xaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  
  
  if(row == 1 & column == 1) grid.text('Jan 18 - Feb 29', y=unit(9.5,'lines'),gp = gpar(fontsize = 12))
  if(row == 1 & column == 2) grid.text('Mar 1 - Apr 6', y=unit(9.5,'lines'),gp = gpar(fontsize = 12))
  if(row == 1 & column == 3) grid.text('Apr 7 - Jun 18', y=unit(9.5,'lines'),gp = gpar(fontsize = 12))
  if(row == 1 & column == 4) grid.text('Jun 19 - Jul 12', y=unit(9.5,'lines'),gp = gpar(fontsize = 12))
  if(row == 1 & column == 5) grid.text('Jul 13 - Dec 31', y=unit(9.5,'lines'),gp = gpar(fontsize = 12))
  
  if(column == 3 & row == 1) grid.text('Avg daily missed importation',y=unit(-2.2,'lines'))
  if(column == 3 & row == 2) grid.text('R',y=unit(-2.2,'lines'))
  if(column == 3 & row == 3) grid.text(bquote(~epsilon[link]~' (%)'),y=unit(-2.2,'lines')) #~epsilon[op]
  if(column == 3 & row == 4) grid.text(bquote(~epsilon[unlink]~' (%)'),y=unit(-2.2,'lines')) #~epsilon[op*'\'']
  if(column == 3 & row == 5) grid.text('k\'',y=unit(-2.2,'lines'))
  
  if(column == 3 & row == 1) grid.text('Lockdown of borders', y=unit(6,'lines'), gp = gpar(fontsize = 8))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}

png('figure/density_wild_noninform_prior_supp.png',height=30,width=30,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=5,ncol=5)))

for(r in 1:5){
  for(c in 1:5){
    paneller(r,c)

    print(c(r,c))
  }
}


grid.text('Density',x=unit(-1,'lines'), rot=90)


popViewport()
popViewport()
dev.off()

rm(paneller)


# wild neg binom total incidence
paneller=function(row = 1,column=1)
{
  
  if(!(row==2 & column==2)) xlm=c(1, 367)
  if(row==2 & column==2) xlm=c(0.5,5.5)
  
  if(row==1 & column==1) ylm=c(0,150)   # total notified cases
  if(row==1 & column==2) ylm=c(0,500)  # missed cases 50CI
  if(row==2 & column==1) ylm=c(0,10000)  # missed cases 95CI
  if(row==2 & column ==2) ylm=c(0,1)    # prop of linked to all cases vs prop of missed to all infections
  
  innermargins = c(2,3,1.5,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  time = 1:366
  if(row==2 & column==2){variant_fit = 'wild_neg_binom_total_mod'; extract_start=time_period[variant == 'W']$doy_start; extract_end=time_period[variant == 'W']$doy_end}
  
  # set colour
  if(row == 1 & column == 1){colMedian =  c('#C3E1EA',CONFIG$cols[6]); colCI = c(CONFIG$colsLight1[4],CONFIG$colsLight1[6])}
  if((row == 1 & column == 2)|(row == 2 & column == 1)){colMedian =  c(CONFIG$cols[5],CONFIG$cols[6]); colCI95 = c(CONFIG$colsLight3[5],CONFIG$colsLight3[6]); colCI50 = c(CONFIG$colsLight2[5],CONFIG$colsLight2[6])}
  
 
  
  # shaded polygon for 2020
  if(row!=2){
    # circuit breaker from 7 Apr to 1 Jun, c(98, 98, 153.5,153.5)
    # phase 1, c(153.5, 153.5, 170.5,170.5)
    # phase 2, c(170.5, 170.5, 362.5,362.5)
    # phase 3, c(362.5,362.5,367,367)
    
    # based on time periods and not policy dates
    grid.polygon(c(98, 98, 153.5,153.5),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#ECECEC',col=NA))
    
    grid.polygon(c(153.5, 153.5, 170.5,170.5),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#F5F5F5',col=NA))
    
    grid.polygon(c(170.5, 170.5, 194.5,194.5),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
    
    grid.polygon(c(194.5,194.5,367,367),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
  }
  
  # grid lines
  if(row==2){
    grid.lines(xlm, c(0.25,0.25), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.5,0.5),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.75,0.75), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  
  # plot total
  if(row==1 & column==1){
    
    # informative and uninformative prior
    for(p in c(2,1)){
      if(p==1) inc_all=inc_inform
      if(p==2) inc_all=inc_uninform
      
      data_model = inc_all$wild_neg_binom_total_mod_daily_local_N_linked_by_doy_isolate + 
        inc_all$wild_neg_binom_total_mod_daily_local_N_unlinked_by_doy_isolate
      
      data_model = data.table(median = apply(data_model, 2, median),
                              upper = apply(data_model, 2, quantile, 0.975),
                              lower = apply(data_model, 2, quantile, 0.025))
      
      grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian[p])) 
      
      grid.polygon(c(time,rev(time)),
                   c(data_model$lower,rev(data_model$upper)),
                   default.units = 'native',gp=gpar(col=NA,fill=colCI[p]))
      
    }
    
  }
  
  # plot missed cases
  if((row==1 & column==2)|(row==2 & column==1)){
    
    for(p in c(2,1)){
      if(p==1) inc_all=inc_inform
      if(p==2) inc_all=inc_uninform
      data_model = inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked + 
        inc_all$wild_neg_binom_total_mod_daily_local_M_linked
      
      data_model = data.table(median = apply(data_model, 2, median),
                              upper95 =  apply(data_model,2, quantile, 0.975),
                              lower95 = apply(data_model,2, quantile, 0.025),
                              upper50 =  apply(data_model,2, quantile, 0.75),
                              lower50 = apply(data_model,2, quantile, 0.25))
      if(row==1 & column==2){
        grid.polygon(c(time,rev(time)),
                     c(data_model$lower50,rev(data_model$upper50)),
                     default.units = 'native',gp=gpar(col=NA,fill=colCI50[p]))
      }
      
      if(row==2 & column==1){
        grid.polygon(c(time,rev(time)),
                     c(data_model$lower95,rev(data_model$upper95)),
                     default.units = 'native',gp=gpar(col=NA,fill=colCI95[p]))
      }
      
      grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian[p]))
      
    }
  }
  
  # plot proportions
  if(row==2 & column==2){
    
    grid.text('Wild-type',x=unit(0.025,'npc'),unit(0.95,'npc'),just = 'left')
    if(variant_fit =='wild_neg_binom_total_mod') period=1:5
    
    for(period_num in period){
      
      # uninformed prior
      inc_all=inc_uninform
      extract_fit = grep(variant_fit, names(inc_all))
      data = inc_all[extract_fit]
      
      for(i in 1:length(data)){data[[i]] = data[[i]][, extract_start[period_num]:extract_end[period_num]]}
      
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
      
      colMedian=CONFIG$cols[6]; colCI50=CONFIG$colsLight2[6]; colCI95=CONFIG$colsLight3[6]
      
      grid.lines(c(period_num+0.125,period_num+0.125), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      grid.lines(c(period_num+0.125,period_num+0.125), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      grid.points(period_num+0.125,ratio$median,default.units = 'native',pch=15,gp=gpar(cex=0.8,col=colMedian))
      
      
      
      # informed prior
      inc_all=inc_inform
      extract_fit = grep(variant_fit, names(inc_all))
      data = inc_all[extract_fit]
      
      for(i in 1:length(data)){data[[i]] = data[[i]][, extract_start[period_num]:extract_end[period_num]]}
      
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
      
      colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]
      
      grid.lines(c(period_num,period_num)-0.125, c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      grid.lines(c(period_num,period_num)-0.125, c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      grid.points(period_num-0.125,ratio$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian))
      
      
      # observed 
      # unlinked = apply(data[[1]],1,sum)
      # linked   = apply(data[[2]],1,sum)
      # # ratio = unlinked/linked
      # ratio = unlinked/(linked+unlinked)
      # 
      # ratio = data.table(median = median(ratio),
      #                    upper95 = quantile(ratio, probs = 0.975),
      #                    lower95 = quantile(ratio, probs = 0.025),
      #                    upper50 = quantile(ratio, probs = 0.75),
      #                    lower50 = quantile(ratio, probs = 0.25))
      # 
      # # ratio = log(ratio)
      # 
      # colMedian=CONFIG$cols[4]; colCI50=CONFIG$colsLight2[4]; colCI95=CONFIG$colsLight3[5]
      # 
      # grid.lines(c(period_num-0.25,period_num-0.25), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      # grid.lines(c(period_num-0.25,period_num-0.25), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      # grid.points(period_num-0.25,ratio$median,default.units = 'native',pch=17,gp=gpar(cex=0.8,col=colMedian))
      
      
    }
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  # plot legend
  space.x=0.15
  if(!(row==2 & column==2)){
    
    if(row==1 & column==2) colCI=colCI50
    if(row==2 & column==1) colCI=colCI95
    
    grid.lines(c(0.75,0.775)-space.x, c(0.9575,0.9575), default.units = 'npc',gp=gpar(col=colMedian[1],lwd=2))
    grid.polygon(c(0.75,0.775,0.775,0.75)-space.x, c(0.885,0.885,0.91,0.91), default.units = 'npc',gp=gpar(col=NA,fill=colCI[1]))
    grid.text('Median (inform)',x=unit(0.785-space.x,'npc'),unit(0.9575,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(!(row==1 & column==2))grid.text('95% CI (inform)',x=unit(0.785-space.x,'npc'),unit(0.8975,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(row==1 & column==2)grid.text('50% CI (inform)',x=unit(0.785-space.x,'npc'),unit(0.8975,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    
    grid.lines(c(0.75,0.775)-space.x, c(0.8375,0.8375), default.units = 'npc',gp=gpar(col=colMedian[2],lwd=2))
    grid.polygon(c(0.75,0.775,0.775,0.75)-space.x, c(0.765,0.765,0.79,0.79), default.units = 'npc',gp=gpar(col=NA,fill=colCI[2]))
    grid.text('Median (non-inform)',x=unit(0.785-space.x,'npc'),unit(0.8375,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(!(row==1 & column==2))grid.text('95% CI (non-inform)',x=unit(0.785-space.x,'npc'),unit(0.7775,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(row==1 & column==2)grid.text('50% CI (non-inform)',x=unit(0.785-space.x,'npc'),unit(0.7775,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    if(row==1 & column==1){
      x_bottom_left = 0.60; y_bottom_left = 0.6375; y_length = 0.03; x_length = 0.3
      colour_hex = c('white', '#FFFFFF', '#FAFAFA', '#F5F5F5', '#ECECEC')
      cols <- colorRampPalette(colour_hex)
      val=c(11,8,5,2,0)
      
      x_start = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-126]
      x_end = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-1]
      
      for(i in 1:125){
        grid.polygon(c(x_start[i],x_end[i],x_end[i],x_start[i]),
                     c(y_bottom_left,y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length),
                     default.units = 'npc',gp=gpar(fill=cols(125)[i], col=cols(125)[i]))
      }
      
      x_start <-  seq(x_bottom_left, x_bottom_left+x_length,length.out=5)
      label <- c('>10', '8', '5', '2', '0')
      
      for(breaks in 1:length(x_start)){
        
        if(breaks %in% c(2:(length(x_start)-1))){
          grid.lines(c(x_start[breaks], x_start[breaks]),
                     c(y_bottom_left, y_bottom_left+(y_length/5)),default.units = 'npc',gp=gpar(col='black'))
          grid.lines(c(x_start[breaks], x_start[breaks]),
                     c(y_bottom_left+y_length, y_bottom_left+y_length-(y_length/5)),default.units = 'npc',gp=gpar(col='black'))
        }
        
        grid.text('Visitor restrictions (pax/day)',x=x_bottom_left,y=0.7275,gp = gpar(fontsize = 6), default.units = 'npc',just = 'left')
        
        grid.text(label[breaks],x=x_start[breaks],y=0.69,gp = gpar(fontsize = 6), default.units = 'npc')
        
      }
      
      grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_length,x_bottom_left+x_length),
                   c(y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length,y_bottom_left),
                   default.units = 'npc',gp=gpar(fill=NA, col='black'))
      
    }
    
  }
  
  if(row==2 & column==2){
    
    colMedian.1=CONFIG$cols[4]; colCI50.1=CONFIG$colsLight2[4]; colCI95.1=CONFIG$colsLight3[4]
    colMedian.2=CONFIG$cols[5]; colCI50.2=CONFIG$colsLight2[5]; colCI95.2=CONFIG$colsLight3[5]
    colMedian.3=CONFIG$cols[6]; colCI50.3=CONFIG$colsLight2[6]; colCI95.3=CONFIG$colsLight3[6]
    
    space.x=0.02
    # grid.lines(c(0.6,0.6)+space.x, c(0.17,0.23), default.units = 'npc',gp=gpar(col=colCI95.1,lwd=2))
    # grid.lines(c(0.6,0.6)+space.x, c(0.185,0.215), default.units = 'npc',gp=gpar(col=colCI50.1,lwd=2))
    # grid.points(0.6+space.x,0.2,default.units = 'npc',pch=15,gp=gpar(cex=0.5,col=colMedian.1))
    
    grid.lines(c(0.6,0.6)+space.x, c(0.095,0.155), default.units = 'npc',gp=gpar(col=colCI95.2,lwd=2))
    grid.lines(c(0.6,0.6)+space.x, c(0.11,0.14), default.units = 'npc',gp=gpar(col=colCI50.2,lwd=2))
    grid.points(0.6+space.x,0.125,default.units = 'npc',pch=16,gp=gpar(cex=0.5,col=colMedian.2))
    
    grid.lines(c(0.6,0.6)+space.x, c(0.02,0.08), default.units = 'npc',gp=gpar(col=colCI95.3,lwd=2))
    grid.lines(c(0.6,0.6)+space.x, c(0.035,0.065), default.units = 'npc',gp=gpar(col=colCI50.3,lwd=2))
    grid.points(0.6+space.x,0.05,default.units = 'npc',pch=17,gp=gpar(cex=0.5,col=colMedian.3))
    
    # grid.text('Unlinked / All cases',x=unit(0.635+space.x,'npc'),unit(0.2,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    grid.text('Missed / All infections ',x=unit(0.635+space.x,'npc'),unit(0.14,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('(Inform)',x=unit(0.635+space.x,'npc'),unit(0.11,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    grid.text('Missed / All infections',x=unit(0.635+space.x,'npc'),unit(0.065,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('(Non-inform)',x=unit(0.635+space.x,'npc'),unit(0.035,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    
    
  }
  
  
  # secondary axis  
  if(!(row==2 & column==2)){
    # grid.xaxis(at=c(60.5, 97.5, 170.5, 194.5),label=FALSE,main=FALSE)
    
    period=c(60.5, 97.5, 170.5, 194.5)
    for(p in 1:length(period)){
      grid.lines(c((period[p]-xlm[1])/(xlm[2]-xlm[1]),(period[p]-xlm[1])/(xlm[2]-xlm[1])), c(0.975,1),default.units = 'npc')
    }
    
    if(row==1){
      grid.text('Wuhan travellers',x=unit(10,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('Global',x=unit(70,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('outbreak',x=unit(85,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      
      grid.text('Partial lockdown',x=unit(107,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('Gradual reopen',x=unit(180,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('Further reopen',x=unit(204,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('Time periods for wild-type',y=unit(1,'npc')+unit(0.5,'lines'))
    }
  }
  
  if(row==2 & column==2){
    grid.text('Outbreak metrics',y=unit(1,'npc')+unit(0.5,'lines'))
  }
  
  
  # main axis ticks
  if(!(row==2 & column==2)) grid.xaxis(at=c(1, 92, 183, 275, 367), label = c('2020\nJan 1', 'Apr 1', 'Jul 1', 'Oct 1', '2021\nJan 1'))
  
  if(row==2 & column==2) grid.xaxis(at=seq(1,5,by=1),label=c('2020 Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
                        gp = gpar(fontsize=8))
  
  # axis label
  if(!row==2 & column==2) grid.text('Incidence',x=unit(-3.5,'lines'),rot=90)
  if(row==2 & column==2) grid.text('Proportion (%)',x=unit(-3.5,'lines'),rot=90)
  
  # plot labels
  if(row==1 & column==1) {grid.yaxis(at=seq(0,150,by=50),label=seq(0,150,by=50))
    grid.text('A',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row==1 & column==2) {grid.yaxis(at=seq(0,500,by=100),label=seq(0,500,by=100))
    grid.text('B',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row==2 & column==1) {grid.yaxis(at=seq(0,10000,by=2500),label=seq(0,10000,by=2500))
    grid.text('C',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row==2 & column==2) {grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
    grid.text('D',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))} 
 
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}

png('figure/incidence_wild_total_noninform_prior_supp.png',height=16,width=20,units='cm',res=300,pointsize=10)
# pdf('figure/incidence_wild_and_delta.pdf',height=9.6,width=8)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2)))

paneller(1,1)
paneller(1,2)
paneller(2,1)
paneller(2,2)


grid.text('Time',y=unit(-1,'lines'),gp = gpar(fontsize = 13))


popViewport()
popViewport()
dev.off()

rm(paneller)

# wild neg binom total density plots
paneller=function(row = 1,column=1)
{
  
  if(row==1) xlm=c(0,5)
  if(row==2) xlm=c(0,2)
  if(row==3) xlm=c(0,5)
  
  if(row == 1 & column == 2) xlm=c(0,200)
  if(row == 1 & column == 4) xlm=c(0,50)
  if(row == 1 & column == 5) xlm=c(0,50)
  
  ylm=c(0,5)
  if(row==1 & column==2) ylm=c(0,0.15)
  if(row==1 & column==4) ylm=c(0,0.5)
  if(row==1 & column==5) ylm=c(0,2)
  if(row==2 & column==3) ylm=c(0,20)
  if(row==2 & column==5) ylm=c(0,20)
  if(row==3 & column==2) ylm=c(0,10)
  if(row==3 & column==3) ylm=c(0,10)
  
  
  innermargins = c(2,2,1,1)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(row==1){param = 1; colLines=c(CONFIG$cols[4],CONFIG$cols[6])}
  if(row==2){param = 2; colLines=c(CONFIG$cols[5],CONFIG$cols[6])}
  if(row==3){param = 5; colLines=c(CONFIG$cols[1],CONFIG$cols[6])}
  
  variant_fit = 'wild_neg_binom_total_period_'
  period_num=column
  
  extract = paste(variant_fit, period_num, sep ='')
  
  for(p in 1:2){
    
    if(p==1) theta_all = theta_inform
    if(p==2) theta_all = theta_uninform
    
    data = theta_all[[extract]][[param]]
    
    if(row==1 & period_num!=3){
      # convert factor of missed imported cases to average missed imported cases

      mean_import = mean(obs_data$daily_arrival_import_N[time_period[period == period_num & variant == 'W']$doy_start:time_period[period == period_num & variant == 'W']$doy_end])
      data = data * mean_import

    } else if(row==1 & period_num==3){
      data[!is.na(data)]=NA
    }
    
    if(!(row==1 & period_num==3)){
      data = density(data)
      grid.lines(data$x, data$y, default.units = 'native',gp=gpar(col=colLines[p],lwd=2))
    }
    
  }
  
  
  if(column==1){
    space.x=0.7
    grid.lines(c(0.75,0.775)-space.x, c(0.95,0.95), default.units = 'npc',gp=gpar(col=colLines[1],lwd=2))
    grid.text('Inform',x=unit(0.785-space.x,'npc'),unit(0.95,'npc'),just = 'left',gp=gpar(fontsize=unit(10,'pt')))
    
    grid.lines(c(0.75,0.775)-space.x, c(0.85,0.85), default.units = 'npc',gp=gpar(col=colLines[2],lwd=2))
    grid.text('Non-inform',x=unit(0.785-space.x,'npc'),unit(0.85,'npc'),just = 'left',gp=gpar(fontsize=unit(10,'pt')))
    
  }
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  if(row==1 & column%in%c(1,3)) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1))
  if(row==1 & column==4) grid.yaxis(at=seq(0,0.5,0.1),label = seq(0,0.5,0.1))
  if(row==1 & column==2) grid.yaxis(at=seq(0,0.15,0.05),label = seq(0,0.15,0.05))
  if(row==1 & column==5) grid.yaxis(at=seq(0,2,0.5),label = seq(0,2,0.5))
  
  if(row==2 & column %in%c(1,2,4)) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1), gp = gpar(fontsize = 12))
  if(row==2 & column==3) grid.yaxis(at=seq(0,20,5),label = seq(0,20,5), gp = gpar(fontsize = 12))
  if(row==2 & column==5) grid.yaxis(at=seq(0,20,5),label = seq(0,20,5), gp = gpar(fontsize = 12))
  
  if(row==3 & column%in%c(1,4,5)) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1))
  if(row==3 & column%in%c(2,3)) grid.yaxis(at=seq(0,10,2),label = seq(0,10,2))
  
  
  if(row == 1 & column %in% c(1,3)) grid.xaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  if(row == 1 & column == 2) grid.xaxis(at=seq(0,200,by=50),label=seq(0,200,by=50))
  if(row == 1 & column == 4) grid.xaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))
  if(row == 1 & column == 5) grid.xaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))
  
  if(row == 2) grid.xaxis(at=seq(0,2,by=0.5),label=seq(0,2,by=0.5), gp = gpar(fontsize = 12))
  if(row == 3) grid.xaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  
  
  if(row == 1 & column == 1) grid.text('Jan 18 - Feb 29', y=unit(9.5,'lines'), gp = gpar(fontsize = 12))
  if(row == 1 & column == 2) grid.text('Mar 1 - Apr 6', y=unit(9.5,'lines'), gp = gpar(fontsize = 12))
  if(row == 1 & column == 3) grid.text('Apr 7 - Jun 18', y=unit(9.5,'lines'), gp = gpar(fontsize = 12))
  if(row == 1 & column == 4) grid.text('Jun 19 - Jul 12', y=unit(9.5,'lines'), gp = gpar(fontsize = 12))
  if(row == 1 & column == 5) grid.text('Jul 13 - Dec 31', y=unit(9.5,'lines'), gp = gpar(fontsize = 12))
  
  if(column == 3 & row == 1) grid.text('Avg daily missed importation',y=unit(-2.2,'lines'))
  if(column == 3 & row == 2) grid.text('R',y=unit(-2.2,'lines'), gp = gpar(fontsize = 12))
  if(column == 3 & row == 3) grid.text('k\'',y=unit(-2.2,'lines'))
  
  if(column == 3 & row == 1) grid.text('Lockdown of borders', y=unit(6,'lines'), gp = gpar(fontsize = 8))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}

png('figure/density_wild_total_noninformative_prior_supp.png',height=18,width=30,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=5)))

for(r in 1:3){
  for(c in 1:5){
    paneller(r,c)
    
    print(c(r,c))
  }
}


grid.text('Density',x=unit(-1,'lines'), rot=90)


popViewport()
popViewport()
dev.off()

rm(paneller)


# delta neg binom total incidence
paneller=function(row = 1,column=1)
{
  
  if(!(row==2 & column==2)) xlm=c(457, 596)
  if(row==2 & column==2) xlm=c(0.5,4.5)
  
  if(row==1 & column==1) ylm=c(0,750)   # total notified cases
  if(row==1 & column==2) ylm=c(0,20000)  # missed cases 50CI
  if(row==2 & column==1) ylm=c(0,20000)  # missed cases 95CI
  if(row==2 & column ==2) ylm=c(0,1)    # prop of linked to all cases vs prop of missed to all infections
  
  innermargins = c(2,3,1.5,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  time = 457:596
  if(row==2 & column==2){variant_fit = 'delta_neg_binom_total_mod'; extract_start=time_period[variant == 'D']$doy_start - 456; extract_end=time_period[variant == 'D']$doy_end - 456}
  
  # set colour
  if(row == 1 & column == 1){colMedian =  c('#C3E1EA',CONFIG$cols[6]); colCI = c(CONFIG$colsLight1[4],CONFIG$colsLight1[6])}
  if((row == 1 & column == 2)|(row == 2 & column == 1)){colMedian =  c(CONFIG$cols[5],CONFIG$cols[6]); colCI95 = c(CONFIG$colsLight3[5],CONFIG$colsLight3[6]); colCI50 = c(CONFIG$colsLight2[5],CONFIG$colsLight2[6])}
  
  
  
  # shaded polygon for 2021
  if(!(row==2 & column==2)){
    
    # phase 3, c(457,457,502.5,502.5)
    # phase 2 heighten alert 16 May to 13 Jun, c(502.5,502.5, 530.5, 530.5)
    # phase 3 heighten alert 14 Jun to 21 July, c(530.5, 530.5, 568.5, 568.5)
    # phase 2 heighten alert 22 Jul to 18 Aug, c(568.5, 568.5, 596, 596)
    
    # based on time periods and not policy dates
    grid.polygon(c(457,457,498.5,498.5),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#FFFFFF',col=NA))
    
    grid.polygon(c(498.5,498.5, 531.5, 531.5),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#F5F5F5',col=NA))
    
    grid.polygon(c(531.5, 531.5, 564.5, 564.5),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#FAFAFA',col=NA))
    
    grid.polygon(c(564.5, 564.5, 596, 596),
                 c(0,ylm[2],ylm[2],0),
                 default.units = 'native',gp=gpar(fill='#F5F5F5',col=NA))
  }
  
  
  # grid lines
  if(row==2 & column==2){
    grid.lines(xlm, c(0.25,0.25), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.5,0.5),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.75,0.75), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
  
  
  # plot total
  if(row==1 & column==1){
    
    # informative and uninformative prior
    for(p in c(2,1)){
      if(p==1) inc_all=inc_inform
      if(p==2) inc_all=inc_uninform
      
      data_model = inc_all$delta_neg_binom_total_mod_daily_local_N_linked_by_doy_isolate + 
        inc_all$delta_neg_binom_total_mod_daily_local_N_unlinked_by_doy_isolate
      
      data_model = data.table(median = apply(data_model, 2, median),
                              upper = apply(data_model, 2, quantile, 0.975),
                              lower = apply(data_model, 2, quantile, 0.025))
      
    
      
      grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian[p])) 
      
      grid.polygon(c(time,rev(time)),
                   c(data_model$lower,rev(data_model$upper)),
                   default.units = 'native',gp=gpar(col=NA,fill=colCI[p]))
      
    }
    
  }
  
  # plot missed cases
  if((row==1 & column==2)|(row==2 & column==1)){
    
    for(p in c(2,1)){
      if(p==1) inc_all=inc_inform
      if(p==2) inc_all=inc_uninform
      
      data_model = inc_all$delta_neg_binom_total_mod_daily_local_M_unlinked + 
        inc_all$delta_neg_binom_total_mod_daily_local_M_linked
      
      data_model = data.table(median = apply(data_model, 2, median),
                              upper95 =  apply(data_model,2, quantile, 0.975),
                              lower95 = apply(data_model,2, quantile, 0.025),
                              upper50 =  apply(data_model,2, quantile, 0.75),
                              lower50 = apply(data_model,2, quantile, 0.25))
      
    
      
      if(row==1 & column==2){
        
        data.upp=rev(data_model$upper50)
        data.upp[data.upp>ylm[2]]=ylm[2]
        
        grid.polygon(c(time,rev(time)),
                     c(data_model$lower50,data.upp),
                     default.units = 'native',gp=gpar(col=NA,fill=colCI50[p]))
      }
      
      if(row==2 & column==1){
        
        data.upp=rev(data_model$upper95)
        data.upp[data.upp>ylm[2]]=ylm[2]
        
        grid.polygon(c(time,rev(time)),
                     c(data_model$lower95,data.upp),
                     default.units = 'native',gp=gpar(col=NA,fill=colCI95[p]))
        
      }
      
      grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian[p], lwd=1))

      if(p==2){
        grid.lines(c(0.875,0.875),c(0.875,0.975),default.units = 'npc',gp=gpar(col=colMedian[p]))
        grid.lines(c(0.86,0.875),c(0.95,0.975),default.units = 'npc',gp=gpar(col=colMedian[p]))
        grid.lines(c(0.875,0.89),c(0.975,0.95),default.units = 'npc',gp=gpar(col=colMedian[p]))
        grid.text('Extend',x=unit(0.9,'npc'),unit(0.955,'npc'),just = 'left',gp=gpar(fontsize=unit(5,'pt')))
        grid.text('beyond',x=unit(0.9,'npc'),unit(0.925,'npc'),just = 'left',gp=gpar(fontsize=unit(5,'pt')))
        grid.text('limit',x=unit(0.9,'npc'),unit(0.895,'npc'),just = 'left',gp=gpar(fontsize=unit(5,'pt')))
      }
      
    }
    
    
  }
  
  # plot proportions
  if(row==2 & column==2){
    
    grid.text('Delta variant',x=unit(0.025,'npc'),unit(0.95,'npc'),just = 'left')
    if(variant_fit =='delta_neg_binom_total_mod') period=1:4
    
    for(period_num in period){
      
      # uninformed prior
      inc_all=inc_uninform
      extract_fit = grep(variant_fit, names(inc_all))
      data = inc_all[extract_fit]
      
      for(i in 1:length(data)){data[[i]] = data[[i]][, extract_start[period_num]:extract_end[period_num]]}
      
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
      
      colMedian=CONFIG$cols[6]; colCI50=CONFIG$colsLight2[6]; colCI95=CONFIG$colsLight3[6]
      
      grid.lines(c(period_num+0.125,period_num+0.125), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      grid.lines(c(period_num+0.125,period_num+0.125), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      grid.points(period_num+0.125,ratio$median,default.units = 'native',pch=15,gp=gpar(cex=0.8,col=colMedian))
      
      
      
      # informed prior
      inc_all=inc_inform
      extract_fit = grep(variant_fit, names(inc_all))
      data = inc_all[extract_fit]
      
      for(i in 1:length(data)){data[[i]] = data[[i]][, extract_start[period_num]:extract_end[period_num]]}
      
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
      
      colMedian=CONFIG$cols[5]; colCI50=CONFIG$colsLight2[5]; colCI95=CONFIG$colsLight3[5]
      
      grid.lines(c(period_num,period_num)-0.125, c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      grid.lines(c(period_num,period_num)-0.125, c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      grid.points(period_num-0.125,ratio$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian))
      
      
      # observed 
      # unlinked = apply(data[[1]],1,sum)
      # linked   = apply(data[[2]],1,sum)
      # # ratio = unlinked/linked
      # ratio = unlinked/(linked+unlinked)
      # 
      # ratio = data.table(median = median(ratio),
      #                    upper95 = quantile(ratio, probs = 0.975),
      #                    lower95 = quantile(ratio, probs = 0.025),
      #                    upper50 = quantile(ratio, probs = 0.75),
      #                    lower50 = quantile(ratio, probs = 0.25))
      # 
      # # ratio = log(ratio)
      # 
      # colMedian=CONFIG$cols[4]; colCI50=CONFIG$colsLight2[4]; colCI95=CONFIG$colsLight3[5]
      # 
      # grid.lines(c(period_num-0.25,period_num-0.25), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      # grid.lines(c(period_num-0.25,period_num-0.25), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      # grid.points(period_num-0.25,ratio$median,default.units = 'native',pch=17,gp=gpar(cex=0.8,col=colMedian))
      
      
    }
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  # plot legend
  if(!(row==2 & column==2)){
    
    space.x=0.7
    if(row==1 & column==1) space.y=0.5525
    if(row==1 & column==2) space.y=0.7
    if(row==2 & column==1) space.y=0
    
    
    if(row==1 & column==2) colCI=colCI50
    if(row==2 & column==1) colCI=colCI95
    
    grid.lines(c(0.75,0.775)-space.x, c(0.9575,0.9575)-space.y, default.units = 'npc',gp=gpar(col=colMedian[1],lwd=2))
    grid.polygon(c(0.75,0.775,0.775,0.75)-space.x, c(0.885,0.885,0.91,0.91)-space.y, default.units = 'npc',gp=gpar(col=NA,fill=colCI[1]))
    grid.text('Median (inform)',x=unit(0.785-space.x,'npc'),unit(0.9575-space.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(!(row==1 & column==2))grid.text('95% CI (inform)',x=unit(0.785-space.x,'npc'),unit(0.8975-space.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(row==1 & column==2)grid.text('50% CI (inform)',x=unit(0.785-space.x,'npc'),unit(0.8975-space.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    
    grid.lines(c(0.75,0.775)-space.x, c(0.8375,0.8375)-space.y, default.units = 'npc',gp=gpar(col=colMedian[2],lwd=2))
    grid.polygon(c(0.75,0.775,0.775,0.75)-space.x, c(0.765,0.765,0.79,0.79)-space.y, default.units = 'npc',gp=gpar(col=NA,fill=colCI[2]))
    grid.text('Median (non-inform)',x=unit(0.785-space.x,'npc'),unit(0.8375-space.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(!(row==1 & column==2))grid.text('95% CI (non-inform)',x=unit(0.785-space.x,'npc'),unit(0.7775-space.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(row==1 & column==2)grid.text('50% CI (non-inform)',x=unit(0.785-space.x,'npc'),unit(0.7775-space.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    if(row==1 & column==1){
      x_bottom_left = 0.05; y_bottom_left = 0.085; y_length = 0.03; x_length = 0.3
      colour_hex = c('white', '#FFFFFF', '#FAFAFA', '#F5F5F5', '#ECECEC')
      cols <- colorRampPalette(colour_hex)
      val=c(11,8,5,2,0)
      
      x_start = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-126]
      x_end = seq(x_bottom_left, x_bottom_left+x_length,length.out = 126)[-1]
      
      for(i in 1:125){
        grid.polygon(c(x_start[i],x_end[i],x_end[i],x_start[i]),
                     c(y_bottom_left,y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length),
                     default.units = 'npc',gp=gpar(fill=cols(125)[i], col=cols(125)[i]))
      }
      
      x_start <-  seq(x_bottom_left, x_bottom_left+x_length,length.out=5)
      label <- c('>10', '8', '5', '2', '0')
      
      for(breaks in 1:length(x_start)){
        
        if(breaks %in% c(2:(length(x_start)-1))){
          grid.lines(c(x_start[breaks], x_start[breaks]),
                     c(y_bottom_left, y_bottom_left+(y_length/5)),default.units = 'npc',gp=gpar(col='black'))
          grid.lines(c(x_start[breaks], x_start[breaks]),
                     c(y_bottom_left+y_length, y_bottom_left+y_length-(y_length/5)),default.units = 'npc',gp=gpar(col='black'))
        }
        
        grid.text('Visitor restrictions (pax/day)',x=x_bottom_left,y=0.175,gp = gpar(fontsize = 6), default.units = 'npc',just = 'left')
        
        grid.text(label[breaks],x=x_start[breaks],y=0.1375,gp = gpar(fontsize = 6), default.units = 'npc')
        
      }
      
      grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_length,x_bottom_left+x_length),
                   c(y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length,y_bottom_left),
                   default.units = 'npc',gp=gpar(fill=NA, col='black'))
      
    }
    
    
  }
  
  if(row==2 & column==2){
    
    colMedian.1=CONFIG$cols[4]; colCI50.1=CONFIG$colsLight2[4]; colCI95.1=CONFIG$colsLight3[4]
    colMedian.2=CONFIG$cols[5]; colCI50.2=CONFIG$colsLight2[5]; colCI95.2=CONFIG$colsLight3[5]
    colMedian.3=CONFIG$cols[6]; colCI50.3=CONFIG$colsLight2[6]; colCI95.3=CONFIG$colsLight3[6]
    
    space.x=0.02
    # grid.lines(c(0.6,0.6)+space.x, c(0.17,0.23), default.units = 'npc',gp=gpar(col=colCI95.1,lwd=2))
    # grid.lines(c(0.6,0.6)+space.x, c(0.185,0.215), default.units = 'npc',gp=gpar(col=colCI50.1,lwd=2))
    # grid.points(0.6+space.x,0.2,default.units = 'npc',pch=15,gp=gpar(cex=0.5,col=colMedian.1))
    
    grid.lines(c(0.6,0.6)+space.x, c(0.095,0.155), default.units = 'npc',gp=gpar(col=colCI95.2,lwd=2))
    grid.lines(c(0.6,0.6)+space.x, c(0.11,0.14), default.units = 'npc',gp=gpar(col=colCI50.2,lwd=2))
    grid.points(0.6+space.x,0.125,default.units = 'npc',pch=16,gp=gpar(cex=0.5,col=colMedian.2))
    
    grid.lines(c(0.6,0.6)+space.x, c(0.02,0.08), default.units = 'npc',gp=gpar(col=colCI95.3,lwd=2))
    grid.lines(c(0.6,0.6)+space.x, c(0.035,0.065), default.units = 'npc',gp=gpar(col=colCI50.3,lwd=2))
    grid.points(0.6+space.x,0.05,default.units = 'npc',pch=17,gp=gpar(cex=0.5,col=colMedian.3))
    
    # grid.text('Unlinked / All cases',x=unit(0.625+space.x,'npc'),unit(0.2,'npc'),just = 'left',gp=gpar(fontsize=unit(6,'pt')))
    
    grid.text('Missed / All infections ',x=unit(0.635+space.x,'npc'),unit(0.14,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('(Inform)',x=unit(0.635+space.x,'npc'),unit(0.11,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    grid.text('Missed / All infections',x=unit(0.635+space.x,'npc'),unit(0.065,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('(Non-inform)',x=unit(0.635+space.x,'npc'),unit(0.035,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    
    
  }
  
  
  # secondary axis  
  if(!(row==2 & column==2)){
    # grid.xaxis(at=c(501.5, 530.5, 558.5),label=FALSE,main=FALSE)
    period=c(498.5, 547.5, 564.5)
    for(p in 1:length(period)){
      grid.lines(c((period[p]-xlm[1])/(xlm[2]-xlm[1]),(period[p]-xlm[1])/(xlm[2]-xlm[1])), c(0.975,1),default.units = 'npc')
    }
    
    
    if(row==1){
      grid.text('India travellers',x=unit(460,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      grid.text('Tighten measures',x=unit(502,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      
      grid.text('Relax measures',x=unit(535,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      
      grid.text('Night club & fishery outbreak',x=unit(551,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      
      grid.text('Tighten measures',x=unit(568,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(10,'pt')), rot=90)
      
      grid.text('Time periods for Delta variant',y=unit(1,'npc')+unit(0.5,'lines'))
    }
  }
  
  if(row==2 & column==2){
    grid.text('Outbreak metrics',y=unit(1,'npc')+unit(0.5,'lines'))
  }
  
  
  # main axis ticks
  if(!(row==2 & column==2)) grid.xaxis(at=c(457, 487, 518, 548, 579), label = c('2021\nApr 1', 'May 1', 'Jun 1', 'Jul 1', 'Aug 1'))
  
  if(row==2 & column==2) grid.xaxis(at=seq(1,4,by=1),label=c('2021 Apr 1-\nMay 12', 'May 13-\nJun 30', 'Ju1 1-\nJul 17', 
                                                             'Jul 18-\nAug 18'),
                                    gp = gpar(fontsize=8))
  
  # axis label
  if(!row==2 & column==2) grid.text('Incidence',x=unit(-3.5,'lines'),rot=90)
  if(row==2 & column==2) grid.text('Proportion (%)',x=unit(-3.5,'lines'),rot=90)
  
  # plot labels
  if(row==1 & column==1) {grid.yaxis(at=seq(0,750,by=250),label=seq(0,750,by=250))
    grid.text('A',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row==1 & column==2) {grid.yaxis(at=seq(0,20000,by=5000),label=seq(0,20000,by=5000))
    grid.text('B',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row==2 & column==1) {grid.yaxis(at=seq(0,20000,by=5000),label=seq(0,20000,by=5000))
    grid.text('C',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row==2 & column==2) {grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25))
    grid.text('D',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))} 
  
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}

png('figure/incidence_delta_total_noninformative_prior_supp.png',height=16,width=20,units='cm',res=300,pointsize=10)
# pdf('figure/incidence_wild_and_delta.pdf',height=9.6,width=8)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2)))

paneller(1,1)
paneller(1,2)
paneller(2,1)
paneller(2,2)


grid.text('Time',y=unit(-1,'lines'),gp = gpar(fontsize = 13))


popViewport()
popViewport()
dev.off()

rm(paneller)

# delta neg binom total density plots
paneller=function(row = 1,column=1)
{
  
  if(row==1) xlm=c(0,50)
  if(row==2 & column!=3) xlm=c(0,5)
  if(row==2 & column==3) xlm=c(0,10)
  if(row==3) xlm=c(0,5)
  
  ylm=c(0,5)
  if(row==2 & column==2) ylm=c(0,15)
  if(row==2 & column==3) ylm=c(0,1)
  if(row==2 & column==4) ylm=c(0,10)
  if(row==3 & column==4) ylm=c(0,10)

  innermargins = c(2,2,1,1)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  if(row==1){param = 1; colLines=c(CONFIG$cols[4],CONFIG$cols[6])}
  if(row==2){param = 2; colLines=c(CONFIG$cols[5],CONFIG$cols[6])}
  if(row==3){param = 5; colLines=c(CONFIG$cols[1],CONFIG$cols[6])}
  
  variant_fit = 'delta_neg_binom_total_period_'
  period_num=column
  
  extract = paste(variant_fit, period_num, sep ='')
  
  for(p in 1:2){
    
    if(p==1) theta_all = theta_inform
    if(p==2) theta_all = theta_uninform
    
    data = theta_all[[extract]][[param]]
    
    if(row==1 & period_num==1){
      # convert factor of missed imported cases to average missed imported cases
      
      mean_import = mean(obs_data$daily_arrival_import_N[time_period[period == period_num & variant == 'W']$doy_start:time_period[period == period_num & variant == 'W']$doy_end])
      data = data * mean_import
      
    } else if(row==1 & period_num!=1){
      data[!is.na(data)]=NA
    }
    
    if(!(row==1 & period_num!=1)){
      data = density(data)
      grid.lines(data$x, data$y, default.units = 'native',gp=gpar(col=colLines[p],lwd=2))
    }
    
  }
  
  
  if(column==1){
    space.x=0.7
    grid.lines(c(0.75,0.775)-space.x, c(0.95,0.95), default.units = 'npc',gp=gpar(col=colLines[1],lwd=2))
    grid.text('Inform',x=unit(0.785-space.x,'npc'),unit(0.95,'npc'),just = 'left',gp=gpar(fontsize=unit(10,'pt')))
    
    grid.lines(c(0.75,0.775)-space.x, c(0.85,0.85), default.units = 'npc',gp=gpar(col=colLines[2],lwd=2))
    grid.text('Non-inform',x=unit(0.785-space.x,'npc'),unit(0.85,'npc'),just = 'left',gp=gpar(fontsize=unit(10,'pt')))
    
  }
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  if(row==1) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1))
  if(row==2 & column==1) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1), gp = gpar(fontsize = 12))
  if(row==2 & column==2) grid.yaxis(at=seq(0,15,5),label = seq(0,15,5), gp = gpar(fontsize = 12))
  if(row==2 & column==3) grid.yaxis(at=seq(0,1,0.25),label = seq(0,1,0.25), gp = gpar(fontsize = 11))
  if(row==2 & column==4) grid.yaxis(at=seq(0,10,2),label = seq(0,10,2), gp = gpar(fontsize = 12))
  if(row==3 & column!=4) grid.yaxis(at=seq(0,5,1),label = seq(0,5,1))
  if(row==3 & column==4) grid.yaxis(at=seq(0,10,2),label = seq(0,10,2))
  
  if(row==1) grid.xaxis(at=seq(0,50,by=10),label=seq(0,50,by=10))
  if(row==2 & column!=3) grid.xaxis(at=seq(0,5,by=1),label=seq(0,5,by=1), gp = gpar(fontsize = 12))
  if(row==2 & column==3) grid.xaxis(at=seq(0,10,by=2),label=seq(0,10,by=2), gp = gpar(fontsize = 12))
  if(row==3) grid.xaxis(at=seq(0,5,by=1),label=seq(0,5,by=1))
  
  
  if(row==1 & column==1) grid.text('Apr 1 - May 12', y=unit(9.5,'lines'), gp = gpar(fontsize = 12))
  if(row==1 & column==2) grid.text('May 13 - Jun 30', y=unit(9.5,'lines'), gp = gpar(fontsize = 12))
  if(row==1 & column==3) grid.text('Jul 1 - Jul 17', y=unit(9.5,'lines'), gp = gpar(fontsize = 12))
  if(row==1 & column==4) grid.text('Jul 18 - Aug 18', y=unit(9.5,'lines'), gp = gpar(fontsize = 12))
  
  if(column==3 & row==1) grid.text('Avg daily missed importation',y=unit(-2.2,'lines'),x=unit(-2.2,'lines'))
  if(column==3 & row==2) grid.text('R',y=unit(-2.2,'lines'),x=unit(-2.2,'lines'))
  if(column==3 & row==3) grid.text('k\'',y=unit(-2.2,'lines'),x=unit(-2.2,'lines'))
  
  if(column!=1 & row == 1) grid.text('Assume no missed\nimported cases due\nto strict border controls', y=unit(6,'lines'), gp = gpar(fontsize = 8))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}

png('figure/density_delta_total_noninformative_prior_supp.png',height=18,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=4)))

for(r in 1:3){
  for(c in 1:4){
    paneller(r,c)
    
    print(c(r,c))
  }
}


grid.text('Density',x=unit(-1,'lines'), rot=90)


popViewport()
popViewport()
dev.off()

rm(paneller)



# zoomed in on time period 3 delta surge
paneller=function(row = 1,column=1)
{
  
  xlm=c(548, 596)
  ylm=c(0,500)  # total notified cases
  
  innermargins = c(2,3,1.5,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  time = 548:596
  
  colObs = '#C3E1EA'; colMedian =  c('#006273',CONFIG$cols[6]); colCI = c(CONFIG$colsLight1[4],CONFIG$colsLight1[6])

  
  # plot all cases
  data_obs = obs_data$daily_local_N[time]
    
  for(t in 1:length(time)){
      
    if(data_obs[t]>0)
    {
      grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                   c(0,data_obs[t],data_obs[t],0),
                   default.units = 'native',gp=gpar(fill=colObs,col=NA)) 
    }
  }
    
  
  for(p in c(2,1)){
    if(p==1) inc_all=inc_inform
    if(p==2) inc_all=inc_uninform

  
    data_model = inc_all$delta_neg_binom_total_mod_daily_local_N_linked_by_doy_isolate + 
      inc_all$delta_neg_binom_total_mod_daily_local_N_unlinked_by_doy_isolate
    
    data_model = data.table(median = apply(data_model, 2, median),
                            upper95 =  apply(data_model,2, quantile, 0.975),
                            lower95 = apply(data_model,2, quantile, 0.025),
                            upper50 =  apply(data_model,2, quantile, 0.75),
                            lower50 = apply(data_model,2, quantile, 0.25))
    
    grid.lines(time,data_model$median[92:140],default.units = 'native',gp=gpar(col=colMedian[p], lwd=1))
    
   
    
  }
  
  
  

  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  # plot legend
  shift.x=0.225; shift.y= 0

  grid.polygon(c(0.75,0.775,0.775,0.75)-shift.x, c(0.945,0.945,0.97,0.97)-shift.y, default.units = 'npc',gp=gpar(col=NA,fill=colObs))

  grid.text('Obs total',x=unit(0.785-shift.x,'npc'),unit(0.9575-shift.y,'npc'),just = 'left',gp=gpar(fontsize=unit(6,'pt')))


  grid.lines(c(0.75,0.775)-shift.x, c(0.8975,0.8975)-shift.y, default.units = 'npc',gp=gpar(col=colMedian[1],lwd=1))
  # grid.polygon(c(0.75,0.775,0.775,0.75)-shift.x, c(0.825,0.825,0.85,0.85)-shift.y, default.units = 'npc',gp=gpar(col=NA,fill=colCI))
  grid.text('Median (inform)',x=unit(0.785-shift.x,'npc'),unit(0.8975-shift.y,'npc'),just = 'left',gp=gpar(fontsize=unit(6,'pt')))
    # grid.text('95% CI',x=unit(0.785-shift.x,'npc'),unit(0.8375-shift.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))

  grid.lines(c(0.75,0.775)-shift.x, c(0.8375,0.8375)-shift.y, default.units = 'npc',gp=gpar(col=colMedian[2],lwd=1))
  # grid.polygon(c(0.75,0.775,0.775,0.75)-shift.x, c(0.825,0.825,0.85,0.85)-shift.y, default.units = 'npc',gp=gpar(col=NA,fill=colCI))
  grid.text('Median (non-inform)',x=unit(0.785-shift.x,'npc'),unit(0.8375-shift.y,'npc'),just = 'left',gp=gpar(fontsize=unit(6,'pt')))
  # grid.text('95% CI',x=unit(0.785-shift.x,'npc'),unit(0.8375-shift.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))

  
  
  # period=c(498.5, 547.5, 564.5)
  # for(p in 1:length(period)){
  #   grid.lines(c((period[p]-xlm[1])/(xlm[2]-xlm[1]),(period[p]-xlm[1])/(xlm[2]-xlm[1])), c(0.975,1),default.units = 'npc')
  # }
  # 
  # grid.text('India travellers',x=unit(460,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(5,'pt')), rot=90)
  # grid.text('Tighten measures',x=unit(502,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(5,'pt')), rot=90)
  # 
  # grid.text('Relax measures',x=unit(535,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(5,'pt')), rot=90)
  # 
  # grid.text('Night club & fishery outbreak',x=unit(551,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(5,'pt')), rot=90)
  # 
  # grid.text('Tighten measures',x=unit(568,'native'),unit(0.95,'npc')+unit(0.5,'lines'),just = 'right',gp=gpar(fontsize=unit(5,'pt')), rot=90)
  # 
  # grid.text('Time periods for Delta variant',y=unit(1,'npc')+unit(0.5,'lines'),gp=gpar(fontsize=unit(5,'pt')))
  
  # main axis ticks
  grid.xaxis(at=c(548, 579), label = c('2021\nJul 1', 'Aug 1'),gp=gpar(fontsize=unit(7,'pt')))
  
  
  
  # axis label
  grid.text('Incidence',x=unit(-3.5,'lines'),rot=90,gp=gpar(fontsize=unit(7,'pt')))
  
  # plot labels
  grid.yaxis(at=seq(0,500,by=100),label=seq(0,500,by=100),gp=gpar(fontsize=unit(7,'pt')))
  # grid.text('A',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))
  
  
  
  
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}

png('figure/test.png',height=8,width=8,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=1)))

paneller()

popViewport()
popViewport()
dev.off()

