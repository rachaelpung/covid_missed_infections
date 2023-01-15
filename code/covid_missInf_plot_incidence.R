source('codes/covid_missInf_load_library.R')
source('codes/covid_missInf_plot_colour.R')


load('input/obs_data.RData')
load('input/time_period.RData')
load('output processed/20220521/inc_all.RData')

# generate incidence plot for wild type fitted using linked and unlinked cases 
# and for delta using total cases and proportion of cases/infection
paneller=function(row = 1,column=1)
{
  
  if(column==1) xlm=c(1, 367)
  if(column==2 & row %in% c(1,2)) xlm=c(457, 596)
  if(column==2 & row == 3) xlm=c(0.5,5.5)
  if(column==2 & row == 4) xlm=c(0.5,4.5)
  
  if(row == 1 & column == 1) ylm=c(0,100)  # imported cases
  if(row == 2 & column == 1) ylm=c(0,100)  # linked cases
  if(row == 3 & column == 1) ylm=c(0,50)   # unlinked cases
  if(row == 4 & column == 1) ylm=c(0,300)  # missed cases
  if(row == 1 & column == 2) ylm=c(0,300)  # total notified cases
  if(row == 2 & column == 2) ylm=c(0,8000) # missed cases
  if(row %in% c(3,4) & column == 2) ylm=c(0,1)
  
  innermargins = c(2,3,1.5,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(column == 1) time = 1:366
  if(column == 2) time = 457:596
  
  if(row == 3 & column == 2){variant_fit = 'wild_neg_binom_mod';        extract_start=time_period[variant == 'W']$doy_start; extract_end=time_period[variant == 'W']$doy_end}
  if(row == 4 & column == 2){variant_fit = 'delta_neg_binom_total_mod'; extract_start=time_period[variant == 'D']$doy_start - 456; extract_end=time_period[variant == 'D']$doy_end - 456}
  
  
  if(row == 1 & column == 1){colObs_all =  CONFIG$colsLight2[1]; colObs_no_shn =  CONFIG$cols[1]}
  if(row == 2 & column == 1){colObs = '#CEE9CD'; colMedian =  '#1C4E1B'; colCI = CONFIG$colsLight1[3]}
  if(row == 3 & column == 1){colObs = '#FCE8DF'; colMedian =  '#A73203'; colCI = CONFIG$colsLight1[2]}
  if(row == 4 & column == 1){colMedian = '#3F2845'; colCI95 = CONFIG$colsLight3[5]; colCI50 = CONFIG$colsLight2[5]}
  if(row == 1 & column == 2){colObs = '#C3E1EA'; colMedian =  '#006273'; colCI = CONFIG$colsLight1[4]}
  if(row == 2 & column == 2){colMedian =  CONFIG$cols[5]; colCI95 = CONFIG$colsLight3[5]; colCI50 = CONFIG$colsLight2[5]}
  
  # shaded polygon for 2020
  if(column==1){
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
  
  # shaded polygon for 2021
  if(column==2 & row %in% c(1,2)){
    
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
  if(row %in% c(3,4) & column == 2){
    grid.lines(xlm, c(0.25,0.25), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.5,0.5),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
    grid.lines(xlm, c(0.75,0.75), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  }
 

 
  # plot imported cases
  if(row == 1 & column == 1){
    
    data_obs_all = obs_data$daily_arrival_import_N
    data_obs_no_shn = obs_data$daily_arrival_import_N_no_shn
    
    for(t in 1:length(time)){
      
      if(data_obs_all[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data_obs_all[t],data_obs_all[t],0),
                     default.units = 'native',gp=gpar(fill=colObs_all,col=NA)) 
      }
      
      if(data_obs_no_shn[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data_obs_no_shn[t],data_obs_no_shn[t],0),
                     default.units = 'native',gp=gpar(fill=colObs_no_shn,col=NA))
      }
      
    }
  }
  
  # plot linked cases
  if(row == 2 & column == 1){
    
    data_model = data.table(median = apply(inc_all$wild_neg_binom_mod_daily_local_N_linked_by_doy_isolate, 2, median),
                            upper = apply(inc_all$wild_neg_binom_mod_daily_local_N_linked_by_doy_isolate,2, quantile, 0.975),
                            lower = apply(inc_all$wild_neg_binom_mod_daily_local_N_linked_by_doy_isolate,2, quantile, 0.025))
    
    data_obs = obs_data$daily_local_N_linked
    
    for(t in 1:length(time)){
      
      if(data_obs[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data_obs[t],data_obs[t],0),
                     default.units = 'native',gp=gpar(fill=colObs,col=NA)) 
      }
    }
    
     
   
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower,rev(data_model$upper)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI))
    
    grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian, lwd=1.5))
  }
  
  # plot unlinked cases
  if(row == 3 & column == 1){
    
    data_model = data.table(median = apply(inc_all$wild_neg_binom_mod_daily_local_N_unlinked_by_doy_isolate, 2, median),
                            upper =  apply(inc_all$wild_neg_binom_mod_daily_local_N_unlinked_by_doy_isolate,2, quantile, 0.975),
                            lower = apply(inc_all$wild_neg_binom_mod_daily_local_N_unlinked_by_doy_isolate,2, quantile, 0.025))
    
    data_obs = obs_data$daily_local_N_unlinked
    
    for(t in 1:length(time)){
      
      if(data_obs[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data_obs[t],data_obs[t],0),
                     default.units = 'native',gp=gpar(fill=colObs,col=NA)) 
      }
    }
    
    grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian, lwd=1.5)) 
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower,rev(data_model$upper)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI))
  }
  
  # plot missed cases
  if(row == 4 & column == 1){
    
    data_model = inc_all$wild_neg_binom_mod_daily_local_M_unlinked + inc_all$wild_neg_binom_mod_daily_local_M_linked
    
    data_model = data.table(median = apply(data_model, 2, median),
                            upper95 =  apply(data_model,2, quantile, 0.975),
                            lower95 = apply(data_model,2, quantile, 0.025),
                            upper50 =  apply(data_model,2, quantile, 0.75),
                            lower50 = apply(data_model,2, quantile, 0.25))
    
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower50,rev(data_model$upper50)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI50))
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower95,rev(data_model$upper95)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI95))
    grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian, lwd=1.5))
    
  }
  
  
  # plot all cases
  if(row == 1 & column == 2){
    
    data_model = inc_all$delta_neg_binom_total_mod_daily_local_N_linked_by_doy_isolate + 
      inc_all$delta_neg_binom_total_mod_daily_local_N_unlinked_by_doy_isolate
    
    data_model = data.table(median = apply(data_model, 2, median),
                            upper = apply(data_model, 2, quantile, 0.975),
                            lower = apply(data_model, 2, quantile, 0.025))
    
    data_obs = obs_data$daily_local_N[time]
    
    for(t in 1:length(time)){
      
      if(data_obs[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data_obs[t],data_obs[t],0),
                     default.units = 'native',gp=gpar(fill=colObs,col=NA)) 
      }
    }
    
    grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian, lwd=1.5)) 
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower,rev(data_model$upper)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI))
  }
  
  
  # plot missed cases
  if(row == 2 & column == 2){
    
    data_model = inc_all$delta_neg_binom_total_mod_daily_local_M_unlinked + 
      inc_all$delta_neg_binom_total_mod_daily_local_M_linked
    
    data_model = data.table(median = apply(data_model, 2, median),
                            upper95 =  apply(data_model,2, quantile, 0.975),
                            lower95 = apply(data_model,2, quantile, 0.025),
                            upper50 =  apply(data_model,2, quantile, 0.75),
                            lower50 = apply(data_model,2, quantile, 0.25))
    
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower50,rev(data_model$upper50)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI50))
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower95,rev(data_model$upper95)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI95))
    grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian, lwd=1.5))
    
  }
  
  if(row %in% c(3,4) & column == 2){
    
    if(row==3) grid.text('Wild-type',x=unit(0.025,'npc'),unit(0.95,'npc'),just = 'left')
    if(row==4) grid.text('Delta variant',x=unit(0.025,'npc'),unit(0.05,'npc'),just = 'left')
 
    if(variant_fit =='wild_neg_binom_mod') period=1:5
    if(variant_fit =='delta_neg_binom_total_mod') period=1:4
    
    for(period_num in period){
      
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
      
      colMedian=CONFIG$cols[4]; colCI50=CONFIG$colsLight2[4]; colCI95=CONFIG$colsLight3[4]
      
      grid.lines(c(period_num-0.2,period_num-0.2), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      grid.lines(c(period_num-0.2,period_num-0.2), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      grid.points(period_num-0.2,ratio$median,default.units = 'native',pch=16,gp=gpar(cex=0.8,col=colMedian))
      
      
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
      
      grid.lines(c(period_num+0.2,period_num+0.2), c(ratio$lower95,ratio$upper95), default.units = 'native',gp=gpar(col=colCI95,lwd=2))
      grid.lines(c(period_num+0.2,period_num+0.2), c(ratio$lower50,ratio$upper50), default.units = 'native',gp=gpar(col=colCI50,lwd=2))
      grid.points(period_num+0.2,ratio$median,default.units = 'native',pch=17,gp=gpar(cex=0.8,col=colMedian))
      
      
    }
    
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  # plot legend
  if(row == 1 & column == 1){
    grid.polygon(c(0.60,0.625,0.625,0.60), c(0.945,0.945,0.97,0.97), default.units = 'npc',gp=gpar(col=NA,fill=colObs_no_shn))
    grid.polygon(c(0.60,0.625,0.625,0.60), c(0.885,0.885,0.91,0.91), default.units = 'npc',gp=gpar(col=NA,fill=colObs_all))
    
    grid.text('Isolated after positive',x=unit(0.635,'npc'),unit(0.9575,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('Qurantined upon arrival',x=unit(0.635,'npc'),unit(0.8975,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    x_bottom_left = 0.60; y_bottom_left = 0.7475-0.01; y_length = 0.03; x_length = 0.3
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
      
      grid.text('Visitor restrictions (pax/day)',x=x_bottom_left,y=0.8375-0.01,gp = gpar(fontsize = 6), default.units = 'npc',just = 'left')
      
      grid.text(label[breaks],x=x_start[breaks],y=0.80-0.01,gp = gpar(fontsize = 6), default.units = 'npc')
      
    }
    
    grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_length,x_bottom_left+x_length),
                 c(y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length,y_bottom_left),
                 default.units = 'npc',gp=gpar(fill=NA, col='black'))
    
    
  }
  
  if(row %in% c(2,3) & column == 1){
    
    grid.polygon(c(0.75,0.775,0.775,0.75), c(0.945,0.945,0.97,0.97), default.units = 'npc',gp=gpar(col=NA,fill=colObs))
    
    if(row==2) grid.text('Obs linked',x=unit(0.785,'npc'),unit(0.9575,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    if(row==3) grid.text('Obs unlinked',x=unit(0.785,'npc'),unit(0.9575,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    grid.lines(c(0.75,0.775), c(0.8975,0.8975), default.units = 'npc',gp=gpar(col=colMedian,lwd=2))
    grid.polygon(c(0.75,0.775,0.775,0.75), c(0.825,0.825,0.85,0.85), default.units = 'npc',gp=gpar(col=NA,fill=colCI))
    
    
    grid.text('Median',x=unit(0.785,'npc'),unit(0.8975,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('95% CI',x=unit(0.785,'npc'),unit(0.8375,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
  }
  
  if(row==1 & column ==2){
    
    
    shift.x=0.735; shift.y=0.4
    
    grid.polygon(c(0.75,0.775,0.775,0.75)-shift.x, c(0.945,0.945,0.97,0.97)-shift.y, default.units = 'npc',gp=gpar(col=NA,fill=colObs))
    
    if(row==1 & column == 2) grid.text('Obs total',x=unit(0.785-shift.x,'npc'),unit(0.9575-shift.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    
    grid.lines(c(0.75,0.775)-shift.x, c(0.8975,0.8975)-shift.y, default.units = 'npc',gp=gpar(col=colMedian,lwd=2))
    grid.polygon(c(0.75,0.775,0.775,0.75)-shift.x, c(0.825,0.825,0.85,0.85)-shift.y, default.units = 'npc',gp=gpar(col=NA,fill=colCI))
    
    
    grid.text('Median',x=unit(0.785-shift.x,'npc'),unit(0.8975-shift.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('95% CI',x=unit(0.785-shift.x,'npc'),unit(0.8375-shift.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
  }
  
  if((row == 4 & column == 1) | (row == 2 & column == 2)){
    
    grid.lines(c(0.75,0.775), c(0.9575,0.9575), default.units = 'npc',gp=gpar(col=colMedian,lwd=2))
    grid.text('Median',x=unit(0.785,'npc'),unit(0.9575,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
     
    grid.polygon(c(0.75,0.775,0.775,0.75), c(0.885,0.885,0.91,0.91), default.units = 'npc',gp=gpar(col=NA,fill=colCI50))
    grid.polygon(c(0.75,0.775,0.775,0.75), c(0.825,0.825,0.85,0.85), default.units = 'npc',gp=gpar(col=NA,fill=colCI95))
    
    grid.text('50% CI',x=unit(0.785,'npc'),unit(0.8975,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('95% CI',x=unit(0.785,'npc'),unit(0.8375,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
  }
  
 
  if((row == 3 & column == 2) | (row == 4 & column == 2)){
    
    colMedian.1=CONFIG$cols[4]; colCI50.1=CONFIG$colsLight2[4]; colCI95.1=CONFIG$colsLight3[4]
    colMedian.2=CONFIG$cols[5]; colCI50.2=CONFIG$colsLight2[5]; colCI95.2=CONFIG$colsLight3[5]
    
    grid.lines(c(0.6,0.6), c(0.13,0.23), default.units = 'npc',gp=gpar(col=colCI95.1,lwd=2))
    grid.lines(c(0.6,0.6), c(0.155,0.205), default.units = 'npc',gp=gpar(col=colCI50.1,lwd=2))
    grid.points(0.6,0.18,default.units = 'npc',pch=16,gp=gpar(cex=0.6,col=colMedian.1))
    
    grid.lines(c(0.6,0.6), c(0.015,0.115), default.units = 'npc',gp=gpar(col=colCI95.2,lwd=2))
    grid.lines(c(0.6,0.6), c(0.04,0.09), default.units = 'npc',gp=gpar(col=colCI50.2,lwd=2))
    grid.points(0.6,0.065,default.units = 'npc',pch=17,gp=gpar(cex=0.6,col=colMedian.2))
    
    
    grid.text('Missed / All infections',x=unit(0.625,'npc'),unit(0.18,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('Unlinked / All cases',x=unit(0.625,'npc'),unit(0.065,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
   
  }
  
 
  # secondary axis  
  if(column == 1){
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
  
  if(column == 2 & row %in% c(1,2)){
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
  
  if(column==2 & row==3){
    grid.text('Outbreak metrics',y=unit(1,'npc')+unit(0.5,'lines'))
  }


  # main axis ticks
  if(column == 1) grid.xaxis(at=c(1, 92, 183, 275, 367), label = c('2020\nJan 1', 'Apr 1', 'Jul 1', 'Oct 1', '2021\nJan 1'))
  if(column == 2 & row%in%c(1,2)) grid.xaxis(at=c(457, 487, 518, 548, 579), label = c('2021\nApr 1', 'May 1', 'Jun 1', 'Jul 1', 'Aug 1'))
  
  if(row == 3 & column == 2) grid.xaxis(at=seq(1,5,by=1),label=c('2020 Jan 18-\nFeb 29', 'Mar 1-\nApr 6', 'Apr 7-\nJun 18', 'Jun 19-\nJul 12', 'Jul 13-\nDec 31'),
                             gp = gpar(fontsize=8))
  if(row == 4 & column == 2) grid.xaxis(at=seq(1,4,by=1),label=c('2021 Apr 1-\nMay 12', 'May 13-\nJun 30', 'Ju1 1-\nJul 17', 
                                                                 'Jul 18-\nAug 18'),
                             gp = gpar(fontsize=8))
  
  # axis label
  if(column == 1 | (column == 2 & row%in%c(1,2))) grid.text('Incidence',x=unit(-3.5,'lines'),rot=90)
  if(row %in% c(3,4) & column == 2) grid.text('Proportion (%)',x=unit(-3.5,'lines'),rot=90)
  
  # plot labels
  if(row == 1 & column == 1) {grid.yaxis(at=seq(0,100,by=10),label=seq(0,100,by=10))
    grid.text('A',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row == 2 & column == 1) {grid.yaxis(at=seq(0,100,by=10),label=seq(0,100,by=10))
    grid.text('B',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row == 3 & column == 1) {grid.yaxis(at=seq(0,50,by=5),label=seq(0,50,by=5))
    grid.text('C',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))} 
  if(row == 4 & column == 1) {grid.yaxis(at=seq(0,300,by=50),label=seq(0,300,by=50))
    grid.text('D',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row == 1 & column == 2) {grid.yaxis(at=seq(0,300,by=50),label=seq(0,300,by=50))
    grid.text('E',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row == 2 & column == 2) {grid.yaxis(at=seq(0,8000,by=2000),label=seq(0,8000,by=2000))
    grid.text('F',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row == 3 &column == 2) {grid.yaxis(at=seq(0,1,0.25),label=seq(0,100,25))
    grid.text('G',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row == 4 &column == 2) {grid.yaxis(at=seq(0,1,0.25),label=seq(0,100,25))
    grid.text('H',x=unit(-2.7,'lines'),y=unit(11,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  
 
 
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}



png('figure/incidence_fig_1.png',height=32,width=20,units='cm',res=300,pointsize=10)
# pdf('figure/incidence_eLife_fig_1.pdf',height=9.6,width=8)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=4,ncol=2)))

# paneller(1,1)
# paneller(2,1)
# paneller(3,1)
# paneller(4,1)
# 
# paneller(1,2)
# paneller(2,2)
# paneller(3,2)
# paneller(4,2)

for(c in 1:2){
  for(r in 1:4){
    paneller(r,c)
  }
}

grid.text('Time',y=unit(-1,'lines'),gp = gpar(fontsize = 13))


popViewport()
popViewport()
dev.off()

rm(paneller)




# generate incidence plot for wild type fitted using total cases
paneller=function(row = 1,column=1)
{
  
  xlm=c(1, 367)
  if(row == 1 & column == 1) ylm=c(0,150)  # total notified cases
  if(row == 1 & column == 2) ylm=c(0,1500) # missed cases
  
  innermargins = c(2,3,1.5,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  time = 1:366
  if(row == 1 & column == 1){colObs = '#C3E1EA'; colMedian =  '#006273'; colCI = CONFIG$colsLight1[4]}
  if(row == 1 & column == 2){colMedian =  CONFIG$cols[5]; colCI95 = CONFIG$colsLight3[5]; colCI50 = CONFIG$colsLight2[5]}
  
  # shaded polygon for 2020
  
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
  
  
  
  # plot all cases
  if(row == 1 & column == 1){
    
    data_model = inc_all$wild_neg_binom_total_mod_daily_local_N_linked_by_doy_isolate + 
      inc_all$wild_neg_binom_total_mod_daily_local_N_unlinked_by_doy_isolate
    
    data_model = data.table(median = apply(data_model, 2, median),
                            upper = apply(data_model, 2, quantile, 0.975),
                            lower = apply(data_model, 2, quantile, 0.025))
    
    data_obs = obs_data$daily_local_N
    
    for(t in 1:length(time)){
      
      if(data_obs[t]>0)
      {
        grid.polygon(time[t]+0.5*c(-1,-1,1,1),
                     c(0,data_obs[t],data_obs[t],0),
                     default.units = 'native',gp=gpar(fill=colObs,col=NA)) 
      }
    }
    
    grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian, lwd=1.5)) 
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower,rev(data_model$upper)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI))
  }
  
  
  # plot missed cases
  if(row == 1 & column == 2){
    
    data_model = inc_all$wild_neg_binom_total_mod_daily_local_M_unlinked + 
      inc_all$wild_neg_binom_total_mod_daily_local_M_linked
    
    data_model = data.table(median = apply(data_model, 2, median),
                            upper95 =  apply(data_model,2, quantile, 0.975),
                            lower95 = apply(data_model,2, quantile, 0.025),
                            upper50 =  apply(data_model,2, quantile, 0.75),
                            lower50 = apply(data_model,2, quantile, 0.25))
    
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower50,rev(data_model$upper50)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI50))
    grid.polygon(c(time,rev(time)),
                 c(data_model$lower95,rev(data_model$upper95)),
                 default.units = 'native',gp=gpar(col=NA,fill=colCI95))
    grid.lines(time,data_model$median,default.units = 'native',gp=gpar(col=colMedian, lwd=1.5))
    
  }
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  
  # plot legend
  if(column == 1){
    space.y=0.13
    space.x=0.15
    
    grid.polygon(c(0.75,0.775,0.775,0.75)-space.x, c(0.945,0.945,0.97,0.97)-space.y, default.units = 'npc',gp=gpar(col=NA,fill=colObs))
    
    grid.text('Obs',x=unit(0.785-space.x,'npc'),unit(0.9575-space.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    grid.lines(c(0.75,0.775)-space.x, c(0.8975,0.8975)-space.y, default.units = 'npc',gp=gpar(col=colMedian,lwd=2))
    grid.polygon(c(0.75,0.775,0.775,0.75)-space.x, c(0.825,0.825,0.85,0.85)-space.y, default.units = 'npc',gp=gpar(col=NA,fill=colCI))
    
    grid.text('Median',x=unit(0.785-space.x,'npc'),unit(0.8975-space.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('95% CI',x=unit(0.785-space.x,'npc'),unit(0.8375-space.y,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    x_bottom_left = 0.60; y_bottom_left = 0.7475-0.01+space.y; y_length = 0.03; x_length = 0.3
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
      
      grid.text('Visitor restrictions (pax/day)',x=x_bottom_left,y=0.8375-0.01+space.y,gp = gpar(fontsize = 6), default.units = 'npc',just = 'left')
      
      grid.text(label[breaks],x=x_start[breaks],y=0.80-0.01+space.y,gp = gpar(fontsize = 6), default.units = 'npc')
      
    }
    
    grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_length,x_bottom_left+x_length),
                 c(y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length,y_bottom_left),
                 default.units = 'npc',gp=gpar(fill=NA, col='black'))
    
  }
  
  if(column == 2){
    
    grid.lines(c(0.75,0.775), c(0.9575,0.9575), default.units = 'npc',gp=gpar(col=colMedian,lwd=2))
    grid.text('Median',x=unit(0.785,'npc'),unit(0.9575,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
    grid.polygon(c(0.75,0.775,0.775,0.75), c(0.885,0.885,0.91,0.91), default.units = 'npc',gp=gpar(col=NA,fill=colCI50))
    grid.polygon(c(0.75,0.775,0.775,0.75), c(0.825,0.825,0.85,0.85), default.units = 'npc',gp=gpar(col=NA,fill=colCI95))
    
    grid.text('50% CI',x=unit(0.785,'npc'),unit(0.8975,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    grid.text('95% CI',x=unit(0.785,'npc'),unit(0.8375,'npc'),just = 'left',gp=gpar(fontsize=unit(7,'pt')))
    
  }
  
  # secondary axis  
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
  
  
  
  grid.xaxis(at=c(1, 92, 183, 275, 367), label = c('2020\nJan 1', 'Apr 1', 'Jul 1', 'Oct 1', '2021\nJan 1'))
  
  if(row == 1 & column == 1) {grid.yaxis(at=seq(0,150,by=50),label=seq(0,150,by=50))
    grid.text('A',x=unit(-2.5,'lines'),y=unit(9,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  if(row == 1 & column == 2) {grid.yaxis(at=seq(0,1500,by=500),label=seq(0,1500,by=500))
    grid.text('B',x=unit(-2.5,'lines'),y=unit(9,'lines'),gp=gpar(fontsize=unit(13,'pt')))}
  
  # grid.text('Partial lockdown',x=unit(12.5,'lines'),y=unit(15,'lines'), gp = gpar(fontsize = 9),rot=90)
  # grid.text('Phase 1',x=unit(19.5,'lines'),y=unit(16.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  # grid.text('Phase 2',x=unit(21.5,'lines'),y=unit(16.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  # grid.text('Phase 3',x=unit(44.9,'lines'),y=unit(16.5,'lines'), gp = gpar(fontsize = 9),rot=90)
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/incidence_supp.png',height=8,width=20,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))

paneller(1,1)
paneller(1,2)

grid.text('Time',y=unit(-1,'lines'))
grid.text('Incidence',x=unit(-1,'lines'),rot=90)

popViewport()
popViewport()
dev.off()

rm(paneller)




