source('codes/covid_missInf_load_library.R')
source('codes/covid_missInf_plot_colour.R')

load('input/obs_data.RData')
load('input/time_period.RData')
load('output processed/20220521/theta_all.RData')

param_time_period = time_period[variant=='W']

# trace plot for wild type
paneller=function(row = 1,column=1)
{
  niter = xmax = 5800
  xlm=c(0,xmax)
  if(row == 1) ylm=c(0,5)
  if(row == 2) ylm=c(0,5)
  if(row == 3) ylm=c(0,1)
  if(row == 4) ylm=c(0,1)
  if(row == 5) ylm=c(0,5)
  
  if(row == 1 & column == 2) ylm=c(0,200)
  if(row == 1 & column == 4) ylm=c(0,50)
  if(row == 1 & column == 5) ylm=c(0,100)
  
  
  innermargins = c(2,3,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  
  iteration = seq_len(niter)
  chain.colour = c(CONFIG$colsLight2[4],CONFIG$colsLight2[3],CONFIG$colsLight2[1],CONFIG$colsLight2[5])
  
  for(chain_num in 1:4){
    data = theta_all[[column]][chain==chain_num][[row]]
    
    # convert to average daily missed imports
    if(row == 1){
      data = data * mean(obs_data$daily_arrival_import_N[param_time_period$doy_start[column]:param_time_period$doy_end[column]])
    }
    
    grid.lines(iteration, data, default.units = 'native',gp=gpar(col=chain.colour[chain_num]))
  }
  

  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.xaxis(at=c(0,2500,5000),label = c(0,2500,5000),gp = gpar(fontsize = 13))
 
  
  if(row == 1 & column %in% c(1,3)) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1),gp = gpar(fontsize = 13))
  if(row == 1 & column == 2) grid.yaxis(at=seq(0,200,by=50),label=seq(0,200,by=50),gp = gpar(fontsize = 13))
  if(row == 1 & column == 4) grid.yaxis(at=seq(0,50,by=10),label=seq(0,50,by=10),gp = gpar(fontsize = 13))
  if(row == 1 & column == 5) grid.yaxis(at=seq(0,100,by=25),label=seq(0,100,by=25),gp = gpar(fontsize = 13))

  if(row == 2) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1),gp = gpar(fontsize = 13))

  
  if(row == 3) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25),gp = gpar(fontsize = 13))
  if(row == 4) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25),gp = gpar(fontsize = 13))
  if(row == 5) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1),gp = gpar(fontsize = 13))
  
  if(row == 1 & column == 1) grid.text('Jan 18 - Feb 29', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 2) grid.text('Mar 1 - Apr 6', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 3) grid.text('Apr 7 - Jun 18', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 4) grid.text('Jun 19 - Jul 12', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 5) grid.text('Jul 13 - Dec 31', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  
  if(column == 1 & row == 1) grid.text('Avg daily missed\nimportation',x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13))
  if(column == 1 & row == 2) grid.text('R',x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13))
  if(column == 1 & row == 3) grid.text(bquote(~epsilon[link]~' (%)'),x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13)) #~epsilon[op]
  if(column == 1 & row == 4) grid.text(bquote(~epsilon[unlink]~' (%)'),x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13)) #~epsilon[op*'\'']
  if(column == 1 & row == 5) grid.text('k\'',x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13))
  
  if(column == 3 & row == 1) grid.text('Lockdown of borders', y=unit(3.5,'lines'), gp = gpar(fontsize = 10))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/trace_wild_supp.png',height=25,width=30,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=5,ncol=5)))

for(r in 1:5){
  for(c in 1:5){
    paneller(r,c)

    print(c(r,c))
  }
}

# paneller(1,1)
# paneller(1,3)

grid.text('iteration',y=unit(-3,'lines'))


popViewport()
popViewport()
dev.off()

rm(paneller)





param_time_period = time_period[variant=='D']

# trace plot for delta
paneller=function(row = 1,column=1)
{
  niter = xmax = 5800
  xlm=c(0,xmax)
  if(row == 1) ylm=c(0,50)
  if(row == 2) ylm=c(0,5)
  # if(row == 3) ylm=c(0,1)
  # if(row == 4) ylm=c(0,1)
  if(row == 3) ylm=c(0,5)
  
  if(row == 2 & column ==3) ylm=c(0,7)
  
  innermargins = c(2,3,2,2)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  
  iteration = seq_len(niter)
  chain.colour = c(CONFIG$colsLight2[4],CONFIG$colsLight2[3],CONFIG$colsLight2[1],CONFIG$colsLight2[5])
  
  if(row==1) extract_row=1
  if(row==2) extract_row=2
  if(row==3) extract_row=5
  
  for(chain_num in 1:4){
    data = theta_all[[column+10]][chain==chain_num][[extract_row]]
    
    # convert to average daily missed imports
    if(row == 1){
      data = data * mean(obs_data$daily_arrival_import_N[param_time_period$doy_start[column]:param_time_period$doy_end[column]])
    }
    
    grid.lines(iteration, data, default.units = 'native',gp=gpar(col=chain.colour[chain_num]))
  }
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  grid.xaxis(at=c(0,2500,5000),label = c(0,2500,5000),gp = gpar(fontsize = 13))
  
  
  if(row == 1) grid.yaxis(at=seq(0,50,by=10),label=seq(0,50,by=10),gp = gpar(fontsize = 13))
  if(row == 2 & column != 3) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1),gp = gpar(fontsize = 13))
  if(row == 2 & column == 3) grid.yaxis(at=seq(0,7,by=1),label=seq(0,7,by=1),gp = gpar(fontsize = 13))
 
  # if(row == 3) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25),gp = gpar(fontsize = 13))
  # if(row == 4) grid.yaxis(at=seq(0,1,by=0.25),label=seq(0,100,by=25),gp = gpar(fontsize = 13))
  if(row == 3) grid.yaxis(at=seq(0,5,by=1),label=seq(0,5,by=1),gp = gpar(fontsize = 13))
  
  if(row == 1 & column == 1) grid.text('Apr 1 - May 12', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 2) grid.text('May 13 - Jun 30', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 3) grid.text('Jul 1 - Jul 17', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 4) grid.text('Jul 18 - Aug 18', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  # if(row == 1 & column == 5) grid.text('Jul 1 - Aug 18', y=unit(6,'lines'),gp = gpar(fontsize = 13))
  
  if(column == 1 & row == 1) grid.text('Avg daily missed\nimportation',x=unit(-3,'lines'),rot=90,gp = gpar(fontsize = 13))
  if(column == 1 & row == 2) grid.text('R',x=unit(-3,'lines'),rot=90,gp = gpar(fontsize = 13))
  # if(column == 1 & row == 3) grid.text(bquote(~epsilon[link]~' (%)'),x=unit(-3,'lines'),rot=90,gp = gpar(fontsize = 13)) #~epsilon[op]
  # if(column == 1 & row == 4) grid.text(bquote(~epsilon[unlink]~' (%)'),x=unit(-3,'lines'),rot=90,gp = gpar(fontsize = 13)) #~epsilon[op*'\'']
  if(column == 1 & row == 3) grid.text('k\'',x=unit(-3,'lines'),rot=90,gp = gpar(fontsize = 13))
  
  if(column != 1 & row == 1) grid.text('Assume no missed\nimported cases due to\nstrict border controls', y=unit(3.5,'lines'), gp = gpar(fontsize = 10))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('figure/trace_delta_supp.png',height=15,width=24,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,2,1)))
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=4)))

for(r in 1:3){
  for(c in 1:4){
    paneller(r,c)
    
    print(c(r,c))
  }
}

grid.text('iteration',y=unit(-3,'lines'))


popViewport()
popViewport()
dev.off()

rm(paneller)
