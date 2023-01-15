source('codes/covid_missInf_load_library.R')
source('codes/covid_missInf_plot_colour.R')

load('input/obs_data.RData')
load('input/time_period.RData')
load('output processed/20220521/theta_all.RData')

decimalnumcount<-function(x){stopifnot(class(x)=="character")
  x<-gsub("(.*)(\\.)|([0]*$)","",x)
  nchar(x)
} 

paneller=function(row = 1,column=1)
{
  # x vs y
  if(row==1){xlm=c(0,1);ylm=c(0,10)}   # e_cf vs rho
  if(row==2){xlm=c(0,1);ylm=c(0,10)}   # e_ct vs rho
  if(row==3){xlm=c(0,2.5);ylm=c(0,10)} # R vs rho
  if(row==4){xlm=c(0,1);ylm=c(0,2.5)} # e_cf vs R
  if(row==5){xlm=c(0,1);ylm=c(0,2.5)} # e_ct vs R
  if(row==6){xlm=c(0,1);ylm=c(0,1)}   # e_cf vs e_ct
  
  innermargins = c(3,4,3,5)
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm,clip=TRUE))
  
  variant_fit = 'wild_neg_binom_period_'
  # param 1=rho; param 2=R; param 3=e_ct; param 4=e_cf
  if(row==1) {param_x = 4; param_y = 1}
  if(row==2) {param_x = 3; param_y = 1}
  if(row==3) {param_x = 2; param_y = 1}
  if(row==4) {param_x = 4; param_y = 2}
  if(row==5) {param_x = 3; param_y = 2}
  if(row==6) {param_x = 4; param_y = 3}
  
  # col_pts_hue = c('#440154FF', '#443A83FF', '#31688EFF', '#21908CFF', '#35B779FF', '#8FD744FF', '#FDE725FF')
  # col_pts_hue = c('#440154FF', '#481567FF', '#482677FF', '#453781FF', '#404788FF',
  #                 '#39568CFF', '#33638DFF', '#2D708EFF', '#287D8EFF', '#238A8DFF',
  #                 '#1F968BFF', '#20A387FF', '#29AF7FFF', '#3CBB75FF', '#55C667FF',
  #                 '#73D055FF', '#95D840FF', '#B8DE29FF', '#DCE319FF', '#FDE725FF')
  col_pts_hue = c('#a6bddb','#0570b0')
  
  
  # grid lines
  # grid.lines(xlm, c(0.25,0.25), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  # grid.lines(xlm, c(0.5,0.5),   default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  # grid.lines(xlm, c(0.75,0.75), default.units = 'native',gp=gpar(col='#F0F0F0',lwd=1))
  
  
  
  if(!(row%in%c(1:3) & column==3)){
    extract = paste(variant_fit, column, sep ='')
    
    x=theta_all[[extract]][[param_x]]
    y=theta_all[[extract]][[param_y]]
    
    xmin=xmax=seq(xlm[1],xlm[2],length.out = 50)
    ymin=ymax=seq(ylm[1],ylm[2],length.out = 50)
    ymin=ymin[1:49];ymax=ymax[2:50]
    
    data=data.table(xmin=xmin[1:49], xmax=xmin[2:50])
    data=data[rep(seq_len(nrow(data)), each=nrow(data))]
    data[,`:=` (ymin=rep(ymin, length.out=nrow(data)), 
                ymax=rep(ymax, length.out=nrow(data)))]
    
    count=sapply(1:nrow(data), function(r){
      
      index=which(x>=data[r,xmin] & x<data[r,xmax])
      length(which(y[index]>=data[r,ymin] & y[index]<data[r,ymax]))
      
    })
    
    data[,count:=count]
    data[,density:=count/length(x)]
    data=data[count!=0]
    
    max_val=signif(max(data$density),1)
    if(max_val<max(data$density)) {
      
      d=decimalnumcount(as.character(max_val))
      max_val=max_val+1*10^(-d)
    }
    
    val=seq(0,max_val,length.out = 2)
   
    pal=gradient_n_pal(colours = col_pts_hue, values = val)
    data[, hex:=pal(density)]
    data[density<0.1*max_val, hex:=lightup(hex,density/(0.1*max_val))]
    
    for(r in 1:nrow(data)){
      
      grid.polygon(c(data[r,xmin],data[r,xmin],data[r,xmax],data[r,xmax]),
                   c(data[r,ymin],data[r,ymax],data[r,ymax],data[r,ymin]),
                   default.units = 'native',gp=gpar(fill=data[r,hex],col=NA))
    }
    
    
    # grid.points(data.x,data.y,default.units = 'native',pch=16,gp=gpar(cex=0.5,col=colMedian)) 
  }
  
  
  
  popViewport()
  pushViewport(plotViewport(innermargins,xscale=xlm,yscale=ylm))
  
  if(!(row%in%c(1:3) & column==3))colourbar(x_bottom_left = 1.09, y_bottom_left = 0.15, y_length = 0.7, x_length = 0.05, palette=col_pts_hue, row=row, data=data, max_val=max_val)
  
  # grid.yaxis(at=seq(0,5,1),label=seq(0,5,1))
  # grid.xaxis(at=seq(0,1,0.25),label=seq(0,1,0.25))
  
  
  if(row==1){
    grid.yaxis(at=seq(0,10,2),label=seq(0,10,2),gp = gpar(fontsize = 13))
    grid.xaxis(at=seq(0,1,0.25),label=seq(0,1,0.25),gp = gpar(fontsize = 13))
    grid.text(bquote(~rho),x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13))
    grid.text(bquote(~epsilon[cf] ~ ' (%)'),y=unit(-2.5,'lines'),gp = gpar(fontsize = 13))
    if(column==3) grid.text('Lockdown of borders',gp = gpar(fontsize = 13))
  }
  if(row==2){
    grid.yaxis(at=seq(0,10,2),label=seq(0,10,2),gp = gpar(fontsize = 13))
    grid.xaxis(at=seq(0,1,0.25),label=seq(0,1,0.25),gp = gpar(fontsize = 13))
    grid.text(bquote(~rho),x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13))
    grid.text(bquote(~epsilon[ct] ~ ' (%)'),y=unit(-2.5,'lines'),gp = gpar(fontsize = 13))
    if(column==3) grid.text('Lockdown of borders',gp = gpar(fontsize = 13))
  }
  if(row==3){
    grid.yaxis(at=seq(0,10,2),label=seq(0,10,2),gp = gpar(fontsize = 13))
    grid.xaxis(at=seq(0,2.5,0.5),label=seq(0,2.5,0.5),gp = gpar(fontsize = 13))
    grid.text(bquote(~rho),x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13))
    grid.text('R',y=unit(-2.5,'lines'),gp = gpar(fontsize = 13))
    if(column==3) grid.text('Lockdown of borders',gp = gpar(fontsize = 13))
  }
  if(row==4){
    grid.yaxis(at=seq(0,2.5,0.5),label=seq(0,2.5,0.5),gp = gpar(fontsize = 13))
    grid.xaxis(at=seq(0,1,0.25),label=seq(0,1,0.25),gp = gpar(fontsize = 13))
    grid.text('R',x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13))
    grid.text(bquote(~epsilon[cf] ~ ' (%)'),y=unit(-2.5,'lines'),gp = gpar(fontsize = 13))
  }
  if(row==5){
    grid.yaxis(at=seq(0,2.5,0.5),label=seq(0,2.5,0.5),gp = gpar(fontsize = 13))
    grid.xaxis(at=seq(0,1,0.25),label=seq(0,1,0.25),gp = gpar(fontsize = 13))
    grid.text('R',x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13))
    grid.text(bquote(~epsilon[ct] ~ ' (%)'),y=unit(-2.5,'lines'),gp = gpar(fontsize = 13))
  }
  if(row==6){
    grid.yaxis(at=seq(0,1,0.25),label=seq(0,1,0.25),gp = gpar(fontsize = 13))
    grid.xaxis(at=seq(0,1,0.25),label=seq(0,1,0.25),gp = gpar(fontsize = 13))
    grid.text(bquote(~epsilon[cf] ~ ' (%)'),x=unit(-2.5,'lines'),rot=90,gp = gpar(fontsize = 13))
    grid.text(bquote(~epsilon[ct] ~ ' (%)'),y=unit(-2.5,'lines'),gp = gpar(fontsize = 13))
  }
  
  # if(row == 1 & column ==1){grid.text('A',x=unit(-3,'lines'),y=unit(13,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  # if(row == 1 & column ==2){grid.text('B',x=unit(-3,'lines'),y=unit(13,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  # if(row == 1 & column ==3){grid.text('C',x=unit(-3,'lines'),y=unit(13,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  # if(row == 1 & column ==4){grid.text('D',x=unit(-3,'lines'),y=unit(13,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  # if(row == 1 & column ==5){grid.text('E',x=unit(-3,'lines'),y=unit(13,'lines'),gp=gpar(fontsize=unit(10,'pt')))}
  
  if(row == 1 & column == 1) grid.text('Jan 18 - Feb 29', y=unit(11,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 2) grid.text('Mar 1 - Apr 6', y=unit(11,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 3) grid.text('Apr 7 - Jun 18', y=unit(11,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 4) grid.text('Jun 19 - Jul 12', y=unit(11,'lines'),gp = gpar(fontsize = 13))
  if(row == 1 & column == 5) grid.text('Jul 13 - Dec 31', y=unit(11,'lines'),gp = gpar(fontsize = 13))
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  
  popViewport()
  popViewport()
  
}

colourbar <- function(x_bottom_left = 1.1, y_bottom_left = 0.4, y_length = 0.2, x_length = 0.25, palette=NULL, row=NULL, data=NULL, max_val=NULL){
  
  
  colour_hex = palette
  cols <- colorRampPalette(colour_hex)
  
  y_start = seq(y_bottom_left, y_bottom_left+y_length,length.out = 126)[-126]
  y_end = seq(y_bottom_left, y_bottom_left+y_length,length.out = 126)[-1]
  
  for(i in 1:125){
    

      grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_length,x_bottom_left+x_length),
                   c(y_start[i],y_end[i],y_end[i],y_start[i]),
                   default.units = 'npc',gp=gpar(fill=cols(125)[i], col=cols(125)[i]))

    
    
    
  }
  
  
  y_start <-  seq(y_bottom_left, y_bottom_left+y_length,length.out=3)
  label <- c('0', max_val/2, max_val)
  
  
  
  
  for(breaks in 1:length(y_start)){
    
    if(breaks %in% c(2:(length(y_start)-1))){
      grid.lines(c(x_bottom_left, x_bottom_left+(x_length/5)),
                 c(y_start[breaks], y_start[breaks]),default.units = 'npc',gp=gpar(col='black'))
      grid.lines(c(x_bottom_left+x_length, x_bottom_left+x_length-(x_length/5)),
                 c(y_start[breaks], y_start[breaks]),default.units = 'npc',gp=gpar(col='black'))
    }
    
    grid.text(label[breaks],x=x_bottom_left+x_length+0.01,y=y_start[breaks],gp = gpar(fontsize = 12), 
              default.units = 'npc', just='left')
    
  }
  
  grid.polygon(c(x_bottom_left,x_bottom_left,x_bottom_left+x_length,x_bottom_left+x_length),
               c(y_bottom_left,y_bottom_left+y_length,y_bottom_left+y_length,y_bottom_left),
               default.units = 'npc',gp=gpar(fill=NA, col='black'))
  
  
  grid.text('Proportion',x=x_bottom_left-0.05, y=0.5,
            default.units = 'npc', gp = gpar(fontsize = 12), rot=90)
  
  
  
}


png('figure/correlation_supp.png',height=48,width=45,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(2,2,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=6,ncol=5)))

for(c in 1:5){
  for(r in 1:6){
    paneller(r,c)
  }
}

# paneller(1,1)
# paneller(1,2)
# paneller(1,3)

popViewport()
popViewport()
dev.off()
