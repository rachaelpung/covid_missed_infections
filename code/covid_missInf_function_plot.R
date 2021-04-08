# colour
CONFIG = list() 
# CONFIG$cols = c('#00468B','#ED0000','#42B540','#0099B4','#925E9F') # BMJ
CONFIG$colsDark = c('#00468BFF','#FDAF91FF','#42B540FF','#0099B4FF','#925E9FFF') # Lancet
CONFIG$colsPastel = c('#9CD6E9', '#E0A6A1', '#C9E0B4', '#039FC6', '#9FADD4')
# blue, pink, green, turquiose, purple


lightup = function(c, alpha)
{
  z=col2rgb(c)/255
  return(rgb(z[1],z[2],z[3],alpha))  # 0.125 for col3
}
CONFIG$colsDark1 = c();for(i in seq(CONFIG$colsDark))CONFIG$colsDark1[i] = lightup(CONFIG$colsDark[i], alpha = 0.6)
CONFIG$colsDark2 = c();for(i in seq(CONFIG$colsDark))CONFIG$colsDark2[i] = lightup(CONFIG$colsDark[i], alpha = 0.4)
CONFIG$colsDark3 = c();for(i in seq(CONFIG$colsDark))CONFIG$colsDark3[i] = lightup(CONFIG$colsDark[i], alpha = 0.25)

rm(lightup,i)
