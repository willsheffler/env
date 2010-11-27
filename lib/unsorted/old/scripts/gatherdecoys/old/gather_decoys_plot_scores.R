all <- read.table('all.fasc',header=T,comment.char='?')
lowscore <- read.table('lowscore.fasc',header=T,comment.char='?')
lowrms <- read.table('lowrms.fasc',header=T,comment.char='?')
gsbr<- read.table('gsbr.fasc',header=T,comment.char='?')
#nat <- read.table('relax_native.fasc',header=T)

xlim = c( 0.0, quantile(all$rms,0.90))
ylim = c( min(all$score) , quantile(all$score,0.90))

postscript(file='score_rms_plot.ps',horizontal=T)

plot(all$rms, all$score, xlim=xlim, ylim=ylim)
points(lowrms$rms , lowrms$score , col=3)
points(lowscore$rms , lowscore$score , col=4)
#points(nat$rms , nat$score , col=5)
points(gsbr$rms , gsbr$score , col=2)


