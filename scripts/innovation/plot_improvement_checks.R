args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
d <- read.delim(args[1], check.names=FALSE)
metrics <- c('single_species_fraction','one_copy_every_species_fraction',
             'mean_species_present','differential_gt20_fraction')
titles <- c('Present in one species','One copy in every species',
            'Mean number of species occupied','Count differential >20')
models <- unique(d$model)
labels <- c(asymmetric='Unequal rates',gamma12='12 gamma categories',
 gamma24='24 gamma categories',zero_error0='Zero error = 0',
 zero_error001='Zero error = 0.001',zero_error_free='Separate zero error',
 hurdle_poisson='Hurdle Poisson root',hurdle_nb='Hurdle NB root',
 combined_poisson_zero0='Combined / Poisson / zero=0',
 combined_nb_zero_free='Combined / NB / free zero error')
if(grepl('\\.png$',args[2])) png(args[2],width=1900,height=1400,res=150) else pdf(args[2],width=13,height=10)
par(mfrow=c(2,2),mar=c(4,12,3,1),oma=c(1,0,2,0))
for(j in seq_along(metrics)) {
 x <- d[d$metric==metrics[j],];x <- x[match(models,x$model),]
 yl <- rev(seq_along(models));limits <- range(c(x$observed,x$simulated_min,x$simulated_max),na.rm=TRUE)
 pad <- max(diff(limits)*.08, .001);limits <- limits+c(-pad,pad)
 plot(x$simulated_mean,yl,xlim=limits,ylim=c(.5,length(models)+.5),yaxt='n',
      xlab=if(j==3) 'Species' else 'Fraction of families',ylab='',pch=19,col='#286387',main=titles[j])
 segments(x$simulated_min,yl,x$simulated_max,yl,col='#286387',lwd=2)
 abline(v=unique(x$observed),col='#b23b2c',lty=2,lwd=2)
 lab <- labels[models];lab[is.na(lab)] <- models[is.na(lab)]
 axis(2,at=yl,labels=lab,las=1,cex.axis=.75)
}
mtext('N11 predictive screening: blue = simulations; red dashed = observed',outer=TRUE,side=3,line=.5,cex=1.1)
mtext('Ranges show replicated datasets at fitted parameters; these are not biological confidence intervals.',outer=TRUE,side=1,line=0,cex=.8)
dev.off()
