# Standalone predictive-check figure from report_n0.py's exported table.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
x <- read.delim(args[1], check.names=FALSE)
keep <- c('single_species_fraction','all_species_present_fraction',
          'mean_species_present','differential_gt20_fraction')
labels <- c('Fraction observed in one species','Fraction observed in all 17 species',
            'Mean number of species occupied','Fraction with count differential >20')
x <- x[match(keep,x$feature),]
stopifnot(!anyNA(x$feature))
if(grepl('\\.png$',args[2],ignore.case=TRUE)) {
  png(args[2],width=1350,height=900,res=150)
} else {
  pdf(args[2],width=9,height=6)
}
par(mfrow=c(2,2),mar=c(3,4,3,1))
for(i in seq_len(nrow(x))) {
  z <- x[i,]
  vals <- c(z$observed,z$replicate_2.5_percentile,z$replicate_97.5_percentile)
  span <- max(diff(range(vals)),.01)
  plot(c(1,2),c(z$observed,z$replicate_mean),pch=19,col=c('#b54930','#386a8a'),
       ylim=range(vals)+c(-.15,.15)*span,xlim=c(.5,2.5),xaxt='n',
       xlab='',ylab='',main=labels[i])
  axis(1,at=c(1,2),labels=c('Observed','Model simulations'))
  arrows(2,z$replicate_2.5_percentile,2,z$replicate_97.5_percentile,
         angle=90,code=3,length=.07,col='#386a8a')
}
dev.off()
