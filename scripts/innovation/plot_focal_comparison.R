args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
d <- read.delim(args[1], check.names=FALSE, stringsAsFactors=FALSE)
pdf(args[2], width=10, height=7)
par(mfrow=c(2,2), mar=c(4,4,3,1))
for (node in unique(d$manuscript_node)) {
  z <- d[d$manuscript_node==node, ]
  main <- if(node==20) 'Social Vespidae' else 'Vespinae'
  highlight <- as.logical(z$manuscript_significant)
  plot(z$MAP_change, -log10(z$branch_p), pch=16, cex=.45,
       col=ifelse(highlight, '#b54930', '#74889a55'),
       xlab='Marginal-MAP child minus parent count',
       ylab='-log10 calibrated transition p', main=main)
  abline(h=2, lty=2, col='#555555')
  legend('topright', 'Manuscript significant HOGs', col='#b54930', pch=16, bty='n', cex=.7)
  old <- z[highlight, ]
  plot(-log10(pmax(as.numeric(old$manuscript_branch_p),1e-8)),
       -log10(old$branch_p), pch=16, col='#b54930',
       xlab='-log10 manuscript Viterbi probability',
       ylab='-log10 new calibrated transition p', main=paste(main, 'published events'))
  abline(h=2, v=2, lty=2, col='#555555')
  # The two axes are different statistics; no equality/reference diagonal is implied.
}
dev.off()
