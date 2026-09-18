args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
d <- read.delim(args[1], check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(length(unique(d$branch_test_statistic))==1)
is_mean <- unique(d$branch_test_statistic)=='mean'
effect_label <- if(is_mean) 'Posterior mean child minus parent count' else 'Marginal-MAP child minus parent count'
p_label <- if(is_mean) '-log10 calibrated mean-change p' else '-log10 calibrated transition p'
branch_threshold <- if('branch_threshold' %in% names(d)) unique(d$branch_threshold) else .01
stopifnot(length(branch_threshold)==1)
pdf(args[2], width=10, height=7)
par(mfrow=c(2,2), mar=c(4,4,3,1))
for (node in unique(d$manuscript_node)) {
  z <- d[d$manuscript_node==node, ]
  main <- if(node==20) 'Social Vespidae' else 'Vespinae'
  highlight <- as.logical(z$manuscript_significant)
  effect <- if(is_mean) z$posterior_mean_change else z$MAP_change
  plot(effect, -log10(z$branch_p), pch=16, cex=.45,
       col=ifelse(highlight, '#b54930', '#74889a55'),
       xlab=effect_label,
       ylab=p_label, main=main)
  abline(h=-log10(branch_threshold), lty=2, col='#555555')
  legend('topright', 'Manuscript significant HOGs', col='#b54930', pch=16, bty='n', cex=.7)
  old <- z[highlight, ]
  plot(-log10(pmax(as.numeric(old$manuscript_branch_p),1e-8)),
       -log10(old$branch_p), pch=16, col='#b54930',
       xlab='-log10 manuscript Viterbi probability',
       ylab=p_label, main=paste(main, 'published events'))
  abline(h=-log10(branch_threshold), v=2, lty=2, col='#555555')
  # The two axes are different statistics; no equality/reference diagonal is implied.
}
dev.off()
