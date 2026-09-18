# Independent verification of thresholds and manuscript overlap from the annotated output.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
x <- read.delim(args[1], check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(length(unique(x$branch_test_statistic))==1)
effect <- if(unique(x$branch_test_statistic)=='mean') x$posterior_mean_change else x$MAP_change
family_threshold <- if('family_threshold' %in% names(x)) x$family_threshold else .05
branch_threshold <- if('branch_threshold' %in% names(x)) x$branch_threshold else .01
call <- x$family_p < family_threshold & x$branch_p < branch_threshold & effect != 0
stopifnot(identical(call, as.logical(x$selected_raw_thresholds)))
old <- as.logical(x$manuscript_significant)
te <- as.logical(x$TE_related)
result <- do.call(rbind, lapply(unique(x$manuscript_node), function(node) {
  take <- x$manuscript_node == node
  new_set <- paste(x[['Family ID']][take & call & !te], x$direction[take & call & !te])
  old_set <- paste(x[['Family ID']][take & old], x$manuscript_direction[take & old])
  shared <- intersect(new_set, old_set)
  data.frame(node=node, statistic=unique(x$branch_test_statistic),
             family_threshold=unique(family_threshold), branch_threshold=unique(branch_threshold),
             manuscript_events=length(old_set), new_nonTE_events=length(new_set),
             shared_same_direction=length(shared),
             manuscript_recall=length(shared)/length(old_set),
             new_set_overlap=if(length(new_set)) length(shared)/length(new_set) else NA_real_,
             Jaccard=length(shared)/length(union(old_set,new_set)))
}))
write.table(result, args[2], sep='\t', quote=FALSE, row.names=FALSE)
print(result)
