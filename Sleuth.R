#load package
library (sleuth)

#read table describing samples and kallisto files
stab = read.table("sleuth_table.txt", header=TRUE)

#make sleuth object
so = sleuth_prep(stab)

#fit models to compare conditions and perform likelihood ratio test between conditions
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

#extract results
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#filter most significant results and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval)

#write values to file
write.table(sleuth_significant, file="sleuth_results.txt",quote = FALSE,row.names = FALSE)
