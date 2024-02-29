library (sleuth)
library (dplyr)

stab = read.table("full_table.txt",header=TRUE)
so = sleuth_prep(stab)
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval)
write.table(sleuth_significant, file="sleuth_results.txt",quote = FALSE,row.names = FALSE)