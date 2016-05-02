
clusters = sapply(readLines("data/WORM/complexes_ratio5.txt"), strsplit, split="\t")

getGOs <- function(x){if(length(x)<4){return("")}else{return(paste(apply(gprofiler(x, organism="celegans", exclude_iea=T, max_set_size=500, correction_method="fdr")[,c(3,9,12)], 1, paste, collapse=","), collapse=";"))}}

gos = sapply(clusters[1:10], getGOs)

out = unlist(lapply(names(gos), function(x){paste(c(gsub("\t", ",", x), gos[[x]]), collapse="\t")}))

write.table(out, "test/Final_worm_complexes.go.txt", quote=F, row.names=F, col.names=F)
