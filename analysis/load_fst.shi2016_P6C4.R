library(data.table)
library(fst)
library(doParallel)

fai.col.names <- c("name", "length", "offset", "linebases", "linewidth")
reference_metadata <- fread("/glusterfs/hisakatha/ensembl_GRCh38.p13/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai", col.names = fai.col.names)

chrs <- reference_metadata[, name]

cat("Start loading fst data\n")
shi2016_P6C4_data_list <- foreach(chr = chrs) %do% {
    read_fst(sprintf("load_and_save.shi2016_P6C4.%s.fst", chr), as.data.table = TRUE)
}
names(shi2016_P6C4_data_list) <- chrs
cat("Loaded shi2016_P6C4_data_list from load_and_save.shi2016_P6C4.*.fst\n")
