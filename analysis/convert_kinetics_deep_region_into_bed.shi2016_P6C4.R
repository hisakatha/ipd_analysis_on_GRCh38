library(data.table)
library(doParallel)
library(fst)
kinetics_deep_region_to_bed <- function (fst_prefix, chrs, coverage_thres_high, coverage_thres_low, out) {
    file.remove(out)
    foreach_ret <- foreach(chr = chrs) %do% {
        gc()
        kinetics_data <- read_fst(sprintf("%s.%s.fst", fst_prefix, chr), as.data.table = TRUE)
        coverage_thres <- ifelse(grepl("^[XY]$", chr, perl = TRUE), coverage_thres_low, coverage_thres_high)
        cat(sprintf("Start processing: %s (coverage_thres = %d)\n", chr, coverage_thres))
        fwrite(kinetics_data[coverage >= coverage_thres, .(refName, tpl - 1, tpl, "deep_kinetics", coverage, ifelse(strand == 0, "+", "-"))] , sep = "\t", file = out, col.names = FALSE, append = TRUE)
        TRUE
    }
    unlist(foreach_ret)
}

fai.col.names <- c("name", "length", "offset", "linebases", "linewidth")
reference_metadata <- fread("/glusterfs/hisakatha/ensembl_GRCh38.p13/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai", col.names = fai.col.names)
chrs <- reference_metadata[, name]
fst_prefix <- "load_and_save.shi2016_P6C4"
coverage_thres_low <- 12
coverage_thres_high <- 25

kinetics_deep_region_to_bed(fst_prefix, chrs, coverage_thres_high, coverage_thres_low, "deep_kinetics_region.shi2016_P6C4.bed")
