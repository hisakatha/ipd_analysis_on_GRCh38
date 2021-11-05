library(data.table)
library(Biostrings)
library(hdf5r)
library(doParallel)
library(fst)

reference <- readDNAStringSet("/glusterfs/hisakatha/ensembl_GRCh38.p13/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa")

chrs <- names(reference)
chrs <- gsub("([^ ]+).*", "\\1", chrs, perl = TRUE)

parse_kinetics_h5 <- function(file, chromosomes) {
    data_hdf5 <- H5File$new(file, mode = "r")
    foreach_ret <- foreach(chr = chromosomes) %do% {
        cat(sprintf("Loading chr: %s\n", chr))
        kinetics_dt <- data.table(
            refName = chr,
            tpl = data_hdf5[[paste0(chr, "/tpl")]][],
            strand = data_hdf5[[paste0(chr, "/strand")]][],
            base = data_hdf5[[paste0(chr, "/base")]][],
            score = data_hdf5[[paste0(chr, "/score")]][],
            tMean = data_hdf5[[paste0(chr, "/tMean")]][],
            tErr = data_hdf5[[paste0(chr, "/tErr")]][],
            modelPrediction = data_hdf5[[paste0(chr, "/modelPrediction")]][],
            ipdRatio = data_hdf5[[paste0(chr, "/ipdRatio")]][],
            coverage = data_hdf5[[paste0(chr, "/coverage")]][],
            frac = data_hdf5[[paste0(chr, "/frac")]][],
            fracLow = data_hdf5[[paste0(chr, "/fracLow")]][],
            fracUp = data_hdf5[[paste0(chr, "/fracUp")]][])
        write_fst(kinetics_dt, sprintf("load_and_save.shi2016_P6C4.%s.fst", chr), compress = 100)
        TRUE
    }
    data_hdf5$close_all()
    stopifnot(unlist(foreach_ret))
}

data_path <- "../shi2016_P6C4/output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5"
parse_kinetics_h5(data_path, chrs)
cat("Loaded data shi2016_P6C4\n")
