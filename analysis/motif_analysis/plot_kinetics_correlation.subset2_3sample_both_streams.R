library(data.table)
#args <- commandArgs(trailingOnly = TRUE)
#motif_dir <- args[1]
cwd_path <- strsplit(getwd(), "/")[[1]]
motif_dir <- cwd_path[length(cwd_path)]

#mean_log2value <- function(kinetics) {
#    occ_threshold <- 100
#    if (kinetics[, .N] == 0 || kinetics[strand == "+", length(unique(src))] < occ_threshold) {
#        data.table(position = numeric(0), strand = character(0), m = numeric(0))
#    } else {
#        # Use only the motif region
#        kinetics[value > 0 & is.finite(value) & strand == "+" & substr(label, 1, 1) == "m"][, .(position, strand, log2value = log2(value))][, .(m = mean(log2value)), by = .(position, strand)]
#    }
#}

occ_threshold <- 100
shi2016_mean_tmp <- fread("motif_ipd.shi2016_P6C4.summary2.csv")
position_motif_max <- shi2016_mean_tmp[region == "Motif", max(position)]
position_upstream_max <- shi2016_mean_tmp[, max(position)]
if(is.finite(position_motif_max)){ stopifnot(position_upstream_max - position_motif_max == 20) }
shi2016_mean <- shi2016_mean_tmp[strand == "+" & (-9 <= position & position <= position_motif_max + 10) & motif_occ >= occ_threshold, .(position, strand, shi2016 = log2ipd_mean)]
#fwrite(shi2016_mean)
#q()

celegans_dir <- "../../../../vs_ce11rel606/jun2018_analysis_high_mapq/motif_analysis/"
ab_mean <- fread(paste0(celegans_dir, motif_dir, "/motif_ipd.ab.c_elegans.summary3.csv"))[strand == "+" & (-9 <= position & position <= position_upstream_max - 10) & motif_occ >= occ_threshold, .(position, strand, ab = log2value_mean)]
cd_mean <- fread(paste0(celegans_dir, motif_dir, "/motif_ipd.cd.c_elegans.summary3.csv"))[strand == "+" & (-9 <= position & position <= position_upstream_max - 10) & motif_occ >= occ_threshold, .(position, strand, cd = log2value_mean)]

# Rename columns
setnames(shi2016_mean, "shi2016", "Human\n/native")
setnames(ab_mean, "ab", "Replicate 1\n/WGA\n/C. elegans")
setnames(cd_mean, "cd", "Replicate 2\n/WGA\n/C. elegans")

ipd_table <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(shi2016_mean, ab_mean, cd_mean))
ipd_table[, c("position", "strand") := NULL]

# function: corrplot.na
source("../plot_kinetics_correlation.corrplot.na.R")

cairo_pdf("plot_kinetics_correlation.subset2_3sample_both_streams.pdf", onefile = TRUE, width = 7, height = 7)
corrplot.na(ipd_table)
invisible(dev.off())
