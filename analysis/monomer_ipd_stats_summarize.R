library(data.table)
convert_sum_data <- function(data) {
    # data: a data.table with columns:
    # kmer_string,kmer_number,position,chromosome,file_index,ipd_sum,ipd_sq_sum,log2_ipd_sum,log2_ipd_sq_sum,prediction_sum,prediction_sq_sum,log2_prediction_sum,log2_prediction_sq_sum,count
    # merge by chromosome and file_index
    data[, lapply(.SD, sum), keyby = .(kmer_string, kmer_number, position), .SDcols = !c("chromosome", "file_index")][,
        .(kmer_string, kmer_number, position, strand = "+",
          ipd_mean = ipd_sum / count, ipd_var = ipd_sq_sum / (count - 1) - ipd_sum ^ 2 / count / (count - 1),
          log2_ipd_mean = log2_ipd_sum / count, log2_ipd_var = log2_ipd_sq_sum / (count - 1) - log2_ipd_sum ^ 2 / count / (count - 1),
          prediction_mean = prediction_sum / count, prediction_var = prediction_sq_sum / (count - 1) - prediction_sum ^ 2 / count / (count - 1),
          log2_prediction_mean = log2_prediction_sum / count, log2_prediction_var = log2_prediction_sq_sum / count - log2_prediction_sum ^ 2 / count / (count - 1), count)]
}

stats_shi2016 <- fread("monomer_ipd_stats.shi2016_P6C4.csv")
summary_shi2016 <- convert_sum_data(stats_shi2016)
fwrite(summary_shi2016, file = "monomer_ipd_stats_summarize.csv")

