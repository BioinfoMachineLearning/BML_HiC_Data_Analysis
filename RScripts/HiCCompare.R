## Biological group Old vs Young

library(multiHiCcompare)
library(BiocParallel)
library(data.table)
numCores <- 20
register(SnowParam(workers = numCores), default = TRUE)

chrs=c("chr1",  "chr2",  "chr3", "chr4", "chr5",
       "chr6",  "chr7",  "chr8", "chr9", "chrX",
       "chr10", "chr11", "chr12", "chr13", "chr14",
       "chr15", "chr16", "chr17", "chr18", "chr19")
wwug
resolutions = c("1000000", "500000", "100000", "50000")
samples <- list('of', 'om', 'yf', 'ym')

ROOT = "D:/Workspace/python-3/BML_HiC_Data_Analysis"
exp_cond = "Old_Young"

for (resolution in resolutions){

  mybiglist <- list()
   for (sample in samples){
        print(sprintf("Reading Sample %s", sample))
        in_table <- read.table(paste(ROOT, sprintf("MuSC_HiC_files/HiC_CPB_normalized/normalized_%s_%s.txt", sample, resolution), sep = "/"), header = FALSE)

        sparse <- HiCcompare::cooler2sparse(in_table)

     for (chr in chrs){
        chrom <- cbind(chr=chr, sparse[chr][[1]])
        mybiglist[[paste(chr, sample, sep="_")]] <- chrom
     }
   }

  for (chr in chrs){
    print(sprintf("Processing Chromosome %s", chr))

    print(sprintf("Creating plot & results directory for chromosome %s", chr))

    plt_dir <- paste(paste(ROOT, sprintf("Old_Young/Plots/HiCCompare_%s/%s", exp_cond, resolution), sep = "/"), chr, sep="/")
    res_dir <- paste(paste(ROOT, sprintf("Old_Young/Results/HiCCompare_%s/%s", exp_cond, resolution), sep = "/"), chr, sep="/")


    dir.create(file.path(plt_dir), recursive = TRUE)
    dir.create(file.path(res_dir), recursive = TRUE)

    hicexp  <- make_hicexp(mybiglist[[paste(chr, "om", sep="_")]],
                           mybiglist[[paste(chr, "of", sep="_")]],
                           mybiglist[[paste(chr, "ym", sep="_")]],
                           mybiglist[[paste(chr, "yf", sep="_")]],
                         groups = c(0, 0, 1, 1),
                         zero.p = 0.8, A.min = 5, filter = TRUE)
    hicexp <- cyclic_loess(hicexp, verbose = FALSE, parallel = TRUE, span = NA)

    write.csv(hicexp@hic_table, paste(res_dir, sprintf("%s_normalized.csv", chr), sep = "/"), row.names = FALSE)


    # make MD plot
    pdf(paste(plt_dir, sprintf("%s_MD_plot.pdf", chr), sep = "/"))
      MD_hicexp(hicexp)
    dev.off()

    hicexp <- hic_exactTest(hicexp, p.method = 'fdr', parallel = FALSE)

    # plot results
    jpeg(paste(plt_dir, sprintf("%s_composite_plot.jpg", chr), sep = "/"))
      MD_composite(hicexp)
    dev.off()

    # plot manhattan
    jpeg(paste(plt_dir, sprintf("%s_manhattan_plot.jpg", chr), sep = "/"))
      manhattan_hicexp(hicexp)
    dev.off()

    td <- topDirs(hicexp, logfc_cutoff = 0.5, logcpm_cutoff = 0.5, p.adj_cutoff = 0.2, return_df = 'pairedbed')
    write.csv(td, paste(res_dir, sprintf("%s_top_pair_bed.csv", chr), sep = "/"), row.names = FALSE)

    if(nrow(td) > 0){

      counts <- topDirs(hicexp, logfc_cutoff = 0.5, logcpm_cutoff = 0.5, p.adj_cutoff = 0.2, return_df = 'bed', pval_aggregate = "max")
      write.csv(counts, paste(res_dir, sprintf("%s_top_bed.csv", chr), sep = "/"), row.names = FALSE)

      # plot pvalue of counts
      jpeg(paste(plt_dir, sprintf("%s_pvalue_counts_plot.jpg", chr), sep = "/"))
        plot_pvals(counts)
      dev.off()

      # plot counts
      jpeg(paste(plt_dir, sprintf("%s__counts_plot.jpg", chr), sep = "/"))
        plot_counts(counts)
      dev.off()

    }
  }
}