## analyze_factorbook_overlap_function.R
## alex amlie-wolf 02/02/16
## a script to analyze the factorbook results

analyze_factorbook_overlap <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh) {
    dir.create(paste0(outdir, "/tables/"), F, T)
    dir.create(paste0(outdir, "/plots/"), F, T)

    ## read in the factorbook file
    factorbook_file <- paste0(datadir, '/factorbook_overlap/', prefix, "_", r2_thresh,
                              "_ld_cutoff_snps_within_", dist_thresh, "_tfbs_overlaps.txt")
    factorbook_df <- read.table(factorbook_file, header=T, sep="\t", quote="", as.is=T)

    if(nrow(factorbook_df)==0) {
        cat("No FactorBook TFBS overlaps found in this dataset!\n")
        ## write a table
        return("No FactorBook TFBS overlaps found in this dataset!")
    } 

    unique_snp_overlap_df <- factorbook_df[!(duplicated(factorbook_df[,c("rsID", "tag_name", "tfbs_chr", "tfbs_start", "tfbs_end", "tf", "score", "cells")])),]
    uniq_tab_outf <- paste0(outdir, "/tables/", prefix, "_", r2_thresh,
                             "_ld_cutoff_snps_within_", dist_thresh, 
                               "_unique_tfbs_overlaps.txt")
    write.table(unique_snp_overlap_df, uniq_tab_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## count how many unique SNPs we get for each tag snp
    snp_counts <- ddply(unique_snp_overlap_df, .(tag_name), summarize, uniq_snp_count = length(unique(rsID)))

    snp_counts$tag_no_rsid <- gsub(":rs.*", "", snp_counts$tag_name)

    make_graphic(paste0(outdir, 'plots/', prefix, '_tfbs_overlapping_snps_per_region_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(snp_counts, aes_string(x=TAG_VAR, y="uniq_snp_count", fill=TAG_VAR)) +
          scale_fill_hue(h=c(180, 270)) + 
          xlab(TAG_LAB) + ylab("Number of unique SNPs overlapping TFBSs") +
          theme_bw() + geom_bar(position="dodge", stat="identity") +
          scale_y_continuous(breaks=seq(0, max(snp_counts$uniq_snp_count)+(5 - max(snp_counts$uniq_snp_count)%%5), by=5)) +
          theme(legend.position="none",
                axis.text.x=element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Number of SNPs overlapping TFBSs per region", r2_thresh, dist_thresh, out_subtitle))
    dev.off()
    
    ## for now, split by TF (need to do cell type too)
    regionwise_overlap_counts <- ddply(unique_snp_overlap_df, .(tag_name), function(x) {
        counts <- count(x$tf)
        colnames(counts) <- c("Transcription_Factor", "Num_Overlaps")
        return(counts)
    })

    regionwise_overlap_counts$tag_no_rsid <- gsub(":rs.*", "", regionwise_overlap_counts$tag_name)
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_tf_counts_by_region_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"), width_ratio=1.5)
    print(ggplot(regionwise_overlap_counts, aes_string(x="Transcription_Factor", y=TAG_VAR)) +
          geom_tile(aes(fill=log10(Num_Overlaps)), colour="white") +
          scale_fill_gradient(low="white", high="steelblue", na.value="white", guide="colorbar") +
          ylab(TAG_LAB) + xlab("Transcription Factor") + theme_bw() + 
          plot_title("Number of SNPs overlapping TF binding sites", r2_thresh, dist_thresh, out_subtitle) +
          theme(axis.text.x=element_text(angle=90, hjust=1, size=AXIS_TEXT_X_SIZE*0.3),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    
}
