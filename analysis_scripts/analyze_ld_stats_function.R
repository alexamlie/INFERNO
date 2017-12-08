## analyze_ld_stats_function.R
## alex amlie-wolf 01/26/2016
## a script that looks at the number of SNPs gained by LD expansion
## part of the automatic analysis pipeline, so this is just a function

## prefix <- param_ref[['outprefix']]
## datadir <- paste0(param_ref[['outdir']], "/ld_expansion/")
## outdir <- paste0(result_outdir, "/ld_stats/")

analyze_ld_stats <- function(prefix, datadir, outdir, out_subttitle, r2_thresh, dist_thresh) {
    dir.create(paste0(outdir, 'plots/'), F, T)
    dir.create(paste0(outdir, 'tables/'), F, T)

    ## read in the data file:
    ld_stats_file <- paste0(datadir, prefix, "_", r2_thresh, "_ld_cutoff_snps_within_",
                            dist_thresh, ".txt")
    ld_stats_df <- read.table(ld_stats_file, header=T, sep="\t", quote="", as.is=T)

    ## add a column to the LD stats dataframe
    ld_stats_df$tag_no_rsid <- gsub(":rs.*", "", ld_stats_df$tag_name)    
    
    ## get unique entries per tag region, write out tables
    ## use the 'global' tagsnp column variable    
    uniq_snp_cols <- colnames(ld_stats_df)[!(colnames(ld_stats_df) %in% tagsnp_cols)]
    ## get all the unique values of the info we care about, and just paste the others together
    unique_snp_df <- ddply(ld_stats_df, uniq_snp_cols, function(x) {
        apply(x[,!(colnames(x) %in% uniq_snp_cols)], 2, paste, collapse=",")
    })

    ## write the file out
    uniq_snp_outf <- paste0(outdir, 'tables/', prefix, '_', r2_thresh,
                            '_ld_cutoff_snps_within_', dist_thresh,
                            '_uniq_snps.txt')
    write.table(unique_snp_df, uniq_snp_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## Analyze LD SNPs
    ## first, just count how many unique SNPs we got from each tag SNP region
    snp_counts <- ddply(ld_stats_df, .(tag_name, tag_no_rsid), summarize, uniq_snp_count = length(unique(rsID)))
    
    max_snp_count <- round_any(max(snp_counts$uniq_snp_count)+1, accuracy=10, f=ceiling)
    if(max_snp_count > 100) {
        y_tick_freq <- 50
    } else {
        y_tick_freq <- 10
    }
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_snp_counts_per_region_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(snp_counts, aes_string(x=TAG_VAR, y="uniq_snp_count",
                                        fill=TAG_VAR)) +
          scale_fill_hue(h=c(180, 270)) + 
          xlab(TAG_LAB) + ylab("Number of unique SNPs in LD expanded set") +
          theme_bw() + geom_bar(position="dodge", stat="identity") +
          scale_y_continuous(breaks=seq(0, max_snp_count, by=y_tick_freq),
                             limits=c(0, max_snp_count), expand=c(0, 0)) + 
          theme(legend.position="none",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Number of SNPs gained by LD expansion", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## look at the distributions of distances from tag SNPs, by region:
    make_graphic(paste0(outdir, 'plots/', prefix, '_dist_to_tagsnp_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(ld_stats_df, aes_string(x=TAG_VAR, fill=TAG_VAR) + aes(y=tag_pos - pos)) +
          scale_fill_hue(h=c(180, 270)) + 
          xlab(TAG_LAB) + ylab("Distance from LD SNP to SNP in LD with tag SNP(bp)") +
          theme_bw() + geom_boxplot() +
          scale_y_continuous(labels=comma) +
          theme(legend.position="none",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of distances from linked SNPs to tag SNP", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## Compare LD stats and investigate distributions
    ## distributions of R^2 values in linked SNPs
    make_graphic(paste0(outdir, 'plots/', prefix, '_R2_distributions_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(ld_stats_df, aes_string(x=TAG_VAR, y="R2", fill=TAG_VAR)) +
          scale_fill_hue(h=c(180, 270)) + 
          xlab(TAG_LAB) + ylab("R^2 Statistic") +
          theme_bw() + geom_boxplot() + 
          theme(legend.position="none",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distribution of R^2 values across tag SNPs", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## also look at D'
    make_graphic(paste0(outdir, 'plots/', prefix, '_dprime_distributions_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(ld_stats_df, aes_string(x=TAG_VAR, y="Dprime", fill=TAG_VAR)) +
          scale_fill_hue(h=c(180, 270)) + 
          xlab(TAG_LAB) + ylab("D\' Statistic") +
          theme_bw() + geom_boxplot() +
          theme(legend.position="none",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distribution of D\' values across tag SNPs", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## next, look at the MAF distributions (only if we have this data)
    if(sum(!is.na(unique_snp_df$MAF)) > 0) {
        make_graphic(paste0(outdir, 'plots/', prefix, '_MAF_distributions_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"))
        print(ggplot(unique_snp_df, aes_string(x=TAG_VAR, y="MAF", fill=TAG_VAR)) +
              scale_fill_hue(h=c(180, 270)) +
              xlab(TAG_LAB) + ylab("Minor Allele Frequency") +
              theme_bw() + geom_boxplot() +
              theme(legend.position="none",
                    axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
              plot_title("Distribution of minor allele frequency values of SNPs linked with tag SNPs", r2_thresh, dist_thresh, out_subtitle))
        dev.off()
    }

    ## also look at the allele counts
    ## we want to count the number of combinations of reference and alternative alleles
    allele_count_tab <- count(unique_snp_df, vars=c("ref", "alt"))
    allele_count_tab$alt <- factor(allele_count_tab$alt, ordered=T,
                                   levels=sort(unique(allele_count_tab$alt), decreasing=T))
    
    ## light yellow to red
    heatmap_col_lims <- c("#FFFFDD", "#73000D")
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_full_allele_rawcount_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)                   
    print(ggplot(allele_count_tab, aes(x=ref, y=alt)) +
          geom_tile(aes(fill=freq), colour=heatmap_col_lims[1]) +
          scale_fill_gradient("Number of SNPs", low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
          theme_bw() +
          scale_x_discrete(position="top") + 
          xlab("Reference Allele") + ylab("Alternative Allele") +
          plot_title("Heatmap of numbers of combinations of reference and alternative alleles",
                     r2_thresh, dist_thresh, out_subtitle, strwrap_width=60) + 
          theme(axis.text.x = element_text(angle=45, hjust=0, vjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                panel.grid.major=element_blank(), panel.grid.minor=element_blank()))
    dev.off()

    make_graphic(paste0(outdir, 'plots/', prefix, '_full_allele_log10_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)
    print(ggplot(allele_count_tab, aes(x=ref, y=alt)) +
          geom_tile(aes(fill=log10(freq)), colour=heatmap_col_lims[1]) +
          scale_fill_gradient("Log10(Number of SNPs)", low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
          theme_bw() +
          scale_x_discrete(position="top") + 
          xlab("Reference Allele") + ylab("Alternative Allele") +
          plot_title("Heatmap of numbers of combinations of reference and alternative alleles",
                     r2_thresh, dist_thresh, out_subtitle, strwrap_width=60) + 
          theme(axis.text.x = element_text(angle=45, hjust=0, vjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                panel.grid.major=element_blank(), panel.grid.minor=element_blank()))   
    dev.off()

}

