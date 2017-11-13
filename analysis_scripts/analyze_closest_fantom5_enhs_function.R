## analyze_closest_fantom5_enhs_function.R
## alex amlie-wolf 02/09/16
## a script that looks at the closest FANTOM5 enhancers for a single parameter setting
## part of automatic analysis

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/closest_fantom5_enhs/")
## tagsnp_cols <- c("tag_rsID", "tag_pos", "tag_MAF", "R2", "Dprime")

analyze_closest_fantom5_enhs <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh) {
    ## -----------------------------------------------------------------------------
    ## read in data
    ## -----------------------------------------------------------------------------
    dir.create(paste0(outdir, 'plots/'), F, T)
    dir.create(paste0(outdir, 'tables/'), F, T)
    
    ## read in the data file:
    closest_enh_file <- paste0(datadir, '/closest_fantom5_enh/', prefix, "_", r2_thresh,
                               "_ld_cutoff_snps_within_", dist_thresh,
                               "_closest_enhancers.txt")
    closest_enh_df <- read.table(closest_enh_file, header=T, sep="\t", quote="", as.is=T)

    if(nrow(closest_enh_df)==0) {
        cat("No closest enhancer data found in this dataset!\n")
        return("No closest enhancer data found in this dataset!")
    } 

    
    ## count how many tissues each closest enhancer comes from
    closest_enh_df <- cbind(closest_enh_df, num_tissues = unlist(lapply(
                                        strsplit(closest_enh_df$enh_tissue, ","), length)))

    ## add a nicely formatted FANTOM5 name description
    closest_enh_df$enh_tissue <- gsub("CL:[0-9]*_|UBERON:[0-9]*_", "", closest_enh_df$enh_tissue)

    ## use the tagsnp_cols 'global' variable   
    uniq_snp_cols <- colnames(closest_enh_df)[!(colnames(closest_enh_df) %in% tagsnp_cols)]
    unique_snp_closest_enh_df <- ddply(closest_enh_df, uniq_snp_cols, function(x) {
        apply(x[,!(colnames(x) %in% uniq_snp_cols)], 2, paste, collapse=",")
    })

    unique_snp_closest_enh_df$tag_no_rsid <- gsub(":rs.*", "", unique_snp_closest_enh_df$tag_name)

    uniq_snp_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                            "_ld_cutoff_snps_within_", dist_thresh,
                            "_uniq_snps_closest_enhs.txt")
    write.table(unique_snp_closest_enh_df, uniq_snp_outf, quote=F, sep="\t", col.names=T, row.names=F)
    
    ## -----------------------------------------------------------------------------
    ## analyze distances of unique snps to enhancers
    ## -----------------------------------------------------------------------------
    make_graphic(paste0(outdir, 'plots/', prefix, '_uniq_snp_dists_to_enhs_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))    
    print(ggplot(unique_snp_closest_enh_df, aes(x="", y=dist_to_enh_midpoint, colour=abs(dist_to_enh_midpoint))) +
          ylab("Distance to midpoint of enhancer transcription") +
          theme_bw() + geom_point(position=position_jitter(width=0.2), alpha=0.4) +
          scale_colour_gradient(low="#56B1F7", high="navy") + 
          geom_boxplot(fill=NA, outlier.size=0) +
          scale_y_continuous(labels=comma, breaks=distance_breaks) +
          theme(legend.position="none", axis.text.x = element_blank(),
                axis.text.y=element_text(size=AXIS_TEXT_Y_SIZE*0.75), title=element_text(size=TITLE_SIZE),
                plot.title = element_text(hjust = 0.5)) +
          xlab("") +
          ## geom_hline(yintercept = enh_window, colour=muted("red"), linetype=2, alpha=0.5) +
          ## geom_hline(yintercept = -enh_window, colour=muted("red"), linetype=2, alpha=0.5) +
          plot_title("Distributions of distances of SNPs from midpoint of closest enhancer transcription", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    make_graphic(paste0(outdir, 'plots/', prefix, '_uniq_snp_dists_to_enhs_by_tagregion_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))    
    print(ggplot(unique_snp_closest_enh_df, aes_string(x=TAG_VAR, y="dist_to_enh_midpoint") +
                 aes(colour=abs(dist_to_enh_midpoint))) +
          xlab(TAG_LAB) + 
          ylab("Distance to midpoint of enhancer transcription") +
          theme_bw() + geom_point(position=position_jitter(width=0.2), alpha=0.4) +
          scale_colour_gradient(low="#56B1F7", high="navy") + 
          geom_boxplot(fill=NA, outlier.size=0) +
          scale_y_continuous(labels=comma, breaks=distance_breaks) +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE*0.75),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          ## geom_hline(yintercept = enh_window, colour=muted("red"), linetype=2, alpha=0.5) +
          ## geom_hline(yintercept = -enh_window, colour=muted("red"), linetype=2, alpha=0.5) +
          plot_title("Distributions of distances of SNPs from midpoint of closest enhancer transcription", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

}
