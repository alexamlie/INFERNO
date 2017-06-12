## analyze_homer_motif_overlaps_function.R
## alex amlie-wolf 07/20/16
## looks at the HOMER motif overlap data

analyze_homer_motif_overlaps <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh, pwm_calc=FALSE) {
    ## -----------------------------------------------------------------------------
    ## read in data
    ## -----------------------------------------------------------------------------
    dir.create(paste0(outdir, 'plots/'), F, T)
    dir.create(paste0(outdir, 'tables/'), F, T)

    ## do the appropriate analysis, depending on whether the PWM changes were calculated
    if (pwm_calc) {
        motif_file <- paste0(datadir, '/homer_motif_overlap/', prefix, "_", r2_thresh,
                                "_ld_cutoff_snps_within_", dist_thresh, "_motif_overlap_and_disruption.txt")
        ## read in this data
        motif_overlap_df <- read.table(motif_file, header=T, sep="\t", quote="", as.is=T)
        if(nrow(motif_overlap_df)==0) {
            ## still write out the table so that other stuff doesn't break
            motif_overlap_df <- cbind(motif_overlap_df, tag_no_rsid=numeric(0))
            ## write this table out
            motif_outf <- paste0(outdir, 'tables/', prefix, '_', r2_thresh, '_ld_cutoff_snps_within_',
                                 dist_thresh, '_homer_motif_overlaps.txt')
            write.table(motif_overlap_df, motif_outf, quote=F, sep="\t", col.names=T, row.names=F)
            cat("No HOMER motif overlaps found in this dataset!\n")
            return("No HOMER motif overlaps found in this dataset!")
        }
        
        ## add annotation for the tag regions without tag rsIDs
        motif_overlap_df$tag_no_rsid <- gsub(":rs.*", "", motif_overlap_df$tag_name)

        ## write this table out
        motif_outf <- paste0(outdir, 'tables/', prefix, '_', r2_thresh, '_ld_cutoff_snps_within_',
                             dist_thresh, '_homer_motif_overlaps.txt')
        write.table(motif_overlap_df, motif_outf, quote=F, sep="\t", col.names=T, row.names=F)

        ## also make this into a bed file of disrupted motifs
        ## to do this, get all the unique motif occurrences
        uniq_motif_overlaps <- unique(motif_overlap_df[,c("motif_chr", "motif_start", "motif_end",
                                                          "tf_name", "log_odds_score", "strand")])
        ## map the log odds score to the bed range (0-1000)
        ## https://stackoverflow.com/questions/18303420/how-to-map-a-vector-to-a-different-range-in-r
        uniq_motif_overlaps$bed_score <- with(uniq_motif_overlaps,
                                              (log_odds_score-min(log_odds_score)) /
                                              max(log_odds_score-min(log_odds_score)) * 999+1)
        ## now output this table in bed format
        motif_bed_outf <- paste0(outdir, 'tables/', prefix, '_', r2_thresh, '_ld_cutoff_snps_within_',
                             dist_thresh, '_homer_motif_overlaps.bed')
        ## first write the track description to this file:
        cat(paste0("track type=bed name=\"Motif hits, ", out_subtitle, "\" description=\"HOMER-identified TFBS disruptions by INFERNO SNPs\" visibility=3 useScore=1\n"), file=motif_bed_outf)
        ## now write the data
        write.table(uniq_motif_overlaps[,c("motif_chr", "motif_start", "motif_end", "tf_name", "bed_score", "strand")], motif_bed_outf, quote=F, sep="\t", row.names=F, col.names=F, append=T)
        
        ## make a figure showing the delta PWM scores
        make_graphic(paste0(outdir, 'plots/', prefix, '_delta_pwm_distributions_by_tagregion_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
        print(ggplot(motif_overlap_df, aes_string(x=TAG_VAR, y="delta_pwm", fill=TAG_VAR)) +
              xlab(TAG_LAB) +
              ylab("Delta PWM score") +
              theme_bw() + geom_violin() +
              geom_jitter(aes_string(fill=TAG_VAR), colour='black', shape=23, alpha=0.5, height=0) +
              scale_fill_hue(h=c(180, 270)) +
              scale_colour_hue(h=c(180, 270)) +
              theme(legend.position="none",
                    axis.text.x=element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
              plot_title("Distributions of PWM changes for motif-overlapping variants", r2_thresh, dist_thresh, out_subtitle))
        dev.off()
        
    } else {
        motif_file <- paste0(datadir, '/homer_motif_overlap/', prefix, "_", r2_thresh,
                                "_ld_cutoff_snps_within_", dist_thresh, "_motif_overlap.txt")
    }

}

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/homer_motif_overlap/")
