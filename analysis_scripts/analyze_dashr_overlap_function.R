## analyze_dashr_overlap_function.R
## alex amlie-wolf 06/22/17
## analyzes overlap with DASHR ncRNA Loci
## part of automatic analysis

analyze_dashr_overlaps <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh) {
    ## -----------------------------------------------------------------------------
    ## read in data
    ## -----------------------------------------------------------------------------
    dir.create(paste0(outdir, 'plots/'), F, T)
    dir.create(paste0(outdir, 'tables/'), F, T)

    dashr_file <- paste0(datadir, '/dashr_ncrna_loci_overlap/', prefix, "_", r2_thresh,
                            "_ld_cutoff_snps_within_", dist_thresh, "_dashr_locus_overlap.txt")
    ## read in this data
    dashr_overlap_df <- read.table(dashr_file, header=T, sep="\t", quote="", as.is=T)
    
    if(nrow(dashr_overlap_df)==0) {
        ## still write out the table so that other stuff doesn't break
        dashr_overlap_df <- cbind(dashr_overlap_df, tag_no_rsid=numeric(0))
        ## write this table out
        dashr_outf <- paste0(outdir, 'tables/', prefix, '_', r2_thresh, '_ld_cutoff_snps_within_',
                             dist_thresh, '_dashr_locus_overlap.txt')
        write.table(dashr_overlap_df, dashr_outf, quote=F, sep="\t", col.names=T, row.names=F)
        cat("No DASHR overlaps found in this dataset!\n")
        return("No DASHR overlaps found in this dataset!")
    }

    ## add annotation for the tag regions without tag rsIDs
    dashr_overlap_df$tag_no_rsid <- gsub(":rs.*", "", dashr_overlap_df$tag_name)

    ## write this table out
    dashr_outf <- paste0(outdir, 'tables/', prefix, '_', r2_thresh, '_ld_cutoff_snps_within_',
                         dist_thresh, '_dashr_locus_overlap.txt')
    write.table(dashr_overlap_df, dashr_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## just make a simple barplot of the types of ncRNA overlaps found
    make_graphic(paste0(outdir, 'plots/', prefix, '_dashr_ncrna_type_counts_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(count(dashr_overlap_df, vars="rna_type"),
                 aes(x=rna_type, y=freq, fill=rna_type)) +
          xlab("Class of small RNA") + ylab("Number of Overlaps") +
          theme_bw() + geom_bar(stat="identity", position="dodge") +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Counts of DASHR overlaps by smRNA type", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    
}

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/dashr_ncrna_loci_overlap/")
