## analyze_fantom5_overlap_function.R
## alex amlie-wolf 02/08/16
## looks at the overlap with FANTOM5 enhancers
## part of automatic analysis

analyze_fantom5_overlap <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh, enh_overlap_type, enh_window, fantom5_class_file, enh_summary_flag=FALSE) {    
    ## ---------------------------------------------------
    ## read in data
    ## ---------------------------------------------------    
    dir.create(paste0(outdir, "/tables/"), F, T)
    dir.create(paste0(outdir, "/plots/"), F, T)

    ## read in the enhancer overlap data file
    if (enh_overlap_type == "midpoint") {
        enh_overlap_file <- paste0(datadir, '/fantom5_overlap/',
                                   prefix, "_", r2_thresh,
                                   "_ld_cutoff_snps_within_", dist_thresh,
                                   "_", enh_window, "bp_around_midpoint_enh_overlaps.txt")
        enh_overlap_df <- read.table(enh_overlap_file, header=T, sep="\t", quote="", as.is=T)
        if (enh_summary_flag) {
            ## read in the data file
            enh_summary_file <- paste0(datadir, "/fantom5_overlap/", prefix, "_", r2_thresh,
                                       "_ld_cutoff_snps_within_", dist_thresh,
                                       "_", enh_window, "bp_around_midpoint_summary.txt")
            enh_summary_df <- read.table(enh_summary_file, header=T, sep="\t", quote="", as.is=T)
        }        
    } else if (enh_overlap_type == "locus") {
        enh_overlap_file <- paste0(datadir, '/fantom5_overlap/',
                                   prefix, "_", r2_thresh,
                                   "_ld_cutoff_snps_within_", dist_thresh,
                                   "_", enh_window, "bp_around_orig_locus_enh_overlaps.txt")
        enh_overlap_df <- read.table(enh_overlap_file, header=T, sep="\t", quote="", as.is=T)
        if (enh_summary_flag) {
            ## read in the data file
            enh_summary_file <- paste0(datadir, "/fantom5_overlap/", prefix, "_", r2_thresh,
                                       "_ld_cutoff_snps_within_", dist_thresh,
                                       "_", enh_window, "bp_around_orig_locus_summary.txt")
            enh_summary_df <- read.table(enh_summary_file, header=T, sep="\t", quote="", as.is=T)
        }        
    }
    
    ## read in the tissue categories
    fantom5_category_df <- read.table(fantom5_class_file, header=T, sep="\t", quote="", as.is=T)
    ## add a column for easier matching
    fantom5_category_df$enh_source <- gsub("_expressed_enhancers.bed", "", fantom5_category_df$FANTOM5.File)
    
    ## add columns to the enhancer data:
    enh_overlap_df$enh_class <- fantom5_category_df$Class[match(enh_overlap_df$enh_source,
                                                               fantom5_category_df$enh_source)]    
    ## clean up the tissue names
    enh_overlap_df$enh_source <- gsub("CL:[0-9]*_|UBERON:[0-9]*_", "", enh_overlap_df$enh_source)
    ## add a column for tag without rsid
    enh_overlap_df$tag_no_rsid <- gsub(":rs.*", "", enh_overlap_df$tag_name)
    ## add a distance column: note that i round down the midpoint, for interpretability
    enh_overlap_df$dist_to_mid <- enh_overlap_df$pos - floor(rowMeans(cbind(enh_overlap_df$enh_start, enh_overlap_df$enh_end-1)))
    ## also add a column for if it overlapped the original (in the case of a 0-window around
    ## the locus, this should always be the case)
    enh_overlap_df$overlap_orig <- ifelse(enh_overlap_df$enh_start <= enh_overlap_df$pos &
                                          enh_overlap_df$pos < enh_overlap_df$enh_end,
                                          "eRNA Overlap",
                                          "No eRNA Overlap")
    
    ## get unique entries per tag region, write out tables
    ## get all the columns except the tagging SNP info (since we want unique per tag region)
    ## use the 'global' tagsnp column variable    
    uniq_snp_overlap_cols <- colnames(enh_overlap_df)[!(colnames(enh_overlap_df) %in% tagsnp_cols)]
    ## get all the unique values of the info we care about, and just paste the others together
    unique_enh_overlap_df <- ddply(enh_overlap_df, uniq_snp_overlap_cols, function(x) {
        apply(x[,!(colnames(x) %in% uniq_snp_overlap_cols)], 2, paste, collapse=",")
    })   
    
    ## write the full file out
    uniq_snp_overlap_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                                    "_ld_cutoff_snps_within_", dist_thresh,
                                    "_", enh_overlap_type, "_", enh_window, "bp_window",
                                    "_uniq_fantom5_overlap_snps.txt")
    write.table(unique_enh_overlap_df, uniq_snp_overlap_outf, quote=F, sep="\t", col.names=T, row.names=F)
    ## also write it without the extra SNP info columns
    uniq_snp_overlap_small_out <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                                         "_ld_cutoff_snps_within_", dist_thresh,
                                         "_", enh_overlap_type, "_", enh_window, "bp_window",
                                         "_uniq_fantom5_overlap_snps_no_tagsnp_info.txt")
    write.table(unique_enh_overlap_df[,uniq_snp_overlap_cols],
                uniq_snp_overlap_small_out, quote=F, sep="\t", col.names=T, row.names=F)    

    ## write a bed file of all the enhancers that we overlap, for browser analysis
    enhancer_bed_data <- unique(enh_overlap_df[,c("enh_source", "enh_chr", "enh_start", "enh_end")])
    enh_bed_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                           "_ld_cutoff_snps_within_", dist_thresh,
                           "_", enh_overlap_type, "_", enh_window, "bp_window",
                           "_fantom5_enh_overlaps.bed")   
    ## first write the track description to this file:
    cat(paste0("track type=bed name=\"FANTOM5 enhancers, ", out_subtitle, "\" description=\"FANTOM5 enhancers overlapping INFERNO SNPs\" visibility=3\n"), file=enh_bed_outf)
    
    write.table(enhancer_bed_data[,c("enh_chr", "enh_start", "enh_end", "enh_source")],
                enh_bed_outf, append=T, quote=F, sep="\t", row.names=F, col.names=F)    
    
    ## ---------------------------------------------------
    ## analysis, generate figures
    ## ---------------------------------------------------    
    ## ---------------------------
    ## look at the individual overlaps    
    ## get distribution of distances from midpoint per tag region
    make_graphic(paste0(outdir, 'plots/', prefix, '_enh_overlap_dists_per_tagregion_', 
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(unique_enh_overlap_df, aes_string(x=TAG_VAR, y="dist_to_mid") +
                 aes(colour=abs(dist_to_mid))) +
          xlab(TAG_LAB) + ylab("Distance to midpoint of enhancer transcription") +
          theme_bw() + geom_boxplot(outlier.size=0) +
          geom_point(position=position_jitter(width=0.2)) +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of distances of SNPs from midpoint of enhancer transcription", r2_thresh, dist_thresh, paste0(enh_window, "bp on either side of ", enh_overlap_type, ", ", out_subtitle)))
    dev.off()

    ## compare the enhancer length distribution between overlapping and non-overlapping
    ## find unique enhancer-SNP combinations
    enh_length_overlap_comparison <- unique_enh_overlap_df[!duplicated(unique_enh_overlap_df[,c("pos", "enh_chr", "enh_start", "enh_end")]), c("pos", "enh_chr", "enh_start", "enh_end", "overlap_orig")]
    enh_length_overlap_comparison$enh_length <- enh_length_overlap_comparison$enh_end - enh_length_overlap_comparison$enh_start
    
    max_length <- round_any(max(enh_length_overlap_comparison$enh_length), 100, f=ceiling)
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_enh_length_dists_', 
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(enh_length_overlap_comparison, aes(x=overlap_orig, y=enh_length, colour=enh_length)) +
          xlab("Overlap with original FANTOM5 eRNA locus") +
          ylab("Length of enhancer") +
          theme_bw() + geom_boxplot(outlier.size=0) +
          geom_point(position=position_jitter(width=0.2)) +
          scale_y_continuous(breaks=seq(0, max_length, by=50), limits=c(0, max_length), expand=c(0, 0)) +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of lengths of enhancers overlapped by variants", r2_thresh, dist_thresh, paste0(enh_window, "bp on either side of ", enh_overlap_type, ", ", out_subtitle)))
    dev.off()
    
    ## also get distribution for all hits
    make_graphic(paste0(outdir, 'plots/', prefix, '_all_enh_overlap_dists_', 
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(unique_enh_overlap_df, aes(x='', y=dist_to_mid, colour=abs(dist_to_mid))) +
          ylab("Distance to midpoint of enhancer transcription") +
          theme_bw() + geom_boxplot(outlier.size=0) +
          geom_point(position=position_jitter(width=0.2)) +          
          theme(legend.position="none", axis.text.x = element_blank(),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of distances of SNPs from midpoint of enhancer transcription", r2_thresh, dist_thresh, paste0(enh_window, "bp on either side of ", enh_overlap_type, ", ", out_subtitle))) 
    dev.off()

    ## look at the number of SNPs that overlapped the original enhancer (looking only at unique
    ## enhancers)
    ## to make sure the table results include both conditions, make this a factor
    unique_enh_overlap_df$overlap_orig <- factor(unique_enh_overlap_df$overlap_orig)

    ## get the number and the proportion of overlapping SNPs per region
    num_erna_overlaps_per_region <- ddply(unique_enh_overlap_df, .(tag_no_rsid, tag_name), function(x) {
        uniq_x <- x[!duplicated(x[,c("rsID", "enh_start", "enh_end")]),]
        ret_table <- data.frame(table(uniq_x$overlap_orig, dnn=c("overlap_orig")))
        ret_table$prop <- ret_table$Freq / sum(ret_table$Freq)
        return(ret_table)
    })

    ## make a color palette vector for this:
    overlap_colors <- c("eRNA Overlap" = "#95001A", "No eRNA Overlap" = "#00674D")
    
    ## plot the numbers of eRNA overlapping and non-overlapping for each tag region
    make_graphic(paste0(outdir, 'plots/', prefix, '_orig_enh_overlap_nums_', 
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(num_erna_overlaps_per_region, aes(x=overlap_orig, y=Freq, fill=overlap_orig)) +
          xlab("") + ylab("Number of SNPs") + 
          theme_bw() + geom_bar(stat="identity", position="stack") +
          facet_grid(reformulate(TAG_VAR)) +
          scale_y_continuous(breaks=seq(0, max(aggregate(Freq~tag_no_rsid+overlap_orig, num_erna_overlaps_per_region, sum)$Freq), by=1)) +
          scale_fill_manual(values=overlap_colors) +
          theme(legend.position="none",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE*0.75),
                strip.text = element_text(size=LEGEND_TEXT_SIZE*0.75), 
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Number of SNPs overlapping original eRNA transcription loci", r2_thresh, dist_thresh, paste0(enh_window, "bp on either side of ", enh_overlap_type, ", ", out_subtitle))) 
    dev.off()

    ## plot the proportions of eRNA overlapping and non-overlapping for each tag region
    make_graphic(paste0(outdir, 'plots/', prefix, '_orig_enh_overlap_props_', 
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(num_erna_overlaps_per_region, aes_string(x=TAG_VAR, y="prop", fill="overlap_orig")) +
          xlab(TAG_LAB) + ylab("Proportion of SNPs") + 
          theme_bw() + geom_bar(stat="identity", position="stack") +
          scale_fill_manual(values=overlap_colors) +
          theme(legend.position="bottom",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE), 
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Proportion of SNPs overlapping original eRNA transcription loci for each tag region", r2_thresh, dist_thresh, paste0(enh_window, "bp on either side of ", enh_overlap_type, ", ", out_subtitle))) 
    dev.off()
    
    ## also look at the distribution of numbers of unique enhancers for each SNP
    num_enhs_per_snp_region <- ddply(unique_enh_overlap_df, .(tag_no_rsid, tag_name), function(x) {
        ## for each unique SNP, count how many enhancers it overlapped
        snp_enh_count <- ddply(x, .(rsID), summarize,
                               num_enhs = length(unique(paste(enh_start, enh_end))))
        return(data.frame(table(snp_enh_count$num_enhs, dnn=c("num_enhs"))))
    })
    ## make it a factor so we get all levels
    num_enhs_per_snp_region$num_enhs <- factor(num_enhs_per_snp_region$num_enhs)
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_enh_overlap_nums_byregion_', 
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(num_enhs_per_snp_region, aes_string(x=TAG_VAR, y="Freq") +
                 aes(fill=factor(num_enhs))) +
          xlab(TAG_LAB) + ylab("Number of enhancer overlaps") + 
          theme_bw() + geom_bar(stat='identity', position='dodge') +
          scale_fill_hue(h=c(60, 200),
                         guide = guide_legend(title="Number of unique enhancers overlapped")) +
          scale_y_continuous(breaks=seq(0, max(num_enhs_per_snp_region$Freq), by=2)) +
          theme(legend.position="bottom",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE), 
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Number of unique enhancers overlapping by each SNP in each tag region", r2_thresh, dist_thresh, paste0(enh_window, "bp on either side of ", enh_overlap_type, ", ", out_subtitle)))
    dev.off()
        
    ## ---------------------------
    ## look at the overlaps using the summary file, if it exists
    if (enh_summary_flag) {
        ## get a cross-reference for the tissue category that each column corresponds to
        summary_sources <- gsub(".", ":",
                                gsub("_num_overlaps|_dists", "", colnames(enh_summary_df)),
                                fixed=T)
        ## this index will allow us to retrieve the tissue categories
        summary_category_idx <- match(paste0(summary_sources, "_expressed_enhancers.bed"),
                                      fantom5_category_df$FANTOM5.File)
        
        colnames(enh_summary_df) <- gsub("CL.[0-9]*_|UBERON.[0-9]*_", "", colnames(enh_summary_df))

        ## make this a named vector so we can use it later (without the CL and UBERON numbers!)
        names(summary_category_idx) <- colnames(enh_summary_df)        
        
        uniq_snp_summ_cols <- colnames(enh_summary_df)[!(colnames(enh_summary_df) %in% tagsnp_cols)]
        unique_enh_summary_df <- ddply(enh_summary_df, uniq_snp_summ_cols, function(x) {
            apply(x[,!(colnames(x) %in% uniq_snp_summ_cols)], 2, paste, collapse=",")
        })
        
        uniq_snp_summ_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                                     "_ld_cutoff_snps_within_", dist_thresh,
                                     "_", enh_overlap_type, "_", enh_window, "bp_window",
                                     "_uniq_fantom5_snps_summary.txt")
        write.table(unique_enh_summary_df, uniq_snp_summ_outf, quote=F, sep="\t", col.names=T, row.names=F)
        
        uniq_snp_summ_small_out <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                                          "_ld_cutoff_snps_within_", dist_thresh,
                                          "_", enh_overlap_type, "_", enh_window, "bp_window",
                                          "_uniq_fantom5_snps_summary_no_tagsnp_info.txt")
        write.table(unique_enh_summary_df[,uniq_snp_summ_cols],
                    uniq_snp_summ_small_out, quote=F, sep="\t", col.names=T, row.names=F)        

        ## ---------------------------------------------------
        ## analysis of summary data (across all tissues)
        ## note that from now on we only use the unique snp:tag region combos
        ## ---------------------------------------------------        
        ## define indexing vectors for the two types of data
        ## first find all the SNP info columns
        snp_info_cols <- grep("overlaps|dists", colnames(unique_enh_summary_df), invert=T)
        ## then find the overlap or distance column                
        num_overlap_idx <- c(snp_info_cols, grep("num_overlaps", colnames(unique_enh_summary_df)))
        dists_idx <- c(snp_info_cols, grep("dists", colnames(unique_enh_summary_df)))

        ## first analyze with respect to FANTOM5 cell type / tissues
        ## look at the distributions across cell types:
        overlap_counts <- melt(unique_enh_summary_df[,num_overlap_idx],
                               id.vars=colnames(unique_enh_summary_df)[snp_info_cols],
                               variable.name="data_source")
        ## add a class column to this
        overlap_counts$enh_class <- fantom5_category_df$Class[summary_category_idx[as.vector(overlap_counts$data_source)]]

        ## ----------------------
        ## plot the number of unique SNPs overlapping enhancers from each tissue / cell type
        overlap_sums <- ddply(overlap_counts, .(data_source, enh_class), summarize,
                              enh_snp_overlaps = sum(!duplicated(rsID[value > 0])))
                
        ## get rid of the 'num_overlaps'
        overlap_sums$data_source <- gsub("_num_overlaps", "", overlap_sums$data_source)

        ## make sure the data sources are sorted by tissue category
        ## order this decreasing because we flip it
        overlap_sum_order <- order(overlap_sums$enh_class, overlap_sums$data_source, decreasing=T)
        overlap_sums$data_source <- factor(overlap_sums$data_source, ordered=T,
                                           levels=overlap_sums$data_source[overlap_sum_order])

        ## get the text color for this plot
        class_text_col <- category_colors[overlap_sums$enh_class[overlap_sum_order]]

        tiss_limit_max <- round_any(max(overlap_sums$enh_snp_overlaps)+1, 5, f=ceiling)
        
        make_graphic(paste0(outdir, 'plots/', prefix, '_enh_overlap_nums_bytiss_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"),
                     height_ratio=1.5)
        print(ggplot(overlap_sums,
                     aes(x=data_source, y=enh_snp_overlaps, fill=enh_class)) +
              cat_col_scale + 
              xlab("Cell or tissue type") + ylab("Number of overlapping variants") +
              coord_flip() + 
              scale_y_continuous(breaks=seq(0, tiss_limit_max, by=5), limits=c(0, tiss_limit_max), expand=c(0, 0)) + 
              theme_bw() + geom_bar(stat="identity", position="stack") +
              guides(fill = guide_legend(title="Tissue category", ncol=1)) +
              theme(legend.position="right",
                    axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE*.5),
                    axis.text.y = element_text(colour=class_text_col, size=AXIS_TEXT_Y_SIZE*.5),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                     
                    legend.text = element_text(size=LEGEND_TEXT_SIZE)) +
              plot_title("Numbers of enhancer overlapping variants across tissue types", r2_thresh, dist_thresh, out_subtitle))
        dev.off()

        ## ----------------------        
        ## look at the distributions of the numbers of SNPs within each data source in each
        ## tissue category

        ## need to make the class factor work right:
        overlap_sums$enh_class <- factor(overlap_sums$enh_class, ordered=T,
                                               levels=sort(unique(overlap_sums$enh_class), decreasing=T))

        make_graphic(paste0(outdir, 'plots/', prefix, '_enh_overlap_num_distribution_byclass_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"),
                     height_ratio=1.5)
        print(ggplot(overlap_sums,
                     aes(x=enh_class, y=enh_snp_overlaps, fill=enh_class)) +
              cat_col_scale + 
              xlab("Tissue category") + ylab("Number of overlapping variants") +
              coord_flip() + 
              scale_y_continuous(breaks=seq(0, tiss_limit_max, by=2), limits=c(0, tiss_limit_max), expand=c(0.01, 0)) + 
              theme_bw() + geom_boxplot() +
              theme(legend.position="none",
                    axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    axis.title = element_text(size=TITLE_SIZE*0.75),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
              plot_title("Distributions of numbers of enhancer-overlapping variants", r2_thresh, dist_thresh, out_subtitle))
        dev.off()        
        
        ## ----------------------        
        ## also make one summed by class
        ## we want to count the unique SNPs that have overlap in any class
        overlap_class_sums <- ddply(overlap_counts, .(enh_class), summarize,
                                    class_snp_overlaps = sum(!duplicated(rsID[value > 0])))

        overlap_class_sums$enh_class <- factor(overlap_class_sums$enh_class, ordered=T,
                                               levels=sort(overlap_class_sums$enh_class, decreasing=T))

        class_limit_max <- round_any(max(overlap_class_sums$class_snp_overlaps), 20, f=ceiling)
        
        make_graphic(paste0(outdir, 'plots/', prefix, '_enh_overlap_nums_bytiss_class_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"),
                     height_ratio=1.5, width_ratio = 1.2)
        print(ggplot(overlap_class_sums,
                        aes(x=enh_class, y=class_snp_overlaps, fill=enh_class)) +
              cat_col_scale + 
              xlab("Tissue category") + ylab("Number of overlapping variants") +
              scale_y_continuous(breaks=seq(0, class_limit_max, by=10),
                                 limits=c(0, class_limit_max), expand=c(0, 0)) +
              coord_flip() + 
              theme_bw() + geom_bar(stat="identity", position="stack") +
              theme(legend.position="none",
                    axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
              plot_title("Numbers of enhancer overlaps across tissue categories", r2_thresh, dist_thresh, out_subtitle))
        dev.off()        
        
        ## proportions of all SNPs that overlap an enhancer in each tissue/cell type
        overlap_props <- melt(ddply(overlap_counts, .(data_source, enh_class), function(x) {
            overlap_prop = sum(!duplicated(x$rsID[x$value > 0])) / nrow(x)
            data.frame(overlap_prop = overlap_prop, non_overlap_prop = 1 - overlap_prop)
        }), id.vars=c('data_source', 'enh_class'))
        
        ## make some nicer looking names
        overlap_props$data_source <- gsub("_num_overlaps", "", overlap_props$data_source)

        ## make sure the data sources are sorted by tissue category
        ## order this decreasing because we flip it
        ## we need the unique tissue and cell types, ordered by class
        ordered_uniq_sources_dec <- unique(overlap_props$data_source[order(overlap_props$enh_class, overlap_props$data_source, decreasing=T)])
        overlap_props$data_source <- factor(overlap_props$data_source, ordered=T,
                                           levels=ordered_uniq_sources_dec)

        ## we want to color the names according to their tissue classes, so we need to match
        ordered_classes <- overlap_props$enh_class[match(ordered_uniq_sources_dec, overlap_props$data_source)]
        ## now we can get the colors by matching the classes
        plot_cols <- category_colors[ordered_classes]

        ## also, we want to have a way other than text color to show the categories!        
        make_graphic(paste0(outdir, 'plots/', prefix, '_enh_overlap_props_bytiss_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"),
                     height_ratio = 1.5)
        print(ggplot(overlap_props,
                     aes(x=data_source, y=value, fill=variable)) +
              xlab("Cell or tissue type") + ylab("Proportion of overlapping variants") +
              theme_bw() + geom_bar(stat="identity", position="stack") + coord_flip() +
              theme(legend.position="bottom",
                    axis.text.x = element_text(size=AXIS_TEXT_X_SIZE*.5),
                    axis.text.y = element_text(colour=plot_cols, size=AXIS_TEXT_Y_SIZE*.5),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                     
                    legend.text = element_text(size=LEGEND_TEXT_SIZE), 
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
              plot_title("Proportions of SNPs with enhancer overlap across tissue types", r2_thresh, dist_thresh, out_subtitle, strwrap_width=30))
        dev.off()

        ## also make a proportion plot by class
        ## i want to know what proportion of the total SNPs had ANY overlap!
        class_overlap_props <- melt(ddply(overlap_counts, .(enh_class), function(x) {
            num_uniq_snps <- nrow(unique(x[,c('chr', 'rsID', 'pos')]))            
            num_uniq_overlap_snps <- nrow(unique(x[x$value > 0, c("chr", "rsID", "pos")]))
            overlap_prop <- num_uniq_overlap_snps / num_uniq_snps
            data.frame(overlap_prop = overlap_prop, non_overlap_prop = 1 - overlap_prop)
        }), id.vars=c('enh_class'))

        class_overlap_props$enh_class <- factor(class_overlap_props$enh_class, ordered=T,
                                                levels=unique(sort(class_overlap_props$enh_class, decreasing=T)))
        
        make_graphic(paste0(outdir, 'plots/', prefix, '_enh_overlap_props_bytiss_class_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"),
                     height_ratio = 1.5)
        print(ggplot(class_overlap_props,
                     aes(x=enh_class, y=value, fill=variable)) +
              xlab("Tissue category") + ylab("Proportion of overlapping variants") +
              theme_bw() + geom_bar(stat="identity", position="stack") + coord_flip() +
              ## we can just use the class color vector
              theme(legend.position="bottom",
                    axis.text.x = element_text(size=AXIS_TEXT_X_SIZE),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE, colour=category_colors[class_overlap_props$enh_class]),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                    legend.title = element_text(size=TITLE_SIZE*0.75), 
                    legend.text = element_text(size=LEGEND_TEXT_SIZE)) +              
              plot_title("Proportions of SNPs with enhancer overlap across tissue categories", r2_thresh, dist_thresh, out_subtitle, strwrap_width=30))
        dev.off()
        
        ## make a heatmap of tag SNP - tissue type enrichment
        tag_snp_tissue_counts <- ddply(overlap_counts, .(tag_name, data_source, enh_class), function(x) {
            num_snp_overlaps <- sum(x$value)
            total_overlap_snps <- sum(overlap_counts$value[overlap_counts$tag_name == unique(x$tag_name)], na.rm=T)
            return(data.frame(num_snp_overlaps = num_snp_overlaps,
                              prop_snp_overlaps = ifelse(total_overlap_snps > 0,
                                  sum(x$value) / total_overlap_snps, 0)))
        })

        ## make some nicer looking names
        tag_snp_tissue_counts$data_source <- gsub("_num_overlaps", "", tag_snp_tissue_counts$data_source)
        tag_snp_tissue_counts$data_source <- gsub("_", " ", tag_snp_tissue_counts$data_source)

        ## first, make sure the data source is ordered:
        ## this time, it is in increasing order
        ordered_uniq_sources_inc <- unique(tag_snp_tissue_counts$data_source[order(tag_snp_tissue_counts$enh_class, tag_snp_tissue_counts$data_source)])

        tag_snp_tissue_counts$data_source <- factor(tag_snp_tissue_counts$data_source, ordered=T,
                                                    levels=ordered_uniq_sources_inc)

        ## get the colors for the axis text labels
        ordered_heatmap_classes <- tag_snp_tissue_counts$enh_class[match(ordered_uniq_sources_inc, tag_snp_tissue_counts$data_source)]
        heatmap_text_cols <- category_colors[ordered_heatmap_classes]
        
        ## raw count heatmap
        ## define a color scheme (from blues)
        #heatmap_col_lims <- colorRampPalette(brewer.pal(9,"Blues"))(2)
        ## white to blue
        heatmap_col_lims <- c("#FFFFFF", "#08306B")

        ## add annotation of tag regions without rsIDs
        tag_snp_tissue_counts$tag_no_rsid <- gsub(":rs.*", "", tag_snp_tissue_counts$tag_name)
        
        make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_enh_count_heatmap_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 2.0)                   
        print(ggplot(tag_snp_tissue_counts, aes_string(x="data_source", y=TAG_VAR)) +
              geom_tile(aes(fill=num_snp_overlaps), colour=heatmap_col_lims[1]) +
              scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
              theme_bw() + 
              xlab("FANTOM5 Tissue/Cell Type") + ylab(TAG_LAB) +
              plot_title("Heatmap of number of linked SNPs overlapping FANTOM5 enhancers by tag region and tissue", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
              theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75, colour=heatmap_text_cols),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                 
                    legend.text = element_text(size=LEGEND_TEXT_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
        dev.off()
        
        ## log10 count heatmap
        make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_enh_log10_count_heatmap_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 2.0)                   
        print(ggplot(tag_snp_tissue_counts, aes_string(x="data_source", y=TAG_VAR)) +
              geom_tile(aes(fill=log10(num_snp_overlaps)), colour=heatmap_col_lims[1]) +
              scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
              #scale_fill_gradient(low="gray80", high="mediumblue", na.value=heatmap_col_lims[1], guide="colorbar") +
              theme_bw() + 
              xlab("FANTOM5 Tissue/Cell Type") + ylab(TAG_LAB) +
              plot_title("Heatmap of number of linked SNPs overlapping FANTOM5 enhancers by tag region and tissue", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
              theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75, colour=heatmap_text_cols),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                 
                    legend.text = element_text(size=LEGEND_TEXT_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
        dev.off()
        
        ## proportions
        make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_enh_prop_heatmap_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 2.0)                   
        print(ggplot(tag_snp_tissue_counts, aes_string(x="data_source", y=TAG_VAR)) +
              theme_bw() +
              geom_tile(aes(fill=prop_snp_overlaps), colour=heatmap_col_lims[1]) +
              scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
              xlab("FANTOM5 Tissue/Cell Type") + ylab(TAG_LAB) +
              plot_title("Heatmap of tissue proportions of enhancer-overlapping SNPs for each tag region", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
              theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75,
                        colour=heatmap_text_cols),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                     
                    legend.text = element_text(size=LEGEND_TEXT_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
        dev.off()

        ## make a heatmap of tag SNP - tissue category enrichment
        tag_snp_tiss_class_counts <- ddply(overlap_counts, .(tag_name, enh_class), function(x) {
            ## only count unique enhancers!
            num_snp_overlaps <- sum(x$value > 0)
            ## find the total number of SNPs that overlap from this tag region
            total_overlap_snps <- sum(overlap_counts$value[overlap_counts$tag_name == unique(x$tag_name)]>0, na.rm=T)
            return(data.frame(num_snp_overlaps = num_snp_overlaps,
                              prop_snp_overlaps = ifelse(total_overlap_snps > 0,
                                  num_snp_overlaps / total_overlap_snps, 0)))
        })

        tag_snp_tiss_class_counts$tag_no_rsid <- gsub(":rs.*", "", tag_snp_tiss_class_counts$tag_name)
        
        ## raw count heatmap
        make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_enh_count_heatmap_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 2.0)                   
        print(ggplot(tag_snp_tiss_class_counts, aes_string(x="enh_class", y=TAG_VAR)) +
              geom_tile(aes(fill=num_snp_overlaps), colour=heatmap_col_lims[1]) +
              scale_fill_gradient(name="# Overlapping SNPs", low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
              theme_bw() + 
              xlab("FANTOM5 Data Source Category") + ylab(TAG_LAB) +
              plot_title("Heatmap of number of linked SNPs overlapping FANTOM5 enhancers by tag region and tissue category", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
              theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE,
                        colour=category_colors[tag_snp_tiss_class_counts$enh_class]),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                     
                    legend.text = element_text(size=LEGEND_TEXT_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
        dev.off()        

        ## make this with black tissue category labels
        make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_enh_count_heatmap_black_labels_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 2.0)
        print(ggplot(tag_snp_tiss_class_counts, aes_string(x="enh_class", y=TAG_VAR)) +
              geom_tile(aes(fill=num_snp_overlaps), colour=heatmap_col_lims[1]) +
              scale_fill_gradient(name="# Overlapping SNPs", low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
              theme_bw() + 
              xlab("FANTOM5 Data Source Category") + ylab(TAG_LAB) +
              plot_title("Heatmap of number of linked SNPs overlapping FANTOM5 enhancers by tag region and tissue category", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
              theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                     
                    legend.text = element_text(size=LEGEND_TEXT_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
        dev.off()        

        ## also make one with greyscale color theme; we set the scale to start at '1' so we can
        ## see even the weak signals
        greyscale_col_lims <- c("#FFFFFF", "gray90", "black")
        
        make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_enh_count_heatmap_greyscale_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 2.0)                   
        print(ggplot(tag_snp_tiss_class_counts, aes_string(x="enh_class", y=TAG_VAR)) +
              geom_tile(aes(fill=num_snp_overlaps), colour=greyscale_col_lims[1]) +
              scale_fill_gradientn(colours=greyscale_col_lims,
                                   values=rescale(c(0, 1, max(tag_snp_tiss_class_counts$num_snp_overlaps))),                                   
                                   name="# Overlapping SNPs", na.value=greyscale_col_lims[1],
                                   guide="colorbar") +
              theme_bw() + 
              xlab("FANTOM5 Data Source Category") + ylab(TAG_LAB) +
              plot_title("Heatmap of number of linked SNPs overlapping FANTOM5 enhancers by tag region and tissue category", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
              theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                     
                    legend.text = element_text(size=LEGEND_TEXT_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
        dev.off()        
        
        ## log10 count heatmap
        make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_enh_log10_count_heatmap_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 2.0)                   
        print(ggplot(tag_snp_tiss_class_counts, aes_string(x="enh_class", y=TAG_VAR)) +
              geom_tile(aes(fill=log10(num_snp_overlaps)), colour=heatmap_col_lims[1]) +
              scale_fill_gradient(name="log10(# overlapping SNPs)", low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
              theme_bw() + 
              xlab("FANTOM5 Data Source Category") + ylab(TAG_LAB) +
              plot_title("Heatmap of number of linked SNPs overlapping FANTOM5 enhancers by tag region and tissue category", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
              theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE,
                        colour=category_colors[tag_snp_tiss_class_counts$enh_class]),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                     
                    legend.text = element_text(size=LEGEND_TEXT_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
        dev.off()        

        ## proportions
        make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_enh_prop_heatmap_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 2.0)                   
        print(ggplot(tag_snp_tiss_class_counts, aes_string(x="enh_class", y=TAG_VAR)) +
              geom_tile(aes(fill=prop_snp_overlaps), colour=heatmap_col_lims[1]) +
              scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
              theme_bw() + 
              xlab("FANTOM5 Data Source Category") + ylab(TAG_LAB) +
              plot_title("Heatmap of tissue class proportions of enhancer-overlapping SNPs for each tag region", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
              theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE,
                        colour=category_colors[tag_snp_tiss_class_counts$enh_class]),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                     
                    legend.text = element_text(size=LEGEND_TEXT_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
        dev.off()                
        
        ## now analyse with respect to SNPs (and across tag SNPs)    
        num_overlap_cols <- colnames(unique_enh_summary_df)[grep("num_overlaps", colnames(unique_enh_summary_df))]

        ## count the number of enhancers each SNP overlaps (across tissues)
        snp_overlap_counts <- cbind(unique_enh_summary_df[,c("chr", "rsID", "pos", "tag_name")],
                                    enh_count = rowSums(unique_enh_summary_df[,num_overlap_cols]))

        ## look at numbers of SNPs with any overlap, by tag region
        num_overlapping_snps <- ddply(snp_overlap_counts, .(tag_name), summarize,
                                      enh_overlapping_snps = sum(enh_count > 0),
                                      total_snps = length(enh_count))
        num_overlapping_snps$tag_no_rsid <- gsub(":rs.*", "", num_overlapping_snps$tag_name)
        
        make_graphic(paste0(outdir, 'plots/', prefix, '_num_overlapping_snps_per_tagregion_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"))                   
        print(ggplot(num_overlapping_snps, aes_string(x=TAG_VAR, y="enh_overlapping_snps", fill=TAG_VAR)) +
              xlab(TAG_LAB) + ylab("Number of overlapping variants") +
              scale_fill_hue(h=c(180, 270)) +
              scale_y_continuous(breaks=seq(0, max(num_overlapping_snps$enh_overlapping_snps)*1.1, by=5)) + 
              theme_bw() + geom_bar(stat="identity", position="stack") +
              theme(legend.position="none",
                    axis.text.x=element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
              plot_title("Numbers of SNPs with any enhancer overlaps per tag SNP region", r2_thresh, dist_thresh, out_subtitle))
        dev.off()

        ## also look at proportions of SNPs with any overlap by tag SNP
        snp_overlap_props <- melt(data.frame(
            tag_name = num_overlapping_snps$tag_name,
            tag_no_rsid = num_overlapping_snps$tag_no_rsid,
            overlap=with(num_overlapping_snps, enh_overlapping_snps / total_snps),
            non_overlap = 1 - with(num_overlapping_snps, enh_overlapping_snps / total_snps)), 
            id.vars=c('tag_name', 'tag_no_rsid'))

        make_graphic(paste0(outdir, 'plots/', prefix, '_prop_overlapping_snps_per_tagregion_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"))                   
        print(ggplot(snp_overlap_props, aes_string(x=TAG_VAR, y="value", fill="variable")) +
              xlab(TAG_LAB) + ylab("Proportion of overlapping variants") +
              theme_bw() + geom_bar(stat="identity", position="stack") +
              scale_y_continuous(breaks=seq(0, 1, by=.1)) + 
              theme(legend.position="bottom",
                    axis.text.x=element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                    legend.title = element_text(size=TITLE_SIZE*0.75),                     
                    legend.text=element_text(size=LEGEND_TEXT_SIZE)) +
              plot_title("Proportions of SNPs with any enhancer overlaps per tag SNP region", r2_thresh, dist_thresh, out_subtitle))
        dev.off()
        
    }

}

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/fantom5_overlap/")
## enh_overlap_type <- "locus"
## enh_window <- 1000
## fantom5_class_file <- "/home/alexaml/data/FANTOM5/Enhancers/fantom5_classes.txt"
## enh_summary_flag <- TRUE
