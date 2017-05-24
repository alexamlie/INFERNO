## analyze_fantom5_eqtl_chromHMM_overlap_function.R
## alex amlie-wolf 03/29/2016
## part of the auto pipeline analysis, this script does integrative analysis of FANTOM5
## enhancer, GTEx eQTL, and Roadmap epigenomic (using chromHMM) data

## a function to create the SNP support heatmap (gets used a lot with different DFs)
## overloads the outdir and prefix arguments for flexibility
## also takes the color schemes that are defined beforehand (so that tag regions have
## consistent color schemes)
snp_support_heatmap <- function(support_df, outdir, prefix, this_tag_name, support_text_cols, tag_name_cols, snp_text_cols, snp_text_size=4) {
    ## ---------------------------
    ## we now have many color possibilities, so we create general categories / palettes to use
    ## this requires a bit of manipulation; we have to find the number of combinations in each
    ## category and make palettes that way. the general idea is that each of the three data
    ## sources is represented by a primary color, and their combinations form spectra of
    ## secondary colors

    ## now do the ones that we need palettes for (because i don't want to have to enumerate
    ## every possible state)
    ## eQTL + eRNA + HMM is range of reddish-brown

    ## count all the unique types of support
    unique_supports <- sort(unique(support_df$support))
    unique_merged_hmm_supports <- sort(unique(support_df$merged_hmm_support))

    eqtl_erna_hmm_supports <- grep("^eQTL\\+eRNA Enh\\+HMM", unique_supports, value=T)
    eqtl_erna_hmm_colors <- colorRampPalette(c("#521C0A", "#8F4E38"))(length(eqtl_erna_hmm_supports))
    names(eqtl_erna_hmm_colors) <- eqtl_erna_hmm_supports

    ## eQTL + HMM is green (yellow + blue)
    eqtl_hmm_supports <- grep("^eQTL\\+HMM", unique_supports, value=T)
    eqtl_hmm_colors <- colorRampPalette(c("#88BDAD", "#1E5645"))(length(eqtl_hmm_supports))
    names(eqtl_hmm_colors) <- eqtl_hmm_supports

    ## eRNA + HMM is purple
    erna_hmm_supports <- grep("^eRNA Enh\\+HMM", unique_supports, value=T)
    erna_hmm_colors <- colorRampPalette(c("#660BAB", "#3D0568"))(length(erna_hmm_supports))
    names(erna_hmm_colors) <- erna_hmm_supports

    ## now we just need the range of HMM colors (range of darker blues)
    hmm_supports <- grep("^HMM .*\\+", unique_supports, value=T)
    hmm_combo_colors <- colorRampPalette(c("#152B55", "#061739"))(length(hmm_supports))
    names(hmm_combo_colors) <- hmm_supports

    support_heatmap_cols <- c("No support"="#FFFFFF", eqtl_color, erna_color, hmm_colors,
                              eqtl_erna_color, eqtl_erna_hmm_colors, eqtl_hmm_colors,
                              erna_hmm_colors, hmm_combo_colors)

    make_graphic(paste0(outdir, prefix, '_roadmap_fantom5_gtex_support_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 3.5, height_ratio = 2.0, type='pdf')
    print(ggplot(support_df, aes(x=rsID, y=tissue_class)) +
          geom_tile(aes(fill=factor(support))) +
          scale_fill_manual(name="Sources of functional support",
                            values = support_heatmap_cols) +
          ## add invisible "points" just to get the legend for the tag SNP region
          geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
          scale_x_discrete(position="top") + 
          scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
          theme_bw() +
          xlab("Linked SNP") + ylab("Tissue Category") +
          plot_title(paste0("Heatmap of sources of functional support for SNPs across ", this_tag_name, " tag region(s) and tissue categories"), r2_thresh, dist_thresh, out_subtitle) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10), ncol=1),
                 fill = guide_legend(override.aes=list(size=10), ncol=1)) +
          theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=snp_text_size, colour=snp_text_cols),
                axis.text.y = element_text(hjust=1, size=30, colour=support_text_cols),
                axis.title = element_text(size=30),
                title = element_text(size=25), legend.text = element_text(size=20),
                legend.position="right", plot.title=element_text(hjust = 0.5)))
    dev.off()

    ## also do the present tissues only
    present_tissues <- as.vector(sort(unique(support_df$tissue_class[support_df$support!="No support"])))
    present_tissue_cols <- category_colors[present_tissues]

    present_tiss_support_df <- support_df[support_df$tissue_class %in% present_tissues,]
    present_tiss_support_df$support <- factor(present_tiss_support_df$support, ordered=T,
                                          levels=c(grep("\\+|No support", unique_supports, val=T, inv=T),
                                              grep("\\+", unique_supports, val=T, inv=F),
                                              "No support"))

    ## only make this plot if there ARE any present tissues
    if (length(present_tissues) > 0) {
        make_graphic(paste0(outdir, prefix, '_roadmap_fantom5_gtex_support_heatmap_present_tissues_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 3.5, height_ratio = 2.0, type='pdf')
        print(ggplot(present_tiss_support_df, aes(x=rsID, y=tissue_class)) +
              geom_tile(aes(fill=factor(support))) +
              scale_fill_manual(name="Sources of functional support",
                                values = support_heatmap_cols) +
              ## add invisible "points" just to get the legend for the tag SNP region
              geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
              scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
              scale_x_discrete(position="top") +               
              theme_bw() +
              xlab("Linked SNP") + ylab("Tissue Category") +
              plot_title(paste0("Heatmap of sources of functional support for SNPs across ", this_tag_name, " tag region(s) and tissue categories, present tissues only"), r2_thresh, dist_thresh, out_subtitle) +
              guides(colour = guide_legend(override.aes=list(shape=15, size=10), ncol=1),
                     fill = guide_legend(override.aes=list(size=10), ncol=1)) +
              theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=snp_text_size, colour=snp_text_cols),
                    axis.text.y = element_text(hjust=1, size=30, colour=support_text_cols),
                    axis.title = element_text(size=30),
                    title = element_text(size=25), legend.text = element_text(size=20),
                    legend.position="right", plot.title=element_text(hjust = 0.5)))
        dev.off()
    }
    
    ## ----------------------
    ## do the same thing but using the merged enhancer states. the color is simpler this time
    support_merged_heatmap_cols <- c("No support"="#FFFFFF", eqtl_color, erna_color,
                                     merge_hmm_color, eqtl_erna_color,
                                     eqtl_merged_hmm_color, erna_merged_hmm_color,
                                     eqtl_erna_merged_hmm_color)

    make_graphic(paste0(outdir, prefix, '_roadmap_fantom5_gtex_support_heatmap_merged_hmm_states_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 3.5, height_ratio = 2.0, type='pdf')
    print(ggplot(support_df, aes(x=rsID, y=tissue_class)) +
          geom_tile(aes(fill=factor(merged_hmm_support))) +
          scale_fill_manual(name="Sources of functional support",
                            values = support_merged_heatmap_cols) +
          ## add invisible "points" just to get the legend for the tag SNP region
          geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
          scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
          scale_x_discrete(position="top") +           
          theme_bw() +
          xlab("Linked SNP") + ylab("Tissue Category") +
          plot_title(paste0("Heatmap of sources of functional support for SNPs across ", this_tag_name, " tag region(s) and tissue categories, with collapsed HMM enhancer states"), r2_thresh, dist_thresh, out_subtitle) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10), ncol=1),
                 fill = guide_legend(override.aes=list(size=10), ncol=1)) +
          ## +1e-09 was the 'hack' number to use for the angle adjustment
          theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=snp_text_size, colour=snp_text_cols),
                axis.text.y = element_text(hjust=1, size=30, colour=support_text_cols),
                axis.title = element_text(size=30),
                title = element_text(size=25), legend.text = element_text(size=20),
                legend.position="right", plot.title=element_text(hjust = 0.5)))
    dev.off()

    ## also make one with black labels instead of the tissue category color scheme, and other changes
    make_graphic(paste0(outdir, prefix, '_roadmap_fantom5_gtex_support_heatmap_merged_hmm_states_black_labels_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 3.5, height_ratio = 2.0, type='pdf')
    print(ggplot(support_df, aes(x=rsID, y=tissue_class)) +
          geom_tile(aes(fill=factor(merged_hmm_support)), color="grey60") +
          scale_fill_manual(name="Sources of functional support",
                            values = support_merged_heatmap_cols) +
          ## add invisible "points" just to get the legend for the tag SNP region
          geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
          scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
          scale_x_discrete(position="top") +           
          theme_bw() +
          xlab("Linked SNP") + ylab("Tissue Category") +
          plot_title(paste0("Heatmap of sources of functional support for SNPs across ", this_tag_name, " tag region(s) and tissue categories, with collapsed HMM enhancer states"), r2_thresh, dist_thresh, out_subtitle) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10), ncol=1),
                 fill = guide_legend(override.aes=list(size=10), ncol=1)) +
          ## +1e-09 was the 'hack' number to use for the angle adjustment
          theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=snp_text_size, colour=snp_text_cols),
                axis.text.y = element_text(hjust=1, size=30),
                axis.title = element_text(size=30),
                title = element_text(size=25), legend.text = element_text(size=20),
                legend.position="right", plot.title=element_text(hjust = 0.5)))
    dev.off()
    
    if (length(present_tissues) > 0) {
        make_graphic(paste0(outdir, prefix, '_roadmap_fantom5_gtex_support_heatmap_merged_hmm_states_present_tissues_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 3.5, height_ratio = 2.0, type='pdf')
        print(ggplot(support_df[support_df$tissue_class %in% present_tissues,],
                     aes(x=rsID, y=tissue_class)) +
              geom_tile(aes(fill=factor(merged_hmm_support))) +
              scale_fill_manual(name="Sources of functional support",
                                values = support_merged_heatmap_cols) +
              ## add invisible "points" just to get the legend for the tag SNP region
              geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
              scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
              scale_x_discrete(position="top") +               
              theme_bw() +
              xlab("Linked SNP") + ylab("Tissue Category") +
              plot_title(paste0("Heatmap of sources of functional support for SNPs across ", this_tag_name, " tag region(s) and tissue categories, with collapsed HMM enhancer states and only present tissues"), r2_thresh, dist_thresh, out_subtitle) +
              guides(colour = guide_legend(override.aes=list(shape=15, size=10), ncol=1),
                     fill = guide_legend(override.aes=list(size=10), ncol=1)) +
              theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=snp_text_size, colour=snp_text_cols),
                    axis.text.y = element_text(hjust=1, size=30, colour=support_text_cols),
                    axis.title = element_text(size=30),
                    title = element_text(size=25), legend.text = element_text(size=20),
                    legend.position="right", plot.title=element_text(hjust = 0.5)))
        dev.off()

        ## also make a plot with black labels for the tissues
        make_graphic(paste0(outdir, prefix, '_roadmap_fantom5_gtex_support_heatmap_merged_hmm_states_present_tissues_black_labels_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 3.5, height_ratio = 2.0, type='pdf')
        print(ggplot(support_df[support_df$tissue_class %in% present_tissues,],
                   aes(x=rsID, y=tissue_class)) +
              geom_tile(aes(fill=factor(merged_hmm_support)), color="grey60") +
              scale_fill_manual(name="Sources of functional support",
                                values = support_merged_heatmap_cols) +
              ## add invisible "points" just to get the legend for the tag SNP region
              geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
              scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
              scale_x_discrete(position="top") +               
              theme_bw() +
              xlab("Linked SNP") + ylab("Tissue Category") +
              plot_title(paste0("Heatmap of sources of functional support for SNPs across ", this_tag_name, " tag region(s) and tissue categories, with collapsed HMM enhancer states and only present tissues"), r2_thresh, dist_thresh, out_subtitle) +
              guides(colour = guide_legend(override.aes=list(shape=15, size=10), ncol=1),
                     fill = guide_legend(override.aes=list(size=10), ncol=1)) +
              theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=snp_text_size, colour=snp_text_cols),
                    axis.text.y = element_text(hjust=1, size=30),
                    axis.title = element_text(size=30),
                    title = element_text(size=25), legend.text = element_text(size=20),
                    legend.position="right", plot.title=element_text(hjust = 0.5)))
        dev.off()        
    }
}

analyze_fantom5_eqtl_chromHMM_overlap <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh, enh_overlap_type, enh_window, fantom5_class_file, gtex_class_file, roadmap_class_file) {
    ## ---------------------------------------------------
    ## read in all the original data and simplify it
    ## ---------------------------------------------------
    dir.create(paste0(outdir, "/tables/"), F, T)
    dir.create(paste0(outdir, "/plots/"), F, T)

    ## ----------------------
    ## first read in all the LD SNPs so that we don't miss any
    ld_stats_file <- paste0(datadir, "/ld_expansion/", prefix, "_", r2_thresh,
                            "_ld_cutoff_snps_within_", dist_thresh, ".txt")
    ld_stats_df <- read.table(ld_stats_file, header=T, sep="\t", quote="", as.is=T)

    ## add a column to the LD stats dataframe
    ld_stats_df$tag_no_rsid <- gsub(":rs.*", "", ld_stats_df$tag_name)

    ## ----------------------
    ## read in the FANTOM5 data
    ## read in the enhancer overlap data file
    if (enh_overlap_type == "midpoint") {
        enh_overlap_file <- paste0(datadir, '/fantom5_overlap/',
                                   prefix, "_", r2_thresh,
                                   "_ld_cutoff_snps_within_", dist_thresh,
                                   "_", enh_window, "bp_around_midpoint_enh_overlaps.txt")
        enh_overlap_df <- read.table(enh_overlap_file, header=T, sep="\t", quote="", as.is=T)
    } else if (enh_overlap_type == "locus") {
        enh_overlap_file <- paste0(datadir, '/fantom5_overlap/',
                                   prefix, "_", r2_thresh,
                                   "_ld_cutoff_snps_within_", dist_thresh,
                                   "_", enh_window, "bp_around_orig_locus_enh_overlaps.txt")
        enh_overlap_df <- read.table(enh_overlap_file, header=T, sep="\t", quote="", as.is=T)
    }

    ## read in the tissue categories
    fantom5_category_df <- read.table(fantom5_class_file, header=T, sep="\t", quote="", as.is=T)
    ## add a column for easier matching
    fantom5_category_df$enh_source <- gsub("_expressed_enhancers.bed", "", fantom5_category_df$FANTOM5.File)

    ## add class columns to the enhancer data:
    enh_overlap_df$enh_class <- fantom5_category_df$Class[match(enh_overlap_df$enh_source,
                                                               fantom5_category_df$enh_source)]
    ## clean up the tissue names
    enh_overlap_df$enh_source <- gsub("CL:[0-9]*_|UBERON:[0-9]*_", "", enh_overlap_df$enh_source)
    ## add a column for tag without rsid
    enh_overlap_df$tag_no_rsid <- gsub(":rs.*", "", enh_overlap_df$tag_name)
    ## TODO: ADD MORE COLUMNS HERE?

    ## get the columns that describe each unique SNP: (this applies to each data source)
    uniq_snp_cols <- colnames(enh_overlap_df)[!(colnames(enh_overlap_df) %in% c("enh_source", "enh_class", "enh_chr", "enh_start", "enh_end", tagsnp_cols))]
    ## now summarize the data for each unique SNP
    enh_uniq_snp_df <- ddply(enh_overlap_df, uniq_snp_cols, summarize,
                             num_enh_tissues = length(unique(enh_source)),
                             enh_tissues = paste(sort(unique(enh_source)), collapse=","),
                             num_enh_classes = length(unique(enh_class)),
                             enh_classes = paste(sort(unique(enh_class)), collapse=","))

    ## ----------------------
    ## read in the eQTL data
    eqtl_overlap_file <- paste0(datadir, '/gtex_eqtl_overlap/',
                                prefix, "_", r2_thresh, "_ld_cutoff_snps_within_",
                                dist_thresh, "_eqtl_overlaps.txt")
    eqtl_overlap_df <- read.table(eqtl_overlap_file, header=T, sep="\t", quote="", as.is=T)

    ## read in the tissue categories
    gtex_category_df <- read.table(gtex_class_file, header=T, sep="\t", quote="", as.is=T)
    ## make it easier to match
    gtex_category_df$tissue <- gsub("_Analysis.snpgenes", "", gtex_category_df$GTEx.Data)

    ## add a category column
    eqtl_overlap_df$eqtl_class <- gtex_category_df$Class[match(eqtl_overlap_df$tissue, gtex_category_df$tissue)]
    ## add a column for tag without rsid
    eqtl_overlap_df$tag_no_rsid <- gsub(":rs.*", "", eqtl_overlap_df$tag_name)

    ## get the unique SNPs from the eQTL, with some summaries
    eqtl_uniq_snp_df <- ddply(eqtl_overlap_df, uniq_snp_cols, summarize,
                              num_eqtl_genes = length(unique(gene_name)),
                              eqtl_genes = paste(sort(unique(gene_name)), collapse=","),
                              num_eqtl_tissues = length(unique(tissue)),
                              eqtl_tissues = paste(sort(unique(tissue)), collapse=","),
                              num_eqtl_classes = length(unique(eqtl_class)),
                              eqtl_classes = paste(sort(unique(eqtl_class)), collapse=","))

    ## merge the enhancer and eQTL data (finding overlapping SNPs)
    enh_eqtl_overlap_snp_df <- join(enh_uniq_snp_df, eqtl_uniq_snp_df, type="inner")

    ## ----------------------
    ## read in the roadmap chromHMM data
    roadmap_state_file <- paste0(datadir, '/roadmap_chromhmm_states/', prefix, "_", r2_thresh,
                                "_ld_cutoff_snps_within_", dist_thresh, "_roadmap_chromHMM_states.txt")
    roadmap_state_df <- read.table(roadmap_state_file, header=T, sep="\t", quote="", as.is=T,
                                   stringsAsFactors = F)

    ## read in the roadmap category file
    roadmap_category_df <- read.table(roadmap_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")
    ## rename the EID column to be easier to match
    colnames(roadmap_category_df)[3] <- "EID"

    ## manipulate the roadmap state data to be more easily analyzed
    ## add a column without the rsID
    roadmap_state_df$tag_no_rsid <- gsub(":rs.*", "", roadmap_state_df$tag_name)

    ## get a table with only unique hits, also rename the state columns
    uniq_roadmap_snp_cols <- colnames(roadmap_state_df)[!(colnames(roadmap_state_df) %in% tagsnp_cols)]
    uniq_roadmap_snp_state_df <- ddply(roadmap_state_df, uniq_roadmap_snp_cols, function(x) {
        apply(x[,!(colnames(x) %in% uniq_roadmap_snp_cols)], 2, paste, collapse=",")
    })#, .progress="text")

    ## rename the state columns to reflect the actual data source
    ## first get an index to these columns
    eid_state_col_idx <- grepl("_state", colnames(uniq_roadmap_snp_state_df))
    eid_orig_cols <- gsub("_state", "", colnames(uniq_roadmap_snp_state_df)[eid_state_col_idx])
    ## now grab the real names
    eid_name_cols <- roadmap_category_df$Standardized.Epigenome.name[match(eid_orig_cols, roadmap_category_df$EID)]
    ## replace the IDs with the real names
    colnames(uniq_roadmap_snp_state_df)[eid_state_col_idx] <- eid_name_cols

    ## also save a vector for the non-tag SNP-related columns
    roadmap_notag_snp_cols <- !(colnames(uniq_roadmap_snp_state_df) %in% tagsnp_cols)

    ## melt the table to make it easier to manipulate
    melt_roadmap_cols <- colnames(uniq_roadmap_snp_state_df)[!eid_state_col_idx & roadmap_notag_snp_cols]

    melt_roadmap_df <- melt(uniq_roadmap_snp_state_df, id.vars=melt_roadmap_cols,
                            measure.vars=colnames(uniq_roadmap_snp_state_df)[eid_state_col_idx],
                            variable.name="tissue", value.name="state")
    ## make the tissue and state into character and not factors
    melt_roadmap_df$tissue <- as.character(melt_roadmap_df$tissue)
    melt_roadmap_df$state <- as.character(melt_roadmap_df$state)

    ## add a tissue class column to this data structure
    melt_roadmap_df$class <- roadmap_category_df$Class[match(melt_roadmap_df$tissue, roadmap_category_df$Standardized.Epigenome.name)]

    ## now grab the enhancer states:
    roadmap_enhancer_snps <- melt_roadmap_df[melt_roadmap_df$state %in% c("6_EnhG", "7_Enh", "12_EnhBiv"),]

    ## summarize these per snp:
    roadmap_uniq_enh_snp_df <- ddply(roadmap_enhancer_snps, uniq_snp_cols, summarize,
                                     num_hmm_tissues = length(unique(tissue)),
                                     hmm_tissues = paste(sort(unique(tissue)), collapse=","),
                                     num_hmm_enh_states = length(unique(state)),
                                     hmm_enh_states = paste(sort(unique(state)), collapse=","),
                                     num_hmm_classes = length(unique(class)),
                                     hmm_classes = paste(sort(unique(class)), collapse=","))   
    
    ## ----------------------
    ## find the SNPs that have all 3 data sources:
    roadmap_enh_eqtl_overlap_snps <- join(roadmap_uniq_enh_snp_df, enh_eqtl_overlap_snp_df, type="inner")
    ## add a column for the classes that overlap between all 3 (if any)
    roadmap_enh_eqtl_overlap_snps$overlap_classes <- apply(roadmap_enh_eqtl_overlap_snps, 1, function(x) {
        hmm_classes <- strsplit(x["hmm_classes"], ",")[[1]]
        enh_classes <- strsplit(x["enh_classes"], ",")[[1]]
        eqtl_classes <- strsplit(x["eqtl_classes"], ",")[[1]]
        common_classes <- intersect(hmm_classes, intersect(enh_classes, eqtl_classes))
        if(length(common_classes) > 0) {
            return(paste(common_classes, collapse=","))
        } else {
            return("None")
        }
    })
        
    ## write out the full overlap file
    overlap_snp_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                               "_ld_cutoff_snps_within_", dist_thresh,
                               "_", enh_overlap_type,
                               "_fantom5_eqtl_roadmap_enh_states_overlap_gene_targets.txt")
    write.table(roadmap_enh_eqtl_overlap_snps, overlap_snp_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## ----------------------
    ## figure out the average number of tissue classes affected by each SNP in each data type
    cat("Average number of tissue classes affected by FANTOM5 enhancer overlapping SNPs:\n")
    cat("Mean:", mean(enh_uniq_snp_df$num_enh_classes), "median:", median(enh_uniq_snp_df$num_enh_classes), "\n")
    cat("Average number of tissue classes affected by Roadmap ChromHMM enhancer overlapping SNPs:\n")
    cat("Mean:", mean(roadmap_uniq_enh_snp_df$num_hmm_classes), "median:", median(roadmap_uniq_enh_snp_df$num_hmm_classes), "\n")
    cat("Average number of tissue classes affected by eQTL overlapping SNPs:\n")
    cat("Mean:", mean(eqtl_uniq_snp_df$num_eqtl_classes), "median:", median(eqtl_uniq_snp_df$num_eqtl_classes), "\n")
        
    
    ## ---------------------------------------------------
    ## make analysis plots
    ## ---------------------------------------------------
    ## ----------------------
    ## we want to know how many SNPs in each tag region are supported by each overlap type
    tagregion_support_counts <- melt(ddply(ld_stats_df, .(tag_name), function(x) {
        this_region_snps <- unique(x$rsID)
        enh_count <- sum(this_region_snps %in% enh_overlap_df$rsID)
        hmm_count <- sum(this_region_snps %in% roadmap_uniq_enh_snp_df$rsID)
        eqtl_count <- sum(this_region_snps %in% eqtl_uniq_snp_df$rsID)
        return(data.frame("FANTOM5 Enhancer" = enh_count, "Roadmap HMM" = hmm_count, "GTEx eQTL" = eqtl_count, stringsAsFactors=F))
    }), id.vars=c("tag_name"), variable.name="annotation", value.name="SNP_Count")
    max_snp_count <- round_any(max(tagregion_support_counts$SNP_Count)+1, accuracy=5, f=ceiling)

    ## write a specific break function for this
    snp_count_breaks <- function(x) {
        if(max(x) >= 100) {
            breaks <- c(0, 10, 25, 50, seq(100, max_snp_count, by=100))
        } else {
            breaks <- c(0, 1, 5, 10, 25, 50, 100)
        }
        breaks
    }
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_indiv_annotation_snp_counts_per_region_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"), height_ratio = 2.0)
    print(ggplot(tagregion_support_counts, aes(x=tag_name, y=SNP_Count, fill=tag_name)) +
          scale_fill_hue(h=c(180, 270)) +
          xlab("Tag SNP rsID and region") + ylab("Number of unique overlapping SNPs") +
          facet_grid(annotation ~ ., scales='free_y') + 
          theme_bw() + geom_bar(position="dodge", stat="identity") +
          scale_y_continuous(breaks=snp_count_breaks, trans="sqrt") +
          geom_text(aes(label=SNP_Count), position=position_dodge(0.9)) +
          ## scale_y_continuous(breaks=seq(0, max_snp_count, by=1),
          ##                    limits=c(0, max_snp_count), expand=c(0, 0)) +
          theme(legend.position="none", axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE), axis.title=element_text(size=TITLE_SIZE),
                plot.title=element_text(size=TITLE_SIZE, hjust=0.5), strip.text=element_text(size=15)) + 
          plot_title("Number of SNPs overlapping FANTOM5 Enhancers, Roadmap HMM Enhancers, and GTEx eQTLs", r2_thresh, dist_thresh, out_subtitle))
    dev.off()    
        
    ## ----------------------
    ## just look at the number of SNPs that are supported by all 3
    threeway_snp_count <- as.data.frame(table(roadmap_enh_eqtl_overlap_snps$tag_name[!duplicated(roadmap_enh_eqtl_overlap_snps[,c("rsID", "tag_name")])], dnn="tag_name"))
    max_snp_count <- round_any(max(threeway_snp_count$Freq)+1, accuracy=5, f=ceiling)

    make_graphic(paste0(outdir, 'plots/', prefix, '_hmm_fantom5_gtex_snp_counts_per_region_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(threeway_snp_count, aes(x=tag_name, y=Freq, fill=tag_name)) +
          scale_fill_hue(h=c(180, 270)) +
          xlab("Tag SNP rsID and region") + ylab("Number of unique SNPs") +
          theme_bw() + geom_bar(position="dodge", stat="identity") +
          scale_y_continuous(breaks=seq(0, max_snp_count, by=1),
                             limits=c(0, max_snp_count), expand=c(0, 0)) +
          theme(legend.position="none", axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE), axis.title=element_text(size=TITLE_SIZE),
                plot.title=element_text(size=TITLE_SIZE, hjust=0.5)) + 
          plot_title("Number of SNPs overlapping FANTOM5 Enhancers, Roadmap HMM Enhancers, and GTEx eQTLs", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## ----------------------
    ## get the supporting data sets across all SNPs and tissue categories
    ## for this, we are just looking for the data source - tissue category combinations that
    ## are present
    enh_snp_category_combos <- melt(unique(enh_overlap_df[,c("rsID", "pos", "tag_name", "enh_class")]),
                                    id.vars=c("rsID", "tag_name", "pos"),
                                    variable.name="data_source", value.name="tissue_class")
    eqtl_snp_category_combos <- melt(unique(eqtl_overlap_df[,c("rsID", "pos", "tag_name", "eqtl_class")]),
                                     id.vars=c("rsID", "tag_name", "pos"),
                                     variable.name="data_source", value.name="tissue_class")
    ## don't actually need to melt this one since we will use the state to index
    roadmap_enh_snp_category_combos <- unique(roadmap_enhancer_snps[,c("rsID", "pos", "tag_name", "state", "class")])
    ## rename these columns to match the other data frames
    colnames(roadmap_enh_snp_category_combos) <- c("rsID", "pos", "tag_name", "data_source", "tissue_class")

    ## get the supporting data sets per SNP and tissue category
    ## keep track of both the separate HMM states as well as the merged
    snp_support_df <- ddply(rbind(enh_snp_category_combos, eqtl_snp_category_combos, roadmap_enh_snp_category_combos),
                            .(rsID, pos, tag_name, tissue_class), function(x) {
                                ## return the amount and type of support
                                ## first check for individual data sources, then do combos
                                eqtl_present <- "eqtl_class" %in% x$data_source
                                enh_present <- "enh_class" %in% x$data_source
                                hmm_enh_present <- "7_Enh" %in% x$data_source
                                hmm_enhG_present <- "6_EnhG" %in% x$data_source
                                hmm_enh_biv_present <- "12_EnhBiv" %in% x$data_source

                                ## make a string describing the types of individual support:
                                support <- ""
                                ## also make one for the merged HMM states
                                merged_hmm_support <- ""
                                if (eqtl_present) {
                                    support <- paste(support, "eQTL", sep="+")
                                    merged_hmm_support <- paste(merged_hmm_support, "eQTL", sep="+") }
                                if(enh_present) {
                                    support <- paste(support, "eRNA Enh", sep="+")
                                    merged_hmm_support <- paste(merged_hmm_support, "eRNA Enh", sep="+") }
                                if(hmm_enh_present) {
                                    support <- paste(support, "HMM Enh", sep="+") }
                                if(hmm_enhG_present) {
                                    support <- paste(support, "HMM Genic Enh", sep="+") }
                                if(hmm_enh_biv_present) {
                                    support <- paste(support, "HMM Biv Enh", sep="+") }

                                if(hmm_enh_present || hmm_enhG_present || hmm_enh_biv_present) {
                                    merged_hmm_support <- paste(merged_hmm_support, "HMM Enh", sep="+") }

                                ## remove the separator at the beginning of the support strings
                                support <- gsub("^\\+", "", support)
                                merged_hmm_support <- gsub("^\\+", "", merged_hmm_support)

                                return(data.frame(support, merged_hmm_support, stringsAsFactors = F))
                            })

    rm(enh_snp_category_combos, eqtl_snp_category_combos, roadmap_enh_snp_category_combos)

    ## we also need to add annotations for all SNP-tissue combos that were not observed
    snp_support_df_full <- merge(expand.grid(rsID = unique(ld_stats_df$rsID),
                                             tissue_class=Reduce(union, list(fantom5_category_df$Class, gtex_category_df$Class, roadmap_category_df$Class)),
                                             stringsAsFactors = F),
                                 snp_support_df, all.x=T)

    ## fix all the NAs we generated
    snp_support_df_full$pos <- ld_stats_df$pos[match(snp_support_df_full$rsID, ld_stats_df$rsID)]
    snp_support_df_full$tag_name <- ld_stats_df$tag_name[match(snp_support_df_full$rsID, ld_stats_df$rsID)]
    snp_support_df_full$support[is.na(snp_support_df_full$support)] <- "No support"
    snp_support_df_full$merged_hmm_support[is.na(snp_support_df_full$merged_hmm_support)] <- "No support"

    ## now we have to order all the factors correctly and get the axis colors (which get reused)
    ## make sure the tissues are ordered
    ordered_support_classes <- as.vector(sort(unique(snp_support_df_full$tissue_class), decreasing=T))
    snp_support_df_full$tissue_class <- factor(snp_support_df_full$tissue_class,
                                                       ordered=T, levels=ordered_support_classes)
    support_text_cols <- category_colors[ordered_support_classes]

    ## we also want to color SNPs by their tag region, so order the SNPs by tag region, then
    ## position, then rsID (if there are any at the same position..)
    snp_support_ordered <- order(snp_support_df_full$tag_name, snp_support_df_full$pos, snp_support_df_full$rsID, decreasing=F)

    ordered_rsIDs <- unique(snp_support_df_full$rsID[snp_support_ordered])
    snp_support_df_full$rsID <- factor(snp_support_df_full$rsID, ordered=T,
                                               levels=ordered_rsIDs)

    tag_name_cols <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(snp_support_df_full$tag_name)))
    names(tag_name_cols) <- unique(snp_support_df_full$tag_name)

    snp_text_cols <- tag_name_cols[snp_support_df_full$tag_name[match(ordered_rsIDs, snp_support_df_full$rsID)]]

    ## find the unique combinations of full states and also get them in the desired factor order
    unique_supports <- sort(unique(snp_support_df_full$support))
    snp_support_df_full$support <- factor(snp_support_df_full$support, ordered=T,
                                          levels=c(grep("\\+|No support", unique_supports, val=T, inv=T),
                                              grep("\\+", unique_supports, val=T, inv=F),
                                              "No support"))

    ## same thing for the merged HMM states
    unique_merged_hmm_supports <- sort(unique(snp_support_df_full$merged_hmm_support))
    snp_support_df_full$merged_hmm_support <- factor(snp_support_df_full$merged_hmm_support, ordered=T,
                                          levels=c("eQTL", "eRNA Enh", "HMM Enh", "eQTL+eRNA Enh",
                                              "eQTL+HMM Enh", "eRNA Enh+HMM Enh",
                                              "eQTL+eRNA Enh+HMM Enh", "No support"))

    ## ----------------------
    ## plot the full heatmap of all states
    ## this goes in the main output directory and uses the normal prefix
    snp_support_heatmap(snp_support_df_full, paste0(outdir, 'plots/'), prefix, "all", support_text_cols, tag_name_cols, snp_text_cols)

    ## ----------------------
    ## also plot heatmaps for each tag region
    for(tag_region in unique(snp_support_df_full$tag_name)) {
        tag_region_name <- strsplit(tag_region, ":")[[1]][1]
        tag_region_name <- gsub("/", "_", tag_region_name)
        tag_outdir <- paste0(outdir, "/plots/", tag_region_name, "/")
        dir.create(tag_outdir, F, T)

        this_tag_support_df <- snp_support_df_full[snp_support_df_full$tag_name==tag_region,]
        ## just to make sure this is sorted
        this_tag_ordered_rsIDs <- unique(as.vector(this_tag_support_df$rsID[order(this_tag_support_df$pos, this_tag_support_df$rsID, decreasing=F)]))
        this_tag_support_df$rsID <- factor(this_tag_support_df$rsID, ordered=T,
                                      levels=this_tag_ordered_rsIDs)

        snp_support_heatmap(this_tag_support_df, tag_outdir, paste0(prefix, '_', tag_region_name),
                            tag_region_name, support_text_cols, tag_name_cols,
                            ## set the snp colors to just be the tag region color
                            tag_name_cols[tag_region], snp_text_size=15)
    }

    ## ----------------------
    ## now summarize the overlaps for each tag region. the approach here is to make heatmaps
    ## for each possible state and find the regions containing SNPs meeting those criteria
    tag_region_support_df <- ddply(snp_support_df_full, .(tissue_class, tag_name), function(x) {
        eqtl <- "eQTL" %in% x$merged_hmm_support
        erna <- "eRNA Enh" %in% x$merged_hmm_support
        hmm <- "HMM Enh" %in% x$merged_hmm_support
        eqtl_erna <- "eQTL+eRNA Enh" %in% x$merged_hmm_support
        eqtl_hmm <- "eQTL+HMM Enh" %in% x$merged_hmm_support
        erna_hmm <- "eRNA Enh+HMM Enh" %in% x$merged_hmm_support
        eqtl_erna_hmm <- "eQTL+eRNA Enh+HMM Enh" %in% x$merged_hmm_support

        ## ## now enumerate the types of overlaps that are present
        ## support <- ""
        ## if (eqtl) {
        ##     support <- paste(support, "eQTL", sep=", ") }
        ## if(erna) {
        ##     support <- paste(support, "eRNA Enh", sep=", ") }
        ## if(hmm) {
        ##     support <- paste(support, "HMM Enh", sep=", ") }
        ## if(eqtl_erna) {
        ##     support <- paste(support, "eQTL+eRNA Enh", sep=", ") }
        ## if(eqtl_hmm) {
        ##     support <- paste(support, "eQTL+HMM Enh", sep=", ") }
        ## if(erna_hmm) {
        ##     support <- paste(support, "eRNA Enh+HMM Enh", sep=", ") }
        ## if(eqtl_erna_hmm) {
        ##     support <- paste(support, "eQTL+eRNA Enh+HMM Enh", sep=", ") }
        ## if(!eqtl && !erna && !hmm && !eqtl_erna && !eqtl_hmm && !erna_hmm && !eqtl_erna_hmm) {
        ##     support <- "No support"
        ## }
        ## ## remove the separator at the beginning of the support string
        ## support <- gsub("^, ", "", support)

        ## add columns for each state, noting that variants that overlap more sources of data
        ## should also count towards the individual (or combination) of sources
        return(data.frame(## support,
                          eqtl=ifelse(eqtl | eqtl_erna | eqtl_hmm | eqtl_erna_hmm,
                              "eQTL", "No overlap"),
                          erna=ifelse(erna | eqtl_erna | erna_hmm | eqtl_erna_hmm,
                              "eRNA Enh", "No overlap"),
                          hmm=ifelse(hmm | eqtl_hmm | erna_hmm | eqtl_erna_hmm,
                              "HMM Enh", "No overlap"),
                          eqtl_erna=ifelse(eqtl_erna | eqtl_erna_hmm,
                              "eQTL+eRNA Enh", "No overlap"),
                          eqtl_hmm=ifelse(eqtl_hmm | eqtl_erna_hmm,
                              "eQTL+HMM Enh", "No overlap"),
                          erna_hmm=ifelse(erna_hmm | eqtl_erna_hmm,
                              "eRNA Enh+HMM Enh", "No overlap"),
                          eqtl_erna_hmm=ifelse(eqtl_erna_hmm,
                              "eQTL+eRNA Enh+HMM Enh", "No overlap"),
                          stringsAsFactors = F))
    })

    melt_tag_region_support <- cbind(melt(tag_region_support_df, id.vars=c("tissue_class", "tag_name"),
                                          variable.name="annotation"),
                                     ## add a width variable to be able to show all 7 overlaps
                                     width=1 / 7)

    ## make sure this is in order
    melt_tag_region_support$tissue_class <- factor(melt_tag_region_support$tissue_class,
                                                 ordered=T, levels=ordered_support_classes)
    support_text_cols <- category_colors[ordered_support_classes]

    melt_tag_region_support$tag_name <- factor(melt_tag_region_support$tag_name, ordered=T,
                                             levels=sort(unique(melt_tag_region_support$tag_name), decreasing=F))

    ## also add a shifting variable for where to put each subtile. this uses 7 equally spaced
    ## shifts around 0 for the 7 possible annotations, calculated here:
    shift_values <- (1:7) / 7 - (1/14) - 1/2

    shift_vec <- ifelse(melt_tag_region_support$annotation=="eqtl", shift_values[1],
                        ifelse(melt_tag_region_support$annotation=="erna", shift_values[2],
                        ifelse(melt_tag_region_support$annotation=="hmm", shift_values[3],
                        ifelse(melt_tag_region_support$annotation=="eqtl_erna", shift_values[4],
                        ifelse(melt_tag_region_support$annotation=="eqtl_hmm", shift_values[5],
                        ifelse(melt_tag_region_support$annotation=="erna_hmm", shift_values[6],
                               ## this is for eqtl_erna_hmm
                               shift_values[7]))))))

    melt_tag_region_support <- cbind(melt_tag_region_support, shift=shift_vec)

    ## make a factor to order the annotation levels:
    melt_tag_region_support$value <- factor(melt_tag_region_support$value, ordered=T,
                                            levels=c("eQTL", "eRNA Enh", "HMM Enh", "eQTL+eRNA Enh",
                                                "eQTL+HMM Enh", "eRNA Enh+HMM Enh",
                                                "eQTL+eRNA Enh+HMM Enh", "No overlap"))
    melt_tag_region_support$tag_no_rsid <- gsub(":rs.*", "", melt_tag_region_support$tag_name)

    ## finally, we need to map our tag names to numbers so that we can actually use the shift!
    tag_name_numbers <- as.numeric(sort(unique(melt_tag_region_support$tag_name)))
    ## set this for the different settings for tag region naming
    if(TAG_VAR=="tag_no_rsid") {
        names(tag_name_numbers) <- sort(unique(melt_tag_region_support$tag_no_rsid))
    } else if(TAG_VAR=="both") {
        names(tag_name_numbers) <- sort(unique(melt_tag_region_support$tag_name))
    }
    
    ## make the combined tag region heatmap!
    bg_color <- "#e6e6e6"

    ## make two plots, one with color text and one without
    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_fantom5_gtex_tag_region_support_heatmap_color_axis_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 3.0, height_ratio = 2.0)
    print(ggplot(melt_tag_region_support,
                 aes(x=as.numeric(tag_name)+shift, y=tissue_class, width=width, fill=value)) +
          geom_tile(size=1) +
          scale_fill_manual(name="Annotation overlaps",
                            values=c(eqtl_color, erna_color, merge_hmm_color,
                                eqtl_erna_color, eqtl_merged_hmm_color, erna_merged_hmm_color,
                                eqtl_erna_merged_hmm_color, "No overlap"=bg_color)) +
          ## add a tile to outline the categories
          geom_tile(aes(x=as.numeric(tag_name), y=tissue_class, width=1),
                    ## old color: #948300
                    color="grey39", fill=NA, size=1) +
          scale_x_continuous(breaks=tag_name_numbers, labels=names(tag_name_numbers),
                             limits=c(0.5, max(tag_name_numbers)+0.5), expand=c(0, 0)) +
          theme_bw() +
          xlab(TAG_LAB) + ylab("Tissue Category") +
          plot_title("Heatmap of sources of functional support across tag regions and tissue categories", r2_thresh, dist_thresh, out_subtitle, strwrap_width=50) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                 fill = guide_legend(override.aes=list(size=10), nrow=3)) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*1.5),
                axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE*1.5, colour=support_text_cols),
                title = element_text(size=TITLE_SIZE*1.5), plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size=LEGEND_TEXT_SIZE*1.75),
                legend.position="bottom"),
          )
    dev.off()

    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_fantom5_gtex_tag_region_support_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 3.0, height_ratio = 2.0)
    print(ggplot(melt_tag_region_support,
                 aes(x=as.numeric(tag_name)+shift, y=tissue_class, width=width, fill=value)) +
          geom_tile(size=1) +
          scale_fill_manual(name="Annotation overlaps",
                            values=c(eqtl_color, erna_color, merge_hmm_color,
                                eqtl_erna_color, eqtl_merged_hmm_color, erna_merged_hmm_color,
                                eqtl_erna_merged_hmm_color, "No overlap"=bg_color)) +
          ## add a tile to outline the categories
          geom_tile(aes(x=as.numeric(tag_name), y=tissue_class, width=1),
                    color="grey39", fill=NA, size=1) +
          scale_x_continuous(breaks=tag_name_numbers, labels=names(tag_name_numbers),
                             limits=c(0.5, max(tag_name_numbers)+0.5), expand=c(0, 0)) +
          theme_bw() +
          xlab(TAG_LAB) + ylab("Tissue Category") +
          plot_title("Heatmap of sources of functional support across tag regions and tissue categories", r2_thresh, dist_thresh, out_subtitle, strwrap_width=50) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                 fill = guide_legend(override.aes=list(size=10), nrow=3)) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*1.5),
                axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE*1.5),
                title = element_text(size=TITLE_SIZE*1.5), plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size=LEGEND_TEXT_SIZE*1.75),
                legend.position="bottom")
          )                    
    dev.off()
    
    ## make separate heatmaps for each state:
    tag_region_support_df$tag_no_rsid <- gsub(":rs.*", "", tag_region_support_df$tag_name)
    
    ## this uses the color variables from the heatmap function!
    for(this_s in c("eqtl", "erna", "hmm", "eqtl_erna", "eqtl_hmm", "erna_hmm", "eqtl_erna_hmm")) {
        if(this_s=="eqtl") {
            this_col_scale <- scale_fill_manual(name="eQTL Overlap",
                                                values=c(eqtl_color, "No overlap"=bg_color))
        } else if(this_s=="erna") {
            this_col_scale <- scale_fill_manual(name="eRNA Enh Overlap",
                                                values=c(erna_color, "No overlap"=bg_color))
        } else if(this_s=="hmm") {
            this_col_scale <- scale_fill_manual(name="HMM Enh Overlap",
                                                values=c(merge_hmm_color, "No overlap"=bg_color))
        } else if(this_s=="eqtl_erna") {
            this_col_scale <- scale_fill_manual(name="eQTL+eRNA Enh Overlap",
                                                values=c(eqtl_erna_color, "No overlap"=bg_color))
        } else if(this_s=="eqtl_hmm") {
            this_col_scale <- scale_fill_manual(name="eQTL+HMM Enh Overlap",
                                                values=c(eqtl_merged_hmm_color,
                                                    "No overlap"=bg_color))
        } else if(this_s=="erna_hmm") {
            this_col_scale <- scale_fill_manual(name="eRNA Enh+HMM Enh Overlap",
                                                values=c(erna_merged_hmm_color,
                                                    "No overlap"=bg_color))
        } else if(this_s=="eqtl_erna_hmm") {
            this_col_scale <- scale_fill_manual(name="eQTL+eRNA Enh+HMM Enh Overlap",
                                                values=c(eqtl_erna_merged_hmm_color,
                                                    "No overlap"=bg_color))
        }

        make_graphic(paste0(outdir, 'plots/', prefix,
                            '_roadmap_fantom5_gtex_tag_region_support_heatmap_',
                            this_s, "_state_only_",
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 2.5, height_ratio = 2.0)
        print(ggplot(tag_region_support_df, aes_string(x=TAG_VAR, y="tissue_class")) +
              geom_tile(aes(fill=factor(get(this_s)))) +
              this_col_scale +
              theme_bw() +
              xlab(TAG_LAB) + ylab("Tissue Category") +
              plot_title(paste0("Heatmap of sources of functional support across tag regions and tissue categories\nState: ", this_s), r2_thresh, dist_thresh, out_subtitle) +
              guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                     fill = guide_legend(override.aes=list(size=10))) +
              theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*1.5),
                    axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE*1.5),
                    title = element_text(size=TITLE_SIZE*1.5), plot.title = element_text(hjust = 0.5),
                    legend.text = element_text(size=LEGEND_TEXT_SIZE*1.75),
                    legend.position="bottom")
              )
        dev.off()
    }

    ## ----------------------
    ## look at the number of SNPs that are supported by both enhancer-related data sources in
    ## at least one tissue class
    snp_enhancer_support_counts <- ddply(snp_support_df_full, .(tag_name), function(x) {
        ## find all the hits that are supported by both data sources
        support_vec <- x$merged_hmm_support=="eRNA Enh+HMM Enh" | x$merged_hmm_support=="eQTL+eRNA Enh+HMM Enh"
        num_enh_hits <- length(unique(x$rsID[support_vec]))
        return(data.frame(num_enh_hits = num_enh_hits,
                          prop_enh_hits = num_enh_hits / length(unique(x$rsID)),
                          stringsAsFactors = F))
    })
    max_enh_snp <- round_any(max(snp_enhancer_support_counts$num_enh_hits)+1, accuracy=5, f=ceiling)

    snp_enhancer_support_counts$tag_no_rsid <- gsub(":rs.*", "", snp_enhancer_support_counts$tag_name)
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_fantom5_enh_snp_counts_per_region_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(snp_enhancer_support_counts, aes_string(x=TAG_VAR, y="num_enh_hits", fill=TAG_VAR)) +
          scale_fill_hue(h=c(180, 270)) +
          xlab(TAG_LAB) + ylab("Number of unique SNPs") +
          theme_bw() + geom_bar(position="dodge", stat="identity") +
          scale_y_continuous(breaks=seq(0, max_enh_snp, by=1),
                             limits=c(0, max_enh_snp), expand=c(0, 0)) +
          theme(legend.position="none",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title=element_text(hjust = 0.5)) +
          plot_title("Number of SNPs overlapping Roadmap HMM Enhancers and FANTOM5 Enhancers", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## also do the proportions
    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_fantom5_enh_snp_prop_per_region_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(snp_enhancer_support_counts, aes_string(x=TAG_VAR, y="prop_enh_hits", fill=TAG_VAR)) +
          scale_fill_hue(h=c(180, 270)) +
          xlab(TAG_LAB) + ylab("Proportion of unique SNPs") +
          theme_bw() + geom_bar(position="dodge", stat="identity") +
          theme(legend.position="none",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title=element_text(hjust = 0.5)) +
          plot_title("Proportion of SNPs overlapping Roadmap HMM Enhancers and FANTOM5 Enhancers", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

}

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/fantom5_eqtl_chromHMM_overlap/")
## enh_overlap_type <- "locus"
## enh_window <- "1000"
