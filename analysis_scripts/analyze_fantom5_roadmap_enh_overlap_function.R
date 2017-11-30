## analyze_fantom5_roadmap_enh_overlap_function.R
## alex amlie-wolf 01/09/17
## a script to compare FANTOM5 and Roadmap enhancer overlaps, part of the automatic analysis

analyze_fantom5_roadmap_enh_overlap <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh, fantom5_overlap_type, enh_window, fantom5_class_file, roadmap_class_file) {
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
    if (fantom5_overlap_type == "midpoint") {
        fantom5_overlap_file <- paste0(datadir, '/fantom5_overlap/',
                                   prefix, "_", r2_thresh,
                                   "_ld_cutoff_snps_within_", dist_thresh,
                                   "_", enh_window, "bp_around_midpoint_enh_overlaps.txt")
        fantom5_overlap_df <- read.table(fantom5_overlap_file, header=T, sep="\t", quote="", as.is=T)
    } else if (fantom5_overlap_type == "locus") {
        fantom5_overlap_file <- paste0(datadir, '/fantom5_overlap/',
                                   prefix, "_", r2_thresh,
                                   "_ld_cutoff_snps_within_", dist_thresh,
                                   "_", enh_window, "bp_around_orig_locus_enh_overlaps.txt")
        fantom5_overlap_df <- read.table(fantom5_overlap_file, header=T, sep="\t", quote="", as.is=T)
    }

    ## read in the tissue categories
    fantom5_category_df <- read.table(fantom5_class_file, header=T, sep="\t", quote="", as.is=T)
    ## add a column for easier matching
    fantom5_category_df$enh_source <- gsub("_expressed_enhancers.bed", "", fantom5_category_df$FANTOM5.File)

    ## add class columns to the enhancer data:
    fantom5_overlap_df$f5_class <- fantom5_category_df$Class[match(fantom5_overlap_df$enh_source,
                                                               fantom5_category_df$enh_source)]
    ## clean up the tissue names
    fantom5_overlap_df$enh_source <- gsub("CL:[0-9]*_|UBERON:[0-9]*_", "", fantom5_overlap_df$enh_source)
    ## add a column for tag without rsid
    fantom5_overlap_df$tag_no_rsid <- gsub(":rs.*", "", fantom5_overlap_df$tag_name)

    ## get the columns that describe each unique SNP: (this applies to each data source)
    uniq_snp_cols <- colnames(fantom5_overlap_df)[!(colnames(fantom5_overlap_df) %in% c("enh_source", "f5_class", "enh_chr", "enh_start", "enh_end", tagsnp_cols))]
    ## now summarize the data for each unique SNP
    fantom5_uniq_enh_snp_df <- ddply(fantom5_overlap_df, uniq_snp_cols, summarize,
                                     num_f5_tissues = length(unique(enh_source)),
                                     f5_tissues = paste(sort(unique(enh_source)), collapse=","),
                                     num_f5_classes = length(unique(f5_class)),
                                     f5_classes = paste(sort(unique(f5_class)), collapse=","))

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
    ## append the real names to the IDs
    colnames(uniq_roadmap_snp_state_df)[eid_state_col_idx] <- paste0(eid_orig_cols, "_", eid_name_cols)
    
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
    melt_roadmap_df$class <- roadmap_category_df$Class[match(gsub("_.*", "", melt_roadmap_df$tissue), roadmap_category_df$EID)]

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

    ## now find any overlapping SNPs
    roadmap_fantom5_overlap_snp_df <- join(fantom5_uniq_enh_snp_df, roadmap_uniq_enh_snp_df, type="inner")
    if(nrow(roadmap_fantom5_overlap_snp_df) > 0) {
        ## add a column for overlapping classes
        roadmap_fantom5_overlap_snp_df$overlap_classes <- apply(roadmap_fantom5_overlap_snp_df, 1, function(x) {
            hmm_classes <- strsplit(x["hmm_classes"], ",")[[1]]
            f5_classes <- strsplit(x["f5_classes"], ",")[[1]]
            common_classes <- intersect(hmm_classes, f5_classes)
            if(length(common_classes) > 0) {
                return(paste(common_classes, collapse=","))
            } else {
                return("None")
            }
        })
        
        ## write out the full overlap file
        overlap_snp_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                                   "_ld_cutoff_snps_within_", dist_thresh,
                                   "_", fantom5_overlap_type,
                                   "_fantom5_roadmap_overlap.txt")
        write.table(roadmap_fantom5_overlap_snp_df, overlap_snp_outf, quote=F, sep="\t", col.names=T, row.names=F)
    }
        
    ## also make a table of all the SNPs with their support from the two different sources
    roadmap_fantom5_union_snp_df <- join(fantom5_uniq_enh_snp_df, roadmap_uniq_enh_snp_df, type="full")
    if(nrow(roadmap_fantom5_union_snp_df) > 0) {
        roadmap_fantom5_union_snp_df$overlap_classes <- apply(roadmap_fantom5_union_snp_df, 1, function(x) {
            if(is.na(x["hmm_classes"]) | is.na(x["f5_classes"])) {
                return("None")
            }
            hmm_classes <- strsplit(x["hmm_classes"], ",")[[1]]
            f5_classes <- strsplit(x["f5_classes"], ",")[[1]]
            common_classes <- intersect(hmm_classes, f5_classes)
            if(length(common_classes) > 0) {
                return(paste(common_classes, collapse=","))
            } else {
                return("None")
            }
        })

        ## write out the union file
        union_snp_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                                 "_ld_cutoff_snps_within_", dist_thresh,
                                 "_", fantom5_overlap_type,
                                 "_fantom5_roadmap_union.txt")
        write.table(roadmap_fantom5_union_snp_df, union_snp_outf, quote=F, sep="\t", col.names=T, row.names=F)
    }
    
    ## ----------------------
    ## analysis plots
    ## want a heatmap for tissue class x tag region support from these two data sources
    f5_snp_category_combos <- melt(unique(fantom5_overlap_df[,c("rsID", "pos", "tag_name", "tag_no_rsid", "f5_class")]),
                                   id.vars=c("rsID", "tag_name", "tag_no_rsid", "pos"),
                                   variable.name="data_source", value.name="tissue_class")
    roadmap_enh_snp_category_combos <- unique(roadmap_enhancer_snps[,c("rsID", "pos", "tag_name", "tag_no_rsid", "state", "class")])
    ## rename these columns to match the other data frames
    colnames(roadmap_enh_snp_category_combos) <- c("rsID", "pos", "tag_name", "tag_no_rsid", "data_source", "tissue_class")

    ## get the supporting data sets per SNP and tissue category
    ## keep track of both the separate HMM states as well as the merged
    snp_support_df <- ddply(rbind(f5_snp_category_combos, roadmap_enh_snp_category_combos),
                            .(rsID, pos, tag_name, tag_no_rsid, tissue_class), function(x) {
                                ## return the amount and type of support
                                ## first check for individual data sources, then do combos
                                f5_present <- "f5_class" %in% x$data_source
                                hmm_enh_present <- "7_Enh" %in% x$data_source
                                hmm_enhG_present <- "6_EnhG" %in% x$data_source
                                hmm_enh_biv_present <- "12_EnhBiv" %in% x$data_source

                                ## make a string describing the types of individual support:
                                support <- ""
                                ## also make one for the merged HMM states
                                merged_hmm_support <- ""
                                if(f5_present) {
                                    support <- paste(support, "FANTOM5 Enh", sep="+")
                                    merged_hmm_support <- paste(merged_hmm_support, "FANTOM5 Enh", sep="+") }
                                if(hmm_enh_present) {
                                    support <- paste(support, "Roadmap Enh", sep="+") }
                                if(hmm_enhG_present) {
                                    support <- paste(support, "Roadmap Genic Enh", sep="+") }
                                if(hmm_enh_biv_present) {
                                    support <- paste(support, "Roadmap Biv Enh", sep="+") }

                                if(hmm_enh_present || hmm_enhG_present || hmm_enh_biv_present) {
                                    merged_hmm_support <- paste(merged_hmm_support, "Roadmap Enh", sep="+") }

                                ## remove the separator at the beginning of the support strings
                                support <- gsub("^\\+", "", support)
                                merged_hmm_support <- gsub("^\\+", "", merged_hmm_support)

                                return(data.frame(support, merged_hmm_support, stringsAsFactors = F))
                            })
    
    ## we also need to add annotations for all SNP-tissue combos that were not observed
    snp_support_df_full <- merge(expand.grid(rsID = unique(ld_stats_df$rsID),
                                             tissue_class=union(fantom5_category_df$Class, roadmap_category_df$Class),
                                             stringsAsFactors = F),
                                 snp_support_df, all.x=T)

    ## fix all the NAs we generated
    snp_support_df_full$pos <- ld_stats_df$pos[match(snp_support_df_full$rsID, ld_stats_df$rsID)]
    snp_support_df_full$tag_name <- ld_stats_df$tag_name[match(snp_support_df_full$rsID, ld_stats_df$rsID)]
    snp_support_df_full$tag_no_rsid <- ld_stats_df$tag_no_rsid[match(snp_support_df_full$rsID, ld_stats_df$rsID)]    
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
                                                     levels=c("FANTOM5 Enh", "Roadmap Enh",
                                                         "FANTOM5 Enh+Roadmap Enh", "No support"))

    ## use this data to summarize the overlaps in each tag region
    tag_region_support_df <- ddply(snp_support_df_full, .(tissue_class, tag_name, tag_no_rsid), function(x) {
        erna <- "FANTOM5 Enh" %in% x$merged_hmm_support
        hmm <- "Roadmap Enh" %in% x$merged_hmm_support
        erna_hmm <- "FANTOM5 Enh+Roadmap Enh" %in% x$merged_hmm_support

        ## add columns for each state, noting that variants that overlap more sources of data
        ## should also count towards the individual (or combination) of sources
        return(data.frame(erna=ifelse(erna | erna_hmm,
                              "FANTOM5 Enh", "No overlap"),
                          hmm=ifelse(hmm | erna_hmm,
                              "Roadmap Enh", "No overlap"),
                          erna_hmm=ifelse(erna_hmm,
                              "FANTOM5 Enh+Roadmap Enh", "No overlap"),
                          stringsAsFactors = F))
    })

    melt_tag_region_support <- cbind(melt(tag_region_support_df, id.vars=c("tissue_class", "tag_name", "tag_no_rsid"),
                                          variable.name="annotation"),
                                     ## add a width variable to be able to show all overlaps
                                     width=1 / 3)

    ## make sure this is in order
    melt_tag_region_support$tissue_class <- factor(melt_tag_region_support$tissue_class,
                                                 ordered=T, levels=ordered_support_classes)
    support_text_cols <- category_colors[ordered_support_classes]

    melt_tag_region_support$tag_name <- factor(melt_tag_region_support$tag_name, ordered=T,
                                             levels=sort(unique(melt_tag_region_support$tag_name), decreasing=F))

    ## also add a shifting variable for where to put each subtile. this uses 3 equally spaced
    ## shifts around 0 for the 3 possible annotations, calculated here:
    ## http://stackoverflow.com/questions/22107666/generating-split-color-rectangles-from-ggplot2-geom-raster
    shift_values <- (1:3) / 3 - (1/6) - 1/2

    shift_vec <- ifelse(melt_tag_region_support$annotation=="erna", shift_values[1],
                        ifelse(melt_tag_region_support$annotation=="hmm", shift_values[2],
                        shift_values[3]))

    melt_tag_region_support <- cbind(melt_tag_region_support, shift=shift_vec)
    
    ## make a factor to order the annotation levels:
    melt_tag_region_support$value <- factor(melt_tag_region_support$value, ordered=T,
                                            levels=c("FANTOM5 Enh", "Roadmap Enh", "FANTOM5 Enh+Roadmap Enh",
                                                "No overlap"))

    ## finally, we need to map our tag names to numbers so that we can actually use the shift!
    tag_name_numbers <- seq(1, length(unique(melt_tag_region_support$tag_name)))
    ## set this for the different settings for tag region naming
    if(TAG_VAR=="tag_no_rsid") {
        names(tag_name_numbers) <- sort(unique(melt_tag_region_support$tag_no_rsid))
    } else if(TAG_VAR=="both") {
        names(tag_name_numbers) <- sort(unique(melt_tag_region_support$tag_name))
    }
    
    ## make the combined tag region heatmap!
    bg_color <- "#e6e6e6"
    
    ## make two plots, one with color text and one without
    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_fantom5_tag_region_support_heatmap_color_axis_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 3.0, height_ratio = 2.0)
    print(ggplot(melt_tag_region_support,
                aes(x=tag_name_numbers[tag_name]+shift, y=tissue_class, width=width, fill=value)) +
          geom_tile(size=1) +
          scale_fill_manual(name="Annotation overlaps",
                            values=c(erna_color, merge_hmm_color, erna_merged_hmm_color,
                                "No overlap"=bg_color)) +
          ## add a tile to outline the categories
          geom_tile(aes(x=tag_name_numbers[tag_name], y=tissue_class, width=1),
                    ## old color: #948300
                    color="grey39", fill=NA, size=1) +
          scale_x_continuous(breaks=tag_name_numbers, labels=names(tag_name_numbers),
                             limits=c(0.5, max(tag_name_numbers)+ifelse(length(tag_name_numbers)==1, 0.51, 0.5)), expand=c(0, 0)) +
          theme_bw() +
          xlab(TAG_LAB) + ylab("Tissue Category") +
          plot_title("Heatmap of sources of functional support across tag regions and tissue categories", r2_thresh, dist_thresh, out_subtitle, strwrap_width=50) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                 fill = guide_legend(override.aes=list(size=10))) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*1.5),
                axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE*1.5, colour=support_text_cols),
                title = element_text(size=TITLE_SIZE*1.5), plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size=LEGEND_TEXT_SIZE*1.75),
                legend.position="bottom"),
          )
    dev.off()

    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_fantom5_tag_region_support_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 3.0, height_ratio = 2.0)
    print(ggplot(melt_tag_region_support,
                 aes(x=tag_name_numbers[tag_name]+shift, y=tissue_class, width=width, fill=value)) +
          geom_tile(size=1) +
          scale_fill_manual(name="Annotation overlaps",
                            values=c(erna_color, merge_hmm_color, erna_merged_hmm_color,
                                "No overlap"=bg_color)) +
          ## add a tile to outline the categories
          geom_tile(aes(x=tag_name_numbers[tag_name], y=tissue_class, width=1),
                    color="grey39", fill=NA, size=1) +
          scale_x_continuous(breaks=tag_name_numbers, labels=names(tag_name_numbers),
                             limits=c(0.5, max(tag_name_numbers)+ifelse(length(tag_name_numbers)==1, 0.51, 0.5)), expand=c(0, 0)) +
          theme_bw() +
          xlab(TAG_LAB) + ylab("Tissue Category") +
          plot_title("Heatmap of sources of functional support across tag regions and tissue categories", r2_thresh, dist_thresh, out_subtitle, strwrap_width=50) +
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

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/fantom5_roadmap_overlap/")
## fantom5_overlap_type <- "locus"
## enh_window <- "1000"
