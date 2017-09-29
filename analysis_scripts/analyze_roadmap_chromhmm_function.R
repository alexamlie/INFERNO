## analyze_roadmap_chromhmm_function.R
## alex amlie-wolf 03/18/16
## analyzes the epigenetic states
## part of automatic analysis

analyze_roadmap_chromHMM <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh, roadmap_class_file) {
    ## -----------------------------------------------------------------------------
    ## read in data
    ## -----------------------------------------------------------------------------
    dir.create(paste0(outdir, 'plots/'), F, T)
    dir.create(paste0(outdir, 'tables/'), F, T)

    ## read in the roadmap file
    roadmap_state_file <- paste0(datadir, '/roadmap_chromhmm_states/', prefix, "_", r2_thresh,
                                "_ld_cutoff_snps_within_", dist_thresh, "_roadmap_chromHMM_states.txt")
    roadmap_state_df <- read.table(roadmap_state_file, header=T, sep="\t", quote="", as.is=T)
    roadmap_state_df$tag_no_rsid <- gsub(":rs.*", "", roadmap_state_df$tag_name)

    ## read in the roadmap category file
    roadmap_category_df <- read.table(roadmap_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")
    ## rename the EID column to be easier to match
    colnames(roadmap_category_df)[3] <- "EID"

    ## get a table with only unique hits, also rename the state columns
    uniq_snp_cols <- colnames(roadmap_state_df)[!(colnames(roadmap_state_df) %in% tagsnp_cols)]
    uniq_snp_state_df <- ddply(roadmap_state_df, uniq_snp_cols, function(x) {
        apply(x[,!(colnames(x) %in% uniq_snp_cols)], 2, paste, collapse=",")
    })#, .progress="text")

    ## rename the state columns to reflect the actual data source
    ## first get an index to these columns
    eid_state_col_idx <- grepl("_state", colnames(uniq_snp_state_df))
    eid_orig_cols <- gsub("_state", "", colnames(uniq_snp_state_df)[eid_state_col_idx])
    ## now grab the real names
    eid_name_cols <- roadmap_category_df$Standardized.Epigenome.name[match(eid_orig_cols, roadmap_category_df$EID)]
    ## append the real names to the IDs
    colnames(uniq_snp_state_df)[eid_state_col_idx] <- paste0(eid_orig_cols, "_", eid_name_cols)
    
    ## write this table out
    uniq_tab_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                            "_ld_cutoff_snps_within_", dist_thresh,
                            "_uniq_snp_roadmap_chromHMM_states.txt")
    write.table(uniq_snp_state_df, uniq_tab_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## also write it without the tag SNP info
    uniq_tab_small_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                                  "_ld_cutoff_snps_within_", dist_thresh,
                                  "_uniq_snp_roadmap_chromHMM_states_no_tagsnp_info.txt")
    ## store a logical vector because we re-use this
    tissname_uniq_snp_cols <- !(colnames(uniq_snp_state_df) %in% tagsnp_cols)
    write.table(uniq_snp_state_df[,colnames(uniq_snp_state_df)[tissname_uniq_snp_cols]],
                uniq_tab_small_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## now melt the table
    ## first have to get the names of columns we want: no tissue names and no tag SNP info
    melt_cols <- colnames(uniq_snp_state_df)[!eid_state_col_idx & tissname_uniq_snp_cols]

    melt_uniq_state_df <- melt(uniq_snp_state_df, id.vars=melt_cols,
                               measure.vars=colnames(uniq_snp_state_df)[eid_state_col_idx],
                               variable.name="tissue", value.name="state")
    melt_uniq_state_df$tissue <- as.character(melt_uniq_state_df$tissue)
    melt_uniq_state_df$state <- as.character(melt_uniq_state_df$state)

    ## add a tissue class column to this data structure
    ## CHANGED 07/18/17 (need to match on EID)
    melt_uniq_state_df$class <- roadmap_category_df$Class[match(gsub("_.*", "", melt_uniq_state_df$tissue), roadmap_category_df$EID)]
    
    ## now grab the enhancer states:
    roadmap_enhancer_snps <- melt_uniq_state_df[melt_uniq_state_df$state %in% c("6_EnhG", "7_Enh", "12_EnhBiv"),]
    if(nrow(roadmap_enhancer_snps)==0) {
        ## still write a spoofed output file
        roadmap_uniq_enh_snp_df <- data.frame(chr=numeric(0), rsID=numeric(0), pos=numeric(0), ref=numeric(0), alt=numeric(0), MAF=numeric(0), tag_name=numeric(0), tag_no_rsid=numeric(0), num_hmm_tissues=numeric(0), hmm_tissues=numeric(0), num_hmm_enh_states=numeric(0), hmm_enh_states=numeric(0), num_hmm_classes=numeric(0), hmm_classes=numeric(0))
        
        enh_tab_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                               "_ld_cutoff_snps_within_", dist_thresh,
                               "_uniq_chromHMM_enh_snps.txt")
        write.table(roadmap_uniq_enh_snp_df, enh_tab_outf, quote=F, sep="\t", col.names=T, row.names=F)
        cat("No variants overlapping Roadmap ChromHMM enhancer found in this dataset!\n")
        return("No variants overlapping Roadmap ChromHMM enhancer found in this dataset!")
    } 

    ## summarize these per snp:
    roadmap_uniq_enh_snp_df <- ddply(roadmap_enhancer_snps, .(chr, rsID, pos, ref, alt, MAF, tag_name, tag_no_rsid), summarize,
                                     num_hmm_tissues = length(unique(tissue)),
                                     hmm_tissues = paste(sort(unique(tissue)), collapse=","),
                                     num_hmm_enh_states = length(unique(state)),
                                     hmm_enh_states = paste(sort(unique(state)), collapse=","),
                                     num_hmm_classes = length(unique(class)),
                                     hmm_classes = paste(sort(unique(class)), collapse=","))   

    ## write this table out
    enh_tab_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                            "_ld_cutoff_snps_within_", dist_thresh,
                            "_uniq_chromHMM_enh_snps.txt")
    write.table(roadmap_uniq_enh_snp_df, enh_tab_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## some other manipulations that get used in several analyses
    melt_uniq_state_df$state <- factor(melt_uniq_state_df$state, ordered=T,
                                       levels=mixedsort(unique(melt_uniq_state_df$state)))

    ## order the tissues by their class and get the text colors to use
    ordered_tissues <- unique(as.vector(melt_uniq_state_df$tissue[order(melt_uniq_state_df$class, melt_uniq_state_df$tissue)]))
    melt_uniq_state_df$tissue <- factor(melt_uniq_state_df$tissue, ordered=T,
                                        levels=ordered_tissues)
    full_tiss_colors <- category_colors[roadmap_category_df$Class[match(gsub("E[0-9]...", "", ordered_tissues), roadmap_category_df$Standardized.Epigenome.name)]]
    
    ## we also want to color SNPs by their tag region, so start by sorting by tag region,
    ## position, and rsID
    snp_order <- order(melt_uniq_state_df$tag_name, melt_uniq_state_df$pos, melt_uniq_state_df$rsID, decreasing=F)
    ordered_rsIDs <- unique(melt_uniq_state_df$rsID[snp_order])
    melt_uniq_state_df$rsID <- factor(melt_uniq_state_df$rsID, ordered=T,
                                      levels=ordered_rsIDs)

    ## -----------------------------------------------------------------------------
    ## generate figures
    ## -----------------------------------------------------------------------------
    ## -------------------------
    ## look at heatmap of states in individual tissues
    tag_name_cols <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(melt_uniq_state_df$tag_name)))
    names(tag_name_cols) <- unique(melt_uniq_state_df$tag_name)

    snp_text_cols <- tag_name_cols[melt_uniq_state_df$tag_name[match(ordered_rsIDs, melt_uniq_state_df$rsID)]]

    ## define custom color scheme for 15 states
    ## c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")

    ## active TSS, flanking active TSS, transcribed at gene 5' and 3', strong transcription,
    ## weak transcription, genic enhancer, enhancer, znf genes + repeats, heterchromatin,
    ## bivalent/poised TSS, flanking bivalent TSS/Enh, bivalent enhancer, repressed polycomb,
    ## weak repressed polycomb, quiescent

    ## this was my manual color definition
    ## state_heatmap_cols <- c(## green-yellow for transcription sites
    ##                         "#627502", "#9ca365", "#90a302", "#90db8a", "#3c4704",
    ##                         ## red for enhancers
    ##                         "#ed0303", "#d94646",
    ##                         ## greens for repeats, heterochromatin
    ##                         "#A4DEA8", "#309c42",
    ##                         ## brown for bivalent TSS and flanking
    ##                         "#E7B881", "#a17740",
    ##                         ## dark red for bivalent enhancer
    ##                         "#b80202",
    ##                         ## blues for repressed states
    ##                         "#028f8a", "#387a78", "#A8D8D2"
    ## this one is consistent with the roadmap colors (http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html)
    state_heatmap_cols <- c(## transcription states
                            rgb(255,0,0,max=255), rgb(255,69,0,max=255), rgb(50,205,50,max=255),
                            rgb(0,128,0,max=255), rgb(0,100,0,max=255),
                            ## enhancers
                            rgb(194,225,5,max=255), rgb(255,255,0,max=255),
                            ## repeats, heterochromatin
                            rgb(102,205,170,max=255), rgb(138,145,208,max=255),
                            ## bivalent TSS and flanking
                            rgb(205,92,92,max=255), rgb(233,150,122,max=255),
                            ## bivalent enhancer
                            rgb(189,183,107,max=255),
                            ## repressed states
                            rgb(128,128,128,max=255), rgb(192,192,192,max=255), rgb(255,255,255,max=255))

    names(state_heatmap_cols) <- mixedsort(unique(melt_uniq_state_df$state))
    melt_uniq_state_df$tissue <- as.character(melt_uniq_state_df$tissue)
    
    ## TODO: more descriptive state names
    make_graphic(paste0(outdir, 'plots/', prefix, '_tissue_state_by_variant_heatmap_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"), width_ratio = 3.5, height_ratio = 2.0, type="pdf")
    ## switch x axis to be on top
    print(ggplot(melt_uniq_state_df, aes(x=rsID, y=tissue)) +
          geom_tile(aes(fill=state)) +
          scale_fill_manual(name="Epigenetic state", values = state_heatmap_cols) +
          ## add invisible "points" just to get the legend for the tag SNP region
          geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
          scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
          theme_bw() +
          xlab("Linked SNP") + ylab("Roadmap tissue source") +
          scale_x_discrete(position="top") + 
          plot_title("Heatmap of ChromHMM-defined epigenetic states for variants across tissues", r2_thresh, dist_thresh, out_subtitle, strwrap_width=100) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                 fill = guide_legend(override.aes=list(size=10))) +
          theme(axis.text.x = element_text(angle=-90+1e-09, hjust=1, vjust=1, size=4, colour=snp_text_cols),
                axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE*0.75, colour=full_tiss_colors),
                title = element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size=LEGEND_TEXT_SIZE)))
    dev.off()

    ## make plots for individual states
    for (this_s in unique(melt_uniq_state_df$state)) {        
        this_new_df <- melt_uniq_state_df
        ## unfactor
        this_new_df$state <- as.vector(this_new_df$state)
        this_new_df$state[this_new_df$state != this_s] <- "other"

        this_heatmap_cols <- c("other"="#FFFFFF", state_heatmap_cols[this_s])

        if (this_s == "8_ZNF/Rpts") {
            this_s <- "8_ZNF_Rpts"
        }

        make_graphic(paste0(outdir, 'plots/', prefix, '_', this_s, '_state_by_variant_heatmap_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"), width_ratio = 3.5, height_ratio = 2.0)
        print(ggplot(this_new_df, aes(x=rsID, y=tissue)) +
              geom_tile(aes(fill=state)) +
              scale_fill_manual(name="Epigenetic state", values = state_heatmap_cols) +
              ## add invisible "points" just to get the legend for the tag SNP region
              geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
              scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
              theme_bw() +
              xlab("Linked SNP") + ylab("Roadmap tissue source") +
              scale_x_discrete(position="top") + 
              plot_title(paste0("Heatmap of ChromHMM state ", this_s, " for variants across tissues"), r2_thresh, dist_thresh, out_subtitle, strwrap_width=100) +
              guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                     fill = guide_legend(override.aes=list(size=10))) +
              theme(axis.text.x = element_text(angle=-90+1e-09, hjust=1, vjust=1, size=4, colour=snp_text_cols),
                    axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE*0.75, colour=full_tiss_colors),
                    title = element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                    legend.text = element_text(size=LEGEND_TEXT_SIZE)))
        dev.off()

    }

    ## next: majority state by tissue class

    ## -------------------------
    ## enhancer analysis
    ## we want to know how many SNPs are characterized as an enhancer in any tissue
    enh_nums_per_tagregion <- melt(ddply(melt_uniq_state_df, .(tag_name), summarize,
                                         num_any_enh = length(unique(rsID[state=="6_EnhG" | state=="7_Enh" | state=="12_EnhBiv"])),
                                         num_genic_enh = length(unique(rsID[state=="6_EnhG"])),
                                         num_enh = length(unique(rsID[state=="7_Enh"])),
                                         num_biv_enh = length(unique(rsID[state=="12_EnhBiv"]))),
                                   id.vars="tag_name")
    
    enh_nums_per_tagregion$tag_no_rsid <- gsub(":rs.*", "", enh_nums_per_tagregion$tag_name)
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_enh_counts_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"), width_ratio=2.0, height_ratio=1.5)
    print(ggplot(enh_nums_per_tagregion, aes_string(x="variable", y="value", fill=TAG_VAR)) +
          scale_fill_hue(h=c(180, 270)) +
          xlab("Enhancer Type") + ylab("Number of SNPs") +
          theme_bw() +
          facet_wrap(reformulate(TAG_VAR), scales="free") +
          geom_bar(stat="identity", position="stack") +
#          scale_y_continuous(breaks=seq(0, max(enh_nums_per_tagregion$value), by=1)) +
          theme(legend.position="none",
                axis.text.x=element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE*0.75),
                strip.text = element_text(size=LEGEND_TEXT_SIZE*0.75), 
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Number of SNPs annotated in Roadmap ChromHMM enhancer states by tag region", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## -------------------------
    ## also want to look at the tissue distributions of enhancer states:
    enh_nums_per_tissue <- ddply(melt_uniq_state_df, .(tissue, class), summarize,
                                         num_any_enh = length(unique(rsID[state=="6_EnhG" | state=="7_Enh" | state=="12_EnhBiv"])),
                                         num_genic_enh = length(unique(rsID[state=="6_EnhG"])),
                                         num_enh = length(unique(rsID[state=="7_Enh"])),
                                         num_biv_enh = length(unique(rsID[state=="12_EnhBiv"])))

    ## make sure the data sources are sorted by tissue category
    enh_num_order <- order(enh_nums_per_tissue$class, enh_nums_per_tissue$tissue, decreasing=T)
    enh_nums_per_tissue$tissue <- factor(enh_nums_per_tissue$tissue, ordered=T,
                                         levels=enh_nums_per_tissue$tissue[enh_num_order])

    ## first make a plot of the full enhancer counts
    ## get the text color for this plot
    class_text_col <- category_colors[enh_nums_per_tissue$class[enh_num_order]]

    tiss_limit_max <- round_any(max(enh_nums_per_tissue$num_any_enh)+1, 5, f=ceiling)

    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_all_enh_counts_bytiss_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"),
                 height_ratio=1.5, width_ratio=1.3)
    print(ggplot(enh_nums_per_tissue,
                 aes(x=tissue, y=num_any_enh, fill=class)) +
          cat_col_scale +
          xlab("Cell or tissue type") + ylab("Number of variants in enhancer ChromHMM state") +
          coord_flip() +
          scale_y_continuous(breaks=seq(0, tiss_limit_max, by=5), limits=c(0, tiss_limit_max), expand=c(0, 0)) +
          theme_bw() + geom_bar(stat="identity", position="stack") +
          guides(fill = guide_legend(title="Tissue category", ncol=1)) +
          theme(legend.position="right",
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=11),
                axis.text.y = element_text(colour = class_text_col, size=13),
                plot.title=element_text(size=18), axis.title=element_text(size=16),
                legend.text=element_text(size=13)) +
          plot_title("Numbers of variants in enhancer Roadmap ChromHMM states across tissues", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## now look at the different types of enhancers
    ## to do this, we will melt the individual enhancer columns:
    melted_enh_nums_per_tissue <- melt(enh_nums_per_tissue, id.vars=c("tissue", "class", "num_any_enh"),
                                       measure.vars=c("num_enh", "num_biv_enh", "num_genic_enh"),
                                       variable.name="enh_type")
    melted_enh_nums_per_tissue$enh_type <- factor(melted_enh_nums_per_tissue$enh_type,
                                                  levels=c("num_enh", "num_genic_enh", "num_biv_enh"))
    
    enh_type_tiss_max <- round_any(max(melted_enh_nums_per_tissue$value)+1, 10, f=ceiling)

    ## make a labeller function for the enhancer types
    enh_type_labeller <- as_labeller(c("num_enh"="Enhancers", "num_genic_enh"="Genic enhancers", "num_biv_enh"="Bivalent enhancers"))
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_enh_type_counts_bytiss_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"),
                 height_ratio=1.75, width_ratio=2.0)
    print(ggplot(melted_enh_nums_per_tissue,
                 aes(x=tissue, y=value, fill=class)) +
          cat_col_scale +
          xlab("Cell or tissue type") + ylab("Number of variants in enhancer ChromHMM state") +
          coord_flip() +
          scale_y_continuous(breaks=seq(0, enh_type_tiss_max, by=10), limits=c(0, enh_type_tiss_max), expand=c(0, 0)) +
          theme_bw() + geom_bar(stat="identity", position="stack") +
          guides(fill = guide_legend(title="Tissue category")) +
          facet_grid(. ~ enh_type, labeller=enh_type_labeller) +
          theme(legend.position="bottom",
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE*0.75),
                axis.text.y = element_text(colour = class_text_col, size=AXIS_TEXT_Y_SIZE*0.5),
                strip.text.x = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                legend.text=element_text(size=LEGEND_TEXT_SIZE)) +
          plot_title("Numbers of variants annotated in enhancer Roadmap ChromHMM states across tissues", r2_thresh, dist_thresh, out_subtitle, strwrap_width = 50))
    dev.off()

    ## -------------------------
    ## next, look at the distributions of enhancers within each tissue class
    ## we need to factor this correctly:
    enh_nums_per_tissue$class <- factor(enh_nums_per_tissue$class, ordered=T,
                                         levels=sort(unique(enh_nums_per_tissue$class), decreasing=T))

    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_all_enh_num_distribution_byclass_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"),
                 height_ratio=1.5)
    print(ggplot(enh_nums_per_tissue,
                 aes(x=class, y=num_any_enh, fill=class)) +
          cat_col_scale +
          xlab("Tissue category") + ylab("Number of variants") +
          coord_flip() +
          scale_y_continuous(breaks=seq(0, tiss_limit_max, by=5), limits=c(0, tiss_limit_max), expand=c(0.01, 0)) +
          theme_bw() + geom_boxplot() +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE*0.75),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("# of variants in enhancer Roadmap ChromHMM states", r2_thresh, dist_thresh, out_subtitle, strwrap_width = 30))
    dev.off()

    ## also do the distribution analysis for the individual enhancer types:
    melted_enh_nums_per_tissue$class <- factor(melted_enh_nums_per_tissue$class, ordered=T,
                                         levels=sort(unique(melted_enh_nums_per_tissue$class), decreasing=T))

    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_enh_type_num_distribution_byclass_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"),
                 height_ratio=1.5, width_ratio = 1.5)
    print(ggplot(melted_enh_nums_per_tissue,
                 aes(x=class, y=value, fill=class)) +
          cat_col_scale +
          xlab("Tissue category") + ylab("Number of variants in enhancer ChromHMM state") +
          coord_flip() +
          scale_y_continuous(breaks=seq(0, enh_type_tiss_max, by=10), limits=c(0, enh_type_tiss_max), expand=c(0.01, 0)) +
          theme_bw() + geom_boxplot() +
          facet_grid(. ~ enh_type, labeller=enh_type_labeller) +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE*0.75),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                strip.text.x = element_text(size=LEGEND_TEXT_SIZE), 
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of numbers of variants in enhancer Roadmap ChromHMM states across tissue classes", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## -------------------------
    ## also do this by tissue class
    enh_nums_per_class <- ddply(melt_uniq_state_df, .(class), summarize,
                                         num_any_enh = length(unique(rsID[state=="6_EnhG" | state=="7_Enh" | state=="12_EnhBiv"])),
                                         num_genic_enh = length(unique(rsID[state=="6_EnhG"])),
                                         num_enh = length(unique(rsID[state=="7_Enh"])),
                                         num_biv_enh = length(unique(rsID[state=="12_EnhBiv"])))

    ## make sure the data sources are sorted by tissue category
    enh_nums_per_class$class <- factor(enh_nums_per_class$class, ordered=T,
                                         levels=sort(enh_nums_per_class$class, decreasing=T))
    ## now make a plot of the total enhancers by class
    class_limit_max <- round_any(max(enh_nums_per_class$num_any_enh)+1, 10, f=ceiling)

    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_all_enh_counts_bytiss_class_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"),
                 height_ratio=1.5)
    print(ggplot(enh_nums_per_class,
                 aes(x=class, y=num_any_enh, fill=class)) +
          cat_col_scale +
          xlab("Tissue class") + ylab("Number of variants in enhancer ChromHMM state") +
          coord_flip() +
          scale_y_continuous(breaks=seq(0, class_limit_max, by=10), limits=c(0, class_limit_max), expand=c(0, 0)) +
          theme_bw() + geom_bar(stat="identity", position="stack") +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=11),
                axis.text.y = element_text(size=15), plot.title=element_text(size=18),
                axis.title=element_text(size=16)) +
          plot_title("Numbers of variants in Roadmap ChromHMM enhancer states", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## now look at the different types of enhancers
    ## to do this, we will melt the individual enhancer columns:
    melted_enh_nums_per_class <- melt(enh_nums_per_class, id.vars=c("class", "num_any_enh"),
                                       measure.vars=c("num_enh", "num_biv_enh", "num_genic_enh"),
                                       variable.name="enh_type")
    melted_enh_nums_per_class$enh_type <- factor(melted_enh_nums_per_class$enh_type,
                                                  levels=c("num_enh", "num_genic_enh", "num_biv_enh"))

    enh_type_class_max <- round_any(max(melted_enh_nums_per_class$value)+1, 20, f=ceiling)

    make_graphic(paste0(outdir, 'plots/', prefix, '_roadmap_enh_type_counts_bytiss_class_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"),
                 height_ratio=1.5, width_ratio = 1.5)
    print(ggplot(melted_enh_nums_per_class,
                 aes(x=class, y=value, fill=class)) +
          cat_col_scale +
          xlab("Tissue category") + ylab("Number of variants in enhancer ChromHMM state") +
          coord_flip() +
          scale_y_continuous(breaks=seq(0, enh_type_class_max, by=20), limits=c(0, enh_type_class_max), expand=c(0, 0)) +
          theme_bw() + geom_bar(stat="identity", position="stack") +
          facet_grid(. ~ enh_type, labeller=enh_type_labeller) +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE*0.75),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE*0.5),
                strip.text.x = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) + 
          plot_title("Numbers of variants in enhancer Roadmap ChromHMM states across tissue categories", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## -------------------------
    ## also want to check the tissue patterns across tag regions
    ## -------------------------
    ## count the number and proportion of SNPs in each tissue
    enh_nums_per_region_x_tissue <- ddply(melt_uniq_state_df, .(tag_name, tag_no_rsid, tissue, class),
                                          function(x) {
                          enh_vec <- x$state=="6_EnhG" | x$state=="7_Enh" | x$state=="12_EnhBiv"
                          num_tissue_enh_snps <- length(unique(x$rsID[enh_vec]))
                          enh_count_vec <- enh_nums_per_tagregion$tag_no_rsid == unique(x$tag_no_rsid) & enh_nums_per_tagregion$variable=="num_any_enh"
                          num_tag_enh_snps <- enh_nums_per_tagregion$value[enh_count_vec]
                          return(data.frame(num_tissue_enh_snps, num_tag_enh_snps, prop_tissue_enh_snps = num_tissue_enh_snps / num_tag_enh_snps, stringsAsFactors = F))
                      })

    ## order the data sources
    region_x_tissue_order <- unique(enh_nums_per_region_x_tissue$tissue[order(enh_nums_per_region_x_tissue$class, enh_nums_per_region_x_tissue$tissue)])
    enh_nums_per_region_x_tissue$tissue <- factor(enh_nums_per_region_x_tissue$tissue, ordered=T,
                                                  levels=region_x_tissue_order)

    ## get the colors for the axis text labels
    ordered_heatmap_classes <- enh_nums_per_region_x_tissue$class[match(region_x_tissue_order, enh_nums_per_region_x_tissue$tissue)]
    heatmap_text_cols <- category_colors[ordered_heatmap_classes]

    ## raw count heatmap
    ## white to blue color scheme
    heatmap_col_lims <- c("#FFFFFF", "#08306B")

    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_all_enh_count_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0, height_ratio = 1.5)
    print(ggplot(enh_nums_per_region_x_tissue, aes_string(x="tissue", y=TAG_VAR)) +
          geom_tile(aes(fill=num_tissue_enh_snps), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
          theme_bw() +
          xlab("Roadmap Tissue / Cell Type") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of linked SNPs overlapping any Roadmap ChromHMM enhancer states by tag region and tissue", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=75, hjust=1, size=AXIS_TEXT_X_SIZE*0.6, colour=heatmap_text_cols),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_all_enh_prop_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0, height_ratio = 1.5)
    print(ggplot(enh_nums_per_region_x_tissue, aes_string(x="tissue", y=TAG_VAR)) +
          geom_tile(aes(fill=prop_tissue_enh_snps), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
          theme_bw() +
          xlab("Roadmap Tissue / Cell Type") + ylab(TAG_LAB) +
          plot_title("Heatmap of proportions of linked SNPs overlapping any Roadmap ChromHMM enhancer states in each tissue, by tag region", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=75, hjust=1, size=AXIS_TEXT_X_SIZE*0.6, colour=heatmap_text_cols),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    ## -------------------------
    ## count the number and proportion of SNPs in each tissue class
    enh_nums_per_region_x_class <- ddply(melt_uniq_state_df, .(tag_no_rsid, tag_name, class),
                                          function(x) {
                          enh_vec <- x$state=="6_EnhG" | x$state=="7_Enh" | x$state=="12_EnhBiv"
                          num_class_enh_snps <- length(unique(x$rsID[enh_vec]))
                          enh_count_vec <- enh_nums_per_tagregion$tag_no_rsid == unique(x$tag_no_rsid) & enh_nums_per_tagregion$variable=="num_any_enh"
                          num_tag_enh_snps <- enh_nums_per_tagregion$value[enh_count_vec]
                          return(data.frame(num_class_enh_snps, num_tag_enh_snps,
                                            prop_class_enh_snps = num_class_enh_snps / num_tag_enh_snps, stringsAsFactors = F))
                      })

    ## raw count heatmap for class
    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_all_enh_count_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)
    print(ggplot(enh_nums_per_region_x_class, aes_string(x="class", y=TAG_VAR)) +
          geom_tile(aes(fill=num_class_enh_snps), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(name="# Overlapping SNPs", low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
          theme_bw() +
          xlab("Roadmap Data Source Category") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of linked SNPs overlapping any Roadmap ChromHMM enhancer states by tag region and tissue class", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75,
                    colour=category_colors[unique(enh_nums_per_region_x_class$class)]),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    ## make this with black tissue category labels
    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_all_enh_count_heatmap_black_labels_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)
    print(ggplot(enh_nums_per_region_x_class, aes_string(x="class", y=TAG_VAR)) +
          geom_tile(aes(fill=num_class_enh_snps), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(name="# Overlapping SNPs", low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
          theme_bw() +
          xlab("Roadmap Data Source Category") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of linked SNPs overlapping any Roadmap ChromHMM enhancer states by tag region and tissue class", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    ## also make one with greyscale color theme; we set the scale to start at '1' so we can
    ## see even the weak signals
    greyscale_col_lims <- c("#FFFFFF", "gray90", "black")

    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_all_enh_count_heatmap_greyscale_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)                   
    print(ggplot(enh_nums_per_region_x_class, aes_string(x="class", y=TAG_VAR)) +
          geom_tile(aes(fill=num_class_enh_snps), colour=greyscale_col_lims[1]) +
          scale_fill_gradientn(colours=greyscale_col_lims,
                               values=rescale(c(0, 1, max(enh_nums_per_region_x_class$num_class_enh_snps))), 
                               name="# Overlapping SNPs", na.value=greyscale_col_lims[1],
                               guide="colorbar") +
          theme_bw() +
          xlab("Roadmap Data Source Category") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of linked SNPs overlapping any Roadmap ChromHMM enhancer states by tag region and tissue class", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()        

    ## proportion heatmap for class
    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_all_enh_prop_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)
    print(ggplot(enh_nums_per_region_x_class, aes_string(x="class", y=TAG_VAR)) +
          geom_tile(aes(fill=prop_class_enh_snps), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar") +
          theme_bw() +
          xlab("Roadmap Tissue Class") + ylab(TAG_LAB) +
          plot_title("Heatmap of proportions of linked SNPs overlapping any Roadmap ChromHMM enhancer states by tag region and tissue class", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75,
                    colour=category_colors[unique(enh_nums_per_region_x_class$class)]),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()
}

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/roadmap_chromHMM/")
