## analyze_eqtl_overlap_function.R
## alex amlie-wolf 02/09/16
## a script to look at GTEx eQTL overlap with our SNPs
## part of the automatic analysis

analyze_eqtl_overlap <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh, gtex_class_file) {
    ## -----------------------------------------------------------------------------
    ## read in data
    ## -----------------------------------------------------------------------------
    dir.create(paste0(outdir, 'plots/'), F, T)
    dir.create(paste0(outdir, 'tables/'), F, T)

    ## read in the full eqtl overlap file
    ## warning: this may use a lot of memory!
    eqtl_overlap_file <- paste0(datadir, '/gtex_eqtl_overlap/', prefix, "_", r2_thresh,
                                "_ld_cutoff_snps_within_", dist_thresh, "_eqtl_overlaps.txt")
    eqtl_overlap_df <- read.table(eqtl_overlap_file, header=T, sep="\t", quote="", as.is=T)

    ## read in the tissue categories
    gtex_category_df <- read.table(gtex_class_file, header=T, sep="\t", quote="", as.is=T)
    ## make it easier to match
    gtex_category_df$tissue <- gsub("_Analysis.snpgenes", "", gtex_category_df$GTEx.Data)
    
    ## add a category column
    eqtl_overlap_df$eqtl_class <- gtex_category_df$Class[match(eqtl_overlap_df$tissue, gtex_category_df$tissue)]
    eqtl_overlap_df$tag_no_rsid <- gsub(":rs.*", "", eqtl_overlap_df$tag_name)
    
    ## use the 'global' tagsnp column variable
    uniq_snp_cols <- colnames(eqtl_overlap_df)[!(colnames(eqtl_overlap_df) %in% tagsnp_cols)]

    uniq_snp_eqtl_df <- ddply(eqtl_overlap_df, uniq_snp_cols, function(x) {
        apply(x[,!(colnames(x) %in% uniq_snp_cols)], 2, paste, collapse=",")
    }, .progress='text')
    
    ## write out the unique table
    uniq_tab_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                             "_ld_cutoff_snps_within_", dist_thresh,
                            "_uniq_snp_eqtl_overlaps.txt")
    ## skip this since its so large
    #write.table(uniq_snp_eqtl_df, uniq_tab_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## also write it without the tag SNP info
    uniq_tab_small_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                             "_ld_cutoff_snps_within_", dist_thresh,
                            "_uniq_snp_eqtl_overlaps_no_tagsnp_info.txt")
    write.table(uniq_snp_eqtl_df[,uniq_snp_cols], uniq_tab_small_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## write the overlapping eQTLs as a bed file
    eqtl_bed_data <- unique(eqtl_overlap_df[,c("chr", "rsID", "pos", "tissue", "gene_name")])
    
    eqtl_bed_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                           "_ld_cutoff_snps_within_", dist_thresh,
                           "_eqtl_overlaps.bed")   
    ## first write the track description to this file:
    cat(paste0("track type=bed name=\"GTEx eQTLs, ", out_subtitle, "\" description=\"GTEx eQTLs overlapping INFERNO SNPs\" visibility=3\n"), file=eqtl_bed_outf)
    
    write.table(cbind(eqtl_bed_data$chr, eqtl_bed_data$pos-1, eqtl_bed_data$pos,
                      paste(eqtl_bed_data$rsID, eqtl_bed_data$tissue, eqtl_bed_data$gene_name,
                            sep="-")),
                      eqtl_bed_outf, append=T, quote=F, sep="\t", row.names=F, col.names=F)
    
    rm(eqtl_overlap_df)
    
    ## we also want to show the regions that have no eQTL overlap, so we need to find those
    input_snp_list <- read.table(paste0(datadir, "/ld_expansion/", prefix, "_", r2_thresh,
                                        "_ld_cutoff_snps_within_", dist_thresh, ".txt"),
                                 header=T, sep="\t", quote="", as.is=T)
    all_regions <- unique(input_snp_list$tag_name)

    ## -----------------------------------------------------------------------------
    ## generate figures
    ## -----------------------------------------------------------------------------
    ## -------------------------
    ## general summary plots
    ## -------------------------
    ## plot the number of LD SNPs that are also eQTLs
    num_ld_eqtl_snps <- ddply(uniq_snp_eqtl_df, .(tag_name), summarize,
                              num_eqtl_ld_snps = length(unique(rsID)))    
    
    num_ld_eqtl_snps <- rbind(num_ld_eqtl_snps,
                         data.frame(tag_name = all_regions[which(!(all_regions %in% uniq_snp_eqtl_df$tag_name))],
                                    num_eqtl_ld_snps = 0, stringsAsFactors = F))

    ## re-order 
    num_ld_eqtl_snps <- num_ld_eqtl_snps[order(num_ld_eqtl_snps$tag_name),]
    num_ld_eqtl_snps$tag_no_rsid <- gsub(":rs.*", "", num_ld_eqtl_snps$tag_name)
    
    ## -------------------------
    ## look at number of linked SNPs per tag SNP that are also eQTLs
    make_graphic(paste0(outdir, 'plots/', prefix, '_num_eqtl_snps_per_tagregion_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))  
    print(ggplot(num_ld_eqtl_snps,
                 aes_string(x=TAG_VAR, y="num_eqtl_ld_snps", fill=TAG_VAR)) +
          xlab(TAG_LAB) + ylab("Number of eQTLs") +
          scale_fill_hue(h=c(180, 270)) +
#          scale_y_discrete(breaks=seq(0, max(num_ld_eqtl_snps$num_eqtl_ld_snps)*1.1, by=10)) +
          theme_bw() + geom_bar(stat="identity", position="stack") +
          theme(legend.position="none",
                axis.text.x=element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Numbers of linked SNPs that are also eQTLs", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## -------------------------
    ## look at the beta value distributions for each region
    make_graphic(paste0(outdir, 'plots/', prefix, '_eqtl_beta_dists_per_tagregion_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))  
    print(ggplot(uniq_snp_eqtl_df, aes_string(x=TAG_VAR, y="beta", fill = TAG_VAR), alpha=0.8) + 
                                       ## fill=tag_name, colour = abs(beta))) +
          xlab(TAG_LAB) + ylab("eQTL Beta values") +
          geom_point(position=position_jitter(width=0.6), pch=23, colour=alpha("black", 0.3)) + 
          scale_fill_hue(h=c(180, 270)) +
          ## scale_colour_gradient(low="#132B43", high="#56B1F7") + 
          theme_bw() + ## geom_boxplot() +
          theme(legend.position="none",
                axis.text.x=element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distribution of eQTL beta values across tag regions", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## -------------------------
    ## compare my MAFs with the GTEx-defined ones, only if we have them
    if(sum(!is.na(uniq_snp_eqtl_df$MAF)) > 0) {
        make_graphic(paste0(outdir, 'plots/', prefix, '_maf_comparison_',
                            r2_thresh, "_ld_", dist_thresh, "_dist"), width_ratio = 2.0)  
        print(ggplot(uniq_snp_eqtl_df, aes_string(x="MAF", y="maf", colour=TAG_VAR)) + 
              xlab("Pipeline-defined MAF") + ylab("GTEx-defined MAF") +
              scale_colour_hue(h=c(180, 270)) +
              facet_grid(reformulate(TAG_VAR)) + 
              theme_bw() + 
              geom_point() + 
              theme(legend.position="none",
                    axis.text.x=element_text(angle=90, size=AXIS_TEXT_X_SIZE*0.75, vjust=0.5, hjust=0.5),
                    axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                    strip.text = element_text(size=LEGEND_TEXT_SIZE*0.75), 
                    title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) + 
              plot_title("Comparisons of pipeline- and GTEx-defined MAFs", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70))
        dev.off()
    }
    
    ## -------------------------
    ## target gene analyses
    ## -------------------------
    ## -------------------------
    ## first just find the number of target genes per region
    num_target_genes_per_region <- ddply(uniq_snp_eqtl_df, .(tag_name, tag_no_rsid), summarize,
                                         num_target_genes = length(unique(gene_name)))

    ## need to figure out what scale to use
    y_scale_size <- ifelse(max(num_target_genes_per_region$num_target_genes) < 10, 2, 10)
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_num_target_genes_per_tagregion_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))  
    print(ggplot(num_target_genes_per_region,
                 aes_string(x=TAG_VAR, y="num_target_genes", fill=TAG_VAR)) +
          xlab(TAG_LAB) + ylab("Number of unique target genes") +
          scale_fill_hue(h=c(180, 270)) +
#          scale_y_discrete(breaks=seq(0, max(num_target_genes_per_region$num_target_genes)*1.1, by=y_scale_size)) +
          theme_bw() + geom_bar(stat="identity", position="stack") +
          theme(legend.position="none",
                axis.text.x=element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Numbers of unique eQTL target genes per tag region", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## now find the number of target genes for each unique SNP in each tag region
    num_target_genes_per_snp <- ddply(uniq_snp_eqtl_df, .(tag_name, tag_no_rsid, rsID), summarize,
                                      num_target_genes = length(unique(gene_name)))

    make_graphic(paste0(outdir, 'plots/', prefix, '_num_target_genes_per_snp_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))  
    print(ggplot(num_target_genes_per_snp,
                 aes_string(x=TAG_VAR, y="num_target_genes", fill=TAG_VAR)) +
          xlab(TAG_LAB) + ylab("Number of unique target genes") +
          scale_fill_hue(h=c(180, 270)) +
          scale_y_continuous(breaks=seq(0, max(num_target_genes_per_snp$num_target_genes), by=1)) +
          theme_bw() + geom_boxplot(alpha=0.7) +
          theme(legend.position="none",
                axis.text.x=element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of numbers of eQTL target genes", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## look at the distributions of distances to the TSSs
    dist_to_tss_per_snp <- ddply(uniq_snp_eqtl_df, .(tag_name, tag_no_rsid, rsID, gene_name), summarize,
                                 tss_distance = unique(tss_distance))
    step_size <- 100000
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_dist_to_tss_per_snp_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))  
    print(ggplot(dist_to_tss_per_snp,
                 aes_string(x=TAG_VAR, y="tss_distance", fill=TAG_VAR)) +
          xlab(TAG_LAB) + ylab("Distance (bp)") +
          scale_fill_hue(h=c(180, 270)) +
          scale_y_continuous(breaks=seq(round_any(min(dist_to_tss_per_snp$tss_distance), accuracy=step_size, f=floor),
                               round_any(max(dist_to_tss_per_snp$tss_distance), accuracy=step_size, f=ceiling)
                                 , by=step_size), labels=comma) + 
          theme_bw() + geom_boxplot() +
          theme(legend.position="none",
                axis.text.x=element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of distances to target gene TSS per tag region", r2_thresh, dist_thresh, out_subtitle))
    dev.off()    

    make_graphic(paste0(outdir, 'plots/', prefix, '_abs_dist_to_tss_per_snp_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))  
    print(ggplot(dist_to_tss_per_snp, aes_string(x=TAG_VAR, fill=TAG_VAR) +
                 aes(y=abs(tss_distance))) +
          xlab(TAG_LAB) + ylab("Absolute value of distance (bp)") +
          scale_fill_hue(h=c(180, 270)) +
          scale_y_continuous(breaks=seq(0,
                                 round_any(max(dist_to_tss_per_snp$tss_distance), accuracy=step_size, f=ceiling)
                                 , by=step_size), labels=comma) + 
          theme_bw() + geom_boxplot() +
          theme(legend.position="none",
                axis.text.x=element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of absolute value of distances to target gene TSS per tag region", r2_thresh, dist_thresh, out_subtitle))
    dev.off()    
    
    ## -------------------------
    ## tissue analyses
    ## -------------------------
    ## -------------------------
    ## look at the number of eQTLs contained in each tissue
    num_eqtls_per_tissue <- ddply(uniq_snp_eqtl_df, .(tissue, eqtl_class), summarize,
                                  num_tissue_eqtls = length(unique(rsID)))
    
    ## add 0-counts for tissues with no eQTLs
    ## get the index
    zero_eqtl_tiss <- which(!(gtex_category_df$tissue %in% num_eqtls_per_tissue$tissue))
    num_eqtls_per_tissue <- rbind(num_eqtls_per_tissue,
                                  data.frame(tissue = gtex_category_df$tissue[zero_eqtl_tiss],
                                             eqtl_class = gtex_category_df$Class[zero_eqtl_tiss],
                                             num_tissue_eqtls = 0))
    
    ## make sure the tissues are sorted by class, decreasing because we flip it
    tissue_order <- order(num_eqtls_per_tissue$eqtl_class, num_eqtls_per_tissue$tissue, decreasing=T)
    num_eqtls_per_tissue$tissue <- factor(num_eqtls_per_tissue$tissue, ordered=T,
                                          levels=num_eqtls_per_tissue$tissue[tissue_order])

    max_eqtl_round <- round_any(max(num_eqtls_per_tissue$num_tissue_eqtls), 100, f=ceiling)
    ## get the text colors:
    class_text_col <- category_colors[num_eqtls_per_tissue$eqtl_class[tissue_order]]

    ## plot the numbers per tissue
    make_graphic(paste0(outdir, 'plots/', prefix, '_eqtl_overlap_nums_bytiss_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"),
                 height_ratio=1.5)
    print(ggplot(num_eqtls_per_tissue,
                 aes(x=tissue, y=num_tissue_eqtls, fill=eqtl_class)) +
          cat_col_scale + 
          xlab("Cell or tissue type") + ylab("Number of SNPs that are also eQTLs") +
          coord_flip() + 
          scale_y_continuous(breaks=seq(0, max_eqtl_round, by=50), limits=c(0, max_eqtl_round), expand=c(0, 0)) + 
          theme_bw() + geom_bar(stat="identity", position="stack") +
          guides(fill = guide_legend(title="Tissue category", ncol=1)) +
          theme(legend.position="right",
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE*.5),
                axis.text.y = element_text(colour=class_text_col, size=AXIS_TEXT_Y_SIZE*.5),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                legend.title = element_text(size=TITLE_SIZE*0.75),                     
                legend.text = element_text(size=LEGEND_TEXT_SIZE)) +
          plot_title("Numbers of eQTLs across tissue amd cell types", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## -------------------------
    ## also look at the distribution of numbers within each tissue class
    ## need to factor correctly:
    num_eqtls_per_tissue$eqtl_class <- factor(num_eqtls_per_tissue$eqtl_class, ordered=T,
                                                  levels=sort(unique(num_eqtls_per_tissue$eqtl_class), decreasing=T))    
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_eqtl_overlap_num_distribution_byclass_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"),
                 height_ratio=1.5)
    print(ggplot(num_eqtls_per_tissue,
                 aes(x=eqtl_class, y=num_tissue_eqtls, fill=eqtl_class)) +
          cat_col_scale + 
          xlab("Tissue category") + ylab("Number of SNPs that are also eQTLs") +
          coord_flip() + 
          scale_y_continuous(breaks=seq(0, max_eqtl_round, by=50), limits=c(0, max_eqtl_round), expand=c(0.01, 0)) + 
          theme_bw() + geom_boxplot() +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                axis.title = element_text(size=TITLE_SIZE*0.75),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of numbers of eQTLs across tissue categories", r2_thresh, dist_thresh, out_subtitle, strwrap_width=30))
    dev.off()
    
    ## -------------------------
    ## similar analysis for tissue class
    num_eqtls_per_tiss_class <- ddply(uniq_snp_eqtl_df, .(eqtl_class), summarize,
                                      num_tiss_class_eqtls = length(unique(rsID)))
    ## add 0-counts for tissue classes with no eQTLs
    zero_eqtl_tiss_class <- which(!(gtex_category_df$Class %in% num_eqtls_per_tiss_class$eqtl_class))
    ## only do this if we need to!
    if (length(zero_eqtl_tiss_class) > 0) {
        num_eqtls_per_tiss_class <- rbind(num_eqtls_per_tiss_class,
                                          unique(data.frame(eqtl_class = gtex_category_df$Class[zero_eqtl_tiss_class],
                                                            num_tiss_class_eqtls = 0)))
    }

    num_eqtls_per_tiss_class$eqtl_class <- factor(num_eqtls_per_tiss_class$eqtl_class, ordered=T,
                                                  levels=sort(num_eqtls_per_tiss_class$eqtl_class, decreasing=T))

    max_eqtl_class_round <- round_any(max(num_eqtls_per_tiss_class$num_tiss_class_eqtls), 100, f=ceiling)

    make_graphic(paste0(outdir, 'plots/', prefix, '_eqtl_overlap_nums_bytiss_class_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"),
                 height_ratio=1.5, width_ratio = 1.2)
    print(ggplot(num_eqtls_per_tiss_class,
                 aes(x=eqtl_class, y=num_tiss_class_eqtls, fill=eqtl_class)) +
          cat_col_scale +
          geom_bar(stat="identity", position="stack") +
          xlab("Tissue category") + ylab("Number of unique SNPs that are also eQTLs") +
          coord_flip() + 
          scale_y_continuous(breaks=seq(0, max_eqtl_class_round, by=50), limits=c(0, max_eqtl_class_round), expand=c(0, 0)) + 
          theme_bw() +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Numbers of eQTLs across tissue categories", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## -------------------------
    ## make a heatmap of tag SNP - tissue enrichment
    ## note that this counts each SNP uniquely for each tissue
    tag_region_tissue_counts <- ddply(uniq_snp_eqtl_df, .(tag_name, tissue, eqtl_class), summarize,
                                      num_eqtl_snps = length(unique(rsID)))
    ## now we have to find all the tissue - tag SNP combinations with no counts and add them
    tag_region_tissue_counts <- merge(expand.grid(tag_name=all_regions, tissue=gtex_category_df$tissue, stringsAsFactors = F),
                               tag_region_tissue_counts, 
                               by=c("tag_name", "tissue"), all.x = T)
    tag_region_tissue_counts$eqtl_class <- gtex_category_df$Class[match(tag_region_tissue_counts$tissue, gtex_category_df$tissue)]
    tag_region_tissue_counts$num_eqtl_snps[is.na(tag_region_tissue_counts$num_eqtl_snps)] <- 0
    
    tag_region_tissue_counts$tissue <- gsub("_", " ", tag_region_tissue_counts$tissue)

    ## make sure that the tissues are ordered correctly
    ordered_uniq_tiss <- unique(tag_region_tissue_counts$tissue[order(tag_region_tissue_counts$eqtl_class, tag_region_tissue_counts$tissue)])

    tag_region_tissue_counts$tissue <- factor(tag_region_tissue_counts$tissue, ordered=T,
                                              levels=ordered_uniq_tiss)

    ## now get the axis text colors
    ordered_heatmap_classes <- tag_region_tissue_counts$eqtl_class[match(ordered_uniq_tiss, tag_region_tissue_counts$tissue)]
    heatmap_text_cols <- category_colors[ordered_heatmap_classes]

    ## raw count heatmap
    ## white to blue
    heatmap_col_lims <- c("#FFFFFF", "#08306B")

    ## add annotation of tag regions without rsIDs
    tag_region_tissue_counts$tag_no_rsid <- gsub(":rs.*", "", tag_region_tissue_counts$tag_name)
    
    ## -------------------------
    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_eqtl_count_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)                   
    print(ggplot(tag_region_tissue_counts, aes_string(x="tissue", y=TAG_VAR)) +
          geom_tile(aes(fill=num_eqtl_snps), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar", name="# eQTL SNPs") +
          theme_bw() + 
          xlab("GTEx eQTL Data Source") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of eQTL SNPs by tag region and GTEx tissue", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75, colour=heatmap_text_cols),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    ## -------------------------
    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_eqtl_log10_count_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)                   
    print(ggplot(tag_region_tissue_counts, aes_string(x="tissue", y=TAG_VAR)) +
          geom_tile(aes(fill=log10(num_eqtl_snps)), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar", name="log10(# eQTL SNPs)") +
          theme_bw() + 
          xlab("GTEx eQTL Data Source") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of eQTL SNPs by tag region and GTEx tissue", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*0.75, colour=heatmap_text_cols),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                 
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    ## now do this counting by eQTL class
    tag_region_class_counts <- ddply(uniq_snp_eqtl_df, .(tag_name, eqtl_class), summarize,
                                      num_eqtl_snps = length(unique(rsID)))
    ## get the class-region combos with no counts
    tag_region_class_counts <- merge(expand.grid(tag_name=all_regions, eqtl_class=unique(gtex_category_df$Class), stringsAsFactors = F),
                                     tag_region_class_counts, by=c("tag_name", "eqtl_class"),
                                     all.x=T)
    tag_region_class_counts$num_eqtl_snps[is.na(tag_region_class_counts$num_eqtl_snps)] <- 0

    tag_region_class_counts$tag_no_rsid <- gsub(":rs.*", "", tag_region_class_counts$tag_name)
    
    ## don't need to re-order, just get the text colors
    heatmap_text_cols <- category_colors[sort(unique(tag_region_class_counts$eqtl_class))]

    ## -------------------------
    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_eqtl_count_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)                   
    print(ggplot(tag_region_class_counts, aes_string(x="eqtl_class", y=TAG_VAR)) +
          geom_tile(aes(fill=num_eqtl_snps), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar", name="# eQTL SNPs") +
          theme_bw() + 
          xlab("GTEx eQTL Data Source Category") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of eQTL SNPs by tag region and GTEx tissue category", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE,
                    colour=heatmap_text_cols),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                     
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    ## make this with black tissue category labels
    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_eqtl_count_heatmap_black_labels_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)                   
    print(ggplot(tag_region_class_counts, aes_string(x="eqtl_class", y=TAG_VAR)) +
          geom_tile(aes(fill=num_eqtl_snps), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar", name="# eQTL SNPs") +
          theme_bw() + 
          xlab("GTEx eQTL Data Source Category") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of eQTL SNPs by tag region and GTEx tissue category", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                     
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    ## also make one with greyscale color theme; we set the scale to start at '1' so we can
    ## see even the weak signals
    greyscale_col_lims <- c("#FFFFFF", "gray90", "black")

    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_eqtl_count_heatmap_greyscale_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)                   
    print(ggplot(tag_region_class_counts, aes_string(x="eqtl_class", y=TAG_VAR)) +
          geom_tile(aes(fill=num_eqtl_snps), colour=greyscale_col_lims[1]) +
          scale_fill_gradientn(colours=greyscale_col_lims,
                               values=rescale(c(0, 1, max(tag_region_class_counts$num_eqtl_snps))),                                   
                               name="# eQTL SNPs", na.value=greyscale_col_lims[1],
                               guide="colorbar") +
          theme_bw() +
          xlab("GTEx eQTL Data Source Category") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of eQTL SNPs by tag region and GTEx tissue category", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                     
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()        
    
    ## -------------------------
    make_graphic(paste0(outdir, 'plots/', prefix, '_tagregion_by_tissue_class_eqtl_log10_count_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)                   
    print(ggplot(tag_region_class_counts, aes_string(x="eqtl_class", y=TAG_VAR)) +
          geom_tile(aes(fill=log10(num_eqtl_snps)), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(low=heatmap_col_lims[1], high=heatmap_col_lims[2], na.value=heatmap_col_lims[1], guide="colorbar", name="log10(# eQTL SNPs)") +
          theme_bw() + 
          xlab("GTEx eQTL Data Source Category") + ylab(TAG_LAB) +
          plot_title("Heatmap of number of eQTL SNPs by tag region and GTEx tissue category", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE,
                    colour=heatmap_text_cols),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.title = element_text(size=TITLE_SIZE*0.75),                     
                legend.text = element_text(size=LEGEND_TEXT_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)))
    dev.off()

    ## -------------------------
    ## now look at the distributions of the number of tissues and tissue classes for each eQTL SNP
    num_tissues_per_eqtl <- melt(ddply(uniq_snp_eqtl_df, .(tag_name, rsID), summarize,
                                       "# of tissues" = length(unique(tissue)),
                                       "# of tissue classes" = length(unique(eqtl_class))),
                                 id.vars=c("tag_name", "rsID"))
    num_tissues_per_eqtl$tag_no_rsid <- gsub(":rs.*", "", num_tissues_per_eqtl$tag_name)
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_tissues_per_eqtl_snp_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"), width_ratio = 2.0)  
    print(ggplot(num_tissues_per_eqtl, aes(x=variable, y=value, fill=variable)) +
          xlab("Tag SNP region") + ylab("Number of Tissues / Tissue Classes") +
          theme_bw() + geom_boxplot() +
          facet_grid(reformulate(TAG_VAR)) + 
          scale_y_continuous(breaks=seq(0, max(num_tissues_per_eqtl$value), by=1)) + 
          theme(legend.position="none",
                axis.text.x=element_text(angle=75, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE*0.5),
                strip.text = element_text(size=LEGEND_TEXT_SIZE*0.75), 
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) + 
          plot_title("Distributions of numbers of tissues and tissue classes for eQTL SNPs in each tag region", r2_thresh, dist_thresh, out_subtitle, strwrap_width=70))
    dev.off()    
    
}

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/gtex_eqtl_overlap/")
## tagsnp_cols <- c("tag_rsID", "tag_pos", "tag_MAF", "R2", "Dprime")
## gtex_class_file <- '/home/alexaml/data/GTEx/gtex_classes.txt'
