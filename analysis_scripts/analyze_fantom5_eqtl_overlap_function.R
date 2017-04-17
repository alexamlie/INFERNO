## analyze_fantom5_eqtl_overlap_function.R
## alex amlie-wolf 01/26/16
## a script to compare eQTL results with FANTOM5 overlaps
## part of the automatic analysis

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, "/fantom5_eqtl_overlap/")
## enh_overlap_type <- "locus"

class_combo_heatmap <- function(input_uniq_snp_df, prefix, outdir, r2_thresh, dist_thresh, extra_title="") {
    ## we need to count how many times each combination occurs:
    eqtl_class_counts <- cbind(type="eqtl", melt(strsplit(input_uniq_snp_df$eqtl_classes, ",")))
    enh_class_counts <- cbind(type="enh", melt(strsplit(input_uniq_snp_df$enh_classes, ",")))

    combo_classes <- ddply(rbind(eqtl_class_counts, enh_class_counts), .(L1), function(x) {
        ## return the data frame containing all the class combinations
        expand.grid(eqtl_class=x$value[x$type=='eqtl'], enh_class=x$value[x$type=='enh'],
                    stringsAsFactors = FALSE)
    })

    combo_class_counts <- melt(table(combo_classes[,c("eqtl_class", "enh_class")]))

    ## sort them into the right order
    ## they need to be opposite in order to get the correct heatmap orientation
    ordered_combo_classes <- sort(unique(as.character(combo_class_counts$enh_class)), decreasing=T)
    combo_class_counts$enh_class <- factor(combo_class_counts$enh_class, ordered=T,
                                           levels=rev(ordered_combo_classes))
    combo_class_counts$eqtl_class <- factor(combo_class_counts$eqtl_class, ordered=T,
                                           levels=ordered_combo_classes)

    ## get the axis colors (using the 'global' variable)
    eqtl_text_cols <- category_colors[ordered_combo_classes]
    enh_text_cols <- category_colors[rev(ordered_combo_classes)]

    ## raw count heatmap
    ## white to blue
    heatmap_col_lims <- c("#FFFFFF", "#08306B")
    label_max <- round_any(max(combo_class_counts$value), 5, ceiling)

    make_graphic(paste0(outdir, prefix, '_enh_eqtl_category_combo_count_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)
    print(ggplot(combo_class_counts, aes(x=enh_class, y=eqtl_class)) +
          geom_tile(aes(fill=value), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(name="Number of SNPs", low=heatmap_col_lims[1],
                              high=heatmap_col_lims[2], na.value=heatmap_col_lims[1],
                              limits=c(0, label_max),
                              guide="colorbar", breaks=seq(0, label_max, by=5)) +
          scale_x_discrete(position="top") + 
          theme_bw() +
          xlab("FANTOM5 Enhancer Tissue Category") + ylab("GTEx eQTL Tissue Category") +
          plot_title(paste0("Heatmap of number of SNPs overlapping both enhancers and eQTLs across tissue category combinations ", extra_title), r2_thresh, dist_thresh, out_subtitle) + 
          theme(axis.text.x = element_text(angle=-45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE, colour=enh_text_cols),
                axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE, colour=eqtl_text_cols),
                title=element_text(size=TITLE_SIZE), legend.text=element_text(size=LEGEND_TEXT_SIZE),
                plot.title=element_text(hjust=0.5)))
    dev.off()

    make_graphic(paste0(outdir, prefix, '_enh_eqtl_category_combo_count_heatmap_withtext_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)
    print(ggplot(combo_class_counts, aes(x=enh_class, y=eqtl_class)) +
          geom_tile(aes(fill=value), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(name="Number of SNPs", low=heatmap_col_lims[1],
                              high=heatmap_col_lims[2], na.value=heatmap_col_lims[1],
                              limits=c(0, label_max),
                              guide="colorbar", breaks=seq(0, label_max, by=5)) +
          theme_bw() +
          scale_x_discrete(position="top") + 
          xlab("FANTOM5 Enhancer Tissue Category") + ylab("GTEx eQTL Tissue Category") +
          plot_title(paste0("Heatmap of number of SNPs overlapping both enhancers and eQTLs across tissue category combinations ", extra_title), r2_thresh, dist_thresh, out_subtitle) +
          theme(axis.text.x = element_text(angle=-45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE, colour=enh_text_cols),
                axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE, colour=eqtl_text_cols),
                title=element_text(size=TITLE_SIZE), legend.text=element_text(size=LEGEND_TEXT_SIZE),
                plot.title=element_text(hjust=0.5)) + 
          ## add labels of the SNP counts
          geom_text(aes(x=enh_class, y=eqtl_class, label=value), color="black", size=4))
    dev.off()

    ## we also want to cluster this heatmap:
    combo_dist <- dist(table(combo_classes[,c("eqtl_class", "enh_class")]), method="manhattan")
    combo_clust <- hclust(combo_dist)

    clustered_order <- labels(combo_dist)[combo_clust$order]
    ## just change the factor level to use the clustering:
    combo_class_counts$enh_class <- factor(combo_class_counts$enh_class, ordered=T,
                                           levels=rev(clustered_order))
    combo_class_counts$eqtl_class <- factor(combo_class_counts$eqtl_class, ordered=T,
                                           levels=clustered_order)

    ## get the axis colors
    eqtl_text_cols <- category_colors[clustered_order]
    enh_text_cols <- category_colors[rev(clustered_order)]

    ## clustered heatmap (this is the same data, we just changed the order)
    make_graphic(paste0(outdir, prefix, '_enh_eqtl_category_clustered_combo_count_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.0)
    print(ggplot(combo_class_counts, aes(x=enh_class, y=eqtl_class)) +
          geom_tile(aes(fill=value), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(name="Number of SNPs", low=heatmap_col_lims[1],
                              high=heatmap_col_lims[2], na.value=heatmap_col_lims[1],
                              limits=c(0, label_max),
                              guide="colorbar", breaks=seq(0, label_max, by=5)) +
          scale_x_discrete(position="top") + 
          theme_bw() +
          xlab("FANTOM5 Enhancer Tissue Category") + ylab("GTEx eQTL Tissue Category") +
          plot_title(paste0("Heatmap of number of SNPs overlapping both enhancers and eQTLs across tissue category combinations ", extra_title), r2_thresh, dist_thresh, out_subtitle) + 
                                        #use 90+1e-09 for vertical
          theme(axis.text.x = element_text(angle=-45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE, colour=enh_text_cols),
                axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE, colour=eqtl_text_cols),
                title=element_text(size=TITLE_SIZE), legend.text=element_text(size=LEGEND_TEXT_SIZE),
                plot.title=element_text(hjust=0.5)))    
    dev.off()

    make_graphic(paste0(outdir, prefix, '_enh_eqtl_category_clustered_combo_count_heatmap_withtext_',
                        r2_thresh, "_ld_", dist_thresh,                        
                        "_dist"), width_ratio = 2.0)
    print(ggplot(combo_class_counts, aes(x=enh_class, y=eqtl_class)) +
          geom_tile(aes(fill=value), colour=heatmap_col_lims[1]) +
          scale_fill_gradient(name="Number of SNPs", low=heatmap_col_lims[1],
                              high=heatmap_col_lims[2], na.value=heatmap_col_lims[1],
                              limits=c(0, label_max),
                              guide="colorbar", breaks=seq(0, label_max, by=5)) +
          theme_bw() +
          scale_x_discrete(position="top") + 
          xlab("FANTOM5 Enhancer Tissue Category") + ylab("GTEx eQTL Tissue Category") +
          plot_title(paste0("Heatmap of number of SNPs overlapping both enhancers and eQTLs across tissue category combinations ", extra_title), r2_thresh, dist_thresh, out_subtitle) +
          theme(axis.text.x = element_text(angle=-45, vjust=1, hjust=1, size=AXIS_TEXT_X_SIZE, colour=enh_text_cols),
                axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE, colour=eqtl_text_cols),
                title=element_text(size=TITLE_SIZE), legend.text=element_text(size=LEGEND_TEXT_SIZE),
                plot.title=element_text(hjust=0.5)) + 
          ## add labels of the SNP counts
          geom_text(aes(x=enh_class, y=eqtl_class, label=value), color="black", size=4))
    dev.off()
}

analyze_fantom5_eqtl_overlap <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh, enh_overlap_type, fantom5_class_file, gtex_class_file) {
    dir.create(paste0(outdir, "/tables/"), F, T)
    dir.create(paste0(outdir, "/plots/"), F, T)

    ## read in the enhancer overlap data file
    if (enh_overlap_type == "midpoint") {
        enh_overlap_target_file <- paste0(datadir, '/correlation_enh_targets/',
                                          prefix, "_", r2_thresh,
                                          "_ld_cutoff_snps_within_", dist_thresh,
                                          "_fantom5_midpoint_overlap_target_genes.txt")
        enh_overlap_target_df <- read.table(enh_overlap_target_file, header=T, sep="\t", quote="", as.is=T)
    } else if (enh_overlap_type == "locus") {
        enh_overlap_target_file <- paste0(datadir, '/correlation_enh_targets/',
                                          prefix, "_", r2_thresh,
                                          "_ld_cutoff_snps_within_", dist_thresh,
                                          "_fantom5_locus_overlap_target_genes.txt")
        enh_overlap_target_df <- read.table(enh_overlap_target_file, header=T, sep="\t", quote="", as.is=T)
    }

    ## clean up the tissue names:
    enh_overlap_target_df$enh_source <- gsub("CL:[0-9]*_|UBERON:[0-9]*_", "", enh_overlap_target_df$enh_source)

    ## read in the tissue categories
    fantom5_category_df <- read.table(fantom5_class_file, header=T, sep="\t", quote="", as.is=T)
    ## make it easier to match
    fantom5_category_df$enh_tissue <- gsub("CL:[0-9]*_|UBERON:[0-9]*_|_expressed_enhancers.bed", "", fantom5_category_df$FANTOM5.File)

    ## add the tissue class column
    enh_overlap_target_df$enh_class <- fantom5_category_df$Class[match(enh_overlap_target_df$enh_source, fantom5_category_df$enh_tissue)]

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

    ## first get some summary information about each enhancer overlapping SNP
    ## to do this, extract the columns that are related to SNPs (and not to enhancers)
    enh_snp_cols <- colnames(enh_overlap_target_df)[!(colnames(enh_overlap_target_df) %in% c("enh_source", "enh_class", "enh_chr", "enh_start", "enh_end", "num_assocs", "refseqId","symbol", "r", "fdr"))]

    enh_target_df <- ddply(enh_overlap_target_df, enh_snp_cols,
                           summarize, num_enh_targets = length(unique(symbol[!is.na(symbol)])),
                           enh_targets = paste(sort(unique(symbol[!is.na(symbol)])), collapse=","),
                           num_enh_tissues = length(unique(enh_source)),
                           enh_tissues = paste(sort(unique(enh_source)), collapse=","),
                           num_enh_classes = length(unique(enh_class)),
                           enh_classes = paste(sort(unique(enh_class)), collapse=","))

    ## then get the same  summary information for each eQTL
    eqtl_snp_cols <- enh_snp_cols

        #colnames(eqtl_overlap_df)[!(colnames(eqtl_overlap_df) %in% c("tissue", "eqtl_class", "gene", "beta", "t_stat", "se", "p_value", "nom_thresh", "min_p", "gene_emp_p", "maf", "gene_name", "ref.1", "alt.1", "gene_source", "gene_type"))]

    eqtl_target_df <- ddply(eqtl_overlap_df, eqtl_snp_cols,
                            summarize, num_eqtl_genes = length(unique(gene_name)),
                            eqtl_genes = paste(sort(unique(gene_name)), collapse=","),
                            num_eqtl_tissues = length(unique(tissue)),
                            eqtl_tissues = paste(sort(unique(tissue)), collapse=","),
                            num_eqtl_classes = length(unique(eqtl_class)),
                            eqtl_classes = paste(sort(unique(eqtl_class)), collapse=","))

    ## now merge these:
    overlapping_snp_df <- join(enh_target_df, eqtl_target_df, type='inner')

    overlap_snp_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                             "_ld_cutoff_snps_within_", dist_thresh,
                               "_", enh_overlap_type,
                               "_fantom5_eqtl_overlap_gene_targets.txt")
    write.table(overlapping_snp_df, overlap_snp_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## also, get only the unique rsID : tag region combinations, to make it easier to parse
    ## use the 'global' tagsnp column variable
    ## we do this differently from the other scripts because we don't have individual enhs / eqtls
    uniq_snp_cols <- enh_snp_cols[!(enh_snp_cols %in% tagsnp_cols)]

    unique_snp_df <- ddply(overlapping_snp_df, uniq_snp_cols, summarize,
                           num_enh_targets = unique(num_enh_targets),
                           enh_targets=unique(enh_targets),
                           num_enh_tissues = unique(num_enh_tissues),
                           enh_tissues = unique(enh_tissues),
                           num_enh_classes = unique(num_enh_classes),
                           enh_classes = unique(enh_classes),
                           num_eqtl_genes = unique(num_eqtl_genes),
                           eqtl_genes = unique(eqtl_genes),
                           num_eqtl_tissues = unique(num_eqtl_tissues),
                           eqtl_tissues = unique(eqtl_tissues),
                           num_eqtl_classes = unique(num_eqtl_classes),
                           eqtl_classes = unique(eqtl_classes))

    unique_snp_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                             "_ld_cutoff_snps_within_", dist_thresh,
                               "_", enh_overlap_type,
                               "_fantom5_eqtl_unique_hits.txt")
    write.table(unique_snp_df, unique_snp_outf, quote=F, sep="\t", col.names=T, row.names=F)

    ## now write the gene sets for each tag region
    dir.create(paste0(outdir, 'tables/tag_region_gene_lists/'), F, T)

    for(this_tag in unique(c(eqtl_target_df$tag_name, enh_target_df$tag_name))) {
        enh_genes <- unique(unlist(lapply(enh_target_df$enh_targets[enh_target_df$tag_name==this_tag], strsplit, split=",")))
        eqtl_genes <- unique(unlist(lapply(eqtl_target_df$eqtl_genes[eqtl_target_df$tag_name==this_tag], strsplit, split=",")))

        all_target_genes <- unique(c(enh_genes, eqtl_genes))
        outf <- paste0(outdir, 'tables/tag_region_gene_lists/',
                       gsub(":|/", "_", this_tag), "_all_genes_", r2_thresh, "_ld_",
                       dist_thresh, "_dist.txt")
        ## write the file out
        cat(all_target_genes, file=outf, sep="\n")

        ## also do this separately for eqtl and enhancer
        outf <- paste0(outdir, 'tables/tag_region_gene_lists/',
                       gsub(":|/", "_", this_tag), "_eqtl_genes_", r2_thresh, "_ld_",
                       dist_thresh, "_dist.txt")
        ## write the file out
        cat(unique(eqtl_genes), file=outf, sep="\n")

        outf <- paste0(outdir, 'tables/tag_region_gene_lists/',
                       gsub(":|/", "_", this_tag), "_enh_genes_", r2_thresh, "_ld_",
                       dist_thresh, "_dist.txt")
        ## write the file out
        cat(unique(enh_genes), file=outf, sep="\n")        
    }
    
    ## plots:
    ## ---------------
    ## heatmap of combinations of categories
    class_combo_heatmap(unique_snp_df, prefix, paste0(outdir, "plots/"), r2_thresh, dist_thresh)

    ## now make these heatmaps for each tag region
    for(tag_region in unique(unique_snp_df$tag_name)) {
        tag_region_name <- strsplit(tag_region, ":")[[1]][1]
        tag_region_name <- gsub("/", "_", tag_region_name)
        tag_outdir <- paste0(outdir, "/plots/", tag_region_name, "/")
        dir.create(tag_outdir, F, T)
        class_combo_heatmap(subset(unique_snp_df, tag_name==tag_region),
                            paste0(prefix, "_", tag_region_name),
                            tag_outdir, r2_thresh, dist_thresh,
                            extra_title=paste(tag_region, "region"))
    }

    ## ---------------
    ## make a heatmap of the amount of functional support for linked SNPs in each tissue category
    ## just get all the unique combinations of SNPs and tissue classes
    enh_snp_category_combos <- melt(unique(enh_overlap_target_df[,c("rsID", "pos", "tag_name", "enh_class")]),
                                    id.vars=c("rsID", "tag_name", "pos"),
                                    variable.name="data_source", value.name="tissue_class")
    eqtl_snp_category_combos <- melt(unique(eqtl_overlap_df[,c("rsID", "pos", "tag_name", "eqtl_class")]),
                                     id.vars=c("rsID", "tag_name", "pos"),
                                     variable.name="data_source", value.name="tissue_class")

    ## now get the supporting data sets per SNP and tissue category
    snp_support_df <- ddply(rbind(enh_snp_category_combos, eqtl_snp_category_combos),
                            .(rsID, pos, tag_name), function(x) {
                        tissue_counts <- table(x$tissue_class)
                        full_support <- tissue_counts==2
                        ## make nicer names for the support variable
                        support_vec <- ifelse(as.vector(x[match(names(tissue_counts)[!full_support], x$tissue_class),"data_source"])=="eqtl_class", "eQTL Only", "Enhancer Only")
                        ## return the DF
                        data.frame(tissue_class=c(names(tissue_counts)[full_support],
                                       names(tissue_counts[!full_support])),
                                   support=c(rep("Enhancer and eQTL", times=sum(full_support)), support_vec), stringsAsFactors = F)
                    })

    ## next, we need to add annotations for all the tissues that weren't observed for each SNP
    snp_support_df_full <- merge(expand.grid(rsID=unique(snp_support_df$rsID),
                                                     tissue_class=union(fantom5_category_df$Class, gtex_category_df$Class), stringsAsFactors = F),
                                    snp_support_df, by.x=c("rsID", "tissue_class"),
                                    all.x=T)

    ## fix all the NAs we just generated
    snp_support_df_full$pos <- snp_support_df$pos[match(snp_support_df_full$rsID, snp_support_df$rsID)]
    snp_support_df_full$tag_name <- snp_support_df$tag_name[match(snp_support_df_full$rsID, snp_support_df$rsID)]
    snp_support_df_full$support[is.na(snp_support_df_full$support)] <- "Neither"

    ## now, make sure the tissues are in order
    ordered_support_classes <- sort(unique(snp_support_df_full$tissue_class), decreasing=T)
    snp_support_df_full$tissue_class <- factor(snp_support_df_full$tissue_class,
                                                       ordered=T, levels=ordered_support_classes)
    support_text_cols <- category_colors[ordered_support_classes]

    ## we also want to color SNPs by their tag region, so order the SNPs by tag region, then
    ## position, then rsID (if there are any at the same position..)
    snp_support_ordered <- order(snp_support_df_full$tag_name, snp_support_df_full$pos, snp_support_df_full$rsID, decreasing=F)

    ordered_rsIDs <- unique(snp_support_df_full$rsID[snp_support_ordered])
    snp_support_df_full$rsID <- factor(snp_support_df_full$rsID, ordered=T,
                                               levels=ordered_rsIDs)

    tag_name_cols <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(snp_support_df$tag_name)))
    names(tag_name_cols) <- unique(snp_support_df$tag_name)

    snp_text_cols <- tag_name_cols[snp_support_df_full$tag_name[match(ordered_rsIDs, snp_support_df_full$rsID)]]

    ## full tissue heatmap
    ## four possibilities: no overlap, enhancer only (red), eQTL only (green), and both (blue)
    support_heatmap_cols <- c("#FFFFFF", "#95001A", "#88BDAD", "#01256E")
    names(support_heatmap_cols) <- c("Neither", "Enhancer Only", "eQTL Only", "Enhancer and eQTL")

    make_graphic(paste0(outdir, 'plots/', prefix, '_snp_num_support_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 3.5, height_ratio = 2.0)
    print(ggplot(snp_support_df_full, aes(x=rsID, y=tissue_class)) +
          geom_tile(aes(fill=factor(support))) +
          scale_fill_manual(name="Sources of functional support",
                            values = support_heatmap_cols) +
          ## add invisible "points" just to get the legend for the tag SNP region
          geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
          scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
          theme_bw() +
          xlab("Linked SNP") + ylab("Tissue Category") +
          scale_x_discrete(position="top") +
          plot_title("Heatmap of sources of functional support for SNPs across tag regions and tissue categories", r2_thresh, dist_thresh, out_subtitle) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                 fill = guide_legend(override.aes=list(size=10))) +
          theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=4, colour=snp_text_cols),
                axis.text.y = element_text(hjust=1, size=30, colour=support_text_cols),
                axis.title = element_text(size=30), 
                title = element_text(size=25), legend.text = element_text(size=20),
                legend.position="right", plot.title=element_text(hjust = 0.5)))
    dev.off()

    ## plot only the tissues that are in the data
    present_tissues <- as.vector(sort(unique(snp_support_df_full$tissue_class[snp_support_df_full$support!="Neither"])))
    present_tissue_cols <- category_colors[present_tissues]

    make_graphic(paste0(outdir, 'plots/', prefix, '_snp_num_support_heatmap_present_tissues_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 3.5, height_ratio = 2.0)
    print(ggplot(snp_support_df_full[snp_support_df_full$tissue_class %in% present_tissues,],
               aes(x=rsID, y=tissue_class)) +
          geom_tile(aes(fill=factor(support))) +
          scale_fill_manual(name="Sources of functional support",
                            values = support_heatmap_cols) +
          scale_x_discrete(position="top") + 
          ## add invisible "points" just to get the legend for the tag SNP region
          geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
          scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
          theme_bw() +
          xlab("Linked SNP") + ylab("Tissue Category") +
          plot_title("Heatmap of sources of functional support for SNPs across tag regions and tissue categories (only categories with SNP overlap)", r2_thresh, dist_thresh, out_subtitle) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                 fill = guide_legend(override.aes=list(size=10))) +
          theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=4, colour=snp_text_cols),
                axis.text.y = element_text(hjust=1, size=30, colour=present_tissue_cols),
                axis.title = element_text(size=25),
                title=element_text(size=25), plot.title=element_text(hjust=0.5),
                legend.text = element_text(size=20), legend.position="right"))
    dev.off()

    ## now make plots for each tag region
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

        make_graphic(paste0(tag_outdir, prefix, '_', tag_region_name, '_snp_num_support_heatmap_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 3.5, height_ratio = 2.0)
        print(ggplot(this_tag_support_df, aes(x=rsID, y=tissue_class)) +
              geom_tile(aes(fill=factor(support))) +
              scale_fill_manual(name="Sources of functional support",
                                values = support_heatmap_cols) +
              scale_x_discrete(position="top") + 
              ## add invisible "points" just to get the legend for the tag SNP region
              geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
              scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
              theme_bw() +
              ylab("Tissue Category") + xlab("Linked SNP") +
              plot_title(paste0("Heatmap of sources of functional support for SNPs in tag region ", tag_region_name, " across all tissue categories"), r2_thresh, dist_thresh, out_subtitle) +
              guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                     fill = guide_legend(override.aes=list(size=10))) +
              theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=4, colour=snp_text_cols),
                    axis.text.y = element_text(hjust=1, size=30, colour=present_tissue_cols),
                    axis.title = element_text(size=25),
                    title=element_text(size=25), plot.title=element_text(hjust=0.5),
                    legend.text = element_text(size=20), legend.position="right"))
        dev.off()

        ## plot only the tissues that are in the data
        present_tissues <- as.vector(sort(unique(this_tag_support_df$tissue_class[this_tag_support_df$support!="Neither"])))
        present_tissue_cols <- category_colors[present_tissues]

        make_graphic(paste0(tag_outdir, prefix, '_', tag_region_name, '_snp_num_support_heatmap_present_tissues_',
                            r2_thresh, "_ld_", dist_thresh,
                            "_dist"), width_ratio = 3.5, height_ratio = 2.0)
        print(ggplot(this_tag_support_df[this_tag_support_df$tissue_class %in% present_tissues,],
                     aes(x=rsID, y=tissue_class)) +
              geom_tile(aes(fill=factor(support))) +
              scale_fill_manual(name="Sources of functional support",
                                values = support_heatmap_cols) +
              scale_x_discrete(position="top") + 
              ## add invisible "points" just to get the legend for the tag SNP region
              geom_point(aes(shape=NA, colour=factor(tag_name)), na.rm=T) +
              scale_color_manual(name="Tag SNP region (axis text color)", values=tag_name_cols) +
              theme_bw() +
              ylab("Tissue Category") + xlab("Linked SNP") +
              plot_title(paste0("Heatmap of sources of functional support for SNPs in tag region ", tag_region_name, " across tissue categories with any SNP overlap"), r2_thresh, dist_thresh, out_subtitle) +
              guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                     fill = guide_legend(override.aes=list(size=10))) +
              theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=1, size=4, colour=snp_text_cols),
                    axis.text.y = element_text(hjust=1, size=30, colour=present_tissue_cols),
                    axis.title = element_text(size=25),
                    title=element_text(size=25), plot.title=element_text(hjust=0.5),
                    legend.text = element_text(size=20), legend.position="right"))
        dev.off()

    }

    ## make a summary plot with high-level tissue vs tag region support
    tag_region_support_df <- ddply(snp_support_df_full, .(tissue_class, tag_name), function(x) {
        if("Enhancer and eQTL" %in% x$support) {
            support <- "Overlapping enhancer and eQTL"
        } else if("Enhancer Only" %in% x$support & "eQTL Only" %in% x$support) {
            support <- "Enhancer and eQTL, non-overlapping"
        } else if("Enhancer Only" %in% x$support) {
            support <- "Enhancer Only"
        } else if("eQTL Only" %in% x$support) {
            support <- "eQTL Only"
        } else {
            support <- "Neither"
        }
        return(data.frame(support = support))
    })

    ## make sure this is in order
    tag_region_support_df$tissue_class <- factor(tag_region_support_df$tissue_class,
                                                 ordered=T, levels=ordered_support_classes)
    support_text_cols <- category_colors[ordered_support_classes]

    tag_region_support_df$tag_name <- factor(tag_region_support_df$tag_name, ordered=T,
                                             levels=sort(unique(tag_region_support_df$tag_name), decreasing=F))
    tag_region_support_df$support <- factor(tag_region_support_df$support, ordered=T,
                                            levels=c("Neither", "Enhancer Only", "eQTL Only", "Enhancer and eQTL, non-overlapping", "Overlapping enhancer and eQTL"))
    tag_region_support_df$tag_no_rsid <- gsub(":rs.*", "", tag_region_support_df$tag_name)

    ## now we have five possibilities:
    tag_support_cols <- c("#FFFFFF", "#95001A", "#88BDAD", "#0450E8", "#01256E")
    names(tag_support_cols) <- c("Neither", "Enhancer Only", "eQTL Only", "Enhancer and eQTL, non-overlapping", "Overlapping enhancer and eQTL")

    make_graphic(paste0(outdir, 'plots/', prefix, '_tag_region_support_heatmap_',
                        r2_thresh, "_ld_", dist_thresh,
                        "_dist"), width_ratio = 2.5, height_ratio = 2.0)
    print(ggplot(tag_region_support_df, aes_string(x=TAG_VAR, y="tissue_class")) +
          geom_tile(aes(fill=factor(support))) +
          scale_fill_manual(name="Sources of functional support",
                            values = tag_support_cols) +
          theme_bw() +
          xlab(TAG_LAB) + ylab("Tissue Category") +
          plot_title("Heatmap of sources of functional support across tag regions and tissue categories", r2_thresh, dist_thresh, out_subtitle) +
          guides(colour = guide_legend(override.aes=list(shape=15, size=10)),
                 fill = guide_legend(override.aes=list(size=10), ncol=2)) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE*1.5),
                axis.text.y = element_text(hjust=1, size=AXIS_TEXT_Y_SIZE*1.5, colour=support_text_cols),
                title = element_text(size=TITLE_SIZE*1.5), plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size=LEGEND_TEXT_SIZE*1.75),
                legend.position="bottom"))
    dev.off()

}
