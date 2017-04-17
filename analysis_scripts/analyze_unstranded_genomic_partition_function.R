## analyze_unstranded_genomic_partition_function.R
## alex amlie-wolf 01/26/2016
## a script that looks at the unstranded genomic partition results 
## part of the automatic analysis, so this is just a function

analyze_unstranded_genomic_partition <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh) {
    dir.create(paste0(outdir, 'plots/'), F, T)
    dir.create(paste0(outdir, 'tables/'), F, T)
    
    partition_f <- paste0(datadir, prefix, "_", r2_thresh, "_ld_cutoff_snps_within_",
                          dist_thresh, "_unstranded_partition_summary.txt") 

    partition_df <- read.table(partition_f, header=T, sep="\t", quote="", as.is=T)
    ## replace the labels with more informative ones
    partition_df$type[partition_df$type=="fp_exon"] <- "5' UTR Exon"
    partition_df$type[partition_df$type=="fp_intron"] <- "5' UTR Intron"
    partition_df$type[partition_df$type=="tp_exon"] <- "3' UTR Exon"
    partition_df$type[partition_df$type=="tp_intron"] <- "3' UTR Intron"
    partition_df$type[partition_df$type=="promoter"] <- "Promoter"
    partition_df$type[partition_df$type=="exon"] <- "mRNA Exon"
    partition_df$type[partition_df$type=="intron"] <- "mRNA Intron"
    partition_df$type[partition_df$type=="repeat"] <- "Repeat"
    partition_df$type[partition_df$type=="intergenic"] <- "Intergenic"        

    partition_order <- c("5' UTR Exon", "5' UTR Intron", "3' UTR Exon", "3' UTR Intron",
                         "Promoter", "mRNA Exon", "mRNA Intron", "Repeat", "Intergenic")
    
    partition_df$type <- factor(partition_df$type, levels=partition_order)

    partition_outf <- paste0(outdir, 'tables/', prefix, '_unstranded_genomic_partition_', r2_thresh, "_ld_", dist_thresh, "_dist.txt")
    write.table(partition_df, file=partition_outf, quote=F, sep="\t",
                row.names=F, col.names=T)
    
    ## Visualize genomic partition
    make_graphic(paste0(outdir, 'plots/', prefix, '_stacked_partition_props_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))                    
    print(ggplot(partition_df, aes(x="", y=proportion, fill=type)) +
          scale_fill_brewer(palette="Paired", name="Genomic Partition") +
          xlab("") + ylab("Genomic Partition Proportions") +
          theme_bw() + geom_bar(stat="identity", position="stack") +
          guides(fill=guide_legend(nrow=3)) +
          theme(legend.position="bottom", axis.text.x = element_blank(),
                legend.text=element_text(size=LEGEND_TEXT_SIZE),
                legend.title=element_text(size=TITLE_SIZE*0.75),                 
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) + 
          plot_title("Proportions of LD SNPs falling in genomic partitions", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    make_graphic(paste0(outdir, 'plots/', prefix, '_partition_props_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))    
    print(ggplot(partition_df, aes(x=type, y=proportion, fill=type)) +
          scale_fill_brewer(palette="Paired", name="Genomic Partition") +
          xlab("Class of genomic element") + ylab("Genomic Partition Proportions") +
          theme_bw() + geom_bar(stat="identity", position="dodge") +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Proportions of LD SNPs falling in genomic partitions", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    make_graphic(paste0(outdir, 'plots/', prefix, '_partition_nums_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(partition_df, aes(x=type, y=number, fill=type)) +
          scale_fill_brewer(palette="Paired", name="Genomic Partition") +
          xlab("Class of genomic element") + ylab("Number of SNPs") +
          theme_bw() + geom_bar(stat="identity", position="dodge") +
          theme(legend.position="none",
                axis.text.x = element_text(angle=45, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Counts of SNPs falling in genomic partitions", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## Analyze genomic partition by tag SNP (region)
    ## first read in all the entry data
    entry_f <- paste0(datadir, prefix, "_", r2_thresh, "_ld_cutoff_snps_within_", dist_thresh,
                      "_unstranded_entrywise_partition.txt") 
    entry_df <- read.table(entry_f, header=T, sep="\t", quote="", as.is=T)
    ## replace with more informative labels
    entry_df$partition[entry_df$partition=="fp_utr_exon"] <- "5' UTR Exon"
    entry_df$partition[entry_df$partition=="fp_utr_intron"] <- "5' UTR Intron"
    entry_df$partition[entry_df$partition=="tp_utr_exon"] <- "3' UTR Exon"
    entry_df$partition[entry_df$partition=="tp_utr_intron"] <- "3' UTR Intron"
    entry_df$partition[entry_df$partition=="promoter"] <- "Promoter"
    entry_df$partition[entry_df$partition=="exon"] <- "mRNA Exon"
    entry_df$partition[entry_df$partition=="intron"] <- "mRNA Intron"
    entry_df$partition[entry_df$partition=="repeat"] <- "Repeat"
    entry_df$partition[entry_df$partition=="intergenic"] <- "Intergenic"        
    
    partition_by_tagsnp <- ddply(entry_df, .(tag_name), function(x) {
        ## we only want unique SNPs in each region to count here:
        uniq_x <- x[!(duplicated(x$rsID)),]
        uniq_x$partition <- factor(uniq_x$partition, levels = partition_order)
        this_tab <- table(uniq_x$partition)
        data.frame(this_tab / sum(this_tab))
    })
    colnames(partition_by_tagsnp) <- c("tag_name", "Partition", "Proportion")

    ## add a 'total' column
    partition_by_tagsnp <- rbind(partition_by_tagsnp,
                                 data.frame(tag_name="Total", Partition=partition_df$type,
                                            Proportion=partition_df$proportion))
    ## add a column for the tag name (without rsID)
    partition_by_tagsnp$tag_no_rsid <- gsub(":rs.*", "", partition_by_tagsnp$tag_name)
    
    ## make sure the total column is last
    partition_by_tagsnp$tag_no_rsid <- factor(partition_by_tagsnp$tag_no_rsid,
                                           levels=c(sort(unique(gsub(":rs.*", "", entry_df$tag_name))), "Total"))
    partition_by_tagsnp$tag_name <- factor(partition_by_tagsnp$tag_name,
                                           levels=c(sort(unique(entry_df$tag_name)), "Total"))
    
    partition_by_tagsnp$Partition <- factor(partition_by_tagsnp$Partition, levels=partition_order)

    make_graphic(paste0(outdir, 'plots/', prefix, '_partition_prop_by_tagsnp_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(partition_by_tagsnp, aes_string(x=TAG_VAR, y="Proportion", fill="Partition")) +
          scale_fill_brewer(palette="Paired", name="Genomic Partition") +
          xlab(TAG_LAB) + ylab("Proportion of SNPs") +
          theme_bw() + geom_bar(stat="identity", position="stack") +
          guides(fill=guide_legend(nrow=3)) +
          theme(legend.position="bottom",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                legend.text=element_text(size=LEGEND_TEXT_SIZE),
                legend.title=element_text(size=TITLE_SIZE*0.75), 
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Proportions of SNPs falling in genomic partitions by region", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

}

## prefix <- param_ref[['outprefix']]
## datadir <- paste0(param_ref[['outdir']], "/unstranded_genomic_partition/")
## outdir <- paste0(result_outdir, "/unstranded_genomic_partition/")
