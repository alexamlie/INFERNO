## analyze_closest_genes_function.R
## alex amlie-wolf 01/26/2016
## analyzing the closest gene results for SNPs
## part of the automatic analysis

## prefix <- param_ref[['outprefix']]
## datadir <- param_ref[['outdir']]
## outdir <- paste0(result_outdir, '/closest_gene/')

analyze_closest_genes <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh) {
    dir.create(paste0(outdir, 'plots/'), F, T)

    closest_genes_file <- paste0(datadir, '/closest_gene/', prefix, "_", r2_thresh,
                                 "_ld_cutoff_snps_within_", dist_thresh, "_closest_genes.txt")
    closest_genes_df <- read.table(closest_genes_file, header=T, sep="\t",
                                   quote="", as.is=T)

    if(nrow(closest_genes_df)==0) {
        cat("No closest genes found in this dataset!\n")
        return("No closest genes found in this dataset!")
    } 
    
    ## count number of overlaps:
    ## count by ID
    num_gene_ids <- unlist(lapply(strsplit(closest_genes_df$closest_gene_ids, ","),
                                  length))
    ## count number of unique gene symbols
    num_uniq_symbols <- unlist(lapply(strsplit(closest_genes_df$closest_gene_symbols, ","),
                                      function(x) {length(unique(x))}))
    ## also get the unique symbols
    uniq_symbols <- unlist(lapply(lapply(strsplit(closest_genes_df$closest_gene_symbols, ","),
                                         unique), paste0, collapse=","))

    ## get the distances too
    uniq_dist <- as.numeric(unlist(lapply(strsplit(as.character(closest_genes_df$distance), ","),
                                          '[', 1)))

    closest_genes_df <- cbind(closest_genes_df, num_gene_ids = num_gene_ids,
                              num_uniq_symbols=num_uniq_symbols, unique_dist=uniq_dist,
                              unique_symbols=uniq_symbols)

    ## Analyze numbers of genes
    ## first check how many unknown gene symbols there were:
    cat("Frequency of unknown gene symbols: ", length(grep("NotFound", closest_genes_df$unique_symbols)) / nrow(closest_genes_df), "\n")
    ## compare to how many times that was the only symbol
    cat("Frequency of only unknown gene symbols: ", sum(closest_genes_df$unique_symbols=="NotFound") / nrow(closest_genes_df), "\n")

    ## checking just the number of unique refseq IDs by tag SNP
    uniq_tag_region_df <- ddply(closest_genes_df, .(tag_name), function(x) {
        ## only look at unique SNPs per region
        return(x[!(duplicated(x$rsID)),])
    })

    uniq_tag_region_df$tag_no_rsid <- gsub(":rs.*", "", uniq_tag_region_df$tag_name)
    
    make_graphic(paste0(outdir, 'plots/', prefix, '_gene_counts_',
                        r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(uniq_tag_region_df, aes_string(x=TAG_VAR, y="num_gene_ids", fill=TAG_VAR)) +
          scale_fill_hue(h=c(180, 270)) +
          xlab(TAG_LAB) + ylab("Number of closest genes (refseq IDs)") +
          theme_bw() + geom_boxplot() +
          theme(legend.position="none",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Number of closest genes by refseq ID", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    gene_count_props <- ddply(closest_genes_df, .(tag_name), function(x) {
        uniq_x <- x[!(duplicated(x$rsID)),]
        data.frame(table(uniq_x$num_uniq_symbols, dnn=c('num_symbols'))/length(uniq_x$num_uniq_symbols))
        })

    gene_count_props$tag_no_rsid <- gsub(":rs.*", "", gene_count_props$tag_name)
    
    ## checking how many different gene symbols each SNP has (including NotFound)
    make_graphic(paste0(outdir, 'plots/', prefix, '_gene_counts_by_symbols_',
                 r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(gene_count_props, aes_string(x=TAG_VAR, y="Freq") +
                 aes(fill=factor(num_symbols))) +
          scale_fill_brewer(palette='Set1', name="Number of gene symbols") +
          xlab(TAG_LAB) + ylab("Proportions of number of unique closest genes by symbol") +
          theme_bw() + geom_bar(stat='identity', position='dodge') +
          theme(legend.position="bottom",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5),
                legend.text = element_text(size=LEGEND_TEXT_SIZE)) +
          plot_title("Proportions of unique closest genes by gene symbol", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

    ## Analyze distance to genes
    max_dist <- max(uniq_tag_region_df$unique_dist[!is.na(uniq_tag_region_df$unique_dist) &
                                                   uniq_tag_region_df$unique_dist != Inf])
    make_graphic(paste0(outdir, 'plots/', prefix, '_gene_dists_',
                 r2_thresh, "_ld_", dist_thresh, "_dist"))
    print(ggplot(uniq_tag_region_df, aes_string(x=TAG_VAR, y="unique_dist", fill=TAG_VAR, colour="unique_dist")) +
          scale_fill_hue(h=c(180, 270)) +
          scale_y_continuous(labels=comma, breaks=seq(0, round_any(max(max_dist), 50000, f=ceiling), by=25000)) +
          xlab(TAG_LAB) + ylab("Distance to closest gene (bp)") +
          theme_bw() + geom_boxplot() +
          geom_point(position=position_jitter(width=0.4), alpha=0.3) +
          theme(legend.position="none",
                axis.text.x = element_text(angle=60, hjust=1, size=AXIS_TEXT_X_SIZE),
                axis.text.y = element_text(size=AXIS_TEXT_Y_SIZE),
                title=element_text(size=TITLE_SIZE), plot.title = element_text(hjust = 0.5)) +
          plot_title("Distributions of distances to closest genes", r2_thresh, dist_thresh, out_subtitle))
    dev.off()

}
