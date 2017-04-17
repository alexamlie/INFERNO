## analyze_tfbs_fantom5_eqtl_overlap_function.R
## alex amlie-wolf 02/02/2016

analyze_tfbs_fantom5_eqtl_overlap <- function(prefix, datadir, outdir, out_subtitle, r2_thresh, dist_thresh, enh_overlap_type) {
    dir.create(paste0(outdir, "/tables/"), F, T)
    
    ## read in the factorbook data
    factorbook_file <- paste0(datadir, '/factorbook_overlap/', prefix, "_", r2_thresh,
                              "_ld_cutoff_snps_within_", dist_thresh, "_tfbs_overlaps.txt")
    factorbook_df <- read.table(factorbook_file, header=T, sep="\t", quote="", as.is=T)

    unique_snp_overlap_df <- factorbook_df[!(duplicated(factorbook_df[,c("rsID", "tag_name", "tfbs_chr", "tfbs_start", "tfbs_end", "tf", "score", "cells")])),]

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
    
    ## read in the eQTL data
    eqtl_overlap_file <- paste0(datadir, '/gtex_eqtl_overlap/',
                                prefix, "_", r2_thresh, "_ld_cutoff_snps_within_",
                                dist_thresh, "_eqtl_overlaps.txt")
    eqtl_overlap_df <- read.table(eqtl_overlap_file, header=T, sep="\t", quote="", as.is=T)

    enh_snp_cols <- colnames(enh_overlap_target_df)[!(colnames(enh_overlap_target_df) %in% c("enh_source", "enh_chr", "enh_start", "enh_end", "num_assocs", "refseqId","symbol", "r", "fdr"))]
    
    enh_target_df <- ddply(enh_overlap_target_df, enh_snp_cols,
                           summarize, num_enh_targets = length(unique(symbol[!is.na(symbol)])),
                           enh_targets = paste(sort(unique(symbol[!is.na(symbol)])), collapse=","),
                           num_enh_tissues = length(unique(enh_source)), 
                           enh_tissues = paste(sort(unique(enh_source)), collapse=","))

    ## then get summary information for each eQTL
    eqtl_snp_cols <- colnames(eqtl_overlap_df)[!(colnames(eqtl_overlap_df) %in% c("tissue", "gene", "beta", "t_stat", "se", "p_value", "nom_thresh", "min_p", "gene_emp_p", "maf", "gene_name", "gene_source", "gene_type"))]
    eqtl_target_df <- ddply(eqtl_overlap_df, eqtl_snp_cols,
                            summarize, num_eqtl_genes = length(unique(gene_name)),
                            eqtl_genes = paste(sort(unique(gene_name)), collapse=","),
                            num_eqtl_tissues = length(unique(tissue)), 
                            eqtl_tissues = paste(sort(unique(tissue)), collapse=","))

    ## now merge these:
    eqtl_enh_snp_df <- join(enh_target_df, eqtl_target_df, type='inner')

    ## then, merge this with the factorbook information
    ## for each SNP, we want information on all the TFs and cell types it is from
    tfbs_snp_cols <- colnames(factorbook_df)[!(colnames(factorbook_df) %in% c("tfbs_chr", "tfbs_start", "tfbs_end", "tf", "score", "cells"))]
    tfbs_simplified_df <- ddply(factorbook_df, tfbs_snp_cols, summarize,
                                num_tfs = length(unique(tf)),
                                tfs = paste(sort(unique(tf)), collapse=","),
                                num_tf_tiss = length(unique(cells)),
                                tfbs_tissues = paste(sort(unique(cells)), collapse=","))

    ## merge all 3:
    tfbs_enh_eqtl_snp_df <- join(eqtl_enh_snp_df, tfbs_simplified_df, type="inner")
    
    merged_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                          "_ld_cutoff_snps_within_", dist_thresh,
                          "_", enh_overlap_type,
                          "_fantom5_eqtl_tfbs_all_hits.txt")
    write.table(tfbs_enh_eqtl_snp_df, merged_outf, quote=F, sep="\t", col.names=T, row.names=F)

    unique_snp_overlap_df <- ddply(tfbs_enh_eqtl_snp_df, .(chr, rsID, pos, tag_name), summarize,
                                   num_enh_targets = unique(num_enh_targets),
                                   enh_targets=unique(enh_targets),
                                   num_enh_tissues = unique(num_enh_tissues), 
                                   enh_tissues = unique(enh_tissues),                           
                                   num_eqtl_genes = unique(num_eqtl_genes), eqtl_genes = unique(eqtl_genes),
                                   num_eqtl_tissues = unique(num_eqtl_tissues),
                                   eqtl_tissues = unique(eqtl_tissues),
                                   num_tfs = unique(num_tfs), tfs = unique(tfs),
                                   num_tf_tiss = unique(num_tf_tiss), tfbs_tissues = unique(tfbs_tissues))
    
    unique_outf <- paste0(outdir, 'tables/', prefix, "_", r2_thresh,
                          "_ld_cutoff_snps_within_", dist_thresh,
                          "_", enh_overlap_type,
                          "_fantom5_eqtl_tfbs_unique_hits.txt")
    write.table(unique_snp_overlap_df, unique_outf, quote=F, sep="\t", col.names=T, row.names=F)

}
