## webgestalt_pathway_analysis.R
## alex amlie-wolf 03/22/18
## a script to automate the running and analysis of webgestalt pathway analyses
## for the colocalization target genes as well as the lncRNA targets
## this script must be run with internet access!

library(ggplot2)
library(plyr)
library(scales)
library(WebGestaltR)
sessionInfo()

## -----------------------------------------------------------------------------
## 0. Table of Contents
## -----------------------------------------------------------------------------
## 0. Table of Contents
## 1. Function Definitions
## 2. Read in gene lists
## 3. Run colocalization target pathway analysis
## 4. Run and analyze tissue-specific colocalization target pathway analysis
## 5. Run lncRNA correlation target pathway analysis
## 6. Visualize tissue-specific lncRNA target pathway enrichments

## -----------------------------------------------------------------------------
## 1. Function Definitions
## -----------------------------------------------------------------------------
make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='pdf') {
## make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='png') {
    if(type=='pdf') {
        pdf(file=paste0(filename, ".pdf"), width=10*width_ratio, height=10*height_ratio, pointsize=12, onefile=FALSE)
    } else if(type=='png') {
        ## use type='cairo' for when X11 doesn't work
        png(filename=paste0(filename, ".png"), width=10*width_ratio, height=10*height_ratio, res=300, units='in', type='cairo')
    } else {
        cat('filetype not supported\n')
    }
}

## -----------------------------------------------------------------------------
## 2. Read in gene lists, set up parameters for pathway analysis
## -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==1) {
    ## only one parameter: the data directory
    datadir <- args[1]
    ## use defaults for the others
    min_pathway_num <- 5
    max_pathway_num <- 1000
    fdr_thresh <- 0.05
} else if(length(args)==4) {
    datadir <- args[1]
    min_pathway_num <- as.numeric(args[2])
    max_pathway_num <- as.numeric(args[3])
    fdr_thresh <- as.numeric(args[4])
} else {
    stop("Usage: Rscript webgestalt_pathway_analysis.R <data directory> <minimum number of genes in pathway> <maximum number of genes in pathway> <FDR threshold>; note that the last 3 arguments are optional and default to 5, 1000, and 0.05")
    ## these parameters are for manual analysis purposes
    ## datadir <- "~/Dropbox/wang_lab/inferno_landscape_analysis/individual_phenotypes/SCZ2_94_indep_regions/"
    ## datadir <- "~/Dropbox/wang_lab/psp_enhancer_snps/PSP_H1H2_tag_variant/analysis_results/updated_tiss_spec_lncRNA_corr/"
    datadir <- "~/Dropbox/wang_lab/psp_enhancer_snps/PSP_H1H2_tag_variant/analysis_results/updated_brain_blood_lncRNA_040518/"

    ## define parameters for the pathway analysis
    min_pathway_num <- 5
    max_pathway_num <- 1000
    fdr_thresh <- 0.05
}

dir.create(paste0(datadir, '/pathway_analysis/'), F, T)

coloc_gene_list <- read.table(list.files(paste0(datadir, '/gtex_gwas_colocalization_analysis/tables/'), pattern="*all_top_genes*", full.names=T), header=F, sep="\t", col.names=c("gene"))

## for tissue-specific analysis, also read in the prioritized colocalization data
top_coloc_data <- read.table(list.files(paste0(datadir, '/gtex_gwas_colocalization_analysis/tables/'), pattern="*gtex_coloc_top_signals*", full.names=T), header=T, sep="\t", stringsAsFactors = F)

summary_file <- paste0(datadir, "/pathway_analysis/webgestalt_pathway_summary.txt")

## -----------------------------------------------------------------------------
## 3. Run colocalization target pathway analysis
## -----------------------------------------------------------------------------
all_coloc_gene_enrichment <- data.frame(stringsAsFactors = FALSE)

coloc_GO_BP_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                       enrichDatabase="geneontology_Biological_Process_noRedundant",
                                      interestGeneFile=list.files(paste0(datadir, '/gtex_gwas_colocalization_analysis/tables/'), pattern="*all_top_genes*", full.names=T),
                                      interestGeneType="genesymbol",
                                      referenceSet="genome",
                                      sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                      is.output=FALSE, maxNum=max_pathway_num)

if(!is.null(coloc_GO_BP_enrichment)) {
    cat(nrow(coloc_GO_BP_enrichment), "enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO BP\n")
    cat(nrow(coloc_GO_BP_enrichment), "enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO BP\n", file=summary_file)

    all_coloc_gene_enrichment <- rbind(all_coloc_gene_enrichment,
                                       data.frame(pathway_type="GO_BP",
                                                  coloc_GO_BP_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))
} else {
    cat("0 enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO BP\n")
    cat("0 enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO BP\n", file=summary_file)
}

coloc_GO_CC_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                       enrichDatabase="geneontology_Cellular_Component_noRedundant",
                                      interestGeneFile=list.files(paste0(datadir, '/gtex_gwas_colocalization_analysis/tables/'), pattern="*all_top_genes*", full.names=T),
                                      interestGeneType="genesymbol",
                                      referenceSet="genome",
                                      sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                      is.output=FALSE, maxNum=max_pathway_num)

if(!is.null(coloc_GO_CC_enrichment)) {
    cat(nrow(coloc_GO_CC_enrichment), "enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO CC\n")
    cat(nrow(coloc_GO_CC_enrichment), "enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO CC\n", file=summary_file, append=T)

    all_coloc_gene_enrichment <- rbind(all_coloc_gene_enrichment,
                                       data.frame(pathway_type="GO_CC",
                                                  coloc_GO_CC_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))    
} else {
    cat("0 enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO CC\n")
    cat("0 enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO CC\n", file=summary_file, append=T)
}

coloc_GO_MF_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                       enrichDatabase="geneontology_Molecular_Function_noRedundant",
                                      interestGeneFile=list.files(paste0(datadir, '/gtex_gwas_colocalization_analysis/tables/'), pattern="*all_top_genes*", full.names=T),
                                      interestGeneType="genesymbol",
                                      referenceSet="genome",
                                      sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                      is.output=FALSE, maxNum=max_pathway_num)

if(!is.null(coloc_GO_MF_enrichment)) {
    cat(nrow(coloc_GO_MF_enrichment), "enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO MF\n")
    cat(nrow(coloc_GO_MF_enrichment), "enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO MF\n", file=summary_file, append=T)

    all_coloc_gene_enrichment <- rbind(all_coloc_gene_enrichment,
                                       data.frame(pathway_type="GO_MF",
                                                  coloc_GO_MF_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))    
} else {
    cat("0 enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO MF\n")
    cat("0 enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in GO MF\n", file=summary_file, append=T)
}

coloc_KEGG_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                       enrichDatabase="pathway_KEGG",
                                      interestGeneFile=list.files(paste0(datadir, '/gtex_gwas_colocalization_analysis/tables/'), pattern="*all_top_genes*", full.names=T),
                                      interestGeneType="genesymbol",
                                      referenceSet="genome",
                                      sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                      is.output=FALSE, maxNum=max_pathway_num)

if(!is.null(coloc_KEGG_enrichment)) {
    cat(nrow(coloc_KEGG_enrichment), "enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in KEGG\n")
    cat(nrow(coloc_KEGG_enrichment), "enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in KEGG\n", file=summary_file, append=T)

    all_coloc_gene_enrichment <- rbind(all_coloc_gene_enrichment,
                                       data.frame(pathway_type="KEGG",
                                                  coloc_KEGG_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))    
} else {
    cat("0 enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in KEGG\n")
    cat("0 enriched pathways found for", nrow(coloc_gene_list), "colocalization targets in KEGG\n", file=summary_file, append=T)
}

## write this out
if(nrow(all_coloc_gene_enrichment) > 0) {
    write.table(all_coloc_gene_enrichment, paste0(datadir, '/pathway_analysis/all_coloc_target_pathway_enrichments.txt'), quote=F, sep="\t", row.names=F)
}

## -----------------------------------------------------------------------------
## 4. Run and analyze tissue-specific colocalization target pathway analysis
## -----------------------------------------------------------------------------
cat("Running tissue-specific colocalization target pathway analysis\n")
{
all_class_start_time <- proc.time()
all_class_pathway_enrichment <- data.frame(stringsAsFactors = FALSE)
for(this_class in sort(unique(top_coloc_data$gtex_tissue_class))) {

    this_gene_list <- top_coloc_data[top_coloc_data$gtex_tissue_class==this_class,]
    
    target_GO_BP_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                           enrichDatabase="geneontology_Biological_Process_noRedundant",
                                          interestGene=this_gene_list$eqtl_gene_name,
                                          interestGeneType="genesymbol",
                                          referenceSet="genome",
                                          sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                          is.output=FALSE, maxNum=max_pathway_num)

    if(!is.null(target_GO_BP_enrichment) & !any((grepl("ERROR", target_GO_BP_enrichment)))) {
        cat(nrow(target_GO_BP_enrichment), "enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO BP\n")
        cat(nrow(target_GO_BP_enrichment), "enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO BP\n", file=summary_file, append=T)
        all_class_pathway_enrichment <- rbind(all_class_pathway_enrichment,
                                              data.frame(class=this_class, pathway_type="GO_BP",
                                                         target_GO_BP_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))
    } else {
        cat("0 enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO BP\n")
        cat("0 enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO BP\n", file=summary_file, append=T)
    }

    target_GO_CC_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                           enrichDatabase="geneontology_Cellular_Component_noRedundant",
                                          interestGene=this_gene_list$eqtl_gene_name,
                                          interestGeneType="genesymbol",
                                          referenceSet="genome",
                                          sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                          is.output=FALSE, maxNum=max_pathway_num)

    if(!is.null(target_GO_CC_enrichment) & !(any(grepl("ERROR", target_GO_CC_enrichment)))) {
        cat(nrow(target_GO_CC_enrichment), "enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO CC\n")
        cat(nrow(target_GO_CC_enrichment), "enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO CC\n", file=summary_file, append=T)
        all_class_pathway_enrichment <- rbind(all_class_pathway_enrichment,
                                              data.frame(class=this_class, pathway_type="GO_CC",
                                                         target_GO_CC_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))         
    } else {
        cat("0 enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO CC\n")
        cat("0 enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO CC\n", file=summary_file, append=T)
    }
    target_GO_MF_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                           enrichDatabase="geneontology_Molecular_Function_noRedundant",
                                          interestGene=this_gene_list$eqtl_gene_name,
                                          interestGeneType="genesymbol",
                                          referenceSet="genome",
                                          sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                          is.output=FALSE, maxNum=max_pathway_num)

    if(!is.null(target_GO_MF_enrichment) & !(any(grepl("ERROR", target_GO_MF_enrichment)))) {
        cat(nrow(target_GO_MF_enrichment), "enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO MF\n")
        cat(nrow(target_GO_MF_enrichment), "enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO MF\n", file=summary_file, append=T)
        all_class_pathway_enrichment <- rbind(all_class_pathway_enrichment,
                                              data.frame(class=this_class, pathway_type="GO_MF",
                                                         target_GO_MF_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))                 
    } else {
        cat("0 enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO MF\n")
        cat("0 enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in GO MF\n", file=summary_file, append=T)
    }

    target_KEGG_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                           enrichDatabase="pathway_KEGG",
                                          interestGene=this_gene_list$eqtl_gene_name,
                                          interestGeneType="genesymbol",
                                          referenceSet="genome",
                                          sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                          is.output=FALSE, maxNum=max_pathway_num)

    if(!is.null(target_KEGG_enrichment) & !any((grepl("ERROR", target_KEGG_enrichment)))) {
        cat(nrow(target_KEGG_enrichment), "enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in KEGG\n")
        cat(nrow(target_KEGG_enrichment), "enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in KEGG\n", file=summary_file, append=T)
        all_class_pathway_enrichment <- rbind(all_class_pathway_enrichment,
                                              data.frame(class=this_class, pathway_type="KEGG",
                                                         target_KEGG_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))
    } else {
        cat("0 enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in KEGG\n")
        cat("0 enriched pathways found for", nrow(this_gene_list), "COLOC target genes in", this_class, "in KEGG\n", file=summary_file, append=T)
    }
    
}
cat("Tissue-specific colocalization target pathway analysis took", (proc.time() - all_class_start_time)[["elapsed"]], 'seconds\n')
}

## add a nice factor for visualization
## NOTE: this isn't currently used for any visualization
all_class_pathway_enrichment$clean_desc <- gsub(" - Homo sapiens (human)", "", all_class_pathway_enrichment$description, fixed=T)
all_class_pathway_enrichment$clean_desc <- factor(all_class_pathway_enrichment$clean_desc, ordered=T, levels=sort(unique(all_class_pathway_enrichment$clean_desc), decreasing=T))

## also add a level to show all the different pathway types
all_class_pathway_enrichment$desc_with_pathway <- paste0(all_class_pathway_enrichment$pathway_type, ": ", all_class_pathway_enrichment$clean_desc)
all_class_pathway_enrichment$desc_with_pathway <- factor(all_class_pathway_enrichment$desc_with_pathway, ordered=T, levels=sort(unique(all_class_pathway_enrichment$desc_with_pathway), decreasing=T))

if(nrow(all_class_pathway_enrichment) > 0) {
    write.table(all_class_pathway_enrichment, paste0(datadir, '/pathway_analysis/coloc_target_class_specific_pathway_enrichments.txt'), quote=F, sep="\t", row.names=F)
}
    
## -----------------------------------------------------------------------------
## 5. Run lncRNA correlation target pathway analysis
## -----------------------------------------------------------------------------
## start with the full target analysis
num_lncRNA_targets <- nrow(read.table(paste0(datadir, '/gtex_lncRNA_correlation_analysis/tables/all_lncRNA_genes_0.5_correlation_threshold.txt'), sep="\t", col.names=c("gene")))

cross_tissue_lncRNA_target_enrichment <- data.frame(stringsAsFactors = FALSE)

all_lncRNA_target_GO_BP_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                       enrichDatabase="geneontology_Biological_Process_noRedundant",
                                      interestGeneFile=paste0(datadir, '/gtex_lncRNA_correlation_analysis/tables/all_lncRNA_genes_0.5_correlation_threshold.txt'),
                                      interestGeneType="genesymbol",
                                      referenceSet="genome",
                                      sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                      is.output=FALSE, maxNum=max_pathway_num)

if(!is.null(all_lncRNA_target_GO_BP_enrichment)) {
    cat(nrow(all_lncRNA_target_GO_BP_enrichment), "enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO BP\n")
    cat(nrow(all_lncRNA_target_GO_BP_enrichment), "enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO BP\n", file=summary_file, append=T)

    cross_tissue_lncRNA_target_enrichment <- rbind(cross_tissue_lncRNA_target_enrichment,
                                                   data.frame(pathway_type="GO_BP",
                                                              all_lncRNA_target_GO_BP_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))    
} else {
    cat("0 enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO BP\n")
    cat("0 enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO BP\n", file=summary_file, append=T)
}

all_lncRNA_target_GO_CC_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                       enrichDatabase="geneontology_Cellular_Component_noRedundant",
                                      interestGeneFile=paste0(datadir, '/gtex_lncRNA_correlation_analysis/tables/all_lncRNA_genes_0.5_correlation_threshold.txt'),
                                      interestGeneType="genesymbol",
                                      referenceSet="genome",
                                      sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                      is.output=FALSE, maxNum=max_pathway_num)

if(!is.null(all_lncRNA_target_GO_CC_enrichment)) {
    cat(nrow(all_lncRNA_target_GO_CC_enrichment), "enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO CC\n")
    cat(nrow(all_lncRNA_target_GO_CC_enrichment), "enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO CC\n", file=summary_file, append=T)

    cross_tissue_lncRNA_target_enrichment <- rbind(cross_tissue_lncRNA_target_enrichment,
                                                   data.frame(pathway_type="GO_CC",
                                                              all_lncRNA_target_GO_CC_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))    
} else {
    cat("0 enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO CC\n")
    cat("0 enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO CC\n", file=summary_file, append=T)
}

all_lncRNA_target_GO_MF_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                       enrichDatabase="geneontology_Molecular_Function_noRedundant",
                                      interestGeneFile=paste0(datadir, '/gtex_lncRNA_correlation_analysis/tables/all_lncRNA_genes_0.5_correlation_threshold.txt'),
                                      interestGeneType="genesymbol",
                                      referenceSet="genome",
                                      sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                      is.output=FALSE, maxNum=max_pathway_num)

if(!is.null(all_lncRNA_target_GO_MF_enrichment)) {
    cat(nrow(all_lncRNA_target_GO_MF_enrichment), "enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO MF\n")
    cat(nrow(all_lncRNA_target_GO_MF_enrichment), "enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO MF\n", file=summary_file, append=T)

    cross_tissue_lncRNA_target_enrichment <- rbind(cross_tissue_lncRNA_target_enrichment,
                                                   data.frame(pathway_type="GO_MF",
                                                              all_lncRNA_target_GO_MF_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))    
} else {
    cat("0 enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO MF\n")
    cat("0 enriched pathways found for", num_lncRNA_targets, "lncRNA targets in GO MF\n", file=summary_file, append=T)
}

all_lncRNA_target_KEGG_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                       enrichDatabase="pathway_KEGG",
                                      interestGeneFile=paste0(datadir, '/gtex_lncRNA_correlation_analysis/tables/all_lncRNA_genes_0.5_correlation_threshold.txt'),
                                      interestGeneType="genesymbol",
                                      referenceSet="genome",
                                      sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                      is.output=FALSE, maxNum=max_pathway_num)

if(!is.null(all_lncRNA_target_KEGG_enrichment)) {
    cat(nrow(all_lncRNA_target_KEGG_enrichment), "enriched pathways found for", num_lncRNA_targets, "lncRNA targets in KEGG\n")
    cat(nrow(all_lncRNA_target_KEGG_enrichment), "enriched pathways found for", num_lncRNA_targets, "lncRNA targets in KEGG\n", file=summary_file, append=T)

    cross_tissue_lncRNA_target_enrichment <- rbind(cross_tissue_lncRNA_target_enrichment,
                                                   data.frame(pathway_type="KEGG",
                                                              all_lncRNA_target_KEGG_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))    
} else {
    cat("0 enriched pathways found for", num_lncRNA_targets, "lncRNA targets in KEGG\n")
    cat("0 enriched pathways found for", num_lncRNA_targets, "lncRNA targets in KEGG\n", file=summary_file, append=T)
}    

## write out the summary table
if(nrow(cross_tissue_lncRNA_target_enrichment) > 0) {
    write.table(cross_tissue_lncRNA_target_enrichment,
                paste0(datadir, '/pathway_analysis/lncRNA_cross_tissue_target_pathway_enrichments.txt'), quote=F, sep="\t", row.names=F)
}

## now do tissue-specific correlation analysis
{
lncrna_all_class_start_time <- proc.time()
lncrna_tissue_specific_target_enrichment <- data.frame(stringsAsFactors = FALSE)
for(class_file in sort(list.files(paste0(datadir, '/gtex_lncRNA_correlation_analysis/tables/class_specific_gene_lists/'), pattern="*correlation_threshold.txt"))) {
## for PSP analysis:
## for(class_file in list.files(paste0(datadir, '/tables/class_specific_gene_lists/'), pattern="*correlation_threshold.txt")) {
    
    this_class <- gsub(" ", "_", gsub("_lncRNA_genes.*", "", class_file))

    this_gene_list <- read.table(paste0(datadir, '/gtex_lncRNA_correlation_analysis/tables/class_specific_gene_lists/', class_file), header=F, col.names="gene", as.is=T)
    ## for PSP analysis:    
    ## this_gene_list <- read.table(paste0(datadir, '/tables/class_specific_gene_lists/', class_file), header=F, col.names="gene", as.is=T)
    
    lncRNA_target_GO_BP_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                           enrichDatabase="geneontology_Biological_Process_noRedundant",
                                          interestGene=this_gene_list$gene,
                                          interestGeneType="genesymbol",
                                          referenceSet="genome",
                                          sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                          is.output=FALSE, maxNum=max_pathway_num)

    if(!is.null(lncRNA_target_GO_BP_enrichment) & !any((grepl("ERROR", lncRNA_target_GO_BP_enrichment)))) {
        cat(nrow(lncRNA_target_GO_BP_enrichment), "enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO BP\n")
        cat(nrow(lncRNA_target_GO_BP_enrichment), "enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO BP\n", file=summary_file, append=T)
        
        lncrna_tissue_specific_target_enrichment <- rbind(lncrna_tissue_specific_target_enrichment,
                                              data.frame(class=this_class, pathway_type="GO_BP",
                                                         lncRNA_target_GO_BP_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))
    } else {
        cat("0 enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO BP\n")
        cat("0 enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO BP\n", file=summary_file, append=T)
    }

    lncRNA_target_GO_CC_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                           enrichDatabase="geneontology_Cellular_Component_noRedundant",
                                          interestGene=this_gene_list$gene,
                                          interestGeneType="genesymbol",
                                          referenceSet="genome",
                                          sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                          is.output=FALSE, maxNum=max_pathway_num)

    if(!is.null(lncRNA_target_GO_CC_enrichment) & !(any(grepl("ERROR", lncRNA_target_GO_CC_enrichment)))) {
        cat(nrow(lncRNA_target_GO_CC_enrichment), "enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO CC\n")
        cat(nrow(lncRNA_target_GO_CC_enrichment), "enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO CC\n", file=summary_file, append=T)
        
        lncrna_tissue_specific_target_enrichment <- rbind(lncrna_tissue_specific_target_enrichment,
                                              data.frame(class=this_class, pathway_type="GO_CC",
                                                         lncRNA_target_GO_CC_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))         
    } else {
        cat("0 enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO CC\n")
        cat("0 enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO CC\n", file=summary_file, append=T)
    }
    lncRNA_target_GO_MF_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                           enrichDatabase="geneontology_Molecular_Function_noRedundant",
                                          interestGene=this_gene_list$gene,
                                          interestGeneType="genesymbol",
                                          referenceSet="genome",
                                          sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                          is.output=FALSE, maxNum=max_pathway_num)

    if(!is.null(lncRNA_target_GO_MF_enrichment) & !(any(grepl("ERROR", lncRNA_target_GO_MF_enrichment)))) {
        cat(nrow(lncRNA_target_GO_MF_enrichment), "enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO MF\n")
        cat(nrow(lncRNA_target_GO_MF_enrichment), "enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO MF\n", file=summary_file, append=T)
        
        lncrna_tissue_specific_target_enrichment <- rbind(lncrna_tissue_specific_target_enrichment,
                                              data.frame(class=this_class, pathway_type="GO_MF",
                                                         lncRNA_target_GO_MF_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))                 
    } else {
        cat("0 enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO MF\n")
        cat("0 enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in GO MF\n", file=summary_file, append=T)
    }

    lncRNA_target_KEGG_enrichment <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                           enrichDatabase="pathway_KEGG",
                                          interestGene=this_gene_list$gene,
                                          interestGeneType="genesymbol",
                                          referenceSet="genome",
                                          sigMethod="fdr", fdrThr=fdr_thresh, minNum=min_pathway_num,
                                          is.output=FALSE, maxNum=max_pathway_num)

    if(!is.null(lncRNA_target_KEGG_enrichment) & !any((grepl("ERROR", lncRNA_target_KEGG_enrichment)))) {
        cat(nrow(lncRNA_target_KEGG_enrichment), "enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in KEGG\n")
        cat(nrow(lncRNA_target_KEGG_enrichment), "enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in KEGG\n", file=summary_file, append=T)
        
        lncrna_tissue_specific_target_enrichment <- rbind(lncrna_tissue_specific_target_enrichment,
                                              data.frame(class=this_class, pathway_type="KEGG",
                                                         lncRNA_target_KEGG_enrichment[,c("geneset", "description", "C", "O", "E", "R", "PValue", "FDR", "OverlapGene_UserID")]))
    } else {
        cat("0 enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in KEGG\n")
        cat("0 enriched pathways found for", nrow(this_gene_list), "lncRNA target genes in", this_class, "in KEGG\n", file=summary_file, append=T)
    }
    
}
cat("Tissue-specific lncRNA target pathway analysis took", (proc.time() - lncrna_all_class_start_time)[["elapsed"]], 'seconds\n')
}

## add a nice factor for visualization
lncrna_tissue_specific_target_enrichment$clean_desc <- gsub(" - Homo sapiens (human)", "", lncrna_tissue_specific_target_enrichment$description, fixed=T)
lncrna_tissue_specific_target_enrichment$clean_desc <- factor(lncrna_tissue_specific_target_enrichment$clean_desc, ordered=T, levels=sort(unique(lncrna_tissue_specific_target_enrichment$clean_desc), decreasing=T))

## also add a level to show all the different pathway types
lncrna_tissue_specific_target_enrichment$desc_with_pathway <- paste0(lncrna_tissue_specific_target_enrichment$pathway_type, ": ", lncrna_tissue_specific_target_enrichment$clean_desc)
lncrna_tissue_specific_target_enrichment$desc_with_pathway <- factor(lncrna_tissue_specific_target_enrichment$desc_with_pathway, ordered=T, levels=sort(unique(lncrna_tissue_specific_target_enrichment$desc_with_pathway), decreasing=T))

if(nrow(lncrna_tissue_specific_target_enrichment) > 0) {
    write.table(lncrna_tissue_specific_target_enrichment, paste0(datadir, '/pathway_analysis/lncRNA_tissue_specific_pathway_enrichments.txt'), quote=F, sep="\t", row.names=F)
}
## write.table(lncrna_tissue_specific_target_enrichment, paste0(datadir, '/pathway_analysis/lncrna_blood_brain_class_pathway_enrichments.txt'), quote=F, sep="\t", row.names=F)

## -----------------------------------------------------------------------------
## 6. Visualize tissue-specific lncRNA target pathway enrichments
## -----------------------------------------------------------------------------
lncrna_all_pathway_combos <- lncrna_tissue_specific_target_enrichment[,c('class', 'desc_with_pathway', 'FDR')]
lncrna_all_pathway_combos <- merge(expand.grid(class=unique(lncrna_all_pathway_combos$class),
                                        desc_with_pathway=unique(lncrna_all_pathway_combos$desc_with_pathway)),
                            lncrna_all_pathway_combos, all.x=T)
lncrna_all_pathway_combos$pathway_type <- unlist(lapply(strsplit(as.character(lncrna_all_pathway_combos$desc_with_pathway), ":"), "[[", 1))
lncrna_all_pathway_combos$pathway_desc <- unlist(lapply(strsplit(as.character(lncrna_all_pathway_combos$desc_with_pathway), ":"), "[[", 2))

make_graphic(paste0(datadir, '/pathway_analysis/lncrna_tissue_specific_pathways_heatmap'), width_ratio=4, height_ratio=2)
print(ggplot(lncrna_all_pathway_combos, aes(y=factor(class), x=factor(desc_with_pathway))) +
      geom_tile(aes(fill=FDR), colour="black", size=2) +
      scale_fill_gradient(low="navyblue", high="orchid", guide="colorbar", na.value="white", limits=c(0, fdr_thresh)) + 
      theme_bw() + ylab("lncRNA tissue class") + xlab("KEGG or GO pathway") +
      theme(axis.text.x = element_text(angle=45, hjust=1, size=25),
            axis.text.y = element_text(size=18),
            legend.text = element_text(size=25), legend.key.height=unit(0.10, "npc"),
            title=element_text(size=30), plot.title = element_text(hjust = 0.5)))
dev.off()

make_graphic(paste0(datadir, '/pathway_analysis/lncrna_tissue_specific_pathways_heatmap_facet'), width_ratio=5, height_ratio=3)
print(ggplot(lncrna_all_pathway_combos, aes(y=factor(class), x=factor(pathway_desc))) +
      geom_tile(aes(fill=FDR), colour="black", size=2) +
      facet_wrap(~ pathway_type, scales="free", shrink=F) + 
      scale_fill_gradient(low="navyblue", high="orchid", guide="colorbar", na.value="white",limits=c(0, fdr_thresh)) + 
      theme_bw() + ylab("lncRNA tissue class") + xlab("KEGG or GO pathway") +
      theme(axis.text.x = element_text(angle=75, hjust=1, size=20),
            axis.text.y = element_text(size=18),
            strip.text = element_text(size=30), 
            legend.text = element_text(size=25), legend.key.height=unit(0.10, "npc"),
            title=element_text(size=30), plot.title = element_text(hjust = 0.5)))
dev.off()

## make specific plots for the different pathway annotations
for(pathway_type in c("KEGG", "GO_BP", "GO_CC", "GO_MF")) {
    these_pathways <- lncrna_all_pathway_combos[grepl(pathway_type, lncrna_all_pathway_combos$desc_with_pathway),]
    these_pathways$desc_with_pathway <- gsub(paste0(pathway_type, ": "), "", these_pathways$desc_with_pathway)
    these_pathways$desc_with_pathway <- factor(these_pathways$desc_with_pathway, ordered=T, levels=sort(unique(these_pathways$desc_with_pathway), decreasing=T))

    make_graphic(paste0(datadir, '/pathway_analysis/lncrna_tissue_specific_', pathway_type, '_pathways_heatmap'), width_ratio=2, height_ratio=1.75 * round_any(length(unique(these_pathways$pathway_desc)), accuracy=100, f=ceiling) / 100)
    print(ggplot(these_pathways, aes(x=factor(class), y=desc_with_pathway)) +
          geom_tile(aes(fill=FDR), colour="black", size=2) +
          scale_fill_gradient(low="navyblue", high="orchid", guide="colorbar", na.value="white", limits=c(0, fdr_thresh)) + 
          theme_bw() + xlab("lncRNA tissue class") +
          ylab(paste0(pathway_type, " pathway")) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=25),
                axis.text.y = element_text(size=18),
                legend.text = element_text(size=25), legend.key.height=unit(0.10, "npc"),
                title=element_text(size=30), plot.title = element_text(hjust = 0.5)))
    dev.off()

    ## also make one with only tissue categories that showed up
    make_graphic(paste0(datadir, '/pathway_analysis/lncrna_tissue_specific_', pathway_type, '_pathways_filtered_category_heatmap'), width_ratio=2, height_ratio=1.75 * round_any(length(unique(these_pathways$pathway_desc)), accuracy=100, f=ceiling) / 100)
    print(ggplot(these_pathways, aes(x=class, y=desc_with_pathway)) +
          geom_tile(aes(fill=FDR), colour="black", size=2) +
          scale_fill_gradient(low="navyblue", high="orchid", guide="colorbar", na.value="white", limits=c(0, fdr_thresh)) + 
          theme_bw() + xlab("lncRNA tissue class") +
          ylab(paste0(pathway_type, " pathway")) +
          theme(axis.text.x = element_text(angle=45, hjust=1, size=25),
                axis.text.y = element_text(size=18),
                legend.text = element_text(size=25), legend.key.height=unit(0.10, "npc"),
                title=element_text(size=30), plot.title = element_text(hjust = 0.5)))
    dev.off()

}
