## calculate_gtex_tissue_PCs.R
## alex amlie-wolf
## 04/04/2018
## a script to calculate tissue-specific principal components

gtex_expr_dir <- "/home/alexaml/data/GTEx/GTEx_v6p/rnaseq/all_chr_gtex_v6p_expression/"
sample_info_file <- "/home/alexaml/data/GTEx/subject_sample_info/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
gtex_class_file <- "/home/alexaml/data/GTEx/gtex_classes.txt"
## put the results in the same directory, for ease in the lncRNA correlation analysis script
outdir <- "/home/alexaml/data/GTEx/GTEx_v6p/rnaseq/all_chr_gtex_v6p_expression/"

gtex_category_df <- read.table(gtex_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")
gtex_category_df$SMTSD <- gsub("_Analysis.snpgenes", "", gtex_category_df$GTEx.Data)
gtex_category_df$SMTSD <- gsub("_", " ", gtex_category_df$SMTSD)

sample_info_df <- read.table(sample_info_file, header=T, sep="\t", quote="", as.is=T)
sample_info_df$inferno_match <- gsub("- ", "", gsub("\\(|\\)", "", sample_info_df$SMTSD))
sample_info_df$inferno_class <- gtex_category_df$Class[match(sample_info_df$inferno_match, gtex_category_df$SMTSD)]

## now read in the full expression matrix and generate the tissue class-specific PCs
{
sample_pc_start <- proc.time()
full_expr_mat <- NULL
for(chr_file in list.files(gtex_expr_dir, pattern="*gtex_expression.txt", full.names=T)) {
    cat("Reading in full expression levels from", basename(chr_file), "\n")
    
    this_chr_expression <- read.table(chr_file, header=T, sep="\t", as.is=T, stringsAsFactors = F)

    ## add this to the full matrix
    full_expr_mat <- rbind(full_expr_mat,
                           as.matrix(this_chr_expression[,3:ncol(this_chr_expression)]))
}

## get the sample info matrix
sample_mat <- t(full_expr_mat) %*% full_expr_mat
## name the samples, using the most recent file
rownames(sample_mat) <- colnames(this_chr_expression)[3:ncol(this_chr_expression)]

## get the principal components for each tissue class
tissue_class_pcs <- list()
for(tissue_class in sort(unique(sample_info_df$inferno_class[!is.na(sample_info_df$inferno_class)]))) {
    ## pull out the samples in this tissue class
    class_samples <- sample_info_df$SAMPID[!is.na(sample_info_df$inferno_class) &
                                           sample_info_df$inferno_class==tissue_class]
    ## get the square matrix for these samples
    sample_match <- rownames(sample_mat) %in% gsub("-", ".", class_samples, fixed=T)
    cat(sum(sample_match), "samples found for", tissue_class, "\n")

    tissue_mat <- sample_mat[sample_match,sample_match]

    ## now calculate the top 10 PCs in this tissue
    this_pca <- prcomp(tissue_mat, scale. = TRUE)

    write.table(data.frame(sample=rownames(tissue_mat),
                           this_pca$rotation[,1:10]),
                paste0(outdir, '/', gsub(" ", "_", tissue_class), '_top10_PCs.txt'),
                quote=F, sep="\t", row.names=F, col.names=T)
}

cat("Computing tissue-specific PCs took", (proc.time() - sample_pc_start)[['elapsed']], 'seconds\n')
}

