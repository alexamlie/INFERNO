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

## also store the variance explained by each principal component
variance_explained <- data.frame(stringsAsFactors = FALSE)
for(tissue_class in sort(unique(sample_info_df$inferno_class[!is.na(sample_info_df$inferno_class)]))) {
    ## pull out the samples in this tissue class
    class_samples <- sample_info_df$SAMPID[!is.na(sample_info_df$inferno_class) &
                                           sample_info_df$inferno_class==tissue_class]
    ## get the square matrix for these samples
    sample_match <- rownames(sample_mat) %in% gsub("-", ".", class_samples, fixed=T)
    cat(sum(sample_match), "samples found for", tissue_class, "\n")

    tissue_mat <- sample_mat[sample_match,sample_match]

    ## calculate the top 10 PCs in this tissue
    this_pca <- prcomp(tissue_mat, scale. = TRUE)

    write.table(data.frame(sample=rownames(tissue_mat),
                           this_pca$x[,1:10]),
                paste0(outdir, '/', gsub(" ", "_", tissue_class), '_top10_PCs.txt'),
                quote=F, sep="\t", row.names=F, col.names=T)

    this_summary <- summary(this_pca)
    
    ## and save the variance explained
    variance_explained <- rbind(variance_explained,
                                data.frame(category=tissue_class,
                                           PC=paste0("PC", 1:10),
                                           variance_prop=this_summary$importance[2,1:10],
                                           cumulative_prop=this_summary$importance[3,1:10]))
}

cat("Computing tissue-specific PCs took", (proc.time() - sample_pc_start)[['elapsed']], 'seconds\n')
}

## now plot the cumulative variance explained per tissue category:
variance_explained$PC <- factor(variance_explained$PC, levels=paste0("PC", 1:10), ordered=T)

library(ggplot2)
dir.create("/home/alexaml/data/INFERNO/gtex_v6p_PC_adjustment", F, T)

fantom5_category_df <- read.table("/home/alexaml/data/FANTOM5/Enhancers/fantom5_classes.txt", header=T, sep="\t", quote="", as.is=T, comment.char="")
gtex_category_df <- read.table("/home/alexaml/data/GTEx/gtex_classes.txt", header=T, sep="\t", quote="", as.is=T, comment.char="")
roadmap_category_df <- read.table("/home/alexaml/data/roadmap/roadmap_classes.txt", header=T, sep="\t", quote="", as.is=T, comment.char="")
## find the full list of classes
all_classes <- sort(union(fantom5_category_df$Class, union(gtex_category_df$Class, roadmap_category_df$Class)))

## now define the color palette
## generated from http://tools.medialab.sciences-po.fr/iwanthue/
## parameters: H 0-360, C 0.46 - 3, L 0.5-1.5
category_colors <- c(rgb(227,191,35, maxColorValue=255), rgb(224,98,247, maxColorValue=255), rgb(36,153,165, maxColorValue=255), rgb(232,2,58, maxColorValue=255), rgb(30,181,54, maxColorValue=255), rgb(253,164,185, maxColorValue=255), rgb(44,115,232, maxColorValue=255), rgb(83,118,43, maxColorValue=255), rgb(203,50,133, maxColorValue=255), rgb(194,196,252, maxColorValue=255), rgb(157,236,192, maxColorValue=255), rgb(141,86,174, maxColorValue=255), rgb(253,133,114, maxColorValue=255), rgb(175,242,253, maxColorValue=255), rgb(149,87,122, maxColorValue=255), rgb(131,233,24, maxColorValue=255), rgb(182,241,140, maxColorValue=255), rgb(252,50,7, maxColorValue=255), rgb(244,148,221, maxColorValue=255), rgb(40,169,123, maxColorValue=255), rgb(247,183,144, maxColorValue=255), rgb(242,184,94, maxColorValue=255), rgb(53,108,145, maxColorValue=255), rgb(198,8,78, maxColorValue=255), rgb(61,115,80, maxColorValue=255), rgb(41,232,215, maxColorValue=255), rgb(122,107,24, maxColorValue=255), rgb(79,153,241, maxColorValue=255), rgb(100,130,128, maxColorValue=255), rgb(166,169,37, maxColorValue=255), rgb(203,137,237, maxColorValue=255), rgb(178,204,231, maxColorValue=255))
names(category_colors) <- all_classes

pdf("/home/alexaml/data/INFERNO/gtex_v6p_PC_adjustment/cumulative_variance_per_category.pdf", width=10, height=10, pointsize=12)
ggplot(variance_explained, aes(x=PC, y=cumulative_prop, color=category, group=category)) +
    geom_point() + geom_line() +
    scale_color_manual(name="", values=category_colors) + 
    theme_bw() + xlab("Principal component") + ylab("Cumulative variance explained") +
    ggtitle("Cumulative proportion explained") +
    guides(color=guide_legend(nrow=5)) + 
    theme(legend.position="bottom", text = element_text(size=25),
          legend.text=element_text(size=18), 
          plot.title = element_text(hjust = 0.5)) 
dev.off()

pdf("/home/alexaml/data/INFERNO/gtex_v6p_PC_adjustment/variance_explained_per_category.pdf", width=10, height=10, pointsize=12)
ggplot(variance_explained, aes(x=PC, y=variance_prop, color=category, group=category)) +
    geom_point() + geom_line() +
    scale_color_manual(name="", values=category_colors) + 
    theme_bw() + xlab("Principal component") + ylab("Variance explained") +
    ggtitle("Cumulative proportion explained") +
    guides(color=guide_legend(nrow=5)) + 
    theme(legend.position="bottom", text = element_text(size=25),
          legend.text=element_text(size=18), 
          plot.title = element_text(hjust = 0.5))
dev.off()

