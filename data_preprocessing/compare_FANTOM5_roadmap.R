## compare_FANTOM5_roadmap.R
## alex amlie-wolf 03/14/18
## a script to plot the results of the FANTOM5 / Roadmap comparison

library(ggplot2)
library(plyr)
library(reshape)
library(scales)
sessionInfo()

## -----------------------------------------------------------------------------
## 0. Table of Contents
## -----------------------------------------------------------------------------
## 0. Table of Contents
## 1. Function Definitions
## 2. Read in data, set output path
## 3. Analyze shared enhancer counts

## -----------------------------------------------------------------------------
## 1. Function Definitions
## -----------------------------------------------------------------------------
## make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='png') {
make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='pdf') {
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
## 2. Read in data, set output path
## -----------------------------------------------------------------------------
datadir <- "/home/alexaml/data/INFERNO/FANTOM5_roadmap_comparison/"
outdir <- paste0(datadir, "/analysis_results/")
dir.create(paste0(outdir, '/plots/'), F, T)
dir.create(paste0(outdir, '/tables/'), F, T)

fantom5_roadmap_summary <- read.table(paste0(datadir, '/FANTOM5_roadmap_overlap_summary.txt'), header=T, sep="\t", as.is=T)
fantom5_roadmap_summary$FANTOM5_unique <- with(fantom5_roadmap_summary, FANTOM5_total - shared_intervals)
fantom5_roadmap_summary$roadmap_unique <- with(fantom5_roadmap_summary, roadmap_total - shared_intervals)

write.table(fantom5_roadmap_summary, paste0(outdir, '/tables/FANTOM5_roadmap_overlap_summary.txt'), quote=F, sep="\t", row.names=F)

fantom5_category_df <- read.table("/home/alexaml/data/FANTOM5/Enhancers/fantom5_classes.txt", header=T, sep="\t", quote="", as.is=T, comment.char="")
gtex_category_df <- read.table("/home/alexaml/data/GTEx/gtex_classes.txt", header=T, sep="\t", quote="", as.is=T, comment.char="")
roadmap_category_df <- read.table("/home/alexaml/data/roadmap/roadmap_classes.txt", header=T, sep="\t", quote="", as.is=T, comment.char="")

## find the full list of classes
all_classes <- sort(union(fantom5_category_df$Class, union(gtex_category_df$Class, roadmap_category_df$Class)))

## -----------------------------------------------------------------------------
## 3. Analyze shared enhancer counts
## -----------------------------------------------------------------------------
melt_overlap_summary <- melt(fantom5_roadmap_summary, id.vars="class")

category_colors <- c(rgb(227,191,35, maxColorValue=255), rgb(224,98,247, maxColorValue=255), rgb(36,153,165, maxColorValue=255), rgb(232,2,58, maxColorValue=255), rgb(30,181,54, maxColorValue=255), rgb(253,164,185, maxColorValue=255), rgb(44,115,232, maxColorValue=255), rgb(83,118,43, maxColorValue=255), rgb(203,50,133, maxColorValue=255), rgb(194,196,252, maxColorValue=255), rgb(157,236,192, maxColorValue=255), rgb(141,86,174, maxColorValue=255), rgb(253,133,114, maxColorValue=255), rgb(175,242,253, maxColorValue=255), rgb(149,87,122, maxColorValue=255), rgb(131,233,24, maxColorValue=255), rgb(182,241,140, maxColorValue=255), rgb(252,50,7, maxColorValue=255), rgb(244,148,221, maxColorValue=255), rgb(40,169,123, maxColorValue=255), rgb(247,183,144, maxColorValue=255), rgb(242,184,94, maxColorValue=255), rgb(53,108,145, maxColorValue=255), rgb(198,8,78, maxColorValue=255), rgb(61,115,80, maxColorValue=255), rgb(41,232,215, maxColorValue=255), rgb(122,107,24, maxColorValue=255), rgb(79,153,241, maxColorValue=255), rgb(100,130,128, maxColorValue=255), rgb(166,169,37, maxColorValue=255), rgb(203,137,237, maxColorValue=255), rgb(178,204,231, maxColorValue=255))
names(category_colors) <- gsub(" ", "_", all_classes)

relevant_melt <- melt_overlap_summary[!(melt_overlap_summary$variable %in% c("FANTOM5_total", "roadmap_total")),]

make_graphic(paste0(outdir, "/plots/FANTOM5_roadmap_enh_comparison_grid"), width_ratio = 1.5)
ggplot(relevant_melt, aes(x=class, y=value, fill=class)) +
    scale_fill_manual(name="Tissue Category", values=category_colors) +
    ## scale_alpha_manual(name="", values=c("shared_intervals"=1.0, "FANTOM5_unique"=0.75, "roadmap_unique"=0.5)) +
    facet_grid(variable ~ ., scales="free") + 
    scale_y_continuous(## breaks=seq(0, 400000, by=50000), limits=c(0, 400000), expand=c(0,0), 
                       labels=comma) +
    xlab("Tissue class") + ylab("Number of enhancers") +
    theme_bw() + geom_bar(position="dodge", stat="identity") +
    theme(legend.position="none",
          text=element_text(size=20), 
          legend.text = element_text(size=20),
          axis.text.x = element_text(angle=60, hjust=1, size=20),
          axis.text.y = element_text(size=20),
          title = element_text(size=25))
dev.off()

make_graphic(paste0(outdir, "/plots/FANTOM5_roadmap_enh_comparison"), width_ratio = 1.5)
ggplot(relevant_melt, aes(x=class, y=value, fill=class, alpha=variable)) +
    scale_fill_manual(name="Tissue Category", values=category_colors, guide='none') +
    scale_alpha_manual(name="", values=c("shared_intervals"=1.0, "FANTOM5_unique"=0.75, "roadmap_unique"=0.5)) +
    scale_y_continuous(breaks=seq(0, 400000, by=50000), limits=c(0, 400000), expand=c(0,0), labels=comma) +
    xlab("Tissue class") + ylab("Number of enhancers") +
    theme_bw() + geom_bar(position="dodge", stat="identity") +
    theme(legend.position="bottom",
          legend.text = element_text(size=20),
          axis.text.x = element_text(angle=60, hjust=1, size=20),
          axis.text.y = element_text(size=20),
          title = element_text(size=25))
dev.off()
