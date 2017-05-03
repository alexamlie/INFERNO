## Rscript_run_full_analysis.R
## alex amlie-wolf 04/19/2016
## runs all the pipeline analysis scripts, from the command line
## requires R >= 3.2.0

library(ggplot2)
library(plyr)
library(scales)
library(reshape2)
library(gplots)
library(RColorBrewer)
# library(cowplot)
library(gtools)
library(data.table)
sessionInfo()

## -----------------------------------------------------------------------------
## 0. Table of Contents
## -----------------------------------------------------------------------------

## 0. Table of Contents
## 1. Function Definitions / 'global' variables
## 2. Read in parameter file
## 3. Source and run functions

## -----------------------------------------------------------------------------
## 1. Function Definitions / 'global' variables
## -----------------------------------------------------------------------------
## these functions are inherited by all the analysis scripts
## ---------------------------------------------------
## "GLOBAL" VARIABLES
## columns to ignore for unique SNP analysis (many analysis functions use this)
tagsnp_cols <- c("tag_rsID", "tag_pos", "tag_MAF", "R2", "Dprime")
## decide whether to include parameters and subttitle in plot title or not
SKIP_SUBTITLES <- TRUE
## also define some default text sizes for the axes, to make it easier to change. these don't
## get used in every single function but it makes it a bit easier to have a consistent scheme
AXIS_TEXT_X_SIZE <- 25
AXIS_TEXT_Y_SIZE <- 25
## this affects all titles (plot, axes, legends), unless they are specifically overridden
TITLE_SIZE <- 30
LEGEND_TEXT_SIZE <- 20
## decide whether you want to use tag rsID, region name, or both
tag_naming <- "name"
if(tag_naming=="name") {
    TAG_VAR <- "tag_no_rsid"
    TAG_LAB <- "Tag region"
} else if(tag_naming=="rsID") {
    ## for now, skip this..
} else if(tag_naming=="both") {
    TAG_VAR <- "tag_name"
    TAG_LAB <- "Tag SNP rsID and region"        
} else {
    cat("Tag name value not recognized; using rsID and name\n")
    TAG_VAR <- "tag_name"
    TAG_LAB <- "Tag SNP rsID and region"            
}

## some color definitions for integrative functions
## eRNA enhancers are red
erna_color <- c("eRNA Enh"="#B3112E")
## HMM enhancers are blue (there are three types, so I use 3 similar light blues)
hmm_colors <- c("HMM Enh"="#6393F3", "HMM Genic Enh"="#2D70F3", "HMM Biv Enh"="#0353F1")
merge_hmm_color <- c("HMM Enh"="#6393F3")
## eQTLs are yellow
eqtl_color <- c("eQTL"="#BACC1C")

## also do combinations with merged enhancer states
## eQTL + eRNA enhancer is orange
eqtl_erna_color <- c("eQTL+eRNA Enh"="#D06900")
erna_merged_hmm_color <- c("eRNA Enh+HMM Enh"="#660BAB")
eqtl_merged_hmm_color <- c("eQTL+HMM Enh"="#88BDAD")
eqtl_erna_merged_hmm_color <- c("eQTL+eRNA Enh+HMM Enh"="#12C702")

## ---------------------------------------------------
## FUNCTIONS
## this is used for combining aes and aes_string calls
## https://stackoverflow.com/questions/28777626/how-do-i-combine-aes-and-aes-string-options
`+.uneval` <- function(a,b) {
        `class<-`(modifyList(a,b), "uneval")
    }

## here we define what kind of graphics output we want:
## the height_ratio and width_ratio arguments say how we should shape the output
## change the default type argument to change the type of figure generated
## TODO: add parameter to force specific file type generation
## make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='png') {
make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='pdf') {
    if(type=='pdf') {
        pdf(file=paste0(filename, ".pdf"), width=10*width_ratio, height=10*height_ratio, pointsize=12, onefile=FALSE)
    } else if(type=='png') {
        ## use type='cairo' for when X11 doesn't work
        png(filename=paste0(filename, ".png"), width=10*width_ratio, height=10*height_ratio, res=300, units='in', type='cairo')
    } else if(type=="svg") {
        svg(file=paste0(filename, ".svg"), width=10*width_ratio, height=10*height_ratio, pointsize=12, onefile=FALSE)
    } else {
        cat('filetype not supported\n')
    }
}

## labelling function for facets
## (probably don't even need this here)
dist_r2_labeller <- function(variable, value) {
    if (variable=="dist_thresh") {
        return(paste(value, "bp"))
    } else if(variable=="r2_thresh") {
        return(paste("R^2 >=", value))
    }
}

## a convenience function to check if a parameter is in the reference vector. only returns true
## if the parameter is not set to None or False (from Python), i.e. the only time we want a
## true return value is if the parameter is set for real
check_param <- function(param_vec, param) {
    ret_val <- FALSE
    if (param %in% names(param_vec)) {
        if(param_vec[[param]] != "None" & param_vec[[param]] != "False") {
            ret_val <- TRUE
        }
    }
    ret_val
}

## a function for getting good axis breaks
distance_breaks <- function(limit_vec) {
    minval <- limit_vec[1]
    maxval <- limit_vec[2]
    return(seq(floor(minval/10000)*10000, ceiling(maxval/10000)*10000, by=10000))
}

## -----------------------------------------------------------------------------
## 2. Read in parameter file
## -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
    parameter_file <- args[1]
    out_subtitle <- args[2]
    SKIP_SUBTITLES <- FALSE
} else if (length(args==3)) {
    parameter_file <- args[1]
    out_subtitle <- args[2]
    SKIP_SUBTITLES <- as.logical(args[3])
    if(is.na(SKIP_SUBTITLES)) {
        stop("Error: SKIP_SUBTITLES must be a logical argument!")
    }
} else {
    cat("Usage: Rscript Rscript_run_full_analysis.R parameter_file out_subtitle <SKIP_SUBTITLES>\n")
    cat("Example command: Rscript Rscript_run_full_analysis.R ~/INFERNO_output/parameters/08_26_2016_16:39:10_parameters.txt \"My experiment\" TRUE\n")
    stop("Cannot parse command line arguments")
}

## we have to define this function after SKIP_SUBTITLES is defined
## a function to plot the title of a plot, containing a switch for including the parameters or not
plot_title <- function(title_string, r2_thresh, dist_thresh, out_subtitle, skip_subtitles=SKIP_SUBTITLES, strwrap_width=35) {
    if(skip_subtitles) {
        ## wrap the string!
        ggtitle(paste(strwrap(title_string, width=strwrap_width), collapse="\n"))
    } else {
        ggtitle(paste0(paste(strwrap(title_string, width=strwrap_width), collapse="\n"), "\nR2 >= ", r2_thresh, ", Distance <= ", dist_thresh, " bp\n", out_subtitle))
    }
}
          
## read in the parameters and make a named vector for reference
parameter_tab <- read.table(parameter_file, header=F, sep="\t", quote="", as.is=T, col.names=c("param", "value"))
param_ref <- parameter_tab$value
names(param_ref) <- parameter_tab$param

result_outdir <- paste0(param_ref[['outdir']], "/analysis_results/")
dir.create(result_outdir, F, T)
if(check_param(param_ref, "skip_ld_expansion")) {
    r2_thresh <- "1.0"
    dist_thresh <- "0"
} else {
    r2_thresh <- param_ref[['ld_threshold']]
    dist_thresh <- param_ref[['ld_check_area']]
}

fantom5_class_file <- '/home/alexaml/data/FANTOM5/Enhancers/fantom5_classes.txt'
gtex_class_file <- '/home/alexaml/data/GTEx/gtex_classes.txt'
roadmap_class_file <- '/home/alexaml/data/roadmap/roadmap_classes.txt'

fantom5_category_df <- read.table(fantom5_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")
gtex_category_df <- read.table(gtex_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")
roadmap_category_df <- read.table(roadmap_class_file, header=T, sep="\t", quote="", as.is=T, comment.char="")

## find the full list of classes
all_classes <- sort(union(fantom5_category_df$Class, union(gtex_category_df$Class, roadmap_category_df$Class)))

## now define the color palette
## generated from http://tools.medialab.sciences-po.fr/iwanthue/
## parameters: H 0-360, C 0.46 - 3, L 0.5-1.5
category_colors <- c(rgb(227,191,35, maxColorValue=255), rgb(224,98,247, maxColorValue=255), rgb(36,153,165, maxColorValue=255), rgb(232,2,58, maxColorValue=255), rgb(30,181,54, maxColorValue=255), rgb(253,164,185, maxColorValue=255), rgb(44,115,232, maxColorValue=255), rgb(83,118,43, maxColorValue=255), rgb(203,50,133, maxColorValue=255), rgb(194,196,252, maxColorValue=255), rgb(157,236,192, maxColorValue=255), rgb(141,86,174, maxColorValue=255), rgb(253,133,114, maxColorValue=255), rgb(175,242,253, maxColorValue=255), rgb(149,87,122, maxColorValue=255), rgb(131,233,24, maxColorValue=255), rgb(182,241,140, maxColorValue=255), rgb(252,50,7, maxColorValue=255), rgb(244,148,221, maxColorValue=255), rgb(40,169,123, maxColorValue=255), rgb(247,183,144, maxColorValue=255), rgb(242,184,94, maxColorValue=255), rgb(53,108,145, maxColorValue=255), rgb(198,8,78, maxColorValue=255), rgb(61,115,80, maxColorValue=255), rgb(41,232,215, maxColorValue=255), rgb(122,107,24, maxColorValue=255), rgb(79,153,241, maxColorValue=255), rgb(100,130,128, maxColorValue=255), rgb(166,169,37, maxColorValue=255), rgb(203,137,237, maxColorValue=255), rgb(178,204,231, maxColorValue=255))
names(category_colors) <- all_classes

cat_col_scale <- scale_fill_manual(name="Tissue Category", values=category_colors)

## read in the LD stats dataframe and check the number of tag regions, updating make_graphic if
## necessary to be wide enough to show all the tag regions..
ld_stats_file <- paste0(param_ref[['outdir']], '/ld_expansion/', param_ref[['outprefix']],
                        "_", r2_thresh, "_ld_cutoff_snps_within_", dist_thresh, ".txt")
ld_stats_df <- read.table(ld_stats_file, header=T, sep="\t", quote="", as.is=T)

num_tags <- length(unique(ld_stats_df$tag_name))

## arbitrary cutoff makes all figures twice as wide..should make this more adaptive
if(num_tags > 30) {
    make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='png') {
        ## make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='pdf') {
        if(type=='pdf') {
            pdf(file=paste0(filename, ".pdf"), width=20*width_ratio, height=10*height_ratio, pointsize=12, onefile=FALSE)
        } else if(type=='png') {
            ## use type='cairo' for when X11 doesn't work
            png(filename=paste0(filename, ".png"), width=20*width_ratio, height=10*height_ratio, res=300, units='in', type='cairo')
        } else {
            cat('filetype not supported\n')
        }
    }
}

## -----------------------------------------------------------------------------
## 3. Source and run functions
## -----------------------------------------------------------------------------
cat("Analyzing LD statistics\n")
source("analyze_ld_stats_function.R")
analyze_ld_stats(param_ref[['outprefix']], paste0(param_ref[['outdir']], "/ld_expansion/"),
                 paste0(result_outdir, "/ld_stats/"), out_subtitle, r2_thresh, dist_thresh)

source("analyze_unstranded_genomic_partition_function.R")
if (check_param(param_ref, "unstranded_partition_dir")) {
    cat("Analyzing genomic partition\n")
    analyze_unstranded_genomic_partition(param_ref[['outprefix']],
                                         paste0(param_ref[['outdir']], "/unstranded_genomic_partition/"),
                                         paste0(result_outdir, "/unstranded_genomic_partition/"),
                                         out_subtitle, r2_thresh, dist_thresh)
}

## source("analyze_closest_genes_function.R")
## if (check_param(param_ref, 'gene_bed_file')) {
##     cat("Analyzing closest genes\n")
##     analyze_closest_genes(param_ref[['outprefix']], param_ref[['outdir']],
##                           paste0(result_outdir, "/closest_gene/"), out_subtitle, r2_thresh, dist_thresh)
## }

source("analyze_fantom5_overlap_function.R")
source("analyze_closest_fantom5_enhs_function.R")
if (check_param(param_ref, 'fantom5_dir')) {
    cat("Analyzing FANTOM5 enhancers\n")
    if (check_param(param_ref, 'enhancer_midpoint_window')) {
        analyze_fantom5_overlap(param_ref[['outprefix']], param_ref[['outdir']],
                                paste0(result_outdir, "/fantom5_overlap/"),
                                out_subtitle, r2_thresh, dist_thresh, "midpoint",
                                param_ref[['enhancer_midpoint_window']],
                                fantom5_class_file, !check_param(param_ref, "skip_enh_summary"))
    }

    if (check_param(param_ref, 'enhancer_locus_window')) {
        analyze_fantom5_overlap(param_ref[['outprefix']], param_ref[['outdir']],
                                paste0(result_outdir, "/fantom5_overlap/"),
                                out_subtitle, r2_thresh, dist_thresh, "locus",
                                param_ref[['enhancer_locus_window']],
                                fantom5_class_file, !check_param(param_ref, "skip_enh_summary"))
    }

    if (!(check_param(param_ref, "skip_closest_enh"))) {
        analyze_closest_fantom5_enhs(param_ref[['outprefix']], param_ref[['outdir']],
                                     paste0(result_outdir, "/closest_fantom5_enhs/"),
                                     out_subtitle, r2_thresh, dist_thresh)
    }
}

source("analyze_eqtl_overlap_function.R")
if (check_param(param_ref, "gtex_dir")) {
    cat("Analyzing GTEx eQTLs\n")
    analyze_eqtl_overlap(param_ref[['outprefix']], param_ref[['outdir']],
                               paste0(result_outdir, "/gtex_eqtl_overlap/"),
                               out_subtitle, r2_thresh, dist_thresh, gtex_class_file)
}

source("analyze_factorbook_overlap_function.R")
if (check_param(param_ref, "factorbook_file")) {
    cat("Analyzing FactorBook TFBSs\n")
    analyze_factorbook_overlap(param_ref[['outprefix']], param_ref[['outdir']],
                               paste0(result_outdir, "/factorbook_overlap/"),
                               out_subtitle, r2_thresh, dist_thresh)
}

source("analyze_roadmap_chromhmm_function.R")
if (check_param(param_ref, "roadmap_chromhmm_dir")) {
    cat("Analyzing Roadmap ChromHMM\n")    
    analyze_roadmap_chromHMM(param_ref[['outprefix']], param_ref[['outdir']],
                             paste0(result_outdir, "/roadmap_chromHMM/"),
                             out_subtitle, r2_thresh, dist_thresh, roadmap_class_file)
}

source("analyze_homer_motif_overlaps_function.R")
## just check for the main file
if (check_param(param_ref, "homer_motif_bed_file")) {
    cat("Analyzing HOMER motif overlaps\n")
    ## check for the PWM and sequence files in the parameters to see if PWM changes were calculated
    pwm_calc <- check_param(param_ref, "homer_motif_pwm_file") & check_param(param_ref, "homer_motif_seq_file")
    ## give the true or false value for whether the PWM changes were calculated or not
    analyze_homer_motif_overlaps(param_ref[['outprefix']], param_ref[['outdir']],
                                 paste0(result_outdir, "/homer_motif_overlap/"),
                                 out_subtitle, r2_thresh, dist_thresh, pwm_calc)
}

## functions for overlaps of combinations of data sources:
source("analyze_fantom5_roadmap_enh_overlap_function.R")
if (check_param(param_ref, 'fantom5_dir') & check_param(param_ref, "roadmap_chromhmm_dir")) {
    cat("Checking FANTOM5 / Roadmap overlap\n")
    if (check_param(param_ref, 'enhancer_midpoint_window')) {
        analyze_fantom5_roadmap_enh_overlap(param_ref[['outprefix']], param_ref[['outdir']],
                                            paste0(result_outdir, "/fantom5_roadmap_overlap/"),
                                            out_subtitle, r2_thresh, dist_thresh, "midpoint",
                                            param_ref[['enhancer_midpoint_window']],
                                            fantom5_class_file, roadmap_class_file)
    }

    if (check_param(param_ref, 'enhancer_locus_window')) {
        analyze_fantom5_roadmap_enh_overlap(param_ref[['outprefix']], param_ref[['outdir']],
                                            paste0(result_outdir, "/fantom5_roadmap_overlap/"),
                                            out_subtitle, r2_thresh, dist_thresh, "locus",
                                            param_ref[['enhancer_locus_window']],
                                            fantom5_class_file, roadmap_class_file)
    }
}

## source("analyze_fantom5_eqtl_overlap_function.R")
## ## TODO: allow skipping the correlation file..
## if (check_param(param_ref, 'fantom5_dir') & check_param(param_ref, 'fantom5_correlation_file') & check_param(param_ref, 'gtex_dir')) {
##     cat("Analyzing FANTOM5 and eQTL overlap\n")
##     if (check_param(param_ref, 'enhancer_midpoint_window')) {
##         analyze_fantom5_eqtl_overlap(param_ref[['outprefix']], param_ref[['outdir']],
##                                      paste0(result_outdir, "/fantom5_eqtl_overlap/"),
##                                      out_subtitle, r2_thresh, dist_thresh, "midpoint",
##                                      fantom5_class_file, gtex_class_file)
##     }

##     if (check_param(param_ref, 'enhancer_locus_window')) {
##         analyze_fantom5_eqtl_overlap(param_ref[['outprefix']], param_ref[['outdir']],
##                                      paste0(result_outdir, "/fantom5_eqtl_overlap/"),
##                                      out_subtitle, r2_thresh, dist_thresh, "locus",
##                                      fantom5_class_file, gtex_class_file)
##     }
## }

## source("analyze_tfbs_fantom5_eqtl_overlap_function.R")
## if (check_param(param_ref, "fantom5_dir") & check_param(param_ref, "gtex_dir") & check_param(param_ref, "factorbook_file")) {
##     cat("Checking FANTOM5, FactorBook, eQTL  overlap\n")
##     if (check_param(param_ref, 'enhancer_midpoint_window')) {
##         analyze_tfbs_fantom5_eqtl_overlap(param_ref[['outprefix']], param_ref[['outdir']],
##                                           paste0(result_outdir, "/tfbs_fantom5_eqtl_overlap/"),
##                                           out_subtitle, r2_thresh, dist_thresh, "midpoint")
##     }
##     if (check_param(param_ref, 'enhancer_locus_window')) {
##         analyze_tfbs_fantom5_eqtl_overlap(param_ref[['outprefix']], param_ref[['outdir']],
##                                           paste0(result_outdir, "/tfbs_fantom5_eqtl_overlap/"),
##                                           out_subtitle, r2_thresh, dist_thresh, "locus")
##     }
## }

source("analyze_fantom5_eqtl_chromHMM_overlap_function.R")
if (check_param(param_ref, "fantom5_dir") & check_param(param_ref, "gtex_dir") & check_param(param_ref, "roadmap_chromhmm_dir")) {
    cat("Analyzing FANTOM5, eQTL, ChromHMM overlap\n")
    if (check_param(param_ref, 'enhancer_midpoint_window')) {
        analyze_fantom5_eqtl_chromHMM_overlap(param_ref[['outprefix']], param_ref[['outdir']],
                                              paste0(result_outdir, "/fantom5_eqtl_chromHMM_overlap/"),
                                              out_subtitle, r2_thresh, dist_thresh, "midpoint",
                                              param_ref[['enhancer_midpoint_window']],
                                              fantom5_class_file, gtex_class_file, roadmap_class_file)
    }
    if (check_param(param_ref, 'enhancer_locus_window')) {
        analyze_fantom5_eqtl_chromHMM_overlap(param_ref[['outprefix']], param_ref[['outdir']],
                                              paste0(result_outdir, "/fantom5_eqtl_chromHMM_overlap/"),
                                              out_subtitle, r2_thresh, dist_thresh, "locus",
                                              param_ref[['enhancer_locus_window']],
                                              fantom5_class_file, gtex_class_file, roadmap_class_file)
    }
}
