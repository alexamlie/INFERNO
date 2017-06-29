# enhancer_only_sample_and_expand_matched_input_variants.R
## alex amlie-wolf 05/03/17
## a script to sample some number of matched variants by LD, distance to nearest TSS, and MAF
## this script just directly uses the input variant list and matches those, then does LD expansion
## unlike sample_and_expand_matched_input_variants.R, this one only bootstraps against enhancer
## annotations (since direct eQTL overlap is not the right way to get eQTL overlap)
## NOTE: this uses around 30 Gb of memory to achieve fast annotation overlapping
## this script runs automatically with Rscript

library(plyr)
library(ggplot2)
library(scales)
library(reshape2)
library(gtools)
sessionInfo()

## -----------------------------------------------------------------------------
## 0. Table of Contents
## -----------------------------------------------------------------------------
## 1. Function definitions
## 2. Parameter settings, read in data
## 3. Perform sampling against input
## 4. Analysis plots for sampling
## 5. LD expansion of sampled sets and input variants
## 6. Analysis plots for LD expansion
## 7. Compute annotation overlap p-values
## 8. Analysis plots for annotation overlap and p-values of non-LD collapsed analysis
## 9. Analysis plots for annotation overlap and p-values of LD collapsed analysis

## -----------------------------------------------------------------------------
## 1. Function Definitions
## -----------------------------------------------------------------------------
# make_graphic <- function(filename, width_ratio=1, height_ratio=1, type='png') {
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

full_analysis_start_time <- proc.time()    
## -----------------------------------------------------------------------------
## 2. Parameter settings, read in data
## -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
cat(args, "\n")
## stop("Done")
if (length(args)==9) {
    ## sampling parameters
    num_samples <- as.numeric(args[1])
    maf_bin_size <- as.numeric(args[2])
    dist_round <- as.numeric(args[3])
    ## parameters for upper thresholds on distance and number of LD partners, beyond which all SNPs
    ## are lumped together into the 'high' category (set to Inf to not do this)
    dist_threshold <- as.numeric(args[4])
    ld_partner_threshold <- as.numeric(args[5])

    ## file with the precomputed background information (for bin matching)
    bg_snp_info_f <- args[6]
    ## precomputed pairwise LD results
    ld_sets_dir <- args[7]
    ## path to the annotation overlaps for all 1kg SNPs
    ref_summary_dir <- args[8]
    
    ## pipeline parameter file for input
    parameter_file <- args[9]    
} else {
    cat(args, '\n')
    stop("Not enough arguments")    
}

parameter_tab <- read.table(parameter_file, header=F, sep="\t", quote="", as.is=T, col.names=c("param", "value"))
param_ref <- parameter_tab$value
names(param_ref) <- parameter_tab$param
## define some variables using the param references
output_dir <- paste0(param_ref[['outdir']], '/background_enh_sampling_match_and_expand/')
if(check_param(param_ref, "skip_ld_expansion")) {
    r2_thresh <- "1.0"
    dist_thresh <- "0"
} else {
    r2_thresh <- param_ref[['ld_threshold']]
    dist_thresh <- param_ref[['ld_check_area']]
}

## check that the background SNP info is a real file
if (!file.exists(bg_snp_info_f)) {
    stop(paste0("Error: background SNP information does not exist for supplied R^2 threshold:", r2_thresh))
}

dir.create(paste0(output_dir, "/samples/"), F, T)
dir.create(paste0(output_dir, "/plots/"), F, T)
dir.create(paste0(output_dir, "/tables/"), F, T)

{
## read in data:
## all 1kg SNPs
snp_info <- read.table(bg_snp_info_f, header=T, sep="\t", quote="",
                       comment.char="", as.is=T)
## fix '.' rsIDs
snp_info$rsID <- ifelse(snp_info$rsID==".", paste0(snp_info$chr, "-", snp_info$pos), snp_info$rsID)

## for the LD expansion, we need to calculate offsets for each chromosome for indexing purposes
chr_offsets <- match(unique(snp_info$chr), snp_info$chr)
names(chr_offsets) <- unique(snp_info$chr)
}

## input SNP file: notice that this uses the sorted input list, not the LD expanded!
input_snp_f <- paste0(param_ref[['outdir']], '/ld_expansion/', param_ref[['outprefix']],
                        "_sorted_snps_", r2_thresh, "_", dist_thresh, ".txt")
full_input_snps <- read.table(input_snp_f, header=F, sep="\t", quote="",
                         col.names=c("chr", "rsID", "tag_name", "pos"), comment.char="", as.is=T)
## get only unique SNPs (this should almost never be necessary)
input_snps <- unique(full_input_snps)
## count the number of regions (this comes in handy throughout the script)
num_regions <- length(unique(input_snps$tag_name))
rm(full_input_snps)

## annotate the input SNPs with the relevant information
annot_input_snps <- merge(input_snps, snp_info, by=c("chr", "rsID", "pos"),
                          all.x=F, all.y=F)

## to find mismatches by rsID:
mismatch_snps <- input_snps[!(input_snps$rsID %in% annot_input_snps$rsID),]
cat(nrow(mismatch_snps), "SNPs did not match up correctly by rsID:", mismatch_snps$rsID, "\n")

## if there are mismatches, try to match by position
if(nrow(mismatch_snps) > 0) {
    mismatch_merge_snps <- merge(mismatch_snps, snp_info, by=c("chr", "pos"), all=F,
                                 suffixes=c(".x", ""))

    if(nrow(mismatch_merge_snps) > 0) {
        cat(nrow(mismatch_merge_snps), "SNP matches found by position. New rsIDs from 1000 Genomes:",
            mismatch_merge_snps$rsID, "\nInput rsIDs:", mismatch_merge_snps$rsID.x, "\n")

        ## make sure the number of SNPs matches up
        if(nrow(mismatch_merge_snps)==nrow(mismatch_snps)) {
            annot_input_snps <- rbind(annot_input_snps, mismatch_merge_snps[,c("chr", "rsID", "pos", "tag_name", "tss_dist", "MAF", "num_ld_partners")])
        } else {
            stop("Did not match the same number of variants!\n")
        }
    } else {
        stop("Could not match rsIDS! :( \n")
    }
}


## -----------------------------------------------------------------------------
## 3. Perform sampling against input
## -----------------------------------------------------------------------------
start_time <- proc.time()

## now define each variable that we want to bin on to count combinations
## TSS distance, rounded to the distance bin parameter. note that we have to use ceiling here,
## otherwise a TSS distance of <(dist_round / 2) will be called 0
input_dist_bins <- round_any(annot_input_snps$tss_dist, dist_round, f=ceiling)
## apply the distance threshold
input_dist_bins <- ifelse(input_dist_bins >= dist_threshold, dist_threshold, input_dist_bins)

## MAF bins, using a floor binning approach
input_maf_bins <- round_any(annot_input_snps$MAF, maf_bin_size, f=floor)

## we don't have to define an LD partner variable since we directly match
input_ld_partners <- annot_input_snps$num_ld_partners
input_ld_partners <- ifelse(input_ld_partners >= ld_partner_threshold, ld_partner_threshold,
                            input_ld_partners)

## now, we have to count the number of input SNPs with each combination for each tag region
input_freq_table <- count(data.frame(region=annot_input_snps$tag_name, input_dist_bins, input_maf_bins, input_ld_partners, stringsAsFactors = F))

## for each combination of binned variables, we will need to sample the same number of SNPs
## from the background distribution and apply the thresholds
bg_dist_bins <- round_any(snp_info$tss_dist, dist_round, f=ceiling)
bg_dist_bins <- ifelse(bg_dist_bins >= dist_threshold, dist_threshold, bg_dist_bins)

bg_maf_bins <- round_any(snp_info$MAF, maf_bin_size, f=floor)

bg_ld_partners <- snp_info$num_ld_partners
bg_ld_partners <- ifelse(bg_ld_partners >= ld_partner_threshold, ld_partner_threshold,
                         bg_ld_partners)

## to make this quicker, pre-compute the required indexing vectors for the background SNPs
bg_dist_vecs <- lapply(sort(unique(input_dist_bins)), function(x) {
    which(bg_dist_bins == x) })
names(bg_dist_vecs) <- as.character(sort(unique(input_dist_bins)))

bg_maf_vecs <- lapply(sort(unique(input_maf_bins)), function(x) {
    which(bg_maf_bins == x) })
names(bg_maf_vecs) <- as.character(sort(unique(input_maf_bins)))

bg_ld_vecs <- lapply(sort(unique(input_ld_partners)), function(x) {
    which(bg_ld_partners == x) })
names(bg_ld_vecs) <- as.character(sort(unique(input_ld_partners)))

## initialize an n x m matrix that stores all the samples, where n is the number of input SNPs
## and m is the number of samples. this matrix stores the (index of the) SNPs matched against
## the input
sample_mat <- matrix(-1, nrow=nrow(annot_input_snps), ncol=num_samples)
## also store which regions each of these came from. this is the same for all samples, so it's just a vector
sample_regions <- vector("character", length=nrow(annot_input_snps))

## we need to track which row of this matrix we're on since sometimes we sample more than one SNP
cur_row <- 1

for(i in seq(nrow(input_freq_table))) {
    ## get the current bin values, the number of variants, and the tag region 
    this_dist <- input_freq_table$input_dist_bins[i]
    this_maf <- input_freq_table$input_maf_bins[i]
    this_ld <- input_freq_table$input_ld_partners[i]
    this_num <- input_freq_table$freq[i]
    this_tag <- input_freq_table$region[i]
    
    ## find the indices of all background SNPs with this combination of bins
    match_idxs <- Reduce(intersect, list(bg_dist_vecs[[as.character(this_dist)]],
                                         bg_maf_vecs[[as.character(this_maf)]],
                                         bg_ld_vecs[[as.character(this_ld)]]))

    cat(paste("Tag region", this_tag, "Freq row", i, "mat row", cur_row, "dist", this_dist, "MAF",
              this_maf, "LD", this_ld, "input SNP num", this_num, "bg matches", length(match_idxs)), "\n")

    if(this_num > length(match_idxs)) { cat("Too many input SNPs!\n") }

    ## store the tag region for the variant(s)
    sample_regions[cur_row:(cur_row + this_num - 1)] <- this_tag
    
    if(length(match_idxs)==1) {
        ## if we only have one match, we have to explicitly use that one SNP, otherwise the
        ## sample function thinks we are giving it the length of the vector to sample from
        sample_mat[cur_row:(cur_row + this_num - 1),] <- match_idxs[1]
    } else {
        for(n in 1:num_samples) {
            this_samp <- sample(as.vector(match_idxs), this_num)
            sample_mat[cur_row:(cur_row + this_num - 1), n] <- this_samp
        }
    }

    cur_row <- cur_row + this_num
}

## write out a file with the sampled rsID indices
write.table(sample_mat, paste0(output_dir, "/samples/input_sample_mat.txt"), quote=F, sep="\t",
            row.names=F, col.names=F)
## also write the regions
write(sample_regions, paste0(output_dir, "/samples/input_sample_regions.txt"), sep="\n")
## ## to read this sample matrix back in, uncomment these:
## sample_mat <- as.matrix(read.table(paste0(output_dir, "/samples/input_sample_mat.txt"),
##                                    header=F, sep="\t", quote="", as.is=T))
## sample_regions <- read.table(paste0(output_dir, "/samples/input_sample_regions.txt"),
##                              header=F, quote="", as.is=T)$V1
cat("Sampling against input and writing out file took:\n")
cat((proc.time() - start_time)[["elapsed"]], 'seconds\n')

## for further analysis, get the sampled SNPs information: we can just do this by directly
## indexing the info table and then make a huge data frame with all the SNP info (including
## both bins and raw values), retaining tag region information
info_cols <- c("chr","rsID", "pos", "tss_dist", "MAF", "num_ld_partners")
## include a reference column for input vs sample and another column for sample number
## use 0 as the sample num for the input
sampled_snp_info <- rbind(cbind(data="input", samp_num=0, annot_input_snps[,info_cols],
                                region=annot_input_snps$tag_name, ld_bin=input_ld_partners,
                                dist_bin=input_dist_bins, maf_bin=input_maf_bins),
                          cbind(data="sample", samp_num=rep(1:num_samples, each=nrow(annot_input_snps)), snp_info[sample_mat, info_cols], region=rep(sample_regions, times=num_samples),
                                ld_bin=bg_ld_partners[sample_mat],
                                dist_bin=bg_dist_bins[sample_mat], maf_bin=bg_maf_bins[sample_mat]))

## -----------------------------------------------------------------------------
## 4. Analysis plots for sampling
## -----------------------------------------------------------------------------
## plot the distributions of how many samples each SNP is found in
## first we count the number of samples containing each SNP
snp_sample_counts <- count(c(sample_mat))
## then we count how many times any SNP was sampled once, twice, etc..
sample_dist_counts <- count(snp_sample_counts$freq)
## define the plot breaks
x_breaks <- seq(0, num_samples, by=num_samples / 20)
y_breaks <- c(0, 100, 500, 1000, 2000, 3000, 4000, seq(5000, round_any(max(table(snp_sample_counts$freq)), accuracy=5000, f=ceiling), by=5000))

dir.create(paste0(output_dir, '/plots/input_sampling'), F, T)

make_graphic(paste0(output_dir, "/plots/input_sampling/", param_ref[['outprefix']], "_snp_sample_num_hist"),
             width_ratio=1.5)
print(ggplot(sample_dist_counts, aes(x=x, y=freq, colour=x)) +
    geom_bar(stat="identity", position="stack") +
    scale_colour_continuous(low="#132B43", high=muted("red")) +
    theme_bw() + xlab("Number of sampled data sets containing SNP") +
    ylab("Number of SNPs") +
    scale_x_continuous(breaks=x_breaks, limits=c(-1, tail(x_breaks, 1)+1), oob=rescale_none) +
    scale_y_continuous(breaks=y_breaks, labels=comma, trans="sqrt") +
    ggtitle("Histogram of number of samples containing each SNP") +
    theme(legend.position="none", axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=25),
          axis.text.y = element_text(size=25), 
          title=element_text(size=30), plot.title = element_text(hjust = 0.5, size=30)))
dev.off()

## -----------------------
## make proportional histograms for input and sample distributions. to do this, we do our own
## per-group proportion counting!
## start with the TSS distance
dist_breaks <- seq(0, round_any(max(sampled_snp_info$dist_bin), dist_round), dist_round)
input_dist_counts <- as.data.frame(table(factor(sampled_snp_info$dist_bin[sampled_snp_info$data=="input"],
                                                levels=dist_breaks)))
sample_dist_counts <- as.data.frame(table(factor(sampled_snp_info$dist_bin[sampled_snp_info$data=="sample"],
                                                 levels=dist_breaks)))

## now combine these into one data frame with the counts and proportions
tss_dist_prop_df <- rbind(data.frame(data_source="input", dist_bin = dist_breaks,
                                     counts = input_dist_counts$Freq,
                                     props = input_dist_counts$Freq/sum(input_dist_counts$Freq)),
                          data.frame(data_source="sample", dist_bin = dist_breaks,
                                     counts = sample_dist_counts$Freq,
                                     props = sample_dist_counts$Freq/sum(sample_dist_counts$Freq)))
## make sure these are ordered correctly
tss_dist_prop_df$data_source <- factor(tss_dist_prop_df$data_source, ordered=T,
                                       levels=c("input", "sample"))

## start with the full distance histogram
dist_labels <- seq(0, tail(dist_breaks, 1)+1, by=dist_round*10)

make_graphic(paste0(output_dir, "/plots/input_sampling/", param_ref[['outprefix']], "_input_and_sampled_tss_dist_bins_prop_hist"))
print(ggplot(tss_dist_prop_df, aes(x=dist_bin, y=props, fill=data_source)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(breaks=dist_labels, limits=c(-1, tail(dist_breaks, 1)+1), oob=rescale_none,
                       labels=if(dist_threshold==Inf) comma_format()(dist_labels) else
                           c(comma_format()(seq(0, tail(dist_breaks, 2)[1], by=dist_round*10)), paste0(comma_format()(dist_threshold), "+"))) +
    scale_y_continuous(breaks=seq(0, 1, by=0.05)) +
    theme_bw() + xlab("Rounded distance to nearest gene (bp)") +
    ylab("Proportion of variants within each class") +
    theme(legend.position="bottom",
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=25),
          axis.text.y = element_text(size=25), 
          title=element_text(size=30), plot.title = element_text(hjust = 0.5, size=30)))
dev.off()

## also make zoomed plots. first find the 95th quantile cutoff
dist_zoom_cutoff <- tss_dist_prop_df$dist_bin[which(cumsum(tss_dist_prop_df$props)>=.95)[1]]

make_graphic(paste0(output_dir, "/plots/input_sampling/", param_ref[['outprefix']], "_input_and_sampled_tss_dist_bins_zoomed_prop_hist"))
print(ggplot(tss_dist_prop_df, aes(x=dist_bin, y=props, fill=data_source)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(limits=c(-1, dist_zoom_cutoff), labels=comma,
                       breaks=seq(0, dist_zoom_cutoff, by=dist_round*5), oob=rescale_none) +
    scale_y_continuous(breaks=seq(0, 1, by=0.05)) +
    theme_bw() + xlab("Rounded distance to nearest gene (bp)") +
    ylab("Proportion of variants within each class") +
    theme(legend.position="bottom",
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=25),
          axis.text.y = element_text(size=25), 
          title=element_text(size=30), plot.title = element_text(hjust = 0.5, size=30)))
dev.off()

## now do the same thing for the MAF data
maf_breaks <- seq(0, 0.5, maf_bin_size)
input_maf_counts <- as.data.frame(table(factor(sampled_snp_info$maf_bin[sampled_snp_info$data=="input"],
                                                levels=maf_breaks)))
sample_maf_counts <- as.data.frame(table(factor(sampled_snp_info$maf_bin[sampled_snp_info$data=="sample"],
                                                 levels=maf_breaks)))

## now combine these into one data frame with the counts and proportions
maf_prop_df <- rbind(data.frame(data_source="input", maf_bin = maf_breaks,
                                counts = input_maf_counts$Freq,
                                props = input_maf_counts$Freq/sum(input_maf_counts$Freq)),
                     data.frame(data_source="sample", maf_bin = maf_breaks,
                                counts = sample_maf_counts$Freq,
                                props = sample_maf_counts$Freq/sum(sample_maf_counts$Freq)))
## make sure these are ordered correctly
maf_prop_df$data_source <- factor(maf_prop_df$data_source, ordered=T,
                                       levels=c("input", "sample"))

make_graphic(paste0(output_dir, "/plots/input_sampling/", param_ref[['outprefix']], "_input_and_sampled_maf_bins_prop_hist"))
print(ggplot(maf_prop_df, aes(x=maf_bin, y=props, fill=data_source)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(limits=c(-0.01, 0.51), breaks=maf_breaks, oob=rescale_none) +
    scale_y_continuous(breaks=seq(0, 1, by=0.01)) +
    theme_bw() + xlab("Minor allele frequency bin") +
    ylab("Proportion of variants") +
    theme(legend.position="bottom",
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=25),
          axis.text.y = element_text(size=25), 
          title=element_text(size=30), plot.title = element_text(hjust = 0.5, size=30)))
dev.off()

if (ld_partner_threshold==Inf) {
    ld_breaks <- seq(0, max(sampled_snp_info$ld_bin), by=1)
    ld_labels <- comma_format()(ld_breaks)
} else {
    ld_breaks <- seq(0, ld_partner_threshold, by=ld_partner_threshold / 20)
    ld_labels <- c(comma_format()(head(ld_breaks, length(ld_breaks)-1)),
                   paste0(comma_format()(ld_partner_threshold), "+"))
}

input_ld_counts <- as.data.frame(table(factor(sampled_snp_info$ld_bin[sampled_snp_info$data=="input"],
                                              levels=ld_breaks)))
sample_ld_counts <- as.data.frame(table(factor(sampled_snp_info$ld_bin[sampled_snp_info$data=="sample"],
                                                 levels=ld_breaks)))

## now combine these into one data frame with the counts and proportions
ld_prop_df <- rbind(data.frame(data_source="input", ld_bin = ld_breaks,
                               counts = input_ld_counts$Freq,
                               props = input_ld_counts$Freq/sum(input_ld_counts$Freq)),
                    data.frame(data_source="sample", ld_bin = ld_breaks,
                               counts = sample_ld_counts$Freq,
                               props = sample_ld_counts$Freq/sum(sample_ld_counts$Freq)))
## make sure these are ordered correctly
ld_prop_df$data_source <- factor(ld_prop_df$data_source, ordered=T,
                                 levels=c("input", "sample"))

## define a custom label for this
ld_partner_labels <- rep("", length(ld_labels))
ld_partner_labels[seq(1, length(ld_partner_labels), 5)] <- ld_labels[seq(1, length(ld_labels), 5)]

make_graphic(paste0(output_dir, "/plots/input_sampling/", param_ref[['outprefix']], "_input_and_sampled_num_ld_partners_prop_hist"))
print(ggplot(ld_prop_df, aes(x=ld_bin, y=props, fill=data_source)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(limits=c(-1, tail(ld_breaks,1)+1),
                       breaks=ld_breaks, oob=rescale_none, labels=ld_partner_labels) +
    scale_y_continuous(breaks=seq(0, 1, by=0.01)) +
    theme_bw() + xlab("Number of LD partners") +
    ylab("Proportion of variants") +
    theme(legend.position="bottom",
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=25),
          axis.text.y = element_text(size=25), 
          title=element_text(size=30), plot.title = element_text(hjust = 0.5, size=30)))
dev.off()

## make a zoomed ld partner plot
ld_zoom_cutoff <- round_any(quantile(sampled_snp_info$ld_bin, .95), 10, ceiling)

make_graphic(paste0(output_dir, "/plots/input_sampling/", param_ref[['outprefix']], "_input_and_sampled_num_ld_partners_zoomed_prop_hist"))
print(ggplot(ld_prop_df, aes(x=ld_bin, y=props, fill=data_source)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(limits=c(-1, ld_zoom_cutoff),
                       breaks=seq(0, ld_zoom_cutoff, by=10), oob=rescale_none) +
    scale_y_continuous(breaks=seq(0, 1, by=0.01)) +
    theme_bw() + xlab("Number of LD partners") +
    ylab("Proportion of variants") +
    theme(legend.position="bottom",
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=25),
          axis.text.y = element_text(size=25), 
          title=element_text(size=30), plot.title = element_text(hjust = 0.5, size=30)))
dev.off()


## -----------------------------------------------------------------------------
## 5. LD expansion of sampled sets and input variants
## -----------------------------------------------------------------------------
start_time <- proc.time()
## in order to do LD expansion split by chromosome, we have to find out what chromosomes each
## of our sampled SNPs came from
sample_chr_mat <- matrix(snp_info[sample_mat,"chr"], nrow=nrow(sample_mat), ncol=num_samples)

## now we make a matrix, since we know how many SNPs should be expanded from each set because
## we match on number of LD partners. we can just sum this number, but we also have to include
## counts for the input SNPs themselves, since the LD partner count doesn't include the 'tag'
expanded_snp_sets <- matrix(NA, nrow=sum(input_ld_partners)+length(input_ld_partners),
                            ncol=num_samples)
## also keep track of which tag regions each one came from
expanded_snp_regions <- matrix("", nrow=sum(input_ld_partners)+length(input_ld_partners),
                               ncol=num_samples)
## and for each expanded SNP, we want to know which LD-tagging variant it came from
## this is for the LD-collapsed bootstrap analysis
expanded_snp_blocks <- matrix(NA, nrow=sum(input_ld_partners)+length(input_ld_partners),
                               ncol=num_samples)

## we also need to store indices to use for filling in each sample in each chromosome (since
## sampling is not matched by chromosome, each sample will have different expanded numbers at
## each chromosome)
samp_idxs <- rep(1, times=num_samples)

## also do this for the input SNPs, which we use for LD-collapsed counts later
input_snp_idx <- match(annot_input_snps$rsID, snp_info$rsID)
expanded_input_snps <- vector("character", length=sum(input_ld_partners)+length(input_ld_partners))
expanded_input_regions <- vector("character", length=sum(input_ld_partners)+length(input_ld_partners))
expanded_input_blocks <- vector("numeric", length=sum(input_ld_partners)+length(input_ld_partners))
input_idx <- 1

for(this_chr in mixedsort(unique(c(sample_chr_mat)))) {
    cat("Performing LD expansion on", this_chr, "\n")
    ## get the file with the precomputed LD information for this chromosome
    this_chr_ld_file <- paste0(ld_sets_dir, this_chr, "_LD_sets_", r2_thresh, "_threshold.txt")
    if (!file.exists(this_chr_ld_file)) {
        stop(paste0("Error: LD sets were not precomputed for supplied R^2 threshold:", r2_thresh, "Run the compute_ld_sets.py script to generate these files"))
    }
    ## now read in this information
    this_chr_ld_sets <- read.table(this_chr_ld_file, header=T, sep="\t", quote="", as.is=T)

    ## ## calculate the size of the expanded set:
    ## ## (this stuff is a sanity check for my other summarization scripts)
    ## this_chr_ld_sets$set_size <- unlist(lapply(strsplit(this_chr_ld_sets$expanded_set, ","), length))
    ## ## also calculate the number of unique SNPs in each expanded set
    ## this_chr_ld_sets$uniq_snp_count <- unlist(lapply(strsplit(this_chr_ld_sets$expanded_set, ","), function(x) {length(unique(x))}))
    ## cat("Number of sets with duplicated SNPs:",
    ##     sum(this_chr_ld_sets$set_size != this_chr_ld_sets$uniq_snp_count), '\n')

    ## now add an indexing column so that we can build an easily indexable data structure
    this_chr_ld_sets$snp_idx <- match(this_chr_ld_sets$rsID, snp_info$rsID)

    ## check how many don't match up
    num_not_found <- sum(is.na(this_chr_ld_sets$snp_idx))
    cat(num_not_found, "LD SNPs were not found in the SNP information DF\n")

    ## these SNPs are SNPs that are observed multiple times in the MAF data, often
    ## corresponding to indels and SNVs at a single position, so we filter these out
    this_chr_ld_sets <- this_chr_ld_sets[!is.na(this_chr_ld_sets$snp_idx),]

    ## ## check for mismatches between the precomputed set size and the SNP info
    ## ## uncomment this to check it, but it's slow
    ## set_sizes <- unlist(lapply(this_chr_ld_sets$expanded_set, function(x) { length(strsplit(x, ",")[[1]]) }))
    ## set_mismatches <- sum(set_sizes != snp_info$num_ld_partners[this_chr_ld_sets$snp_idx])
    ## cat(set_mismatches, "SNP LD sets did not match up with summary information\n")
    ## if(set_mismatches > 0) {
    ##     ## find the average difference
    ##     set_differences <- abs(set_sizes - snp_info$num_ld_partners[this_chr_ld_sets$snp_idx])
    ##     cat("Mean size difference was", mean(set_differences[set_differences!=0]), "\n")
    ## }
    
    ## next, we use the offset information to map these matching indices to start at 1 for the chromosome we're on
    this_chr_offset <- chr_offsets[this_chr]
    this_chr_ld_sets$norm_snp_idx <- this_chr_ld_sets$snp_idx - (this_chr_offset-1)

    ## initialize a list that we can use to index and get expanded SNP sets
    ## this list initially contains just individual SNPs, then we expand the sets
    snp_set_list <- as.list(snp_info$rsID[snp_info$chr==this_chr])
    ## for all the SNPs that we have expanded sets for, add these to the vectors
    snp_set_list <- replace(snp_set_list, this_chr_ld_sets$norm_snp_idx,
                            mapply(c, snp_set_list[this_chr_ld_sets$norm_snp_idx],
                                   lapply(this_chr_ld_sets$expanded_set, function(x) {
                                       strsplit(x, split=",")[[1]] })))

    ## now go through each sample and sample the LD expanded sets
    for(s in seq(num_samples)) {
        ## index the ('tag') SNPs in this sample that are from this chromosome
        samp_chr_idx <- which(sample_chr_mat[,s]==this_chr)
        this_samp_idxs <- sample_mat[samp_chr_idx,s]
        ## also get the tag regions of these SNPs
        this_samp_regions <- sample_regions[samp_chr_idx]

        ## get the expanded set of SNP rsIDs
        expanded_snps <- snp_set_list[this_samp_idxs-(this_chr_offset-1)]
        ## collapse this into a single vector to put it in the expanded matrix
        this_expanded_set <- unlist(expanded_snps)
        this_set_size <- length(this_expanded_set)
        ## to store the region names, we have to figure out the size of the expanded sets for
        ## each variant
        this_expanded_regions <- rep(this_samp_regions,times=unlist(lapply(expanded_snps, length)))
        ## same for the LD blocks: just use the indices of the SNPs
        this_expanded_blocks <- rep(this_samp_idxs, times=unlist(lapply(expanded_snps, length)))
        
        ## add this to the expanded SNP matrices, using the index offset for this sample
        if(this_set_size > 0) {
            expanded_snp_sets[samp_idxs[s]:(samp_idxs[s]+this_set_size-1),s] <- this_expanded_set
            expanded_snp_regions[samp_idxs[s]:(samp_idxs[s]+this_set_size-1),s] <- this_expanded_regions
            expanded_snp_blocks[samp_idxs[s]:(samp_idxs[s]+this_set_size-1),s] <- this_expanded_blocks
            ## now update the index offset for this sample
            samp_idxs[s] <- samp_idxs[s] + this_set_size
        }
    }

    ## also do this for the input SNPs
    ## index the input SNPs that are from this chromosome
    input_chr_idx <- which(annot_input_snps$chr==this_chr)
    ## get the index into the full SNP information
    this_input_idxs <- input_snp_idx[input_chr_idx]
    ## also get the tag regions of these SNPs
    this_input_regions <- annot_input_snps$tag_name[input_chr_idx]

    ## get the expanded set of SNP rsIDs
    expanded_snps <- snp_set_list[this_input_idxs-(this_chr_offset-1)]
    ## collapse this into a single vector to put it in the expanded matrix
    this_expanded_set <- unlist(expanded_snps)
    this_set_size <- length(this_expanded_set)
    ## to store the region names, we have to figure out the size of the expanded sets for
    ## each variant
    this_expanded_regions <- rep(this_input_regions, times=unlist(lapply(expanded_snps, length)))
    ## same for the LD blocks: just use the indices of the SNPs
    this_expanded_blocks <- rep(this_input_idxs, times=unlist(lapply(expanded_snps, length)))

    ## add this to the expanded SNP matrices, using the index offset for this sample
    if(this_set_size > 0) {
        expanded_input_snps[input_idx:(input_idx+this_set_size-1)] <- this_expanded_set
        expanded_input_regions[input_idx:(input_idx+this_set_size-1)] <- this_expanded_regions
        expanded_input_blocks[input_idx:(input_idx+this_set_size-1)] <- this_expanded_blocks
        ## now update the index offset for the input
        input_idx <- input_idx + this_set_size
    }
    
}

## finally, we need to convert this to indices into the snp_info matrix
expanded_snp_idx_mat <- matrix(match(c(expanded_snp_sets), snp_info$rsID),
                               nrow=nrow(expanded_snp_sets), ncol=num_samples)
## no need to store the rsIDs
rm(expanded_snp_sets)

## same for the input data
expanded_input_idxs <- match(expanded_input_snps, snp_info$rsID)
rm(expanded_input_snps)

## write out the files of the full expanded SNP sets and their tag regions
write.table(expanded_snp_idx_mat, paste0(output_dir, "/samples/expanded_sample_mat.txt"), quote=F,
            sep="\t", row.names=F, col.names=F)
write.table(expanded_snp_regions, paste0(output_dir, "/samples/expanded_sample_regions.txt"), quote=F,
            sep="\t", row.names=F, col.names=F)
write.table(expanded_snp_blocks, paste0(output_dir, "/samples/expanded_sample_blocks.txt"), quote=F,
            sep="\t", row.names=F, col.names=F)
## same for the input SNP info
write.table(expanded_input_idxs, paste0(output_dir, "/samples/expanded_input_idxs.txt"), quote=F,
            sep="\t", row.names=F, col.names=F)
write.table(expanded_input_regions, paste0(output_dir, "/samples/expanded_input_regions.txt"), quote=F,
            sep="\t", row.names=F, col.names=F)
write.table(expanded_input_blocks, paste0(output_dir, "/samples/expanded_input_blocks.txt"), quote=F,
            sep="\t", row.names=F, col.names=F)

## ## to read them back in, uncomment this:
## expanded_snp_idx_mat <- as.matrix(read.table(paste0(output_dir, "/samples/expanded_sample_mat.txt"),
##                                    header=F, sep="\t", quote="", as.is=T))
## expanded_snp_regions <- as.matrix(read.table(paste0(output_dir, "/samples/expanded_sample_regions.txt"),
##                                    header=F, sep="\t", quote="", as.is=T))
## expanded_snp_blocks <- as.matrix(read.table(paste0(output_dir, "/samples/expanded_sample_blocks.txt"),
##                                    header=F, sep="\t", quote="", as.is=T))

## expanded_input_idxs <- as.vector(read.table(paste0(output_dir, "/samples/expanded_input_idxs.txt"),
##                                    header=F, sep="\t", quote="", as.is=T)$V1)
## expanded_input_regions <- as.vector(read.table(paste0(output_dir, "/samples/expanded_input_regions.txt"),
##                                    header=F, sep="\t", quote="", as.is=T)$V1)
## expanded_input_blocks <- as.vector(read.table(paste0(output_dir, "/samples/expanded_input_blocks.txt"),
##                                    header=F, sep="\t", quote="", as.is=T)$V1)

cat("Auto LD expansion took:\n")
cat((proc.time() - start_time)[["elapsed"]], 'seconds\n')

## -----------------------------------------------------------------------------
## 6. Analysis plots for LD expansion
## -----------------------------------------------------------------------------
## look at the number of unique SNPs in each sampled set
uniq_snp_counts <- apply(expanded_snp_idx_mat, 2, function(x) {length(unique(x))})
uniq_snp_counts <- data.frame(count=uniq_snp_counts)

## also just check the number of unique columns (samples) that there are:
num_uniq_samples <- length(unique(apply(expanded_snp_idx_mat, 2, function(x) {sort(unique(x))})))
cat(num_uniq_samples, "different samples after LD expansion\n")

## also check how many unique ones there were in the input
input_ld_snp_f <- paste0(param_ref[['outdir']], '/ld_expansion/', param_ref[['outprefix']],
                        "_", r2_thresh, "_ld_cutoff_snps_within_", dist_thresh, ".txt")
## read in the LD expanded file
input_ld_snps <- read.table(input_ld_snp_f, header=T, sep="\t", quote="",
                            comment.char="", as.is=T)
num_uniq_input_ld_snps <- length(unique(input_ld_snps$rsID))

cat(num_uniq_input_ld_snps, "unique variants in INFERNO LD-expanded set\n")
cat(length(unique(expanded_input_idxs)), "unique variants in manual input expansion set\n")

## define some variables for the plotting range
min_scale <- round_any(min(min(uniq_snp_counts$count), num_uniq_input_ld_snps), accuracy=100, f=floor)
max_scale <- round_any(max(uniq_snp_counts$count), accuracy=100, f=ceiling)

make_graphic(paste0(output_dir, '/plots/', param_ref[['outprefix']], '_uniq_ld_snp_dists'))
print(ggplot(uniq_snp_counts, aes(x=count)) +
      geom_histogram(binwidth=1) +
      ## add a line for the sample count
      geom_vline(xintercept=num_uniq_input_ld_snps, color=muted("red"),
                 size=0.75, linetype=3) +
      scale_x_continuous(breaks=seq(min_scale, max_scale, by=10), limits=c(min_scale, max_scale)) +
      theme_bw() + xlab("Number of unique SNPs") + ylab("Number of samples") +
      ggtitle("Histogram of number of unique LD-expanded SNPs per sampled set") +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=25),
          axis.text.y = element_text(size=25), 
          title=element_text(size=30), plot.title = element_text(hjust = 0.5, size=30)))
dev.off()

## we want to see how many samples each SNP is found in for the LD expanded sets
ld_snp_samp_counts <- count(c(expanded_snp_idx_mat))
## then we count the frequency of these occurrences
ld_snp_dist_counts <- count(ld_snp_samp_counts$freq)

## define the plot breaks
x_breaks <- seq(0, num_samples, by=num_samples / 20)
y_breaks <- c(0, 100, 500, 1000, 2000, 3000, 4000, seq(5000, round_any(max(table(ld_snp_samp_counts$freq)), accuracy=5000, f=ceiling), by=5000))

make_graphic(paste0(output_dir, '/plots/', param_ref[['outprefix']], '_ld_expanded_snp_sample_num_hist'), width_ratio=1.5)
print(ggplot(ld_snp_dist_counts, aes(x=x, y=freq, colour=x)) +
    geom_bar(stat="identity", position="stack") +
    scale_colour_continuous(low="#132B43", high=muted("red")) +
    theme_bw() + xlab("Number of sampled data sets containing SNP") +
    ylab("Number of SNPs") +
    scale_x_continuous(breaks=x_breaks, limits=c(-1, tail(x_breaks, 1)+1), oob=rescale_none) +
    scale_y_continuous(breaks=y_breaks, labels=comma, trans="sqrt") +
    ggtitle("Histogram of number of samples containing each SNP after LD expansion") +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=25),
          axis.text.y = element_text(size=25), 
          title=element_text(size=30), plot.title = element_text(hjust = 0.5, size=30)))
dev.off()

## -----------------------------------------------------------------------------
## 7. Compute annotation overlap p-values
## -----------------------------------------------------------------------------
{
annot_start <- proc.time()

## read in all the different tissue classes
all_classes <- scan(paste0(ref_summary_dir, "/all_classes.txt"), "character", sep="\t")

## make logical matrices for the different annotations
enh_overlap_mat <- matrix(FALSE, nrow(snp_info), length(all_classes),
                          dimnames=list(snp_info$rsID, all_classes))
hmm_overlap_mat <- matrix(FALSE, nrow(snp_info), length(all_classes),
                          dimnames=list(snp_info$rsID, all_classes))

## use the 1000bp locus window for enhancers
enh_locus_overlaps <- read.table(paste0(ref_summary_dir, "/uniq_class_enh_locus_snps.txt"),
                                 header=F, sep="\t", quote="", comment.char="",
                                 as.is=T, col.names=c("chr", "rsID", "pos", "tissue_class"))
## add an indexing column so we don't have to match on row names
enh_locus_overlaps$snp_idx <- match(enh_locus_overlaps$rsID, snp_info$rsID)
## also add a tissue indexing column
enh_locus_overlaps$tiss_idx <- match(enh_locus_overlaps$tissue_class, colnames(enh_overlap_mat))
## now fill in the entries with overlap
enh_overlap_mat[cbind(enh_locus_overlaps$snp_idx, enh_locus_overlaps$tiss_idx)] <- TRUE
rm(enh_locus_overlaps)

## read in roadmap HMM data
hmm_overlaps <- read.table(paste0(ref_summary_dir, "/uniq_class_roadmap_hmm_snps.txt"),
                            header=F, sep="\t", quote="", comment.char="", as.is=T,
                            col.names=c("chr", "rsID", "pos", "tissue_class"))
## add an indexing column
hmm_overlaps$snp_idx <- match(hmm_overlaps$rsID, snp_info$rsID)
## also add a tissue indexing column
hmm_overlaps$tiss_idx <- match(hmm_overlaps$tissue_class, colnames(hmm_overlap_mat))
## now fill in the entries with overlap
hmm_overlap_mat[cbind(hmm_overlaps$snp_idx, hmm_overlaps$tiss_idx)] <- TRUE
rm(hmm_overlaps)

## now read in the annotation overlap counts for the input
## this is NOT collapsed by LD blocks
input_annotation_counts <- read.table(paste0(param_ref[['outdir']],
                                             '/summaries/tissue_class_annotation_counts_',
                                             r2_thresh, '_ld_', dist_thresh, '_dist.txt'),
                                      header=T, sep="\t", quote="", comment.char="", as.is=T)
## parse this down to only include the enhancer annotations for this analysis
input_annotation_counts <- input_annotation_counts[grep("eqtl", input_annotation_counts$annotation, invert=T),]

annot_names <- unique(input_annotation_counts$annotation)
## now we want to convert this into a count matrix so that we can easily compare
## rows for each tissue category, columns for each annotation and combo:
input_annot_mat <- acast(input_annotation_counts, tissue_class ~ annotation, value.var="count")

## in addition to full input counts, also read in the counts split by tag region
region_annotation_counts <- read.table(paste0(param_ref[['outdir']],
                                             '/summaries/tag_region_tissue_class_annotation_counts_',
                                             r2_thresh, '_ld_', dist_thresh, '_dist.txt'),
                                      header=T, sep="\t", quote="", comment.char="", as.is=T)
## also parse this one down
region_annotation_counts <- region_annotation_counts[grep("eqtl", region_annotation_counts$annotation, invert=T),]

## make a 3d array for each tag region, tissue class, and annotation
input_region_annot_arr <- acast(region_annotation_counts, tissue_class ~ annotation ~ tag_region, value.var="count")

## make a matrix to count the number of sampled sets that exceed the counts for each tissue and
## annotation
bootstrap_count_mat <- matrix(data=0, nrow=nrow(input_annot_mat), ncol=ncol(input_annot_mat),
                              dimnames=dimnames(input_annot_mat))
## make a 3d array for the region-specific counting
bootstrap_region_count_arr <- array(data=0, dim=c(nrow(input_region_annot_arr),
                                                  ncol(input_region_annot_arr), num_regions),
                                    dimnames=dimnames(input_region_annot_arr))

## LD-collapsed bootstrap analysis
## make bootstrap matrices / arrays for the LD-collapsed counts
## at this point, we haven't actually found the input variant overlaps yet
collapsed_bootstrap_count_mat <- matrix(data=0, nrow=nrow(input_annot_mat), ncol=ncol(input_annot_mat),
                                        dimnames=dimnames(input_annot_mat))
## make a 3d array for the region-specific counting
collapsed_bootstrap_region_count_arr <- array(data=0, dim=c(nrow(input_region_annot_arr),
                                                  ncol(input_region_annot_arr), num_regions),
                                              dimnames=dimnames(input_region_annot_arr))

## make a large data frame to save the results of each sample count
sample_counts <- data.frame(sample_num=rep(1:num_samples, each=length(all_classes)*length(annot_names)),
                            tissue_class=rep(rownames(input_annot_mat), times=num_samples*length(annot_names)),
                            annotation=rep(colnames(input_annot_mat), each=length(all_classes), times=num_samples),
                            count=rep(0, times=num_samples*length(annot_names)*length(all_classes)))

cat("Reading in annotation took:\n")
cat((proc.time() - annot_start)[['elapsed']], 'seconds\n')
}

## now we need to go through each sampled dataset, grab the annotation overlaps, and
## compare them to the observed count matrices/arrays
start_time <- proc.time()
## first do the annotation overlap on the input variants
## note that these input matrices are both for the LD collapsed analysis
manual_input_annot_mat <- matrix(data=0, nrow=length(all_classes), ncol=length(annot_names),
                                 dimnames=dimnames(input_annot_mat))
## make a 3d array for the tag region analysis
manual_input_region_annot_arr <- array(data=0, dim=c(nrow(input_region_annot_arr),
                                                   ncol(input_region_annot_arr), num_regions),
                                       dimnames=dimnames(input_region_annot_arr))

## we only want to use the unique expanded variants to avoid double counting, so define a
## vector to use for this since we use it to index the blocks and regions too
uniq_input_snps <- !duplicated(expanded_input_idxs)
uniq_input_snp_and_regions <- !duplicated(cbind(expanded_input_idxs, expanded_input_regions))

## use the sample indices to directly index the overlap matrices and convert to numeric by
## multiplying them by 1
input_enh_tiss_overlaps <- enh_overlap_mat[expanded_input_idxs[uniq_input_snps],] * 1
input_hmm_tiss_overlaps <- hmm_overlap_mat[expanded_input_idxs[uniq_input_snps],] * 1
## also get the combination overlaps
input_enh_hmm_tiss_overlaps <- (input_enh_tiss_overlaps & input_hmm_tiss_overlaps) * 1
## now, collapse these counts by LD block!!
## multiple hits in a single LD block only count once (per tag region)
input_ld_collapsed_enh_counts <- aggregate(input_enh_tiss_overlaps,
                                           by=list(expanded_input_blocks[uniq_input_snps],
                                               expanded_input_regions[uniq_input_snps]),
                                           function(x) {ifelse(sum(x) > 0, 1, 0)})
input_ld_collapsed_hmm_counts <- aggregate(input_hmm_tiss_overlaps,
                                           by=list(expanded_input_blocks[uniq_input_snps],
                                               expanded_input_regions[uniq_input_snps]),
                                           function(x) {ifelse(sum(x) > 0, 1, 0)})
input_ld_collapsed_enh_hmm_counts <- aggregate(input_enh_hmm_tiss_overlaps,
                                               by=list(expanded_input_blocks[uniq_input_snps],
                                                   expanded_input_regions[uniq_input_snps]),
                                               function(x) {ifelse(sum(x) > 0, 1, 0)})

## for convenience, create matrices that don't have the grouping column
input_enh_counts <- subset(input_ld_collapsed_enh_counts, select=-c(Group.1, Group.2))
input_hmm_counts <- subset(input_ld_collapsed_hmm_counts, select=-c(Group.1, Group.2))
input_enh_hmm_counts <- subset(input_ld_collapsed_enh_hmm_counts, select=-c(Group.1, Group.2))

## now count the annotation overlaps and fill in the columns
manual_input_annot_mat[,"locus_enh"] <- colSums(input_enh_counts)
manual_input_annot_mat[,"roadmap_hmm_enh"] <- colSums(input_hmm_counts)
manual_input_annot_mat[,"locus_enh+roadmap_hmm_enh"] <- colSums(input_enh_hmm_counts)

## we also need to split the manually expanded input by tag regions
## need to recalculate the annotation overlaps because now we allow the same SNP to count towards multiple regions
input_region_enh_tiss_overlaps <- enh_overlap_mat[expanded_input_idxs[uniq_input_snp_and_regions],] * 1
input_region_hmm_tiss_overlaps <- hmm_overlap_mat[expanded_input_idxs[uniq_input_snp_and_regions],] * 1
## also get the combination overlaps
input_region_enh_hmm_tiss_overlaps <- (input_region_enh_tiss_overlaps & input_region_hmm_tiss_overlaps) * 1

## now, collapse these counts by LD block!!
## multiple hits in a single LD block only count once (per tag region)
input_region_ld_collapsed_enh_counts <- aggregate(input_region_enh_tiss_overlaps,
                                                  by=list(expanded_input_blocks[uniq_input_snp_and_regions],
                                                      expanded_input_regions[uniq_input_snp_and_regions]),
                                                  function(x) {ifelse(sum(x) > 0, 1, 0)})
input_region_ld_collapsed_hmm_counts <- aggregate(input_region_hmm_tiss_overlaps,
                                                  by=list(expanded_input_blocks[uniq_input_snp_and_regions],
                                                      expanded_input_regions[uniq_input_snp_and_regions]),
                                                  function(x) {ifelse(sum(x) > 0, 1, 0)})
input_region_ld_collapsed_enh_hmm_counts <- aggregate(input_region_enh_hmm_tiss_overlaps,
                                                      by=list(expanded_input_blocks[uniq_input_snp_and_regions],
                                                          expanded_input_regions[uniq_input_snp_and_regions]),
                                                      function(x) {ifelse(sum(x) > 0, 1, 0)})

## for convenience, create matrices that don't have the grouping column
input_region_enh_counts <- subset(input_region_ld_collapsed_enh_counts, select=-c(Group.1, Group.2))
input_region_hmm_counts <- subset(input_region_ld_collapsed_hmm_counts, select=-c(Group.1, Group.2))
input_region_enh_hmm_counts <- subset(input_region_ld_collapsed_enh_hmm_counts, select=-c(Group.1, Group.2))

## to get the right regions, just use the second grouping variable of the region specific data
manual_input_region_annot_arr[,"locus_enh",] <- t(rowsum(input_region_enh_counts, input_region_ld_collapsed_enh_counts$Group.2))
manual_input_region_annot_arr[,"roadmap_hmm_enh",] <- t(rowsum(input_region_hmm_counts, input_region_ld_collapsed_hmm_counts$Group.2))
manual_input_region_annot_arr[,"locus_enh+roadmap_hmm_enh",] <- t(rowsum(input_region_enh_hmm_counts, input_region_ld_collapsed_enh_hmm_counts$Group.2))

for(s in seq(num_samples)) {
    ## pull out the sampled snps
    this_samp <- expanded_snp_idx_mat[,s]
    ## also get the regions associated with each
    this_regions <- expanded_snp_regions[,s]
    ## for the LD collapsed analysis, get indices for the LD blocks
    this_blocks <- expanded_snp_blocks[,s]
    
    ## make an annotation count matrix
    this_annot_mat <- matrix(data=0, nrow=length(all_classes), ncol=length(annot_names),
                             dimnames=dimnames(input_annot_mat))
    ## make a 3d array for the tag region analysis
    this_region_annot_arr <- array(data=0, dim=c(nrow(input_region_annot_arr),
                                              ncol(input_region_annot_arr), num_regions),
                                    dimnames=dimnames(input_region_annot_arr))

    ## same as for the input, we only want to use unique variants
    uniq_sample_snps <- !duplicated(this_samp)
    ## we also want to get unique ones per region for the split region counts
    uniq_sample_snp_and_regions <- !duplicated(cbind(this_samp, this_regions))
    
    ## use the sample indices to directly index the overlap matrices and convert to numeric by
    ## multiplying them by 1
    enh_tiss_overlaps <- enh_overlap_mat[this_samp[uniq_sample_snps],] * 1
    hmm_tiss_overlaps <- hmm_overlap_mat[this_samp[uniq_sample_snps],] * 1
    ## also get the combination overlaps
    enh_hmm_tiss_overlaps <- (enh_tiss_overlaps & hmm_tiss_overlaps) * 1

    ## now count the annotation overlaps and fill in the columns
    this_annot_mat[,"locus_enh"] <- colSums(enh_tiss_overlaps)
    this_annot_mat[,"roadmap_hmm_enh"] <- colSums(hmm_tiss_overlaps)
    this_annot_mat[,"locus_enh+roadmap_hmm_enh"] <- colSums(enh_hmm_tiss_overlaps)

    ## we also need to split this by tag regions
    ## because we have a different set of unique SNPs by region, recompute overlaps
    enh_region_tiss_overlaps <- enh_overlap_mat[this_samp[uniq_sample_snp_and_regions],] * 1
    hmm_region_tiss_overlaps <- hmm_overlap_mat[this_samp[uniq_sample_snp_and_regions],] * 1
    ## also get the combination overlaps
    enh_hmm_region_tiss_overlaps <- (enh_region_tiss_overlaps & hmm_region_tiss_overlaps) * 1
    
    this_region_annot_arr[,"locus_enh",] <- t(rowsum(enh_region_tiss_overlaps,
                                                     this_regions[uniq_sample_snp_and_regions]))
    this_region_annot_arr[,"roadmap_hmm_enh",] <- t(rowsum(hmm_region_tiss_overlaps,
                                                           this_regions[uniq_sample_snp_and_regions]))
    this_region_annot_arr[,"locus_enh+roadmap_hmm_enh",] <- t(rowsum(enh_hmm_region_tiss_overlaps,
                                                                     this_regions[uniq_sample_snp_and_regions]))
    
    ## now compare to the input matrix and increment the counts: which annotation-tissue class
    ## combinations are observed the same or more times in this sample as compared to the
    ## input?
    high_count_idx <- which(this_annot_mat >= input_annot_mat)
    bootstrap_count_mat[high_count_idx] <- bootstrap_count_mat[high_count_idx]+1

    ## same procedure for the tag region split analysis
    region_count_idx <- which(this_region_annot_arr >= input_region_annot_arr)
    bootstrap_region_count_arr[region_count_idx] <- bootstrap_region_count_arr[region_count_idx]+1
    
    ## store these counts in the sample count matrix (this relies on the annotation matrix
    ## columns and rows being in the same order as the input matrix, which was used to generate
    ## the sample_counts DF)
    sample_counts$count[sample_counts$sample_num==s] <- melt(this_annot_mat)$value
    ## sample_region_counts$count[sample_region_counts$sample_num==s] <- melt(this_region_annot_arr)$value

    ## also do the LD collapsed analysis
    ## also make count matrices for LD collapsed counting
    this_ld_collapsed_annot_mat <- matrix(data=0, nrow=length(all_classes),
                                          ncol=length(annot_names),
                                          dimnames=dimnames(input_annot_mat))
    ## make a 3d array for the tag region analysis
    this_ld_collapsed_region_annot_arr <- array(data=0, dim=c(nrow(input_region_annot_arr),
                                                  ncol(input_region_annot_arr), num_regions),
                                                dimnames=dimnames(input_region_annot_arr))
    
    ## multiple hits in a single LD block only count once (per region), and only unique ones
    this_ld_collapsed_enh_counts <- aggregate(enh_tiss_overlaps,
                                              by=list(this_blocks[uniq_sample_snps],
                                                  this_regions[uniq_sample_snps]),
                                              function(x) {ifelse(sum(x) > 0, 1, 0)})
    this_ld_collapsed_hmm_counts <- aggregate(hmm_tiss_overlaps,
                                              by=list(this_blocks[uniq_sample_snps],
                                                  this_regions[uniq_sample_snps]),
                                              function(x) {ifelse(sum(x) > 0, 1, 0)})
    this_ld_collapsed_enh_hmm_counts <- aggregate(enh_hmm_tiss_overlaps,
                                                  by=list(this_blocks[uniq_sample_snps],
                                                      this_regions[uniq_sample_snps]),
                                                  function(x) {ifelse(sum(x) > 0, 1, 0)})

    ## for convenience, create matrices that don't have the grouping column
    this_enh_counts <- subset(this_ld_collapsed_enh_counts, select=-c(Group.1, Group.2))
    this_hmm_counts <- subset(this_ld_collapsed_hmm_counts, select=-c(Group.1, Group.2))
    this_enh_hmm_counts <- subset(this_ld_collapsed_enh_hmm_counts, select=-c(Group.1, Group.2))

    ## now count the annotation overlaps and fill in the columns
    this_ld_collapsed_annot_mat[,"locus_enh"] <- colSums(this_enh_counts)
    this_ld_collapsed_annot_mat[,"roadmap_hmm_enh"] <- colSums(this_hmm_counts)
    this_ld_collapsed_annot_mat[,"locus_enh+roadmap_hmm_enh"] <- colSums(this_enh_hmm_counts)

    ## do the same thing for the split tag regions
    this_region_ld_collapsed_enh_counts <- aggregate(enh_region_tiss_overlaps,
                                              by=list(this_blocks[uniq_sample_snp_and_regions],
                                                  this_regions[uniq_sample_snp_and_regions]),
                                              function(x) {ifelse(sum(x) > 0, 1, 0)})
    this_region_ld_collapsed_hmm_counts <- aggregate(hmm_region_tiss_overlaps,
                                              by=list(this_blocks[uniq_sample_snp_and_regions],
                                                  this_regions[uniq_sample_snp_and_regions]),
                                              function(x) {ifelse(sum(x) > 0, 1, 0)})
    this_region_ld_collapsed_enh_hmm_counts <- aggregate(enh_hmm_region_tiss_overlaps,
                                                  by=list(this_blocks[uniq_sample_snp_and_regions],
                                                      this_regions[uniq_sample_snp_and_regions]),
                                                  function(x) {ifelse(sum(x) > 0, 1, 0)})

    ## for convenience, create matrices that don't have the grouping column
    this_region_enh_counts <- subset(this_region_ld_collapsed_enh_counts, select=-c(Group.1, Group.2))
    this_region_hmm_counts <- subset(this_region_ld_collapsed_hmm_counts, select=-c(Group.1, Group.2))
    this_region_enh_hmm_counts <- subset(this_region_ld_collapsed_enh_hmm_counts, select=-c(Group.1, Group.2))
    
    this_ld_collapsed_region_annot_arr[,"locus_enh",] <- t(rowsum(this_region_enh_counts, this_region_ld_collapsed_enh_counts$Group.2))
    this_ld_collapsed_region_annot_arr[,"roadmap_hmm_enh",] <- t(rowsum(this_region_hmm_counts, this_region_ld_collapsed_enh_counts$Group.2))
    this_ld_collapsed_region_annot_arr[,"locus_enh+roadmap_hmm_enh",] <- t(rowsum(this_region_enh_hmm_counts, this_region_ld_collapsed_enh_counts$Group.2))

    ## now compare this to the input LD collapsed matrix
    high_collapsed_count_idx <- which(this_ld_collapsed_annot_mat >= manual_input_annot_mat)
    collapsed_bootstrap_count_mat[high_collapsed_count_idx] <-
        collapsed_bootstrap_count_mat[high_collapsed_count_idx]+1

    ## same procedure for the tag region split analysis
    collapsed_region_count_idx <- which(this_ld_collapsed_region_annot_arr >= manual_input_region_annot_arr)
    collapsed_bootstrap_region_count_arr[collapsed_region_count_idx] <-
        collapsed_bootstrap_region_count_arr[collapsed_region_count_idx]+1    
}
cat("Bootstrap counting for p-value calculation took:\n")
cat((proc.time() - start_time)[['elapsed']], 'seconds\n')

## remove the annotation matrices to save memory
rm(enh_overlap_mat, hmm_overlap_mat)

## ------------------------------
## analysis on fully expanded sets (no LD collapsing)
## --------------
## all region analysis
## generate the (uncorrected) p value matrix
## we add one to all the values to represent the input sample
empirical_pvals <- (bootstrap_count_mat+1) / (num_samples+1)
dimnames(empirical_pvals) <- dimnames(input_annot_mat)
cat(sum(empirical_pvals < 0.05), 'unadjusted tests were significant for combined region analysis\n')

## now correct these for multiple testing
## correct all the p-values using Benjamini-Hochberg correction together
bh_adj_pvals <- matrix(p.adjust(empirical_pvals, "BH"), nrow=nrow(empirical_pvals),
                       ncol=ncol(empirical_pvals), dimnames=dimnames(empirical_pvals))

dimnames(bh_adj_pvals) <- dimnames(empirical_pvals)

cat(sum(bh_adj_pvals < 0.05), "adjusted tests were significant for combined tag region\n")

## --------------
## split tag region analysis
## generate the (uncorrected) p value array
## we add one to all the values to represent the input sample
split_empirical_pvals <- (bootstrap_region_count_arr+1) / (num_samples+1)
dimnames(split_empirical_pvals) <- dimnames(input_region_annot_arr)
cat(sum(split_empirical_pvals < 0.05), 'unadjusted tests were significant for split tag region analysis\n')

## do this split by region
for(i in seq(dim(split_empirical_pvals)[3])) {
    cat("Significant unadjusted tests for", dimnames(split_empirical_pvals)[[3]][i], "\n")
    cat(sum(split_empirical_pvals[,,i] < 0.05), "\n")
}

## use BH correction within each tag region
split_bh_adj_pvals <- aperm(aaply(split_empirical_pvals, 3, function(tag_mat) {
    corr_mat <- matrix(p.adjust(tag_mat, method="BH"), nrow=nrow(tag_mat), ncol=ncol(tag_mat),
                       dimnames=dimnames(tag_mat))
    return(corr_mat)
}), perm=c(2, 3, 1))

cat(sum(split_bh_adj_pvals < 0.05), "adjusted tests were significant for split tag region analysis\n")

## do this split by region
for(i in seq(dim(split_bh_adj_pvals)[3])) {
    cat("Significant adjusted tests for", dimnames(split_bh_adj_pvals)[[3]][i], "\n")
    cat(sum(split_bh_adj_pvals[,,i] < 0.05), "\n")
}

## ------------------------------
## analysis on LD collapsed sets
## --------------
## all region analysis
## generate the (uncorrected) p value matrix
## we add one to all the values to represent the input sample
collapsed_empirical_pvals <- (collapsed_bootstrap_count_mat+1) / (num_samples+1)
dimnames(collapsed_empirical_pvals) <- dimnames(input_annot_mat)
cat(sum(collapsed_empirical_pvals < 0.05), 'unadjusted tests were significant for LD collapsed combined region analysis\n')

## now correct these for multiple testing
## correct all the p-values using Benjamini-Hochberg correction together
collapsed_bh_adj_pvals <- matrix(p.adjust(collapsed_empirical_pvals, "BH"),
                                 nrow=nrow(collapsed_empirical_pvals),
                                 ncol=ncol(collapsed_empirical_pvals),
                                 dimnames=dimnames(collapsed_empirical_pvals))

dimnames(collapsed_bh_adj_pvals) <- dimnames(collapsed_empirical_pvals)

cat(sum(collapsed_bh_adj_pvals < 0.05), "adjusted tests were significant for LD collapsed combined tag region analysis\n")

## --------------
## split tag region analysis
## generate the (uncorrected) p value array
## we add one to all the values to represent the input sample
collapsed_split_empirical_pvals <- (collapsed_bootstrap_region_count_arr+1) / (num_samples+1)
dimnames(collapsed_split_empirical_pvals) <- dimnames(input_region_annot_arr)
cat(sum(collapsed_split_empirical_pvals < 0.05), 'unadjusted tests were significant for LD collapsed split tag region analysis\n')

## do this split by region
for(i in seq(dim(collapsed_split_empirical_pvals)[3])) {
    cat("Significant unadjusted tests for", dimnames(collapsed_split_empirical_pvals)[[3]][i], "\n")
    cat(sum(collapsed_split_empirical_pvals[,,i] < 0.05), "\n")
}

## use BH correction within each tag region
collapsed_split_bh_adj_pvals <- aperm(aaply(collapsed_split_empirical_pvals, 3,
    function(tag_mat) {
    corr_mat <- matrix(p.adjust(tag_mat, method="BH"), nrow=nrow(tag_mat), ncol=ncol(tag_mat),
                       dimnames=dimnames(tag_mat))
    return(corr_mat)
}), perm=c(2, 3, 1))

cat(sum(collapsed_split_bh_adj_pvals < 0.05), " adjusted tests were significant for split tag region analysis\n")

## do this split by region
for(i in seq(dim(collapsed_split_bh_adj_pvals)[3])) {
    cat("Significant adjusted tests for", dimnames(collapsed_split_bh_adj_pvals)[[3]][i], "\n")
    cat(sum(collapsed_split_bh_adj_pvals[,,i] < 0.05), "\n")
}

## -----------------------------------------------------------------------------
## 8. Analysis plots for annotation overlap and p-values of non-LD collapsed analysis
## -----------------------------------------------------------------------------
## -------------------
## combined region analysis
## first, a boxplot of the p-value distributions. to do this, we need to melt the p-value matrices
pval_df <- rbind(cbind(pval="Raw P-value", melt(empirical_pvals, varnames=c("tissue_class", "annotation"))),
                 cbind(pval="BH-adjusted P-value", melt(bh_adj_pvals, varnames=c("tissue_class", "annotation"))))
pval_df$annotation <- factor(pval_df$annotation, ordered=T,
                             levels=c("locus_enh", "roadmap_hmm_enh", "locus_enh+roadmap_hmm_enh"),
                             labels=c("FANTOM5 Enhancer", "Roadmap HMM Enhancer", "FANTOM5 Enh+Roadmap HMM Enh"))
pval_df$tissue_class <- factor(pval_df$tissue_class, ordered=T,
                               levels=sort(unique(pval_df$tissue_class), dec=T))

write.table(pval_df, paste0(output_dir, '/tables/all_region_bootstrap_results.txt'), quote=F, sep="\t", row.names=F, col.names=T)
## ## to read in these results:
## pval_df <- read.table(paste0(output_dir, '/tables/all_region_bootstrap_results.txt'), header=T, sep="\t", quote="", as.is=T)
## pval_df$annotation <- factor(pval_df$annotation, ordered=T,
##                              levels=c("GTEx eQTL", "FANTOM5 Enhancer", "Roadmap HMM Enhancer",
##                                  "FANTOM5 Enh+GTEx eQTL", "GTEx eQTL+Roadmap HMM Enh", "FANTOM5 Enh+Roadmap HMM Enh", "FANTOM5 Enh+GTEx eQTL+Roadmap HMM Enh"))
## pval_df$tissue_class <- factor(pval_df$tissue_class, ordered=T,
##                                levels=sort(unique(pval_df$tissue_class), dec=T))

## -------------------
## make the boxplot
make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_annotation_overlap_pval_boxplots"))
print(ggplot(pval_df, aes(x=annotation, y=value, fill=pval)) +
      geom_boxplot() + scale_y_continuous(breaks=seq(0, 1, by=0.05)) +
      scale_fill_discrete(h=c(100, 200)) +
      theme_bw() + xlab("Annotation") +
      ylab("P-value") + ggtitle("P-value distributions for each annotation") +
      theme(legend.position="right",
            legend.text=element_text(size=20), legend.title=element_text(size=20),
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(vjust=1, hjust=1, size=15),
            axis.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## -------------------
## also make a histogram of the p-values across all annotations
make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_raw_pval_hist"))
print(ggplot(pval_df[pval_df$pval=="Raw P-value",], aes(x=value)) +
      geom_histogram(binwidth=0.01) +
      scale_x_continuous(breaks=seq(0, 1, by=0.05)) +
      theme_bw() + xlab("P-value") +
      ylab("Number of tissue-annotation combinations") + ggtitle("Raw P-value histogram") +
      theme(legend.position="none",
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(vjust=1, hjust=1, size=15),            
            axis.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_bh_pval_hist"))
print(ggplot(pval_df[pval_df$pval=="BH-adjusted P-value",], aes(x=value)) +
      geom_histogram(binwidth=0.01) +
      scale_x_continuous(breaks=seq(0, 1, by=0.05)) +
      theme_bw() + xlab("P-value") +
      ylab("Number of tissue-annotation combinations") + ggtitle("BH-adjusted P-value histogram") +     
      theme(legend.position="none",
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(vjust=1, hjust=1, size=15),
            title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## make histograms split by annotation type
pval_df$annot_wrap <- unlist(lapply(strwrap(pval_df$annotation, width=15, simplify=F),
                                    paste, collapse="\n"))
sorted_wrap_levels <- unlist(lapply(strwrap(levels(pval_df$annotation), width=15, simplify=F),
                                    paste, collapse="\n"))
pval_df$annot_wrap <- factor(pval_df$annot_wrap, ordered=T,
                             levels=sorted_wrap_levels)

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_raw_pval_split_annot_hist"))
print(ggplot(pval_df[pval_df$pval=="Raw P-value",], aes(x=value, fill=annotation)) +
      geom_histogram(binwidth=0.01) +
      scale_x_continuous(breaks=seq(0, 1, by=0.05)) +
      facet_grid(annot_wrap ~ ., scales="free") +
      theme_bw() + xlab("P-value") +
      ylab("Number of tissue-annotation combinations") + ggtitle("Raw P-value histogram") +
      theme(legend.position="none",
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(vjust=1, hjust=1, size=15),
            strip.text = element_text(size=20), 
            axis.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_bh_pval_split_annot_hist"))
print(ggplot(pval_df[pval_df$pval=="BH-adjusted P-value",], aes(x=value, fill=annotation)) +
      geom_histogram(binwidth=0.01) +
      scale_x_continuous(breaks=seq(0, 1, by=0.05)) +
      facet_grid(annot_wrap ~ ., scales="free") +
      theme_bw() + xlab("P-value") +
      ylab("Number of tissue-annotation combinations") + ggtitle("BH-adjusted P-value histogram") +
      theme(legend.position="none",
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(angle=45, vjust=1, hjust=1, size=15),
            strip.text = element_text(size=20),             
            axis.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## -------------------
## make a heatmap of the p-values
## define the breaks
hm_breaks <- c(0, 0.05, 1)

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_annotation_by_tissue_raw_pval_heatmap"))
print(ggplot(pval_df[pval_df$pval=="Raw P-value",],
             aes(x=annotation, y=tissue_class)) +
      geom_tile(aes(fill=value), colour="#000000") +
      ## use 4 colors here so that all significant hits are fairly red
      scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                           values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                           limits=c(0, 1), name="P-value", guide="colorbar") +
      theme_bw() + ylab("Tissue category") +
      xlab("Annotation") + ggtitle("Heatmap of p-values") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## also make one with text in the boxes for the p-values
make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_annotation_by_tissue_raw_pval_heatmap_with_text"))
print(ggplot(pval_df[pval_df$pval=="Raw P-value",],
             aes(x=annotation, y=tissue_class)) +
      geom_tile(aes(fill=value), colour="#000000") +
      scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                           values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                           limits=c(0, 1), name="P-value", guide="colorbar") +
      theme_bw() + ylab("Tissue category") +
      xlab("Annotation") + ggtitle("Heatmap of p-values") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)) +
      geom_text(aes(x=annotation, y=tissue_class,
                    label=format(round(value, digits=4), scientific=F)), color="black", size=4))
dev.off()

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_annotation_by_tissue_bh_pval_heatmap"))
print(ggplot(pval_df[pval_df$pval=="BH-adjusted P-value",],
             aes(x=annotation, y=tissue_class)) +
      geom_tile(aes(fill=value), colour="#000000") +
      ## use 4 colors here so that all significant hits are fairly red
      scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                           values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                           limits=c(0, 1), name="BH-adjusted p-value", guide="colorbar") +
      theme_bw() + ylab("Tissue category") +
      xlab("Annotation") + ggtitle("Heatmap of adjusted p-values") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## also make one with text in the boxes for the p-values
make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_annotation_by_tissue_bh_pval_heatmap_with_text"))
print(ggplot(pval_df[pval_df$pval=="BH-adjusted P-value",],
             aes(x=annotation, y=tissue_class)) +
      geom_tile(aes(fill=value), colour="#000000") +
      scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                           values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                           limits=c(0, 1), name="BH-adjusted p-value", guide="colorbar") +
      theme_bw() + ylab("Tissue category") +
      xlab("Annotation") + ggtitle("Heatmap of adjusted p-values") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)) +
      geom_text(aes(x=annotation, y=tissue_class,
                    label=format(round(value, digits=4), scientific=F)), color="black", size=4))
dev.off()

## -------------------
## split tag region analysis
## melt the results
split_pval_df <- rbind(cbind(pval="Raw P-value", melt(split_empirical_pvals, varnames=c("tissue_class", "annotation", "tag_region"))),
                 cbind(pval="BH-adjusted P-value", melt(split_bh_adj_pvals, varnames=c("tissue_class", "annotation", "tag_region"))))
split_pval_df$annotation <- factor(split_pval_df$annotation, ordered=T,
                             levels=c("gtex_eqtl", "locus_enh", "roadmap_hmm_enh", "locus_enh+gtex_eqtl", "gtex_eqtl+roadmap_hmm_enh", "locus_enh+roadmap_hmm_enh", "locus_enh+gtex_eqtl+roadmap_hmm_enh"),
                             labels=c("GTEx eQTL", "FANTOM5 Enhancer", "Roadmap HMM Enhancer",
                                 "FANTOM5 Enh+GTEx eQTL", "GTEx eQTL+Roadmap HMM Enh", "FANTOM5 Enh+Roadmap HMM Enh", "FANTOM5 Enh+GTEx eQTL+Roadmap HMM Enh"))
split_pval_df$tissue_class <- factor(split_pval_df$tissue_class, ordered=T,
                               levels=sort(unique(split_pval_df$tissue_class), dec=T))

write.table(split_pval_df, paste0(output_dir, '/tables/split_region_bootstrap_results.txt'), quote=F, sep="\t", row.names=F, col.names=T)
## to read in this data
## split_pval_df <- read.table(paste0(output_dir, '/tables/split_region_bootstrap_results.txt'), header=T, sep="\t", quote="", as.is=T)
## split_pval_df$annotation <- factor(split_pval_df$annotation, ordered=T,
##                              levels=c("GTEx eQTL", "FANTOM5 Enhancer", "Roadmap HMM Enhancer",
##                                  "FANTOM5 Enh+GTEx eQTL", "GTEx eQTL+Roadmap HMM Enh", "FANTOM5 Enh+Roadmap HMM Enh", "FANTOM5 Enh+GTEx eQTL+Roadmap HMM Enh"))
## split_pval_df$tissue_class <- factor(split_pval_df$tissue_class, ordered=T,
##                                levels=sort(unique(split_pval_df$tissue_class), dec=T))

dir.create(paste0(output_dir, '/plots/split_tag_regions/'), F, T)
for(this_tag in unique(split_pval_df$tag_region)) {
    tag_out <- gsub("/", "_", strsplit(this_tag, ":")[[1]][1])
    dir.create(paste0(output_dir, "/plots/split_tag_regions/", tag_out), F, T)

    ## define the breaks
    hm_breaks <- c(0, 0.05, 1)

    ## grab the relevant data
    this_pval_df <- subset(split_pval_df, pval=="BH-adjusted P-value" & tag_region==this_tag)
    
    make_graphic(paste0(output_dir, "/plots/split_tag_regions/", tag_out, '/',
                        param_ref[['outprefix']], "_", tag_out, "_annotation_by_tissue_bh_pval_heatmap"))
    print(ggplot(this_pval_df, 
                 aes(x=annotation, y=tissue_class)) +
          geom_tile(aes(fill=value), colour="#000000") +
          ## use 4 colors here so that all significant hits are fairly red
          scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                               values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                               limits=c(0, 1), name="BH-adjusted p-value", guide="colorbar") +
          theme_bw() + ylab("Tissue category") +
          xlab("Annotation") + ggtitle(paste("Heatmap of adjusted p-values,", this_tag)) +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
    dev.off()

    ## also make one with text in the boxes for the p-values
    make_graphic(paste0(output_dir, "/plots/split_tag_regions/", tag_out, '/',
                        param_ref[['outprefix']], "_", tag_out, "_annotation_by_tissue_bh_pval_heatmap_with_text"))
    print(ggplot(this_pval_df,
                 aes(x=annotation, y=tissue_class)) +
          geom_tile(aes(fill=value), colour="#000000") +
          scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                               values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                               limits=c(0, 1), name="BH-adjusted p-value", guide="colorbar") +
          theme_bw() + ylab("Tissue category") +
          xlab("Annotation") + ggtitle(paste("Heatmap of adjusted p-values,", this_tag)) +
          theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
                axis.text.y = element_text(size=15), axis.title=element_text(size=20),
                legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
                legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)) +
          geom_text(aes(x=annotation, y=tissue_class,
                    label=format(round(value, digits=4), scientific=F)), color="black", size=4))
    dev.off()

    ## also make plots for the uncorrected p-values    
    ## grab the relevant data
    this_raw_pval_df <- subset(split_pval_df, pval=="Raw P-value" & tag_region==this_tag)
    
    make_graphic(paste0(output_dir, "/plots/split_tag_regions/", tag_out, '/',
                        param_ref[['outprefix']], "_", tag_out, "_annotation_by_tissue_raw_pval_heatmap"))
    print(ggplot(this_raw_pval_df, 
                 aes(x=annotation, y=tissue_class)) +
          geom_tile(aes(fill=value), colour="#000000") +
          ## use 4 colors here so that all significant hits are fairly red
          scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                               values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                               limits=c(0, 1), name="Raw p-value", guide="colorbar") +
          theme_bw() + ylab("Tissue category") +
          xlab("Annotation") + ggtitle(paste("Heatmap of unadjusted p-values,", this_tag)) +
          theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
                axis.text.y = element_text(size=15), axis.title=element_text(size=20),
                legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
                legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
    dev.off()

    ## also make one with text in the boxes for the p-values
    make_graphic(paste0(output_dir, "/plots/split_tag_regions/", tag_out, '/',
                        param_ref[['outprefix']], "_", tag_out, "_annotation_by_tissue_raw_pval_heatmap_with_text"))
    print(ggplot(this_raw_pval_df,
                 aes(x=annotation, y=tissue_class)) +
          geom_tile(aes(fill=value), colour="#000000") +
          scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                               values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                               limits=c(0, 1), name="Raw p-value", guide="colorbar") +
          theme_bw() + ylab("Tissue category") +
          xlab("Annotation") + ggtitle(paste("Heatmap of unadjusted p-values,", this_tag)) +
          theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
                axis.text.y = element_text(size=15), axis.title=element_text(size=20),
                legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
                legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)) +
          geom_text(aes(x=annotation, y=tissue_class,
                        label=format(round(value, digits=4), scientific=F)), color="black", size=4))
    dev.off()
}

## -------------------
## make histograms of the sample count distributions
## define the color scheme for the tissues (from run_full_analysis.R)
category_colors <- c(rgb(227,191,35, maxColorValue=255), rgb(224,98,247, maxColorValue=255), rgb(36,153,165, maxColorValue=255), rgb(232,2,58, maxColorValue=255), rgb(30,181,54, maxColorValue=255), rgb(253,164,185, maxColorValue=255), rgb(44,115,232, maxColorValue=255), rgb(83,118,43, maxColorValue=255), rgb(203,50,133, maxColorValue=255), rgb(194,196,252, maxColorValue=255), rgb(157,236,192, maxColorValue=255), rgb(141,86,174, maxColorValue=255), rgb(253,133,114, maxColorValue=255), rgb(175,242,253, maxColorValue=255), rgb(149,87,122, maxColorValue=255), rgb(131,233,24, maxColorValue=255), rgb(182,241,140, maxColorValue=255), rgb(252,50,7, maxColorValue=255), rgb(244,148,221, maxColorValue=255), rgb(40,169,123, maxColorValue=255), rgb(247,183,144, maxColorValue=255), rgb(242,184,94, maxColorValue=255), rgb(53,108,145, maxColorValue=255), rgb(198,8,78, maxColorValue=255), rgb(61,115,80, maxColorValue=255), rgb(41,232,215, maxColorValue=255), rgb(122,107,24, maxColorValue=255), rgb(79,153,241, maxColorValue=255), rgb(100,130,128, maxColorValue=255), rgb(166,169,37, maxColorValue=255), rgb(203,137,237, maxColorValue=255), rgb(178,204,231, maxColorValue=255))
names(category_colors) <- all_classes

cat_col_scale <- scale_fill_manual(name="Tissue Category", values=category_colors)

## first replace all the annotation labels
sample_counts$annotation <- ifelse(sample_counts$annotation=="gtex_eqtl", "eqtl",
                            ifelse(sample_counts$annotation=="gtex_eqtl+roadmap_hmm_enh", "eqtl+hmm",
                            ifelse(sample_counts$annotation=="locus_enh", "enh",
                            ifelse(sample_counts$annotation=="locus_enh+gtex_eqtl", "enh+eqtl",
                            ifelse(sample_counts$annotation=="locus_enh+gtex_eqtl+roadmap_hmm_enh", "enh+eqtl+hmm",
                            ifelse(sample_counts$annotation=="locus_enh+roadmap_hmm_enh", "enh+hmm",
                            "hmm"))))))
sample_counts$annotation <- factor(sample_counts$annotation, ordered=T,
                                   levels=c("eqtl", "enh", "hmm", "eqtl+hmm", "enh+eqtl", "enh+hmm", "enh+eqtl+hmm"))

## convert the input into a similar form
input_annot_counts <- melt(input_annot_mat, varnames=c("tissue_class", "annotation"), value.name="count")
input_annot_counts$annotation <- ifelse(input_annot_counts$annotation=="gtex_eqtl", "eqtl",
                            ifelse(input_annot_counts$annotation=="gtex_eqtl+roadmap_hmm_enh", "eqtl+hmm",
                            ifelse(input_annot_counts$annotation=="locus_enh", "enh",
                            ifelse(input_annot_counts$annotation=="locus_enh+gtex_eqtl", "enh+eqtl",
                            ifelse(input_annot_counts$annotation=="locus_enh+gtex_eqtl+roadmap_hmm_enh", "enh+eqtl+hmm",
                            ifelse(input_annot_counts$annotation=="locus_enh+roadmap_hmm_enh", "enh+hmm",
                            "hmm"))))))
input_annot_counts$annotation <- factor(input_annot_counts$annotation, ordered=T,
                                   levels=c("eqtl", "enh", "hmm", "eqtl+hmm", "enh+eqtl", "enh+hmm", "enh+eqtl+hmm"))

## make a summary plot (which isn't all that useful)
make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_sample_annotation_count_distributions"), width_ratio = 2.5)
print(ggplot(sample_counts, aes(x=count, fill=tissue_class)) +
    cat_col_scale +
    geom_histogram(binwidth=1, center=0.5, closed="left") +
    ## add lines for the sample counts
    geom_vline(aes(xintercept=count), data=input_annot_counts, color=muted("red"),
               size=0.75, linetype=3) +
    facet_grid(annotation ~ tissue_class, scales="free") +
    theme_bw() +
    xlab("SNP count") + ggtitle("Sample count distributions") +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
          legend.position="none",
          title=element_text(size=20), plot.title = element_text(hjust = 0.5)))
dev.off()

## now we have to make separate plots for each annotation
for(annot in unique(input_annot_counts$annotation)) {
    make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_sample_", annot, "_count_distributions"), width_ratio=2.0)
    print(ggplot(sample_counts[sample_counts$annotation==annot,], aes(x=count, fill=tissue_class)) +
          geom_histogram(binwidth=1, center=0.5, closed="left") +
          ## add lines for the sample counts
          geom_vline(aes(xintercept=count), color=muted("red"), size=0.75, linetype=3,
                     data=input_annot_counts[input_annot_counts$annotation==annot,]) +
          cat_col_scale +
          facet_wrap(~ tissue_class, scales="free") +
          theme_bw() + xlab("SNP count") + ylab("# of samples") +
          ggtitle(paste0("Sample count distributions for ", annot)) +
          theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
                legend.position="none", axis.title=element_text(size=16),
                strip.text=element_text(size=16),
                title=element_text(size=20), plot.title = element_text(hjust = 0.5)))
    dev.off()
}

## -----------------------------------------------------------------------------
## 9. Analysis plots for annotation overlap and p-values of LD collapsed analysis
## -----------------------------------------------------------------------------
## combined region analysis
## first, a boxplot of the p-value distributions. to do this, we need to melt the p-value matrices
collapsed_pval_df <- rbind(cbind(pval="Raw P-value", melt(collapsed_empirical_pvals, varnames=c("tissue_class", "annotation"))),
                 cbind(pval="BH-adjusted P-value", melt(collapsed_bh_adj_pvals, varnames=c("tissue_class", "annotation"))))
collapsed_pval_df$annotation <- factor(collapsed_pval_df$annotation, ordered=T,
                             levels=c("locus_enh", "roadmap_hmm_enh", "locus_enh+roadmap_hmm_enh"),
                             labels=c("FANTOM5 Enhancer", "Roadmap HMM Enhancer", "FANTOM5 Enh+Roadmap HMM Enh"))
collapsed_pval_df$tissue_class <- factor(collapsed_pval_df$tissue_class, ordered=T,
                               levels=sort(unique(collapsed_pval_df$tissue_class), dec=T))

write.table(collapsed_pval_df, paste0(output_dir, '/tables/all_region_collapsed_bootstrap_results.txt'), quote=F, sep="\t", row.names=F, col.names=T)
## ## to read in these results:
## collapsed_pval_df <- read.table(paste0(output_dir, '/tables/all_region_collapsed_bootstrap_results.txt'), header=T, sep="\t", quote="", as.is=T)
## collapsed_pval_df$annotation <- factor(collapsed_pval_df$annotation, ordered=T,
##                              levels=c("GTEx eQTL", "FANTOM5 Enhancer", "Roadmap HMM Enhancer",
##                                  "FANTOM5 Enh+GTEx eQTL", "GTEx eQTL+Roadmap HMM Enh", "FANTOM5 Enh+Roadmap HMM Enh", "FANTOM5 Enh+GTEx eQTL+Roadmap HMM Enh"))
## collapsed_pval_df$tissue_class <- factor(collapsed_pval_df$tissue_class, ordered=T,
##                                levels=sort(unique(collapsed_pval_df$tissue_class), dec=T))

## -------------------
## make the boxplot
make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_collapsed_annotation_overlap_pval_boxplots"))
print(ggplot(collapsed_pval_df, aes(x=annotation, y=value, fill=pval)) +
      geom_boxplot() + scale_y_continuous(breaks=seq(0, 1, by=0.05)) +
      scale_fill_discrete(h=c(100, 200)) +
      theme_bw() + xlab("Annotation") +
      ylab("P-value") + ggtitle("LD-collapsed P-value distributions for each annotation") +
      theme(legend.position="right",
            legend.text=element_text(size=20), legend.title=element_text(size=20),            
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(vjust=1, hjust=1, size=15),
            axis.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## -------------------
## also make a histogram of the p-values across all annotations
make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_collapsed_raw_pval_hist"))
print(ggplot(collapsed_pval_df[collapsed_pval_df$pval=="Raw P-value",], aes(x=value)) +
      geom_histogram(binwidth=0.01) +
      scale_x_continuous(breaks=seq(0, 1, by=0.05)) +
      theme_bw() + xlab("P-value") +
      ylab("Number of tissue-annotation combinations") +
      ggtitle("Raw P-value histogram, LD collapsed") +
      theme(legend.position="none",
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(vjust=1, hjust=1, size=15),            
            axis.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_collapsed_bh_pval_hist"))
print(ggplot(collapsed_pval_df[collapsed_pval_df$pval=="BH-adjusted P-value",], aes(x=value)) +
      geom_histogram(binwidth=0.01) +
      scale_x_continuous(breaks=seq(0, 1, by=0.05)) +
      theme_bw() + xlab("P-value") +
      ylab("Number of tissue-annotation combinations") +
      ggtitle("BH-adjusted P-value histogram, LD collapsed") +
      theme(legend.position="none",
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(vjust=1, hjust=1, size=15),
            title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## make histograms split by annotation type
collapsed_pval_df$annot_wrap <- unlist(lapply(strwrap(collapsed_pval_df$annotation, width=15, simplify=F),
                                    paste, collapse="\n"))
sorted_wrap_levels <- unlist(lapply(strwrap(levels(collapsed_pval_df$annotation), width=15, simplify=F),
                                    paste, collapse="\n"))
collapsed_pval_df$annot_wrap <- factor(collapsed_pval_df$annot_wrap, ordered=T,
                                       levels=sorted_wrap_levels)

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_collapsed_raw_pval_split_annot_hist"))
print(ggplot(collapsed_pval_df[collapsed_pval_df$pval=="Raw P-value",], aes(x=value, fill=annotation)) +
      geom_histogram(binwidth=0.01) +
      scale_x_continuous(breaks=seq(0, 1, by=0.05)) +
      facet_grid(annot_wrap ~ ., scales="free") +
      theme_bw() + xlab("P-value") +
      ylab("Number of tissue-annotation combinations") +
      ggtitle("Raw P-value histogram, LD collapsed") +
      theme(legend.position="none",
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(vjust=1, hjust=1, size=15),
            strip.text = element_text(size=20),
            axis.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_collapsed_bh_pval_split_annot_hist"))
print(ggplot(collapsed_pval_df[collapsed_pval_df$pval=="BH-adjusted P-value",], aes(x=value, fill=annotation)) +
      geom_histogram(binwidth=0.01) +
      scale_x_continuous(breaks=seq(0, 1, by=0.05)) +
      facet_grid(annot_wrap ~ ., scales="free") +
      theme_bw() + xlab("P-value") +
      ylab("Number of tissue-annotation combinations") +
      ggtitle("BH-adjusted P-value histogram, LD collapsed") +
      theme(legend.position="none",
            axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(angle=45, vjust=1, hjust=1, size=15),
            strip.text = element_text(size=20),
            axis.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## -------------------
## make a heatmap of the p-values
## define the breaks
hm_breaks <- c(0, 0.05, 1)

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_collapsed_annotation_by_tissue_raw_pval_heatmap"))
print(ggplot(collapsed_pval_df[collapsed_pval_df$pval=="Raw P-value",],
             aes(x=annotation, y=tissue_class)) +
      geom_tile(aes(fill=value), colour="#000000") +
      ## use 4 colors here so that all significant hits are fairly red
      scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                           values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                           limits=c(0, 1), name="P-value", guide="colorbar") +
      theme_bw() + ylab("Tissue category") +
      xlab("Annotation") + ggtitle("Heatmap of p-values, LD collapsed") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## also make one with text in the boxes for the p-values
make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_collapsed_annotation_by_tissue_raw_pval_heatmap_with_text"))
print(ggplot(collapsed_pval_df[collapsed_pval_df$pval=="Raw P-value",],
             aes(x=annotation, y=tissue_class)) +
      geom_tile(aes(fill=value), colour="#000000") +
      scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                           values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                           limits=c(0, 1), name="P-value", guide="colorbar") +
      theme_bw() + ylab("Tissue category") +
      xlab("Annotation") + ggtitle("Heatmap of p-values, LD collapsed") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)) +
      geom_text(aes(x=annotation, y=tissue_class,
                    label=format(round(value, digits=4), scientific=F)), color="black", size=4))
dev.off()

make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_collapsed_annotation_by_tissue_bh_pval_heatmap"))
print(ggplot(collapsed_pval_df[collapsed_pval_df$pval=="BH-adjusted P-value",],
             aes(x=annotation, y=tissue_class)) +
      geom_tile(aes(fill=value), colour="#000000") +
      ## use 4 colors here so that all significant hits are fairly red
      scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                           values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                           limits=c(0, 1), name="BH-adjusted p-value", guide="colorbar") +
      theme_bw() + ylab("Tissue category") +
      xlab("Annotation") + ggtitle("Heatmap of adjusted p-values, LD collapsed") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
dev.off()

## also make one with text in the boxes for the p-values
make_graphic(paste0(output_dir, "/plots/", param_ref[['outprefix']], "_collapsed_annotation_by_tissue_bh_pval_heatmap_with_text"))
print(ggplot(collapsed_pval_df[collapsed_pval_df$pval=="BH-adjusted P-value",],
             aes(x=annotation, y=tissue_class)) +
      geom_tile(aes(fill=value), colour="#000000") +
      scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                           values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                           limits=c(0, 1), name="BH-adjusted p-value", guide="colorbar") +
      theme_bw() + ylab("Tissue category") +
      xlab("Annotation") + ggtitle("Heatmap of adjusted p-values, LD collapsed") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)) +
      geom_text(aes(x=annotation, y=tissue_class,
                    label=format(round(value, digits=4), scientific=F)), color="black", size=4))
dev.off()

## -------------------
## split tag region analysis
## melt the results
collapsed_split_pval_df <- rbind(cbind(pval="Raw P-value", melt(collapsed_split_empirical_pvals, varnames=c("tissue_class", "annotation", "tag_region"))),
                 cbind(pval="BH-adjusted P-value", melt(collapsed_split_bh_adj_pvals, varnames=c("tissue_class", "annotation", "tag_region"))))
collapsed_split_pval_df$annotation <- factor(collapsed_split_pval_df$annotation, ordered=T,
                             levels=c("gtex_eqtl", "locus_enh", "roadmap_hmm_enh", "locus_enh+gtex_eqtl", "gtex_eqtl+roadmap_hmm_enh", "locus_enh+roadmap_hmm_enh", "locus_enh+gtex_eqtl+roadmap_hmm_enh"),
                             labels=c("GTEx eQTL", "FANTOM5 Enhancer", "Roadmap HMM Enhancer",
                                 "FANTOM5 Enh+GTEx eQTL", "GTEx eQTL+Roadmap HMM Enh", "FANTOM5 Enh+Roadmap HMM Enh", "FANTOM5 Enh+GTEx eQTL+Roadmap HMM Enh"))
collapsed_split_pval_df$tissue_class <- factor(collapsed_split_pval_df$tissue_class, ordered=T,
                               levels=sort(unique(collapsed_split_pval_df$tissue_class), dec=T))

write.table(collapsed_split_pval_df, paste0(output_dir, '/tables/split_region_collapsed_bootstrap_results.txt'), quote=F, sep="\t", row.names=F, col.names=T)
## ## to read in this data
## collapsed_split_pval_df <- read.table(paste0(output_dir, '/tables/split_region_collapsed_bootstrap_results.txt'), header=T, sep="\t", quote="", as.is=T)
## collapsed_split_pval_df$annotation <- factor(collapsed_split_pval_df$annotation, ordered=T,
##                              levels=c("GTEx eQTL", "FANTOM5 Enhancer", "Roadmap HMM Enhancer",
##                                  "FANTOM5 Enh+GTEx eQTL", "GTEx eQTL+Roadmap HMM Enh", "FANTOM5 Enh+Roadmap HMM Enh", "FANTOM5 Enh+GTEx eQTL+Roadmap HMM Enh"))
## collapsed_split_pval_df$tissue_class <- factor(collapsed_split_pval_df$tissue_class, ordered=T,
##                                levels=sort(unique(collapsed_split_pval_df$tissue_class), dec=T))

dir.create(paste0(output_dir, '/plots/ld_collapsed_split_tag_regions/'), F, T)
for(this_tag in unique(collapsed_split_pval_df$tag_region)) {
    ## TODO: fix this!
    tag_out <- gsub("/", "_", strsplit(this_tag, ":")[[1]][1])
    dir.create(paste0(output_dir, "/plots/ld_collapsed_split_tag_regions/", tag_out), F, T)

    ## define the breaks
    hm_breaks <- c(0, 0.05, 1)

    ## grab the relevant data
    this_pval_df <- subset(collapsed_split_pval_df, pval=="BH-adjusted P-value" & tag_region==this_tag)
    
    make_graphic(paste0(output_dir, "/plots/ld_collapsed_split_tag_regions/", tag_out, '/',
                        param_ref[['outprefix']], "_", tag_out, "_collapsed_annotation_by_tissue_bh_pval_heatmap"))
    print(ggplot(this_pval_df, 
                 aes(x=annotation, y=tissue_class)) +
          geom_tile(aes(fill=value), colour="#000000") +
          ## use 4 colors here so that all significant hits are fairly red
          scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                               values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                               limits=c(0, 1), name="BH-adjusted p-value", guide="colorbar") +
          theme_bw() + ylab("Tissue category") +
          xlab("Annotation") + ggtitle(paste("Heatmap of adjusted p-values,", this_tag,
                                             "\nLD collapsed")) +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
            axis.text.y = element_text(size=15), axis.title=element_text(size=20),
            legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
            legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
    dev.off()

    ## also make one with text in the boxes for the p-values
    make_graphic(paste0(output_dir, "/plots/ld_collapsed_split_tag_regions/", tag_out, '/',
                        param_ref[['outprefix']], "_", tag_out, "_collapsed_annotation_by_tissue_bh_pval_heatmap_with_text"))
    print(ggplot(this_pval_df,
                 aes(x=annotation, y=tissue_class)) +
          geom_tile(aes(fill=value), colour="#000000") +
          scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                               values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                               limits=c(0, 1), name="BH-adjusted p-value", guide="colorbar") +
          theme_bw() + ylab("Tissue category") +
          xlab("Annotation") + ggtitle(paste("Heatmap of adjusted p-values,", this_tag,
                                             "\nLD collapsed")) +
          theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
                axis.text.y = element_text(size=15), axis.title=element_text(size=20),
                legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
                legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)) +
          geom_text(aes(x=annotation, y=tissue_class,
                    label=format(round(value, digits=4), scientific=F)), color="black", size=4))
    dev.off()

    ## also make plots for the uncorrected p-values    
    ## grab the relevant data
    this_raw_pval_df <- subset(collapsed_split_pval_df, pval=="Raw P-value" & tag_region==this_tag)
    
    make_graphic(paste0(output_dir, "/plots/ld_collapsed_split_tag_regions/", tag_out, '/',
                        param_ref[['outprefix']], "_", tag_out, "_collapsed_annotation_by_tissue_raw_pval_heatmap"))
    print(ggplot(this_raw_pval_df, 
                 aes(x=annotation, y=tissue_class)) +
          geom_tile(aes(fill=value), colour="#000000") +
          ## use 4 colors here so that all significant hits are fairly red
          scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                               values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                               limits=c(0, 1), name="Raw p-value", guide="colorbar") +
          theme_bw() + ylab("Tissue category") +
          xlab("Annotation") + ggtitle(paste("Heatmap of unadjusted p-values,", this_tag,
                                             "\nLD collapsed")) +
          theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
                axis.text.y = element_text(size=15), axis.title=element_text(size=20),
                legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
                legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)))
    dev.off()

    ## also make one with text in the boxes for the p-values
    make_graphic(paste0(output_dir, "/plots/ld_collapsed_split_tag_regions/", tag_out, '/',
                        param_ref[['outprefix']], "_", tag_out, "_collapsed_annotation_by_tissue_raw_pval_heatmap_with_text"))
    print(ggplot(this_raw_pval_df,
                 aes(x=annotation, y=tissue_class)) +
          geom_tile(aes(fill=value), colour="#000000") +
          scale_fill_gradientn(colours=c("red", muted("red"), "#08306B", muted("blue")),
                               values=c(0, 0.05, 0.0501, 1), breaks=hm_breaks, labels=format(hm_breaks),
                               limits=c(0, 1), name="Raw p-value", guide="colorbar") +
          theme_bw() + ylab("Tissue category") +
          xlab("Annotation") + ggtitle(paste("Heatmap of unadjusted p-values,", this_tag,
                                             "\nLD collapsed")) +
          theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
                axis.text.y = element_text(size=15), axis.title=element_text(size=20),
                legend.text=element_text(size=20), legend.key.height=unit(0.10, "npc"),
                legend.title=element_text(size=20), plot.title = element_text(size=30, hjust = 0.5)) +
          geom_text(aes(x=annotation, y=tissue_class,
                        label=format(round(value, digits=4), scientific=F)), color="black", size=4))
    dev.off()
}

cat("Entire enhancer bootstrapping analysis took:\n")
cat((proc.time() - full_analysis_start_time)[["elapsed"]], 'seconds\n')
