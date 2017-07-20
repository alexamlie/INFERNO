"""
pval_expand_tagsnp_set.py
Alex Amlie-Wolf, 05/31/17

Takes in a list of tag SNPs and expands using distance and significance thresholds
"""

import argparse, os, sys

# input_snp_list = "/home/alexaml/data/enhancer_snp_pipeline/input_data/IGAP_with_genes.txt"
# gwas_directory = "/home/alexaml/data/IGAP_data/"
# output_file = "/home/alexaml/code/ad_enhancer_snps/input_parsing/temp_test.txt"
# sig_threshold = 10.0
# dist_threshold = 500000

def expand_tagsnp_list(input_snp_list, summary_stats_file, output_file, output_allele_info, rsid_col, pos_col, chr_col, pval_col, sig_multiplier=None, locus_thresh=None, dist_threshold=500000):
    """
    does the actual expansion by p-value
    this method uses a single summary statistics file as input
    """
    ## create the output directory, if needed
    try:
        os.makedirs(os.path.dirname(output_file))
    except OSError:
        pass
    
    ## first read in the tag SNPs, store in dictionary
    tag_snp_dict = {}

    ## we will write the tag SNPs to the output
    with open(output_file, 'wb') as output_snps:
        ## read in the input SNPs
        num_tagsnps = 0
        with open(input_snp_list, 'rb') as input_snps:
            for line in input_snps:
                num_tagsnps += 1
                snp_data = line.strip().split("\t")
                ## store rsID, region name, and position, tagged by chromosome
                if snp_data[0] in tag_snp_dict:
                    tag_snp_dict[snp_data[0]].append([snp_data[1], snp_data[2], snp_data[3]])
                else:
                    tag_snp_dict[snp_data[0]] = [[snp_data[1], snp_data[2], snp_data[3]]]
                ## write to output (if we don't need to get their information)
                if not output_allele_info:
                    output_snps.write("\t".join([snp_data[0], snp_data[1], snp_data[2]+":"+snp_data[1],
                                                snp_data[3]])+'\n')

        ## loop through the summary stats file once to get tag SNP p-values
        ## keep track of any chromosomes not found in the tag SNPs for reporting
        mismatch_chrs = []        
        print "Looping through summary stats file to get tag SNP information"
        num_found_tagsnps = 0
        with open(summary_stats_file, 'rb') as summary_stats:
            for line in summary_stats:
                this_snp_data = line.strip().split("\t")
                chr_orig = this_snp_data[chr_col]
                chr_string = "chr"+str(this_snp_data[chr_col])
                ## figure out which chromosome annotation to use
                if chr_orig in tag_snp_dict.keys():
                    chr_key = chr_orig
                elif chr_string in tag_snp_dict.keys():
                    chr_key = chr_string
                else:
                    if chr_orig not in mismatch_chrs:
                        mismatch_chrs.append(chr_orig)
                        print "Chromosome %s / %s from summary file not found in tag variants file!" % (chr_orig, chr_string)
                    ## keep looping through the file since we don't have tag SNPs on this chromosome
                    continue
                        
                for tagsnp in tag_snp_dict[chr_key]:
                    if this_snp_data[rsid_col]==tagsnp[0]:
                        sys.stdout.write("Found matching entry for "+tagsnp[0]+" ")
                        num_found_tagsnps += 1
                        tagsnp.append(float(this_snp_data[pval_col]))
                        if output_allele_info:
                            ## write the full input line back out
                            output_snps.write("\t".join([tagsnp[0]]+this_snp_data)+"\n")
        sys.stdout.write("\n")
        sys.stdout.flush()

        print "Found "+str(num_found_tagsnps)+" out of "+str(num_tagsnps)+" tag variants"
        if num_found_tagsnps == 0:
            print "No matching variants found in summary statistics!"
            return
        
        print "Looping through summary statistics and performing p-value expansion"
                                                
        ## now loop through again to actually do the expansion
        ## track how many SNPs we expand, including the tag variants
        num_expanded_snps = num_tagsnps
        
        with open(summary_stats_file, 'rb') as summary_stats:
            for line in summary_stats:
                this_snp_data = line.strip().split("\t")
                chr_orig = str(this_snp_data[chr_col])
                chr_string = "chr"+str(this_snp_data[chr_col])
                ## figure out which chromosome annotation to use
                if chr_orig in tag_snp_dict.keys():
                    chr_key = chr_orig
                elif chr_string in tag_snp_dict.keys():
                    chr_key = chr_string
                else:
                    ## we should never see a mismatch chr we didn't already find, but check anyway
                    if chr_orig not in mismatch_chrs:
                        mismatch_chrs.append(chr_orig)
                        print "Chromosome %s / %s from summary file not found in tag variants!" % (chr_orig, chr_string)
                    ## keep looping through the file because we don't have tag SNPs on this chromosome
                    continue
                
                for tagsnp in tag_snp_dict[chr_key]:
                    ## make sure we actually found the p-value for this one
                    if(len(tagsnp) < 4):
                        continue
                    ## first check distance and make sure it's not a tag SNP
                    this_snp_pos = int(this_snp_data[pos_col])
                    if (abs(int(tagsnp[2])-this_snp_pos) <= dist_threshold and tagsnp[0]!=this_snp_data[rsid_col]):
                        ## check whichever p-value approach we are using and write it out if it
                        ## passes:
                        this_snp_pval = float(this_snp_data[pval_col])
                        if sig_multiplier and this_snp_pval <= tagsnp[3] * sig_multiplier:
                            num_expanded_snps += 1
                            if output_allele_info:
                               output_snps.write("\t".join([tagsnp[0]]+this_snp_data)+"\n")
                            else:      
                                output_snps.write("\t".join([chr_key, this_snp_data[rsid_col], tagsnp[1]+":"+tagsnp[0], str(this_snp_pos)])+"\n")
                        if locus_thresh and this_snp_pval <= locus_thresh:
                            num_expanded_snps += 1
                            if output_allele_info:
                               output_snps.write("\t".join([tagsnp[0]]+this_snp_data)+"\n")
                            else:      
                                output_snps.write("\t".join([chr_key, this_snp_data[rsid_col], tagsnp[1]+":"+tagsnp[0], str(this_snp_pos)])+"\n")
                                    
    print "Expanded from "+str(num_tagsnps)+" tagging variants to "+str(num_expanded_snps)+" variants"
    
def expand_tagsnp_list_consistent_effects(input_snp_list, summary_stats_file, output_file, output_allele_info, rsid_col, pos_col, chr_col, pval_col, beta_col, maf_col, sig_multiplier=None, locus_thresh=None, dist_threshold=500000):
    """   
    this method uses a single summary statistics file as input, and that file must include beta
    effect values and minor allele frequencies. this is to make sure that any expanded variants
    are in the same direction as the tag SNPs
    """
    ## create the output directory, if needed
    try:
        os.makedirs(os.path.dirname(output_file))
    except OSError:
        pass
    
    ## first read in the tag SNPs, store in dictionary
    tag_snp_dict = {}

    ## we will write the tag SNPs to the output
    with open(output_file, 'wb') as output_snps:
        ## read in the input SNPs
        num_tagsnps = 0
        with open(input_snp_list, 'rb') as input_snps:
            for line in input_snps:
                num_tagsnps += 1
                snp_data = line.strip().split("\t")
                ## store rsID, region name, and position, tagged by chromosome
                if snp_data[0] in tag_snp_dict:
                    tag_snp_dict[snp_data[0]].append([snp_data[1], snp_data[2], snp_data[3]])
                else:
                    tag_snp_dict[snp_data[0]] = [[snp_data[1], snp_data[2], snp_data[3]]]
                ## write to output (if we don't need to get their information)
                if not output_allele_info:
                    output_snps.write("\t".join([snp_data[0], snp_data[1], snp_data[2]+":"+snp_data[1],
                                                snp_data[3]])+'\n')

        ## loop through the summary stats file once to get tag SNP p-values
        ## keep track of any chromosomes not found in the tag SNPs for reporting
        mismatch_chrs = []
        print "Looping through summary stats file to get tag SNP information"
        num_found_tagsnps = 0
        with open(summary_stats_file, 'rb') as summary_stats:
            for line in summary_stats:
                this_snp_data = line.strip().split("\t")
                chr_orig = this_snp_data[chr_col]
                chr_string = "chr"+str(this_snp_data[chr_col])
                ## figure out which chromosome annotation to use
                if chr_orig in tag_snp_dict.keys():
                    chr_key = chr_orig
                elif chr_string in tag_snp_dict.keys():
                    chr_key = chr_string
                else:
                    if chr_orig not in mismatch_chrs:
                        mismatch_chrs.append(chr_orig)
                        print "Chromosome %s / %s from summary file not found in tag variants file!" % (chr_orig, chr_string)
                    ## keep looping through the file since we don't have tag SNPs on this chromosome
                    continue
                        
                for tagsnp in tag_snp_dict[chr_key]:
                    if this_snp_data[rsid_col]==tagsnp[0]:
                        sys.stdout.write("Found matching entry for "+tagsnp[0]+" ")
                        num_found_tagsnps += 1
                        tagsnp.append(float(this_snp_data[pval_col]))
                        ## also include the direction of effect, relative to minor allele
                        this_effect = float(this_snp_data[beta_col])
                        if float(this_snp_data[maf_col]) > 0.50:
                            this_effect = -this_effect
                        tagsnp.append(this_effect)
                        if output_allele_info:
                            ## write the full input line back out
                            output_snps.write("\t".join([tagsnp[0]]+this_snp_data)+"\n")
                            
        sys.stdout.write("\n")
        sys.stdout.flush()

        print "Found "+str(num_found_tagsnps)+" out of "+str(num_tagsnps)+" tag variants"
        if num_found_tagsnps == 0:
            print "No matching variants found in summary statistics!"
            return
        
        print "Looping through summary statistics and performing p-value expansion"
                                                
        ## now loop through again to actually do the expansion
        ## track how many SNPs we expand, including the tag variants
        num_expanded_snps = num_tagsnps
        
        with open(summary_stats_file, 'rb') as summary_stats:
            for line in summary_stats:
                this_snp_data = line.strip().split("\t")
                chr_orig = str(this_snp_data[chr_col])
                chr_string = "chr"+str(this_snp_data[chr_col])
                ## figure out which chromosome annotation to use
                if chr_orig in tag_snp_dict.keys():
                    chr_key = chr_orig
                elif chr_string in tag_snp_dict.keys():
                    chr_key = chr_string
                else:
                    ## we should never see a mismatch chr we didn't already find, but check anyway
                    if chr_orig not in mismatch_chrs:
                        mismatch_chrs.append(chr_orig)
                        print "Chromosome %s / %s from summary file not found in tag variants!" % (chr_orig, chr_string)
                    ## keep looping through the file because we don't have tag SNPs on this chromosome
                    continue
                
                for tagsnp in tag_snp_dict[chr_key]:
                    ## make sure we actually found data for this tag SNP
                    if(len(tagsnp) < 5):
                        continue
                    ## first check distance and make sure it's not a tag SNP
                    this_snp_pos = int(this_snp_data[pos_col])
                    if (abs(int(tagsnp[2])-this_snp_pos) <= dist_threshold and tagsnp[0]!=this_snp_data[rsid_col]):
                        this_snp_effect = float(this_snp_data[beta_col])
                        if float(this_snp_data[maf_col]) > 0.50:
                            this_snp_effect = -this_snp_effect
                        if tagsnp[4] * this_snp_effect > 0:                                
                            ## now check whichever p-value approach we are using and write it
                            ## out if it passes:
                            this_snp_pval = float(this_snp_data[pval_col])
                            if sig_multiplier and this_snp_pval <= tagsnp[3] * sig_multiplier:
                                num_expanded_snps += 1
                                if output_allele_info:
                                   output_snps.write("\t".join([tagsnp[0]]+this_snp_data)+"\n")
                                else:      
                                    output_snps.write("\t".join([chr_key, this_snp_data[rsid_col], tagsnp[1]+":"+tagsnp[0], str(this_snp_pos)])+"\n")
                            if locus_thresh and this_snp_pval <= locus_thresh:
                                num_expanded_snps += 1
                                if output_allele_info:
                                   output_snps.write("\t".join([tagsnp[0]]+this_snp_data)+"\n")
                                else:      
                                    output_snps.write("\t".join([chr_key, this_snp_data[rsid_col], tagsnp[1]+":"+tagsnp[0], str(this_snp_pos)])+"\n")
                                                                        
    print "Expanded from "+str(num_tagsnps)+" tagging variants to "+str(num_expanded_snps)+" variants"
                                            
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Expand tag SNPs by distance and significance")
    parser.add_argument("--output_allele_info", action='store_true', help="Whether to output all the allele information (the full data from the summary statistics file) for the expanded SNP set")
    parser.add_argument("--sig_multiplier", type=float, help="The multiplier range for significance to use (i.e. a value of 10 means one order of magnitude)")
    parser.add_argument("--locus_thresh", type=float, help="The absolute p-value cutoff to use for locus-wide significance analysis (only one of sig_multiplier and locus_thresh may be specified)")
    parser.add_argument("--dist_threshold", default=500000, type=int, help="How far away from each tag SNP to look, in base pairs.")
    ## TODO: allow definition of header line in summary stats and reference columns by name.
    ## when i do this, just read in the first line of the summary stats file and define the
    ## columns so i can use the same methods
    parser.add_argument("--rsid_col", type=int, required=True, help="The column number containing SNP rsIDs (1-based). This is required.")
    parser.add_argument("--pos_col", type=int, required=True, help="The column number containing SNP positions  (1-based). This is required.")
    parser.add_argument("--pval_col", type=int, required=True, help="The column number containing SNP p-values (1-based). This is required.")
    parser.add_argument("--chr_col", type=int, required=True, help="The column number containing SNP p-values (1-based). This is required.")
    parser.add_argument("--beta_col", type=int, help="The column number containing SNP beta values. Required to define expanded sets that all have consistent minor allele effect directions.") 
    parser.add_argument("--maf_col", type=int, help="The column number containing SNP MAF values. Required to define expanded sets that all have consistent minor allele effect directions.")
    parser.add_argument("input_snp_list", help="The list of tag SNPs (chr, rsID, name, pos)")
    parser.add_argument("summary_stats", help="The file containig summary statistics.")
    parser.add_argument("output_file", help="The path to the desired output file")

    pargs = parser.parse_args()

    if not (pargs.rsid_col and pargs.pos_col and pargs.chr_col and pargs.pval_col):
        print "This script requires rsID, position, chromosome, and p-value column definitions!"
    elif pargs.sig_multiplier and pargs.locus_thresh:
        print "Cannot specify both relative and absolute thresholds!"
    elif pargs.sig_multiplier:
        if pargs.beta_col and pargs.maf_col:
            print "Expanding variants with consistent effect directions"
            expand_tagsnp_list_consistent_effects(pargs.input_snp_list, pargs.summary_stats,
                                                  pargs.output_file, pargs.output_allele_info,
                                                  pargs.rsid_col-1, pargs.pos_col-1,
                                                  pargs.chr_col-1, pargs.pval_col-1,
                                                  pargs.beta_col-1, pargs.maf_col-1,
                                                  sig_multiplier = pargs.sig_multiplier,
                                                  dist_threshold = pargs.dist_threshold)
        else:
            print "Expanding variants with no consideration of effect direction"
            expand_tagsnp_list(pargs.input_snp_list, pargs.summary_stats, pargs.output_file,
                               pargs.output_allele_info, pargs.rsid_col-1,
                               pargs.pos_col-1, pargs.chr_col-1, pargs.pval_col-1,
                               sig_multiplier = pargs.sig_multiplier,
                               dist_threshold = pargs.dist_threshold)
    elif pargs.locus_thresh:
        if 0.0 <= pargs.locus_thresh <= 1.0:
            if pargs.beta_col and pargs.maf_col:
                expand_tagsnp_list_consistent_effects(pargs.input_snp_list, pargs.summary_stats,
                                                    pargs.output_file, pargs.output_allele_info,
                                                    pargs.rsid_col-1, pargs.pos_col-1,
                                                    pargs.chr_col-1, pargs.pval_col-1,
                                                    pargs.beta_col-1, pargs.maf_col-1,
                                                    locus_thresh = pargs.locus_thresh,
                                                    dist_threshold = pargs.dist_threshold) 

            else:
                expand_tagsnp_list(pargs.input_snp_list, pargs.summary_stats, pargs.output_file,
                                   pargs.output_allele_info, pargs.rsid_col-1,
                                   pargs.pos_col-1, pargs.chr_col-1, pargs.pval_col-1,
                                   locus_thresh = pargs.locus_thresh,
                                   dist_threshold = pargs.dist_threshold) 
        else:
            print "Invalid locus-wide significance threshold (should be between 0 and 1)"
    else:
        print "Need either relative or absolute thresholds!"        
