"""
ld_prune_snp_set.py

Alex Amlie-Wolf 05/31/17
Takes in an input file in INFERNO format and the 1,000 genomes directory and generates a vcf file
for the input SNPs, then runs PLINK on that file to do LD pruning, and finally converts this back
to INFERNO input format
"""

import argparse, glob, os, gzip, subprocess, copy

def convert_to_vcf(ld_threshold, input_snp_list, vcf_ref_dir, output_dir):
    try:
        os.makedirs(output_dir)
    except OSError:
        pass
    
    ## get the vcf files
    vcf_files = glob.glob(vcf_ref_dir+"/*.vcf.gz")
    vcf_chrs = [str.split(os.path.basename(x), ".")[0] for x in vcf_files]

    ## read in the input SNP list and store the SNPs by chromosome
    input_snp_dict = {}
    with open(input_snp_list, 'rb') as snpin:
        for line in snpin:
            this_data = line.strip().split("\t")
            this_chr = this_data[0]
            this_rsID = this_data[1]
            this_tag = this_data[2]
            this_pos = this_data[3]
            if this_chr in input_snp_dict:
                input_snp_dict[this_chr].append([this_rsID, this_tag, this_pos])
            else:
                input_snp_dict[this_chr] = [[this_rsID, this_tag, this_pos]]

    ## copy this to another dict that we actually go through and remove stuff from. we don't
    ## want to lose this one because we need it to convert back to INFERNO format
    search_snp_dict = copy.deepcopy(input_snp_dict)
                
    ## now read through the VCF files and find our SNPs, writing them to the output
    with open(output_dir+"/converted_input_snp_list.vcf", "wb") as vcf_out:
        ## we only want to write the header once so just track if we're on the first loop
        firstloop = True
        for this_chr in input_snp_dict:
            print "Analyzing %s" % (this_chr)
            this_vcf_file = filter(lambda x:this_chr+"." in x, vcf_files)[0]
            ## TODO: check that this file exists
            with gzip.open(this_vcf_file, 'rb') as kg_vcf_data:
                for line in kg_vcf_data:
                    if line[0]=="#":
                        if firstloop:
                            vcf_out.write(line)
                    ## otherwise, only check it if it's an SNV
                    else:
                        line_data = line.strip().split("\t")
                        this_pos = line_data[1]
                        this_rsID = line_data[2]
                        this_maj_allele = line_data[3]
                        this_min_allele = line_data[4]
                        if len(this_maj_allele)==1 and len(this_min_allele)==1:
                            ## check the SNPs we have left on this chromosome
                            for input_snp in search_snp_dict[this_chr]:
                                snp_rsID = input_snp[0]
                                snp_pos = input_snp[2]
                                ## just check by position here
                                if snp_pos == this_pos:
                                    ## print a message if we have mismatched rsIDs
                                    if this_rsID != snp_rsID:
                                        print "rsID mismatch at position %s. 1kg: %s, input: %s" % (snp_pos, this_rsID, snp_rsID)
                                        ## set the input rsID to this
                                        for idx in range(len(input_snp_dict[this_chr])):
                                            if input_snp_dict[this_chr][idx][0]==snp_rsID:
                                                input_snp_dict[this_chr][idx][0] = this_rsID
                                                break
                                    ## write out the vcf line
                                    vcf_out.write(line)
                                    ## remove this SNP from the search dictionary and move on
                                    search_snp_dict[this_chr].remove(input_snp)
                                    break
                            ## if we're out of SNPs, stop looping through the 1kg file
                            if len(search_snp_dict[this_chr])==0:
                                break
                ## update the loop track variable
                firstloop = False

    ## write out a list of the SNPs that weren't found
    num_unmatched = 0
    with open(output_dir+"/unmatched_variants.txt", "wb") as unmatched_variants:
        for this_chr in search_snp_dict:
            for input_snp in search_snp_dict[this_chr]:
                num_unmatched += 1
                unmatched_variants.write("\t".join([this_chr, input_snp[0], input_snp[1], input_snp[2]])+"\n")
    print str(num_unmatched)+" variants were not found in 1kg data"                
                
    ## run PLINK on this vcf file, getting a pruned set. this just runs directly, creating a
    ## file called plink.prune.in containing the remaining rsIDs
    ## parameters: use a 500kb window, step by 1 variant at a time, and use the parameter LD thresh
    subprocess.call(["plink", "--vcf", output_dir+"/converted_input_snp_list.vcf",
                     "--indep-pairwise", "500kb", "1", str(ld_threshold)], cwd=output_dir)
    
    ## convert this pruned set to INFERNO-compatible input using stored tag data
    with open(output_dir+"/plink.prune.in", 'rb') as pruned_snp_set, open(output_dir+"/pruned_set_pipeline_input.txt", "wb") as converted_out:
        for line in pruned_snp_set:
            this_rsID = line.strip()
            ## we only have rsIDs so we have to check everything..
            snp_found = False
            for this_chr in input_snp_dict:
                for this_snp in input_snp_dict[this_chr]:
                    if this_rsID==this_snp[0]:
                        converted_out.write("\t".join([this_chr, this_snp[0], this_snp[1], this_snp[2]])+'\n')
                        snp_found = True
                        break
                if snp_found:
                    break

            if not snp_found:
                print "Error: SNP %s was not recovered in input data!" % (this_rsID)
                
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Pulls out SNP data from 1,000 genomes and converts it to vcf", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("ld_threshold", type=float, help="The threshold for LD pruning")    
    parser.add_argument("input_snp_list", help="The file containing the SNPs (in INFERNO format: chr, rsID, tag_name, position")
    parser.add_argument("vcf_ref_dir", help="The directory containing the 1000 genomes vcf files")
    parser.add_argument("output_dir", help="The path to the desired output directory (vcf files will be written as well as input files for the INFERNO pipeline)")

    pargs = parser.parse_args()

    convert_to_vcf(pargs.ld_threshold, pargs.input_snp_list, pargs.vcf_ref_dir, pargs.output_dir)
