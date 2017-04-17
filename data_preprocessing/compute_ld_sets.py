"""
compute_ld_sets.py
Alex Amlie-Wolf, started 06/14/16

A script that looks through pairwise LD files and outputs all expanded sets of variants at a given
threshold. It uses allele information to only look at single nucleotide variants. Also note that
this will only output information about SNPs that have >=1 LD partner, so variants in their own LD
blocks have to be considered outside of this script
"""

import argparse, glob, os, gzip, time

# ld_threshold = 0.7
# ld_dir = '/home/alexaml/data/1000_genomes/phase1_release_v3/EUR/pairwise_ld'
# maf_dir = '/home/alexaml/data/1000_genomes/phase1_release_v3/EUR/MAF_info'
# output_dir = '/home/alexaml/data/1000_genomes/phase1_release_v3/EUR/precomputed_ld_sets'

def compute_ld_sets(ld_threshold, ld_dir, maf_dir, output_dir):
    """
    the main function to figure out what the LD expanded sets are
    """
    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    ## now we have to get all the files for each directory and make sure they match up
    ld_files = glob.glob(ld_dir+"/*.ld.gz")
    maf_files = glob.glob(maf_dir+"/*.frq")
    ld_chrs = [str.split(os.path.basename(x), ".")[0] for x in ld_files]
    maf_chrs = [str.split(os.path.basename(x), ".")[0] for x in maf_files]

    ## get the common chromosomes:
    chr_set = list(set(ld_chrs).intersection(set(maf_chrs)))

    ## check for any chromosomes that are not in all data sets
    ld_maf_non_overlap_chrs = set(ld_chrs).symmetric_difference(set(maf_chrs))
    if len(ld_maf_non_overlap_chrs) > 0:
        print "Found non-overlapping chromosomes between LD and MAF data:"
        for x in ld_maf_non_overlap_chrs:
            print x

    ## for convenience and nice output, process the chromosomes in order
    chr_nums = [int(x.replace("chr", "")) for x in chr_set]
    ## http://stackoverflow.com/questions/6422700/how-to-get-indices-of-a-sorted-array-in-python
    chr_order = [i[0] for i in sorted(enumerate(chr_nums), key=lambda x:x[1])]

    for chr_idx in chr_order:
        this_chr = chr_set[chr_idx]
        start_time = time.clock()
        print "Analyzing %s" % (this_chr)
        ## initialize storage dict for SNPs and their partners:
        snp_partners = {}
        ## define another dict indexed by position that contains the rsIDs of each SNP so that
        ## we can efficiently remove SNPs from consideration
        snp_positions = {}
        ## for this purpose, we also define a variable to track the earliest position we
        ## currently have
        cur_min_pos = -1

        this_ld_file = filter(lambda x:this_chr+"." in x, ld_files)[0]
        this_maf_file = filter(lambda x:this_chr+"." in x, maf_files)[0]

        ## first read through the MAF file to store the alleles for each SNP
        print "Reading through MAF file to save allele information"
        this_chr_alleles = {}
        with open(this_maf_file, 'rb') as this_maf:
            maf_header = next(this_maf).strip().split("\t")
            maf_idx = {maf_header[x]:x for x in range(len(maf_header))}
            for line in this_maf:
                line_data = line.strip().split("\t")
                this_snp_key = line_data[maf_idx['rsID']]+":"+line_data[maf_idx['pos']]
                this_alleles = [line_data[maf_idx['A1']], line_data[maf_idx['A2']]]

                ## only do anything if we have an SNV
                if len(this_alleles[0])==1 and len(this_alleles[1])==1:
                    if this_snp_key not in this_chr_alleles:
                        this_chr_alleles[this_snp_key] = this_alleles
                    else:
                        print "Error: SNV %s observed multiple times in MAF data!" % (this_snp_key)
                        print "Using first observed alleles.."

        print "Performing LD precomputation"
        output_file = output_dir+"/"+this_chr+"_LD_sets_"+str(ld_threshold)+"_threshold.txt"
        with gzip.open(this_ld_file, 'rb') as this_ld, open(output_file, 'wb') as out_partners:
            out_partners.write("\t".join(["rsID", "expanded_set"])+"\n")

            ld_header = next(this_ld).strip().split("\t")
            ld_idx = {ld_header[x]:x for x in range(len(ld_header))}
            for line in this_ld:
                line_data = line.strip().split("\t")
                ## check if this pair meets the LD threshold, and don't process if not
                if float(line_data[ld_idx['R2']]) >= ld_threshold:
                    ## if they match, we have to add both SNPs to the partner dict
                    ## save the positions and rsIDs from this line
                    snp_A_pos = int(line_data[ld_idx['BP_A']])
                    snp_B_pos = int(line_data[ld_idx['BP_B']])
                    snp_A_rsID = line_data[ld_idx['SNP_A']]
                    snp_B_rsID = line_data[ld_idx['SNP_B']]

                    ## the key needs to be defined before replacing "." rsIDs
                    snp_A_key = snp_A_rsID+":"+str(snp_A_pos)
                    snp_B_key = snp_B_rsID+":"+str(snp_B_pos)
                    ## we check our keys; if either are not in the allele dict, we just skip this whole entry because the allele dict only stores SNVs
                    try:
                        snp_A_alleles = this_chr_alleles[snp_A_key]
                    except KeyError:
                        ## don't print a message because this is probably a non-SNV
                        continue
                    try:
                        snp_B_alleles = this_chr_alleles[snp_B_key]
                    except KeyError:
                        continue

                    ## for rsIDs, replace '.' with chr-pos format (for output)
                    snp_A_rsID = snp_A_rsID if snp_A_rsID != "." else this_chr+"-"+str(snp_A_pos)
                    snp_B_rsID = snp_B_rsID if snp_B_rsID != "." else this_chr+"-"+str(snp_B_pos)

                    ## now add the matching SNP rsIDs, indexed by rsID. note that we know they
                    ## are both SNVs because we found both of them in the allele dict
                    if snp_A_rsID in snp_partners:
                        snp_partners[snp_A_rsID].append(snp_B_rsID)
                    else:
                        snp_partners[snp_A_rsID] = [snp_B_rsID]
                        ## if we haven't seen it yet, we have to add it to the position dict
                        if snp_A_pos in snp_positions:
                            ## we might have multiple SNPs at a position
                            snp_positions[snp_A_pos].append(snp_A_rsID)
                        else:
                            snp_positions[snp_A_pos] = [snp_A_rsID]

                    if snp_B_rsID in snp_partners:
                        snp_partners[snp_B_rsID].append(snp_A_rsID)
                    else:
                        snp_partners[snp_B_rsID] = [snp_A_rsID]
                        if snp_B_pos in snp_positions:
                            snp_positions[snp_B_pos].append(snp_B_rsID)
                        else:
                            snp_positions[snp_B_pos] = [snp_B_rsID]

                    ## if the SNP A position is past our current minimum position, write the
                    ## partners of that position (and any others before the SNP A position) out
                    if snp_A_pos > cur_min_pos:
                        try:
                            min_pos_rsIDs = snp_positions[cur_min_pos]
                        except KeyError:
                            ## this is just for the case of the first position we see
                            cur_min_pos = snp_A_pos
                            continue
                        ## write out the LD partners for this minimum position
                        for snp_rsID in min_pos_rsIDs:
                            ## note that we convert this to a set before writing it out so that
                            ## we only have unique partners
                            out_partners.write("\t".join([snp_rsID, ",".join(set(snp_partners[snp_rsID]))])+"\n")
                            ## remove this SNP from the LD partner list
                            del snp_partners[snp_rsID]
                        ## remove the position from the position list
                        del snp_positions[cur_min_pos]

                        ## now loop through the positions, writing out the information until we
                        ## are past the current A SNP position.  we need to store the positions
                        ## we're removing so that we don't edit the dict while looping
                        pos_to_remove = []
                        for snp_pos in sorted(snp_positions.keys()):
                            if snp_pos < snp_A_pos:
                                this_pos_rsIDs = snp_positions[snp_pos]
                                for snp_rsID in this_pos_rsIDs:
                                    ## same thing here: convert to a set to get only unique output
                                    out_partners.write("\t".join([snp_rsID, ",".join(set(snp_partners[snp_rsID]))])+"\n")
                                    del snp_partners[snp_rsID]
                                pos_to_remove.append(snp_pos)
                            else:
                                cur_min_pos = snp_pos
                                break
                        for pos in pos_to_remove:
                            del snp_positions[pos]

            ## now write out whatever partners remain, if any
            for snp_rsID in snp_partners:
                out_partners.write("\t".join([snp_rsID, ",".join(set(snp_partners[snp_rsID]))])+"\n")

        end_time = time.clock()
        print "Time for this chromosome: %s" % (str(end_time - start_time))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Compute sets of LD variants using pairwise LD information. Note that this only considers single nucleotide variants!", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ## TODO: optional argument for output prefix in data directories
    parser.add_argument("ld_threshold", type=float, help="The threshold to use for counting the LD blocks")
    parser.add_argument("ld_dir", help="The directory containing the LD results (as generated by calculate_pairwise_ld.sh script)")
    parser.add_argument("maf_dir", help="The directory containing the MAF and allele information for the SNPs of interest (as generated by calculate_maf.sh script)")
    parser.add_argument("output_dir", help="The path to the desired output directory")

    pargs = parser.parse_args()

    compute_ld_sets(pargs.ld_threshold, pargs.ld_dir, pargs.maf_dir, pargs.output_dir)
