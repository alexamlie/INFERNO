"""
summarize_ld_maf_dist_results.py
Alex Amlie-Wolf, started 04/19/2016

A script that takes LD, MAF, distance to TSS information and summarizes them for each SNP of interest
"""

import argparse, glob, os, resource, gzip

# ld_threshold = 0.7
# ld_dir = '/home/alexaml/data/1000_genomes/phase1_release_v3/EUR/pairwise_ld'
# maf_dir = '/home/alexaml/data/1000_genomes/phase1_release_v3/EUR/MAF_info'
# tss_dist_dir = '/home/alexaml/data/1000_genomes/phase1_release_v3/EUR/dist_to_tss'
# output_dir = '/home/alexaml/data/1000_genomes/phase1_release_v3/EUR/snp_maf_tss_ld_summary'

def summarize_snp_info(ld_threshold, ld_dir, maf_dir, tss_dist_dir, output_dir, store_alleles):
    """
    the main function to summarize information for each SNP
    """
    ## first make the output directory
    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    ## now we have to get all the files for each directory and match them up
    ld_files = glob.glob(ld_dir+"/*.ld.gz")
    maf_files = glob.glob(maf_dir+"/*.frq")
    dist_files  = glob.glob(tss_dist_dir+"/*_closest_tss.txt")

    ## assume that the chromosomes are contained in the first period-separated field of the
    ## filename
    ld_chrs = [str.split(os.path.basename(x), ".")[0] for x in ld_files]
    maf_chrs = [str.split(os.path.basename(x), ".")[0] for x in maf_files]
    dist_chrs = [str.split(os.path.basename(x), ".")[0] for x in dist_files]

    ## get the common chromosomes:
    chr_set = list(set(ld_chrs).intersection(set(maf_chrs)).intersection(set(dist_chrs)))

    ## check for any chromosomes that are not in all data sets
    ld_maf_non_overlap_chrs = set(ld_chrs).symmetric_difference(set(maf_chrs))
    if len(ld_maf_non_overlap_chrs) > 0:
        print "Found non-overlapping chromosomes between LD and MAF data:"
        for x in ld_maf_non_overlap_chrs:
            print x

    ld_dist_non_overlap_chrs = set(ld_chrs).symmetric_difference(set(dist_chrs))
    if len(ld_dist_non_overlap_chrs) > 0:
        print "Found non-overlapping chromosomes:"
        for x in ld_dist_non_overlap_chrs:
            print x

    dist_maf_non_overlap_chrs = set(dist_chrs).symmetric_difference(set(maf_chrs))
    if len(dist_maf_non_overlap_chrs) > 0:
        print "Found non-overlapping chromosomes:"
        for x in dist_maf_non_overlap_chrs:
            print x

    if store_alleles:
        outfile = output_dir+"/snp_maf_tss_dist_"+str(ld_threshold)+"_ld_info_with_alleles.txt"
        outf_header = ['chr', 'rsID', 'pos', 'tss_dist', 'MAF', 'num_ld_partners', 'maj', 'min']
    else:
        outfile = output_dir+"/snp_maf_tss_dist_"+str(ld_threshold)+"_ld_info.txt"
        outf_header = ['chr', 'rsID', 'pos', 'tss_dist', 'MAF', 'num_ld_partners']
    outf_idx = {outf_header[x]:x for x in range(len(outf_header))}

    with open(outfile, 'wb') as outf:
        outf.write("\t".join(outf_header)+"\n")

        ## now go through our chromosomes and make the comparisons
        ## for convenience and nice output, process the chromosomes in order
        chr_nums = [int(x.replace("chr", "")) for x in chr_set]
        chr_order = [i[0] for i in sorted(enumerate(chr_nums), key=lambda x:x[1])]

        for chr_idx in chr_order:
            this_chr = chr_set[chr_idx]
            print "Analyzing LD, MAF, and TSS distance for chromosome %s" % (this_chr)

            this_ld_file = filter(lambda x:this_chr+"." in x, ld_files)[0]
            this_maf_file = filter(lambda x:this_chr+"." in x, maf_files)[0]
            this_dist_file = filter(lambda x:this_chr+"." in x, dist_files)[0]

            ## initialize dict to store various information about each SNP
            this_chr_dict = {}
            ## initialize a separate dict to store the partners of each SNP
            snp_partner_dict = {}

            ## start by reading through the MAF file to get the allele and MAF information, and
            ## also because we only want to store single nucleotide variants in the dict
            with open(this_maf_file, 'rb') as this_maf:
                maf_header = next(this_maf).strip().split("\t")
                maf_idx = {maf_header[x]:x for x in range(len(maf_header))}
                for line in this_maf:
                    line_data = line.strip().split("\t")
                    this_snp_key = line_data[maf_idx['rsID']]+":"+line_data[maf_idx['pos']]
                    ## check that it's a single nucleotide variant
                    this_alleles = [line_data[maf_idx['A1']], line_data[maf_idx['A2']]]
                    ## only process it if it is
                    if len(this_alleles[0])==1 and len(this_alleles[1])==1:
                        ## if it's not in the dict, initialize it
                        if this_snp_key not in this_chr_dict:
                            ## store the info: rsID, position, distance to TSS (initialized to
                            ## -1), MAF, number of LD partners (initialized to 0), major allele,
                            ## and minor allele
                            this_chr_dict[this_snp_key] = [line_data[maf_idx['rsID']], line_data[maf_idx['pos']], -1, line_data[maf_idx['MAF']], 0] + this_alleles
                            ## also initialize the set storage for this SNP
                            snp_partner_dict[this_snp_key] = set()
                        else:
                            print "Error: SNP %s observed multiple times in MAF data!" % (this_snp_key)
                            print "Using first observed alleles.."

            ## now go through the distance file to get the TSS distances
            with open(this_dist_file, 'rb') as this_dist:
                for line in this_dist:
                    line_data = line.strip().split("\t")
                    ## the key is the rsID and the position
                    this_snp_key = line_data[3]+":"+line_data[1]
                    ## if it's in the dict, it's an SNV so we want to store it
                    if this_snp_key in this_chr_dict:
                        ## we also check to see if it's been set before
                        cur_dist = this_chr_dict[this_snp_key][2]
                        if cur_dist==-1:
                            this_chr_dict[this_snp_key][2] = line_data[-1]
                        else:
                            ## check if our distance is different
                            if cur_dist != line_data[-1]:
                                print "SNP %s observed with different distance values! Using first observed one." % (this_snp_key)

            ## now read through the LD data to get the number of LD partners
            with gzip.open(this_ld_file, 'rb') as this_ld:
                ld_header = next(this_ld).strip().split("\t")
                ld_idx = {ld_header[x]:x for x in range(len(ld_header))}
                for line in this_ld:
                    line_data = line.strip().split("\t")

                    snp_a_key = line_data[ld_idx["SNP_A"]]+":"+line_data[ld_idx["BP_A"]]
                    snp_b_key = line_data[ld_idx["SNP_B"]]+":"+line_data[ld_idx["BP_B"]]

                    ## check for meeting LD threshold
                    if float(line_data[ld_idx['R2']]) >= ld_threshold:
                        ## get the data: for these, if either SNP cannot be found, it means
                        ## that one or both is not a single nucleotide variant so we want to
                        ## skip this entry entirely
                        try:
                            snp_a_data = this_chr_dict[snp_a_key]
                        except KeyError:
                            continue
                        try:
                            snp_b_data = this_chr_dict[snp_b_key]
                        except KeyError:
                            continue
                        ## now we know they are both single nucleotide variants:
                        ## we have to check whether the partner was observed already
                        if snp_b_key not in snp_partner_dict[snp_a_key]:
                            ## increment the partner count for the A SNP
                            this_chr_dict[snp_a_key][4] += 1
                            ## add SNP B to the A SNP set
                            snp_partner_dict[snp_a_key].add(snp_b_key)

                        if snp_a_key not in snp_partner_dict[snp_b_key]:
                            ## increment the partner count for the B  SNP
                            this_chr_dict[snp_b_key][4] += 1
                            ## add SNP A to the B SNP set
                            snp_partner_dict[snp_b_key].add(snp_a_key)
                            
            ## write the results out
            for snp in this_chr_dict:
                snp_data = this_chr_dict[snp]
                ## check the various quantities
                if snp_data[3] == "-1":
                    print "Error: MAF not found for snp %s!" % (snp_data[0]+":"+snp_data[1])
                if snp_data[2] == "-1":
                    print "Error: TSS distance not found for snp %s!" % (snp_data[0]+":"+snp_data[1])
                if snp_data[5]=="NA" or snp_data[6]=="NA":
                    print "Error: alleles not found for snp %s!" % (snp_data[0]+":"+snp_data[1])
                if store_alleles:
                    outf.write("\t".join([this_chr]+[str(x) for x in snp_data])+'\n')
                else:
                    outf.write("\t".join([this_chr]+[str(x) for x in snp_data[:5]])+'\n')

            print "Finished processing chromosome %s; current memory usage %s Mb" % (this_chr, str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000.0))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Summarize LD, MAF, and distance to TSS information for SNPs", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ## TODO: optional argument for output prefix in data directories to only analyze some
    parser.add_argument("--ld_threshold", type=float, default=0.7, help="The threshold to use for counting the number of LD partners")
    parser.add_argument("--store_alleles", action="store_true", help="Optional flag to include allele information in the output. This makes the output file larger but can be useful, especially for skipping LD expansion in the main pipeline script.")
    parser.add_argument("ld_dir", help="The directory containing the LD results (as generated by calculate_pairwise_ld.sh script)")
    parser.add_argument("maf_dir", help="The directory containing the MAF results (as generated by calculate_maf.sh script)")
    parser.add_argument("tss_dist_dir", help="The directory containing the distance to TSS results (as generated by calculate_tss_distance.sh script)")
    parser.add_argument("output_dir", help="The path to the desired output directory")

    pargs = parser.parse_args()

    summarize_snp_info(pargs.ld_threshold, pargs.ld_dir, pargs.maf_dir, pargs.tss_dist_dir, pargs.output_dir, pargs.store_alleles)
