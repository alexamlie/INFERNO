"""
find_ensembl_genes_by_chr.py
alex amlie-wolf 09/06/2016
pulls out all ensembl genes and splits them into files for each chromosome
"""

import argparse, os, re

def parse_ensembl_genes(ensembl_ref, outdir):
    try:
        os.makedirs(outdir)
    except OSError:
        pass

    ## we need to create files for each chromosome (that we care about). these will be
    ## initialized as we encounter each of them
    chr_file_dict = {}

    ## a regex to check that we are in a 'real' chromosome
    nonref_chr_pattern = re.compile("^chr.*\_.*$")                        
    
    with open(ensembl_ref, 'rb') as input_genes:
        gene_header = next(input_genes).strip().split("\t")
        gene_idx = {gene_header[x]:x for x in range(len(gene_header))}
        ## now loop through the genes and find the appropriate ones
        for line in input_genes:
            this_data = line.strip().split("\t")
            this_chr = this_data[gene_idx['chrom']]
            
            ## if we don't have a 'normal' chromosome, move on
            if nonref_chr_pattern.match(this_chr) or this_chr=="chrM":
                continue

            ## otherwise, write to the correct file
            ## initialize it if necessary
            if this_chr not in chr_file_dict:
                print "Parsing chromosome %s" % (this_chr)
                chr_file_dict[this_chr] = open(outdir+"/"+this_chr+"_genes.txt", 'wb')
                chr_file_dict[this_chr].write("transcript_id\tgene_id\n")
            ## now write it out to the file
            chr_file_dict[this_chr].write("\t".join([this_data[gene_idx['name']], this_data[gene_idx['name2']]])+"\n")

    ## once we're done reading through everything, close all the files
    for this_chr in chr_file_dict:
        chr_file_dict[this_chr].close()
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Pull out ensembl genes for each chromosome")
    parser.add_argument("ensembl_ref", help="The Ensembl reference file")
    parser.add_argument("outdir", help="The desired output directory where files for each chromosome will be written")

    pargs = parser.parse_args()

    parse_ensembl_genes(pargs.ensembl_ref, pargs.outdir)

