"""
generate_FANTOM5_facet_promoter_expression.py
Alex Amlie-Wolf, 09/06/17

Takes in the FANTOM5 CAGE peak based expression table (RLE normalized) for human with
annotation file (hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz) and facet annotations
and generates facet-level promoter peak annotations
"""

import argparse, os, sys, gzip, subprocess

def generate_facet_expression(tpm_file, output_dir, cell_facet_f, tissue_facet_f):
    ## read in the cell line -> facet mappings
    cell_facets = {}
    with open(cell_facet_f, 'rb') as cell_data:
        header_read = False
        for line in cell_data:
            if not header_read:
                header_data = line.strip().split("\t")
                cell_idx = {header_data[x]:x for x in range(len(header_data))}
                header_read = True
            else:
                line_data = line.strip().split("\t")
                ## match the enhancer approach
                this_facet = line_data[cell_idx['Facet ontology ID']]+"_"+line_data[cell_idx['Facet ontology term']].replace(" ", "_")
                cell_facets[line_data[cell_idx['Sample ID']]] = this_facet

    ## read in the tissue -> facet mappings
    tissue_facets = {}
    with open(tissue_facet_f, 'rb') as tissue_data:
        header_read = False
        for line in tissue_data:
            if not header_read:
                header_data = line.strip().split("\t")
                tissue_idx = {header_data[x]:x for x in range(len(header_data))}
                header_read = True
            else:
                line_data = line.strip().split("\t")
                ## match the enhancer approach
                this_facet = line_data[tissue_idx['Facet ontology ID']]+"_"+line_data[tissue_idx['Facet ontology term']].replace(" ", "_")
                tissue_facets[line_data[tissue_idx['Sample ID']]] = this_facet                

    ## now make a dictionary to store peaks for every facet
    facet_peaks = {}
    facet_counts = {}
    for facet in set(cell_facets.values()):
        facet_peaks[facet] = set()
        facet_counts[facet] = 0
    for facet in set(tissue_facets.values()):
        facet_peaks[facet] = set()
        facet_counts[facet] = 0
                        
    ## skip through the comment lines of the expression file
    header_read = False
    lines_read = 0
    with gzip.open(tpm_file, 'rb') as tpm_data:
        for line in tpm_data:
            if line[0]!="#" and not header_read:
                ## this is the header, read it in
                header_data = line.strip().split("\t")
                header_idx = {x:header_data[x] for x in range(len(header_data))}
                ## make the names nicer, we only need the IDs
                for i in range(7, len(header_data)):
                    header_idx[i] = header_idx[i].split(".")[2]                    
                header_read = True
            elif line[0]!="#":
                lines_read += 1
                if lines_read % 50000==0:
                    print "Read "+str(lines_read)+" lines"
                ## this is a data line, parse the peak information
                line_data = line.strip().split("\t")
                if line_data[0][0:3]=="chr":                                        
                    this_chr = line_data[0].split(":")[0]
                    this_start = line_data[0].split(":")[1].split("..")[0]
                    this_end = line_data[0].split("..")[1].split(",")[0]
                    this_strand = line_data[0].split("..")[1].split(",")[1]
                    ## get the short description for this promoter
                    this_description = line_data[1]
                    ## now loop through all the samples
                    for i in range(7, len(header_idx)):
                        ## use a threshold of 1 TPM
                        if float(line_data[i]) > 1:
                            if header_idx[i] in tissue_facets:
                                this_facet = tissue_facets[header_idx[i]]
                            elif header_idx[i] in cell_facets:
                                this_facet = cell_facets[header_idx[i]]
                            else:
                                continue
                            facet_peaks[this_facet].add("\t".join([this_chr, this_start, this_end, this_strand, this_description]))
                            facet_counts[this_facet] += 1
            ## otherwise, we skip comment lines

    ## output the facet-level files
    for facet in facet_peaks:
        this_set = facet_peaks[facet]
        print str(len(this_set))+" unique peaks found in "+facet+", observed "+str(facet_counts[facet])+" times across data sources in facet"
        peak_outf = output_dir+"/"+facet+"_expressed_promoters.bed"
        with open(peak_outf, 'wb') as peak_out:
            for peak in this_set:
                peak_data = peak.split("\t")
                peak_out.write("\t".join([peak_data[0], peak_data[1], peak_data[2],
                                          peak_data[4], "0", peak_data[3]])+"\n")

        ## now sort this file
        with open(peak_outf+".sorted", 'wb') as sort_out:
            subprocess.call(['sort', '-k1,1V', '-k2,2n', '-k3,3n', '-k6,6', peak_outf], stdout=sort_out)
        subprocess.call(["mv", peak_outf+".sorted", peak_outf])
                                             
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Process FANTOM5 promoter data")
    parser.add_argument("tpm_file", help="The FANTOM5 peak based expression table (hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz")
    parser.add_argument("output_dir", help="The directory to write the facet-level promoter peabed files to")
    parser.add_argument("cell_facet_table_file", help="The FANTOM5 file containing the sample ID -> facet mappings for cell lines")
    parser.add_argument("tissue_facet_table_file", help="The FANTOM5 file containing the sample ID -> facet mappings for tissues")

    pargs = parser.parse_args()

    ## create the output directory, if needed
    try:
        os.makedirs(os.path.dirname(pargs.output_dir))
    except OSError:
        pass
    
    ## run the expression analysis
    generate_facet_expression(pargs.tpm_file, pargs.output_dir, pargs.cell_facet_table_file, pargs.tissue_facet_table_file)
