"""
expand_and_annotate_snps.py
Alex Amlie-Wolf, started 11/10/2015

This is the main method for generating annotation data given a list of tag SNPs.
First, it expands the list of SNPs to include all SNPs within a user-specified base-pair distance
limit and LD threshold (this can be set to 0 to only include the input list).
Then, it annotates this list of SNPs with several genomic features
"""

import argparse, subprocess, sys, os, datetime, commands, gzip, re, glob, pickle, math

def create_logging_function(logfile, loglevel="full"):
    """
    a function to take care of logging as well as printing output to stdout
    """
    if loglevel=="full":
        def logging_function(outstring):
            print outstring
            logfile.write(outstring+"\n")
    elif loglevel=="print":
        def logging_function(outstring):
            print outstring
    elif loglevel=="save":
        def logging_function(outstring):
            logfile.write(outstring+"\n")
    elif loglevel=="quiet":
        def logging_function(outstring):
            pass
    else:
        print "Logging level %s not supported!" % (loglevel)
        sys.exit(1)
    return logging_function            

def expand_ld_snps(input_snp_list, kg_pop, kg_dir, ld_threshold, ld_check_area, outdir, outprefix, logging_function, buffer_size=None, plink_path=None, force_recalc=False, remove_ld_data=False):
    """
    given input SNPs, expand LD region
    returns path to the file
    """
    ## note that we put this analysis into a new folder
    outfile_path = outdir+"/ld_expansion/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)+".txt"
    if not os.path.isfile(outfile_path) or force_recalc:
        ## if we don't have a plink path, check to make sure it's otherwise there
        if not plink_path:
            plink_found = False
            for path in os.environ.get('PATH').split(":"):
                if os.path.exists(path+"/plink"):
                    plink_found = True
            if not plink_found:
                sys.exit("Plink not found in PATH. This program is required (>=v1.9)")

        try:
            os.makedirs(outdir+"/ld_expansion/")
        except OSError:
            pass

        ## sort the input SNP list by chr and position (only if they haven't been already)
        sorted_snp_list = outdir+"/ld_expansion/"+outprefix+"_sorted_snps_"+str(ld_threshold)+"_"+str(ld_check_area)+".txt"
        if not os.path.isfile(sorted_snp_list):
            with open(sorted_snp_list, 'wb') as sort_snps:
                subprocess.call(['sort', '-k1,1V', '-k4,4n', input_snp_list], stdout=sort_snps)

        ## find the number of lines (to know how often we should give output messages)
        num_snps = 0
        with open(sorted_snp_list, 'rb') as sort_snps:
            for line in sort_snps:
                num_snps += 1
        sparse_output = True if num_snps > 10000 else False
                
        ## define our header
        outf_header = ['chr', 'rsID', 'pos', 'ref', 'alt', 'MAF', 'tag_rsID', 'tag_pos', 'tag_MAF', 'tag_name', 'R2', 'Dprime']
        outf_idx = {outf_header[x]:x for x in range(len(outf_header))}
        with open(sorted_snp_list, 'rb') as sort_snps, open(outfile_path, 'wb') as outf:
            outf.write("\t".join(outf_header)+"\n")
            ## store all the SNPs on a given chromosome to find their neighbors
            this_chr = ""
            this_chr_snps = {}
            ## two counting variables: the SNPs being actually analyzed, and the SNPs we've
            ## looped through. the first is just for output purposes while the second is used
            ## for the binning implementation
            analyzed_snps = 0
            buffer_snp_ctr = 0
            for cur_snp in sort_snps: 
                cur_snp_data = cur_snp.strip().split("\t")
                cur_snp_chr = cur_snp_data[0]
                cur_snp_id = cur_snp_data[1]
                cur_snp_name = cur_snp_data[2]
                cur_snp_pos = cur_snp_data[3]
                ## if we have '.' rsID, replace with a more informative name (in the line,
                ## since we save the whole line data)
                if cur_snp_id == '.':
                    cur_snp_data[1] = cur_snp_chr+"-"+cur_snp_pos
                
                cur_snp_key = ":".join([cur_snp_data[1], cur_snp_name, cur_snp_pos])
                ## count towards the loop counter
                buffer_snp_ctr += 1
                ## possible conditions: 1. no buffering, and we're on the same
                ## chromosome. 2. buffering, and we're on the same chromosome and haven't
                ## reached the buffer size. 3. buffering, we're on the same chromosome and have
                ## reached the buffer size. 4. we're not on the same chromosome (buffering or
                ## not)
                if (cur_snp_chr == this_chr and 
                    (not buffer_size or (buffer_size and buffer_snp_ctr % buffer_size != 0))):
                    if cur_snp_key not in this_chr_snps:
                        this_chr_snps[cur_snp_key] = cur_snp_data
                    else:
                        logging_function("Found duplicated SNP ID in input: %s" % (cur_snp_key))
                ## in any other case, we go through the analysis (either different chromosome,
                ## or reached buffer size) 
                else:
                    ## if we had been previously looking at a chromosome:
                    if this_chr:
                        ## TODO: make this less hard-coded / allow for other references
                        kg_file = kg_dir+"/"+kg_pop+"/"+this_chr+".phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel."+kg_pop+".vcf.gz"

                        ## loop through each SNP and open its neighboring SNP file (if we are
                        ## performing LD expansion)
                        if ld_check_area != 0 or ld_threshold != 1.0:
                            neighbor_snp_files = {}
                        ## this contains a bunch of information for each SNP, needed for
                        ## finding neighbors
                        snp_match_info = {}
                        ## we also want a dict to store all the alleles for this chromosome
                        allele_dict = {}

                        ## go through all the unique SNPs on this chromosome (or buffer)
                        for snp_key in this_chr_snps:
                            snp_data = this_chr_snps[snp_key]
                            snp_chr = snp_data[0]
                            snp_id = snp_data[1]
                            snp_name = snp_data[2]
                            snp_pos = int(snp_data[3])
                            
                            analyzed_snps += 1
                            if not sparse_output:
                                logging_function("Analyzing SNP number %s - %s at position %s:%s, region %s" % (analyzed_snps, snp_id, snp_chr, str(snp_pos), snp_name))
                            ## if we have a buffer size, use that
                            elif buffer_size and analyzed_snps % (buffer_size / 10)==0:
                                logging_function("Analyzing SNP number %s - %s at position %s:%s, region %s" % (analyzed_snps, snp_id, snp_chr, str(snp_pos), snp_name))
                            ## otherwise, just use a generic approach
                            elif analyzed_snps % 50000==0:
                                logging_function("Analyzing SNP number %s - %s at position %s:%s, region %s" % (analyzed_snps, snp_id, snp_chr, str(snp_pos), snp_name))
                                
                            ## only analyze the neighboring SNPs if this is a unique key (we
                            ## wouldn't expect duplicates here since the key includes rsID,
                            ## position, and tag region)
                            if snp_key not in snp_match_info:
                                ## each SNP is indexed by ID and position
                                allele_key = ":".join([snp_id, str(snp_pos)])
                                neighboring_snps_file = outdir+"/ld_expansion/snp_data/"+allele_key+"/"+allele_key+"_"+str(ld_check_area)+"_bp_neighboring_snps.vcf"

                                if ld_check_area != 0 or ld_threshold != 1.0:
                                    try:
                                        os.makedirs(outdir+"/ld_expansion/snp_data/"+allele_key)
                                    except OSError:
                                        pass
                                    neighbor_snp_files[snp_key] = open(neighboring_snps_file, 'wb')
                                ## track the SNP ID to use for searching
                                search_rsID = snp_id
                                ## also we need to track whether we find a match for this SNP
                                this_snp_found = False
                                ## store all this information for future reference (including
                                ## the position, for convenience)
                                snp_match_info[snp_key] = [search_rsID, this_snp_found, neighboring_snps_file, snp_pos]
                                
                        ## we also want to store a list of the SNPs that we need to find the
                        ## neighbors for (which gets updated as we go through the file)
                        cur_search_snps = snp_match_info.keys()

                        if not buffer_size:
                            logging_function("Analyzing SNPs from chromosome %s" % (this_chr))
                        else:
                            logging_function("Analyzing SNPs from chromosome %s, buffer number %s" % (this_chr, str(buffer_snp_ctr/buffer_size)))
                        
                        ## now read through the reference file
                        with gzip.open(kg_file, 'r') as kg_data:
                            for line in kg_data:
                                ## figure out if it's a header line
                                if line[0]=="#":
                                    if ld_check_area != 0 or ld_threshold != 1.0:
                                        for snp_key in cur_search_snps:
                                            ## write the header line for all SNPs
                                            neighbor_snp_files[snp_key].write(line)
                                ## if it's not a header, only loop through the SNPs if it's a
                                ## single nucleotide variant (and we have any search SNPs left)
                                else:
                                    line_data = line.strip().split("\t")
                                    this_pos = int(line_data[1])
                                    this_rsID = line_data[2]
                                    this_maj_allele = line_data[3]
                                    this_min_allele = line_data[4]
                                    if len(this_maj_allele)==1 and len(this_min_allele)==1:
                                        ## we need to store the SNPs that we will remove at this 1kg line
                                        removed_snps = []
                                        for snp_key in cur_search_snps:
                                            ## get the original id of the SNP, for comparison purposes
                                            orig_snp_id = snp_key.split(":")[0]
                                            ## pull out the position and ID (to match with) of this SNP
                                            snp_id = snp_match_info[snp_key][0]
                                            snp_pos = snp_match_info[snp_key][3]
                                            start_pos = snp_pos - ld_check_area
                                            end_pos = snp_pos + ld_check_area
                                            ## check the position
                                            if this_pos >= start_pos and this_pos <= end_pos:
                                                ## if the 1kg variant has a "." rsID, replace
                                                ## it in the line to be consistent with our
                                                ## naming:
                                                if this_rsID == ".":
                                                    ## this can't be a ':' or it messes up later
                                                    ## splitting functions
                                                    line_data[2] = this_chr+"-"+str(this_pos)
                                                    this_rsID = this_chr+"-"+str(this_pos)
                                                ## store the alleles of the reference SNP!
                                                allele_key = ":".join([this_rsID, str(this_pos)])
                                                    
                                                ## store the reference and alt alleles
                                                if allele_key not in allele_dict:
                                                    ## store the reference and alt alleles
                                                    ## in the allele dict, using the rsID
                                                    ## and pos as key
                                                    allele_dict[allele_key] = line_data[3:5]
                                                
                                                ## check if it's the search SNP by position
                                                if this_pos == snp_pos:
                                                    ## check if this SNP matches the input SNP ID
                                                    if orig_snp_id != this_rsID:
                                                        logging_function("SNP %s does not match 1000 genomes rsID. Using ID %s for this input SNP" % (orig_snp_id, this_rsID))
                                                        ## set the search rsID for PLINK
                                                        snp_match_info[snp_key][0] = this_rsID
                                                        ## in this case, we also want to store
                                                        ## the alleles using the input SNP id
                                                        orig_allele_key = ":".join([orig_snp_id, str(this_pos)])
                                                        if orig_allele_key not in allele_dict:
                                                            allele_dict[orig_allele_key] = line_data[3:5]
                                                    ## set the boolean that this SNP was found
                                                    snp_match_info[snp_key][1] = True
                                                if ld_check_area != 0 or ld_threshold != 1.0:
                                                    ## write the line to the neighboring SNP file
                                                    neighbor_snp_files[snp_key].write("\t".join(line_data)+"\n")
                                            # if we've passed the end position, stop looking at
                                            # this SNP (mark it for removal and close the file)
                                            elif this_pos > end_pos:
                                                removed_snps.append(snp_key)
                                                if ld_check_area != 0 or ld_threshold != 1.0:
                                                    neighbor_snp_files[snp_key].close()
                                        ## now remove the SNPs that we've passed at this iteration
                                        [cur_search_snps.remove(x) for x in removed_snps]
                                        ## if we closed the last file, stop looking through the
                                        ## reference
                                        if len(cur_search_snps) == 0:
                                            break
                                            
                            ## after looping through the reference, close any files that may
                            ## happen to be open still (should only happen if the window around a
                            ## SNP extends past the chromosome end)
                            if ld_check_area != 0 or ld_threshold != 1.0:
                                for snp_key in cur_search_snps:
                                    neighbor_snp_files[snp_key].close()

                        ## now, we can run plink for each SNP. note that we are now running
                        ## this loop on all combinations of SNPs and regions, but if it was
                        ## already run, we just skip the LD expansion
                        for snp_key in this_chr_snps:
                            ## get the information for this snp
                            snp_data = this_chr_snps[snp_key]
                            snp_chr = snp_data[0]
                            ## this is the input rsID
                            snp_id = snp_data[1]
                            snp_name = snp_data[2]
                            ## this time, save it as a string
                            snp_pos = str(snp_data[3])
                            allele_key = ":".join([snp_id, str(snp_pos)])
                            
                            this_snp_info = snp_match_info[snp_key]
                            search_rsID = this_snp_info[0]
                            this_snp_found = this_snp_info[1]
                            neighboring_snps_file = this_snp_info[2]

                            ## if we are skipping the LD expansion, just write it out directly
                            ## so that we can skip calling PLINK. for now, skip the MAF
                            ## calculation
                            if ld_check_area == 0 and ld_threshold == 1.0:
                                if this_snp_found:
                                    outf.write("\t".join([snp_chr, snp_id, snp_pos]+
                                                         allele_dict[allele_key]+
                                                         ["NA", snp_id, snp_pos, "NA", snp_name, "1.0", "1.0"])+"\n")
                                else:
                                    logging_function("SNP %s not found in 1,000 genomes by ID or by position. No allele information is being written for this SNP!" % (snp_id))
                                    outf.write("\t".join([snp_chr, snp_id, snp_pos, "NA", "NA",
                                                          "NA", snp_id, snp_pos, "NA",
                                                          snp_name, "1.0", "1.0"])+"\n")
                            else:
                                ## convert the LD checking area to kb
                                ld_kb_check_area = ld_check_area / 1000
                                ## output prefix
                                plink_prefix = (outdir+"/ld_expansion/snp_data/"+allele_key+"/"+
                                                "_".join([allele_key, str(ld_check_area), "bp", str(ld_threshold), "cutoff", "snps"]))
                                ## only run analysis if we found the SNP and haven't done this
                                ## analysis yet (even if we are forcing recalculation, this is
                                ## because of potential duplicated IDs and SNPs belonging to more
                                ## than one tag region):
                                if not plink_path and this_snp_found and not os.path.isfile(plink_prefix+".ld"):
                                    with open(plink_prefix+"_pipeline.log", 'wb') as plink_log:
                                        subprocess.call(["plink", "--allow-no-sex", "--vcf", neighboring_snps_file,
                                                        "--r2", "with-freqs", "dprime",
                                                        "--ld-snp", search_rsID,
                                                        "--ld-window", "99999",
                                                        "--ld-window-kb", str(ld_kb_check_area),
                                                        "--ld-window-r2",
                                                        str(ld_threshold),
                                                        "--out", plink_prefix],
                                                        stdout=plink_log)
                                elif plink_path and this_snp_found and not os.path.isfile(plink_prefix+".ld"):
                                    with open(plink_prefix+"_pipeline.log", 'wb') as plink_log:
                                        subprocess.call([plink_path, "--allow-no-sex", "--vcf", neighboring_snps_file,
                                                        "--r2", "with-freqs", "dprime",
                                                        "--ld-snp", search_rsID,
                                                        "--ld-window", "99999",
                                                        "--ld-window-kb", str(ld_kb_check_area),
                                                        "--ld-window-r2",
                                                        str(ld_threshold),
                                                        "--out", plink_prefix],
                                                        stdout=plink_log)
                                elif os.path.isfile(plink_prefix+".ld"):
                                    logging_function("No need to recompute LD SNPs for SNP ID %s" % (snp_id))
                                else:
                                    logging_function("SNP %s not found in 1,000 genomes by ID or by position. Skipping!" % (snp_id))

                                ## now output to the master file, if it exists and we found the SNP
                                if os.path.isfile(plink_prefix+".ld") and this_snp_found:
                                    with open(plink_prefix+".ld", 'rb') as plinkout:
                                        # read in the header, make an index
                                        plink_header = re.sub("\s+", "\t", next(plinkout).strip()).split("\t")
                                        plink_idx = {plink_header[x]:x for x in range(len(plink_header))}
                                        for ld_line in plinkout:
                                            ld_data = re.sub("\s+", "\t", ld_line.strip()).split("\t")
                                            ## check if we found the search rsID (replace with tag snp ID)
                                            if ld_data[plink_idx['SNP_B']]==search_rsID:
                                                outf.write("\t".join(["chr"+ld_data[plink_idx['CHR_A']], snp_id,
                                                                      ld_data[plink_idx['BP_B']]]+
                                                                      allele_dict[allele_key]+
                                                                      [ld_data[plink_idx['MAF_B']],
                                                                       snp_id, ld_data[plink_idx['BP_A']],
                                                                       ld_data[plink_idx['MAF_A']],
                                                                       snp_name, ld_data[plink_idx['R2']],
                                                                       ld_data[plink_idx['DP']]])+'\n')
                                            else:
                                                ld_snp_key = ":".join([ld_data[plink_idx['SNP_B']],
                                                                       ld_data[plink_idx['BP_B']]])
                                                if ld_snp_key in allele_dict:
                                                    outf.write("\t".join(["chr"+ld_data[plink_idx['CHR_A']],
                                                                      ld_data[plink_idx['SNP_B']],
                                                                      ld_data[plink_idx['BP_B']]]+
                                                                      allele_dict[ld_snp_key]+
                                                                      [ld_data[plink_idx['MAF_B']],
                                                                       snp_id, ld_data[plink_idx['BP_A']],
                                                                       ld_data[plink_idx['MAF_A']],
                                                                       snp_name, ld_data[plink_idx['R2']],
                                                                       ld_data[plink_idx['DP']]])+'\n')
                                                else:
                                                    logging_function("SNP %s not in allele dict" % (ld_snp_key))
                                                    outf.write("\t".join(["chr"+ld_data[plink_idx['CHR_A']],
                                                                      ld_data[plink_idx['SNP_B']],
                                                                      ld_data[plink_idx['BP_B']]]+
                                                                      ["NA", "NA"]+
                                                                      [ld_data[plink_idx['MAF_B']],
                                                                       snp_id, ld_data[plink_idx['BP_A']],
                                                                       ld_data[plink_idx['MAF_A']],
                                                                       snp_name, ld_data[plink_idx['R2']],
                                                                       ld_data[plink_idx['DP']]])+'\n')
                                                    
                                elif not os.path.isfile(plink_prefix+".ld"):
                                    logging_function("Plink failed to generate an LD file for snp %s (search ID %s)!" % (snp_id, search_rsID))

                    ## now reset the tracking variables (to the SNP we're actually on!)
                    this_chr = cur_snp_chr
                    this_chr_snps = {}
                    ## add this SNP to the dict
                    this_chr_snps[cur_snp_key] = cur_snp_data

            ## now analyze the last remaining chromosome (or buffer)
            if not buffer_size:
                logging_function("Analyzing SNPs from chromosome %s" % (this_chr))
            else:
                logging_function("Analyzing SNPs from chromosome %s, buffer number %s" % (this_chr, str(buffer_snp_ctr/buffer_size)))

            ## todo: make this less hard-coded / allow for other references
            kg_file = kg_dir+"/"+kg_pop+"/"+this_chr+".phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel."+kg_pop+".vcf.gz"

            ## loop through each SNP and open its neighboring SNP file
            if ld_check_area != 0 or ld_threshold != 1.0:
                neighbor_snp_files = {}
            ## this contains a bunch of information for each SNP for matching neighbors
            snp_match_info = {}
            ## finally, we want a dict to store all the alleles for this chromosome
            allele_dict = {}

            for snp_key in this_chr_snps:
                snp_data = this_chr_snps[snp_key]
                snp_chr = snp_data[0]
                snp_id = snp_data[1]
                snp_name = snp_data[2]
                snp_pos = int(snp_data[3])

                analyzed_snps += 1
                if not sparse_output:
                    logging_function("Analyzing SNP number %s - %s at position %s:%s, region %s" % (analyzed_snps, snp_id, snp_chr, str(snp_pos), snp_name))
                ## if we have a buffer size, use that
                elif buffer_size and analyzed_snps % (buffer_size / 10)==0:
                    logging_function("Analyzing SNP number %s - %s at position %s:%s, region %s" % (analyzed_snps, snp_id, snp_chr, str(snp_pos), snp_name))
                ## otherwise, just use a generic approach
                elif analyzed_snps % 50000==0:
                    logging_function("Analyzing SNP number %s - %s at position %s:%s, region %s" % (analyzed_snps, snp_id, snp_chr, str(snp_pos), snp_name))

                ## only analyze the neighboring SNPs if this is a unique ID (although the ID
                ## contains name, position, and tag region info, so they should always be
                ## unique..)
                if snp_key not in snp_match_info:
                    ## each SNP is indexed by ID and position
                    allele_key = ":".join([snp_id, str(snp_pos)])
                    neighboring_snps_file = outdir+"/ld_expansion/snp_data/"+allele_key+"/"+allele_key+"_"+str(ld_check_area)+"_bp_neighboring_snps.vcf"
                    
                    if ld_check_area != 0 or ld_threshold != 1.0:
                        try:
                            os.makedirs(outdir+"/ld_expansion/snp_data/"+allele_key)
                        except OSError:
                            pass
                        neighbor_snp_files[snp_key] = open(neighboring_snps_file, 'wb')                        
                    ## track the SNP ID to use for searching
                    search_rsID = snp_id
                    ## also we need to track whether we find a match for this SNP
                    this_snp_found = False
                    ## store all this information for future reference (also by rsID)
                    snp_match_info[snp_key] = [search_rsID, this_snp_found, neighboring_snps_file, snp_pos]
                    
            ## we also want to store a list of the SNPs that we need to find the
            ## neighbors for
            cur_search_snps = snp_match_info.keys()
            
            ## now read through the reference file
            with gzip.open(kg_file, 'r') as kg_data:
                for line in kg_data:
                    ## figure out if it's a header line
                    if line[0]=="#":
                        if ld_check_area != 0 or ld_threshold != 1.0:
                            for snp_key in cur_search_snps:
                                ## write the header line for all SNPs
                                neighbor_snp_files[snp_key].write(line)
                    else:
                        ## if it's not a header, only loop through the SNPs if it's
                        ## a single nucleotide variant and we have any search SNPs left
                        line_data = line.strip().split("\t")
                        this_pos = int(line_data[1])
                        this_rsID = line_data[2]
                        this_maj_allele = line_data[3]
                        this_min_allele = line_data[4]
                        if len(this_maj_allele)==1 and len(this_min_allele)==1:
                            ## store the SNPs we will remove at this line
                            removed_snps = []
                            for snp_key in cur_search_snps:
                                ## get the original id of the SNP, for comparison purposes
                                orig_snp_id = snp_key.split(":")[0]
                                ## get the position and matching rsID for this SNP
                                snp_id = snp_match_info[snp_key][0]
                                snp_pos = snp_match_info[snp_key][3]
                                allele_key = ":".join([snp_id, str(snp_pos)])
                                start_pos = snp_pos - ld_check_area
                                end_pos = snp_pos + ld_check_area
                                ## check the position 
                                if this_pos >= start_pos and this_pos <= end_pos:
                                    ## if it has a "." rsID, replace it in the line:
                                    if this_rsID == ".":
                                        ## this can't be a ':' or it messes up later analysis
                                        line_data[2] = this_chr+"-"+str(this_pos)
                                        this_rsID = this_chr+"-"+str(this_pos)
                                    ## store the alleles of the reference SNP!
                                    allele_key = ":".join([this_rsID, str(this_pos)])

                                    ## store the reference and alt alleles
                                    if allele_key not in allele_dict:
                                        ## store the reference and alt alleles in the allele
                                        ## dict, using the rsID and pos as key
                                        allele_dict[allele_key] = line_data[3:5]

                                    ## check if it's the search SNP by position
                                    if this_pos == snp_pos:
                                        ## check if this SNP matches the input SNP ID
                                        if orig_snp_id != this_rsID:
                                            logging_function("SNP %s does not match 1000 genomes rsID. Using ID %s for this input SNP" % (orig_snp_id, this_rsID))
                                            ## set the search rsID for PLINK
                                            snp_match_info[snp_key][0] = this_rsID
                                            ## in this case, we also want to store the alleles
                                            ## using the input SNP id
                                            orig_allele_key = ":".join([orig_snp_id, str(this_pos)])
                                            if orig_allele_key not in allele_dict:
                                                allele_dict[orig_allele_key] = line_data[3:5]
                                            
                                        ## set the boolean that this SNP was found
                                        snp_match_info[snp_key][1] = True
                                    if ld_check_area != 0 or ld_threshold != 1.0:
                                        ## write the line to the neighboring SNP file
                                        neighbor_snp_files[snp_key].write("\t".join(line_data)+"\n")
                                        
                                # if we've passed the end position, stop looking at
                                # this SNP (mark it for removal and close the file)
                                elif this_pos > end_pos:
                                    removed_snps.append(snp_key)
                                    if ld_check_area != 0 or ld_threshold != 1.0:
                                        neighbor_snp_files[snp_key].close()
                                        
                            ## now remove the SNPs we marked as having passed                    
                            [cur_search_snps.remove(x) for x in removed_snps]
                            ## if we closed the last file, stop looking
                            if len(cur_search_snps) == 0:
                                break
                                
                ## after looping through the reference, close any files that may
                ## happen to be open still (will only happen if the window around a
                ## SNP extends past the chromosome end)
                if ld_check_area != 0 or ld_threshold != 1.0:
                    for snp_key in cur_search_snps:
                        neighbor_snp_files[snp_key].close()

            ## now, we can run plink for each SNP. note that we are now running
            ## this loop on all combinations of SNPs and regions, but if it was
            ## already run, we just skip the LD expansion
            for snp_key in this_chr_snps:
                ## get the information for this snp
                snp_data = this_chr_snps[snp_key]
                snp_chr = snp_data[0]
                ## this is the input rsID
                snp_id = snp_data[1]
                snp_name = snp_data[2]
                snp_pos = str(snp_data[3])
                allele_key = ":".join([snp_id, str(snp_pos)])

                this_snp_info = snp_match_info[snp_key]
                search_rsID = this_snp_info[0]
                this_snp_found = this_snp_info[1]
                neighboring_snps_file = this_snp_info[2]

                ## if we are skipping the LD expansion, just write it out directly
                ## for now, skip the MAF calculation
                if ld_check_area == 0 and ld_threshold == 1.0:
                    if this_snp_found:
                        outf.write("\t".join([snp_chr, snp_id, snp_pos]+allele_dict[allele_key]+
                                             ["NA", snp_id, snp_pos, "NA", snp_name, "1.0", "1.0"])+"\n")
                    else:
                        logging_function("SNP %s not found in 1,000 genomes by ID or by position. No allele information is being written for this SNP!" % (snp_id))
                        outf.write("\t".join([snp_chr, snp_id, snp_pos, "NA", "NA",
                                              "NA", snp_id, snp_pos, "NA",
                                              snp_name, "1.0", "1.0"])+"\n")
                else:
                    ## convert the LD checking area to kb
                    ld_kb_check_area = ld_check_area / 1000
                    ## output prefix
                    plink_prefix = (outdir+"/ld_expansion/snp_data/"+allele_key+"/"+
                                    "_".join([allele_key, str(ld_check_area), "bp", str(ld_threshold), "cutoff", "snps"]))
                    ## only run analysis if we found the SNP and haven't done this
                    ## analysis yet (even if we are forcing recalculation, this is
                    ## because of potential duplicated IDs and SNPs belonging to more
                    ## than one tag region):
                    if not plink_path and this_snp_found and not os.path.isfile(plink_prefix+".ld"):
                        with open(plink_prefix+"_pipeline.log", 'wb') as plink_log:
                            subprocess.call(["plink", "--allow-no-sex", "--vcf", neighboring_snps_file,
                                            "--r2", "with-freqs", "dprime",
                                            "--ld-snp", search_rsID,
                                            "--ld-window", "99999",
                                            "--ld-window-kb", str(ld_kb_check_area),
                                            "--ld-window-r2",
                                            str(ld_threshold), 
                                            "--out", plink_prefix],
                                            stdout=plink_log)
                    elif plink_path and this_snp_found and not os.path.isfile(plink_prefix+".ld"):
                        with open(plink_prefix+"_pipeline.log", 'wb') as plink_log:
                            subprocess.call([plink_path, "--allow-no-sex", "--vcf", neighboring_snps_file,
                                            "--r2", "with-freqs", "dprime",
                                            "--ld-snp", search_rsID,
                                            "--ld-window", "99999",
                                            "--ld-window-kb", str(ld_kb_check_area),
                                            "--ld-window-r2",
                                            str(ld_threshold),
                                            "--out", plink_prefix],
                                            stdout=plink_log)
                    elif os.path.isfile(plink_prefix+".ld"):
                        logging_function("No need to recompute LD SNPs for SNP ID %s" % (snp_id))
                    elif not this_snp_found:
                        logging_function("SNP %s not found in 1,000 genomes by ID or by position. Skipping!" % (snp_id))

                    ## now output to the master file, if it exists and we found the SNP
                    if os.path.isfile(plink_prefix+".ld") and this_snp_found:
                        with open(plink_prefix+".ld", 'rb') as plinkout:
                            # read in the header, make an index
                            plink_header = re.sub("\s+", "\t", next(plinkout).strip()).split("\t")
                            plink_idx = {plink_header[x]:x for x in range(len(plink_header))}
                            for ld_line in plinkout:
                                ld_data = re.sub("\s+", "\t", ld_line.strip()).split("\t")
                                ## check if we found the search rsID (replace with tag snp ID)
                                if ld_data[plink_idx['SNP_B']]==search_rsID:
                                    outf.write("\t".join(["chr"+ld_data[plink_idx['CHR_A']], snp_id,
                                                          ld_data[plink_idx['BP_B']]]+
                                                          allele_dict[allele_key]+
                                                          [ld_data[plink_idx['MAF_B']],
                                                           snp_id, ld_data[plink_idx['BP_A']],
                                                           ld_data[plink_idx['MAF_A']],
                                                           snp_name, ld_data[plink_idx['R2']],
                                                           ld_data[plink_idx['DP']]])+'\n')
                                else:
                                    ld_snp_key = ":".join([ld_data[plink_idx['SNP_B']],
                                                           ld_data[plink_idx['BP_B']]])
                                    if ld_snp_key in allele_dict:
                                        outf.write("\t".join(["chr"+ld_data[plink_idx['CHR_A']],
                                                          ld_data[plink_idx['SNP_B']],
                                                          ld_data[plink_idx['BP_B']]]+
                                                          allele_dict[ld_snp_key]+
                                                          [ld_data[plink_idx['MAF_B']],
                                                           snp_id, ld_data[plink_idx['BP_A']],
                                                           ld_data[plink_idx['MAF_A']],
                                                           snp_name, ld_data[plink_idx['R2']],
                                                           ld_data[plink_idx['DP']]])+'\n')
                                    else:
                                        logging_function("SNP %s not in allele dict" % (ld_snp_key))
                                        outf.write("\t".join(["chr"+ld_data[plink_idx['CHR_A']],
                                                          ld_data[plink_idx['SNP_B']],
                                                          ld_data[plink_idx['BP_B']]]+
                                                          ["NA", "NA"]+
                                                          [ld_data[plink_idx['MAF_B']],
                                                           snp_id, ld_data[plink_idx['BP_A']],
                                                           ld_data[plink_idx['MAF_A']],
                                                           snp_name, ld_data[plink_idx['R2']],
                                                           ld_data[plink_idx['DP']]])+'\n')
                                        
                    else:
                        logging_function("Plink failed to generate an LD file for snp %s (search ID %s)!" % (snp_id, search_rsID))

        ## sort the output (chromosome, position, rsID, region name, tag rsID)
        ## first get the sort strings
        chr_sort = "-k"+",".join([str(outf_idx['chr']+1)]*2)+"V"
        pos_sort = "-k"+",".join([str(outf_idx['pos']+1)]*2)+"n"
        rsid_sort = "-k"+",".join([str(outf_idx['rsID']+1)]*2)
        name_sort = "-k"+",".join([str(outf_idx['tag_name']+1)]*2)
        tag_rsid_sort = "-k"+",".join([str(outf_idx['tag_rsID']+1)]*2)

        with open(outfile_path+".sorted", 'wb') as sort_out:
            subprocess.call(["sort", chr_sort, pos_sort, rsid_sort, name_sort, tag_rsid_sort, outfile_path], stdout=sort_out)
        subprocess.call(["mv", outfile_path+".sorted", outfile_path])

        ## generate a bed file
        try:
            os.makedirs(outdir+"/ld_expansion/bed_files/")
        except OSError:
            pass

        ## generate a bed file from the LD SNPs
        snp_bed_prefix = outdir+"/ld_expansion/bed_files/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)
        with open(snp_bed_prefix+".tmp", 'wb') as bedout, open(outfile_path, 'rb') as snpin:
            # get the header
            snp_header_data = next(snpin).strip().split('\t')
            snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

            for snpline in snpin:
                this_data = snpline.strip().split('\t')
                end_pos = int(this_data[snp_idx['pos']])+1
                bedout.write("\t".join([this_data[snp_idx['chr']],
                                        this_data[snp_idx['pos']],
                                        str(end_pos),
                                        ":".join([this_data[snp_idx['rsID']],
                                                  this_data[snp_idx['tag_rsID']],
                                                  this_data[snp_idx['tag_name']]])])+'\n')

        ## sort this file (note that we don't use the 'V' flag for the chromosome sort because this isn't compatible with bedtools..)
        with open(snp_bed_prefix+".bed", "wb") as sort_out:
            subprocess.call(['sort', '-k1,1', '-k2,2n', snp_bed_prefix+".tmp"], stdout=sort_out)
        subprocess.call(['rm', snp_bed_prefix+".tmp"])

        ## if we want to, get rid of all the SNP data
        if remove_ld_data and (ld_check_area != 0 or ld_threshold != 1.0):
            logging_function("Removing LD data files")
            neighbor_snp_files = glob.glob(outdir+"/ld_expansion/snp_data/*/*_"+str(ld_check_area)+"_*")
            for snp_f in neighbor_snp_files:
                subprocess.call(["rm", snp_f])
                ## if the directory is empty, remove it
                try:
                    os.rmdir(os.path.dirname(snp_f))
                except OSError as ex:
                    continue
            ## try deleting the snp data directory, if it's not empty
            try:
                os.rmdir(outdir+"/ld_expansion/snp_data/")
            except OSError as ex:
                logging_function("Could not delete snp data folder!")
    else:
        logging_function("Using existing file: "+outfile_path)
        logging_function("To force recalculation, give the --force_ld_recalc flag to the Python script.")

    return outfile_path

def compute_closest_genes(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, gene_bed_file, kgxref_file, bedtools_bin_dir, code_dir, logging_function):
    """
    uses bedtools to compute closest genes for each SNP
    """
    try:
        os.makedirs(outdir+"/closest_gene/")
    except OSError:
        pass

    ## use the LD snp bed file:
    ld_snp_bed = outdir+"/ld_expansion/bed_files/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)+".bed"

    bedout_prefix = outdir+"/closest_gene/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)

    ## then use the bedtools wrapper script to run bedtools intersect:
    ## parameters: -t all means report all features with the same distance
    ## -d says to report the distance to A
    snp_bedtools_outf = bedout_prefix+"_bedtools_output.txt"
    with open(snp_bedtools_outf, 'wb') as bedtools_out:
        ## TODO: add bedtools version output
        ## convert the bin dir to a string so that if it's NoneType, you don't get an error
        bed_call = subprocess.Popen([code_dir+"/bedtools_wrapper_script.sh", str(bedtools_bin_dir),
                         "closest", "-a", ld_snp_bed, "-b", gene_bed_file,
                         "-t", "all", "-d"], stdout=bedtools_out)
    ## use the wrapper script to parse it
    parse_call = subprocess.Popen([code_dir+"/parse_closest_gene_output.sh", snp_bedtools_outf])
    
    ## if we have cross-reference, read that into a dict
    if kgxref_file:
        with open(kgxref_file, 'rb') as xref:
            xref_header_data = next(xref).strip().split("\t")
            xref_header_dict = {xref_header_data[x]:x for x in range(len(xref_header_data))}
            xref_dict = {}
            for xref_line in xref:
                xref_data = xref_line.strip().split('\t')
                refseq = xref_data[xref_header_dict['refseq']]
                if refseq:
                    xref_dict[refseq] = xref_data[xref_header_dict['geneSymbol']]

    ## now combine this with the LD SNP list (in a new file)
    closest_genes_outf = bedout_prefix+"_closest_genes.txt"
    with open(snp_bedtools_outf, 'rb') as bed_out, open(ld_snp_file, 'rb') as snpin, open(closest_genes_outf, 'wb') as parsed_outf:
        # get the header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

        ## write the header:
        if kgxref_file:
            parsed_outf.write("\t".join(snp_header_data+
                                        ["closest_gene_ids", "closest_gene_symbols", "distance"])+'\n')
        else:
            parsed_outf.write("\t".join(snp_header_data+
                                        ["closest_gene_ids", "distance"])+'\n')

        ## start reading the overlap file
        bedline = bed_out.readline()
        if bedline:
            bed_data = bedline.strip().split('\t')
        else:
            logging_function("Bed file is empty!")
            ## in this case, just return so we don't get an error
            return None

        ## define dictionary for the closest_gene file:
        ## original SNP position is bed file start, SNP:LD SNP:tag region ID is fourth column,
        ## the distance output is the very last column, and gene is 8th column (4th in bed file)
        ## this should be robust to the type of bed file you get for the gene data
        closest_gene_dict = {'pos':1, 'snps':3, 'gene': 7, 'dist':len(bed_data)-1}

        ## track which SNP we're on in the closest gene data
        ## uses all pieces of information (rsID, LD SNP rsID, and tag region)
        tag_snp_data = bed_data[closest_gene_dict['snps']].split(":")
        closest_gene_rsID = tag_snp_data[0]
        closest_gene_tag_rsID = tag_snp_data[1]
        closest_gene_tag_name = ":".join(tag_snp_data[2:])

        closest_gene_pos = bed_data[closest_gene_dict['pos']]

        ## initialize variables to track the genes of interest
        this_snp_genes = []
        ## track the distances (for now, also track all of these, to double check)
        this_snp_dist = []

        ## loop through LD SNP list, collecting all entries for each one
        for snpline in snpin:
            this_snp_data = snpline.strip().split('\t')
            this_snp_rsID = this_snp_data[snp_idx['rsID']]
            this_snp_tag_rsID = this_snp_data[snp_idx['tag_rsID']]
            this_snp_pos = this_snp_data[snp_idx['pos']]
            this_snp_tag_name = this_snp_data[snp_idx['tag_name']]

            ## check for mismatches (this means we've missed something)
            if (this_snp_rsID != closest_gene_rsID or
                this_snp_tag_rsID != closest_gene_tag_rsID or
                this_snp_tag_name != closest_gene_tag_name or
                this_snp_pos != closest_gene_pos):
                logging_function("Closest gene and LD files have mismatch! Potential missed SNP")
                logging_function("LD SNP: %s, %s, %s. Closest_Gene: %s, %s, %s" % (this_snp_rsID, this_snp_tag_rsID, this_snp_tag_name, closest_gene_rsID, closest_gene_tag_rsID, closest_gene_tag_name))
                logging_function("SNP position: %s. Closest_Gene position: %s\n" % (this_snp_pos, closest_gene_pos))

            ## match by position and rsID of the LD SNP as well as rsID of the tagging SNP
            ## as long as we match and have more in the closest gene file, keep reading in genes 
            while (this_snp_pos == closest_gene_pos and
                   this_snp_tag_rsID == closest_gene_tag_rsID and
                   this_snp_rsID == closest_gene_rsID and
                   this_snp_tag_name == closest_gene_tag_name and
                   bedline):
                ## add the gene and distance
                this_snp_genes.append(bed_data[closest_gene_dict['gene']])
                this_snp_dist.append(bed_data[closest_gene_dict['dist']])

                ## read in the next line from the closest gene file
                bedline = bed_out.readline()
                if bedline:
                    bed_data = bedline.strip().split("\t")
                    tag_snp_data = bed_data[closest_gene_dict['snps']].split(":")
                    closest_gene_rsID = tag_snp_data[0]
                    closest_gene_tag_rsID = tag_snp_data[1]
                    closest_gene_tag_name = ":".join(tag_snp_data[2:])
                    closest_gene_pos = bed_data[closest_gene_dict['pos']]
                else:
                    # this shouldn't happen
                    closest_gene_rsID = 'DoneIterating'
                    closest_gene_tag_rsID = 'DoneITerating'
                    closest_gene_pos = -1
                    break

            ## we are now either on a different intersection SNP, or are out of intersects
            ## either way, we want to write this data:
            if kgxref_file:
                ## get the genes                
                this_snp_symbols = [xref_dict[x] if x in xref_dict else 'NotFound'
                                    for x in this_snp_genes]
                parsed_outf.write("\t".join(this_snp_data)+"\t"+",".join(this_snp_genes)+"\t"+
                                  ",".join(this_snp_symbols)+"\t" +
                                  ",".join(this_snp_dist)+"\n")
                
            else:
                parsed_outf.write("\t".join(this_snp_data)+"\t"+",".join(this_snp_genes)+"\t"+
                                ",".join(this_snp_dist)+"\n")
            ## reset the tracking variables
            this_snp_genes = []
            this_snp_dist = []
            
    return closest_genes_outf

# partition_dir = "/home/alexaml/data/refgenomes/hg19/hg19_utr_refseq/resorted_full_utrs/"
def calculate_genomic_partition(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, partition_dir, logging_function):
    """
    determines what kind of genomic element each SNP falls in
    split by positive and negative strands
    """
    try:
        os.makedirs(outdir+"/genomic_partition/")
    except OSError:
      pass

    ## get all the file paths
    pos_5utr_exons_f = partition_dir+"/pos_files/final_files/parsed_pos_5utr_exons.merged.bed"
    pos_5utr_introns_f = partition_dir+"/pos_files/final_files/pos_n5e_5utr_introns.bed"
    pos_3utr_exons_f = partition_dir+"/pos_files/final_files/pos_n5e5i_3utr_exons.bed"
    pos_3utr_introns_f = partition_dir+"/pos_files/final_files/pos_n5e5i3e_3utr_introns.bed"
    pos_promoters_f = partition_dir+"/pos_files/final_files/pos_n5e5i3e3i_promoters.bed"
    pos_exons_f = partition_dir+"/pos_files/final_files/pos_n5e5i3e3ip_exons.bed"
    pos_introns_f = partition_dir+"/pos_files/final_files/pos_n5e5i3e3ipe_introns.bed"
    pos_repeats_f = partition_dir+"/pos_files/final_files/pos_n5e5i3e3ipei_repeats.bed"

    neg_5utr_exons_f = partition_dir+"/neg_files/final_files/parsed_neg_5utr_exons.merged.bed"
    neg_5utr_introns_f = partition_dir+"/neg_files/final_files/neg_n5e_5utr_introns.bed"
    neg_3utr_exons_f = partition_dir+"/neg_files/final_files/neg_n5e5i_3utr_exons.bed"
    neg_3utr_introns_f = partition_dir+"/neg_files/final_files/neg_n5e5i3e_3utr_introns.bed"
    neg_promoters_f = partition_dir+"/neg_files/final_files/neg_n5e5i3e3i_promoters.bed"
    neg_exons_f = partition_dir+"/neg_files/final_files/neg_n5e5i3e3ip_exons.bed"
    neg_introns_f = partition_dir+"/neg_files/final_files/neg_n5e5i3e3ipe_introns.bed"
    neg_repeats_f = partition_dir+"/neg_files/final_files/neg_n5e5i3e3ipei_repeats.bed"

    partition_outprefix = outdir+"/genomic_partition/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)
    
    ## write a single output file
    ## contains each entry with its positive and negative strand partition
    with open(ld_snp_file, 'rb') as snpin, open(pos_5utr_exons_f, 'rb') as pos_fp_exons, open(pos_5utr_introns_f, 'rb') as pos_fp_introns, open(pos_3utr_exons_f, 'rb') as pos_tp_exons, open(pos_3utr_introns_f, 'rb') as pos_tp_introns, open(pos_promoters_f, 'rb') as pos_promoters, open(pos_exons_f, 'rb') as pos_exons, open(pos_introns_f, 'rb') as pos_introns, open(pos_repeats_f, 'rb') as pos_repeats, open(neg_5utr_exons_f, 'rb') as neg_fp_exons, open(neg_5utr_introns_f, 'rb') as neg_fp_introns, open(neg_3utr_exons_f, 'rb') as neg_tp_exons, open(neg_3utr_introns_f, 'rb') as neg_tp_introns, open(neg_promoters_f, 'rb') as neg_promoters, open(neg_exons_f, 'rb') as neg_exons, open(neg_introns_f, 'rb') as neg_introns, open(neg_repeats_f, 'rb') as neg_repeats, open(partition_outprefix+"_entrywise_partition.txt", 'wb') as entry_out:
        # get the LD SNP header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

        # define a convenience dictionary to make subsetting the bed files clearer
        bed_idx = {} 
        bed_idx["chr"] = 0
        bed_idx["start"] = 1
        bed_idx["end"] = 2    
        bed_idx["name"] = 3
        
        ## write the output headers
        entry_out.write("\t".join(snp_header_data + ['pos_partition', 'neg_partition'])+'\n')

        ## store lists of each type of element
        cur_pos_fp_exons = []
        cur_pos_fp_introns = []
        cur_pos_tp_exons = []
        cur_pos_tp_introns = []        
        cur_pos_promoters = []
        cur_pos_exons = []
        cur_pos_introns = []
        cur_pos_repeats = []
        
        cur_neg_fp_exons = []
        cur_neg_fp_introns = []
        cur_neg_tp_exons = []
        cur_neg_tp_introns = []        
        cur_neg_promoters = []
        cur_neg_exons = []
        cur_neg_introns = []
        cur_neg_repeats = []                
        
        ## start reading in these elements
        cur_pos_fp_exon = pos_fp_exons.readline().strip().split('\t')
        cur_pos_fp_intron = pos_fp_introns.readline().strip().split('\t')        
        cur_pos_tp_exon = pos_tp_exons.readline().strip().split('\t')
        cur_pos_tp_intron = pos_tp_introns.readline().strip().split('\t')        
        cur_pos_promoter = pos_promoters.readline().strip().split('\t')
        cur_pos_exon = pos_exons.readline().strip().split('\t')
        cur_pos_intron = pos_introns.readline().strip().split('\t')
        cur_pos_repeat = pos_repeats.readline().strip().split('\t')

        cur_neg_fp_exon = neg_fp_exons.readline().strip().split('\t')
        cur_neg_fp_intron = neg_fp_introns.readline().strip().split('\t')        
        cur_neg_tp_exon = neg_tp_exons.readline().strip().split('\t')
        cur_neg_tp_intron = neg_tp_introns.readline().strip().split('\t')        
        cur_neg_promoter = neg_promoters.readline().strip().split('\t')
        cur_neg_exon = neg_exons.readline().strip().split('\t')
        cur_neg_intron = neg_introns.readline().strip().split('\t')
        cur_neg_repeat = neg_repeats.readline().strip().split('\t')

        ## track which chromosome we're on
        this_chr = ''

        ## a regex to check that we are in a 'real' chromosome
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        

        ## track how many times each type of element was found
        element_list = ['fp_exon', 'fp_intron', 'tp_exon', 'tp_intron',
                        'promoter', 'exon', 'intron', 'repeat', 'intergenic']
        element_counts = {}
        for x in element_list:
            element_counts['pos_'+x] = 0.0
            element_counts['neg_'+x] = 0.0

        entry_ctr = 0.0
        for entry in snpin:
            entry_ctr += 1.0
            entry_data = entry.strip().split('\t')
            entry_pos = int(entry_data[snp_idx['pos']])

            this_rsID = entry_data[snp_idx['rsID']]
            this_tag_rsID = entry_data[snp_idx['tag_rsID']]
            
            ## while we're on a new chromosome, reset all the tracking variables:
            if entry_data[snp_idx['chr']] != this_chr:
                this_chr = entry_data[snp_idx['chr']]
                ## store lists of each type of element
                cur_pos_fp_exons = []
                cur_pos_fp_introns = []
                cur_pos_tp_exons = []
                cur_pos_tp_introns = []
                cur_pos_promoters = []
                cur_pos_exons = []
                cur_pos_introns = []
                cur_pos_repeats = []
                
                cur_neg_fp_exons = []
                cur_neg_fp_introns = []
                cur_neg_tp_exons = []
                cur_neg_tp_introns = []        
                cur_neg_promoters = []
                cur_neg_exons = []
                cur_neg_introns = []
                cur_neg_repeats = []
                
                ## read in data until we match chromosomes
                while (len(cur_pos_fp_exon) > 1
                    and cur_pos_fp_exon[bed_idx['chr']] != this_chr):
                    cur_pos_fp_exon = pos_fp_exons.readline().strip().split('\t')
                while (len(cur_pos_fp_intron) > 1
                       and cur_pos_fp_intron[bed_idx['chr']] != this_chr):
                    cur_pos_fp_intron = pos_fp_introns.readline().strip().split('\t')
                while (len(cur_pos_tp_exon) > 1
                       and cur_pos_tp_exon[bed_idx['chr']] != this_chr):
                    cur_pos_tp_exon = pos_tp_exons.readline().strip().split('\t')
                while (len(cur_pos_tp_intron) > 1
                       and cur_pos_tp_intron[bed_idx['chr']] != this_chr):
                    cur_pos_tp_intron = pos_tp_introns.readline().strip().split('\t')
                while (len(cur_pos_promoter) > 1
                       and cur_pos_promoter[bed_idx['chr']] != this_chr):
                    cur_pos_promoter = pos_promoters.readline().strip().split('\t')
                while (len(cur_pos_exon) > 1
                       and cur_pos_exon[bed_idx['chr']] != this_chr):
                    cur_pos_exon = pos_exons.readline().strip().split('\t')
                while (len(cur_pos_intron) > 1
                       and cur_pos_intron[bed_idx['chr']] != this_chr):
                    cur_pos_intron = pos_introns.readline().strip().split('\t')
                while (len(cur_pos_repeat) > 1
                       and cur_pos_repeat[bed_idx['chr']] != this_chr):
                    cur_pos_repeat = pos_repeats.readline().strip().split('\t')
                    
                while (len(cur_neg_fp_exon) > 1
                    and cur_neg_fp_exon[bed_idx['chr']] != this_chr):
                    cur_neg_fp_exon = neg_fp_exons.readline().strip().split('\t')
                while (len(cur_neg_fp_intron) > 1
                       and cur_neg_fp_intron[bed_idx['chr']] != this_chr):
                    cur_neg_fp_intron = neg_fp_introns.readline().strip().split('\t')
                while (len(cur_neg_tp_exon) > 1
                       and cur_neg_tp_exon[bed_idx['chr']] != this_chr):
                    cur_neg_tp_exon = neg_tp_exons.readline().strip().split('\t')
                while (len(cur_neg_tp_intron) > 1
                       and cur_neg_tp_intron[bed_idx['chr']] != this_chr):
                    cur_neg_tp_intron = neg_tp_introns.readline().strip().split('\t')
                while (len(cur_neg_promoter) > 1
                       and cur_neg_promoter[bed_idx['chr']] != this_chr):
                    cur_neg_promoter = neg_promoters.readline().strip().split('\t')
                while (len(cur_neg_exon) > 1
                       and cur_neg_exon[bed_idx['chr']] != this_chr):
                    cur_neg_exon = neg_exons.readline().strip().split('\t')
                while (len(cur_neg_intron) > 1
                       and cur_neg_intron[bed_idx['chr']] != this_chr):
                    cur_neg_intron = neg_introns.readline().strip().split('\t')
                while (len(cur_neg_repeat) > 1
                       and cur_neg_repeat[bed_idx['chr']] != this_chr):
                    cur_neg_repeat = neg_repeats.readline().strip().split('\t')
                logging_function('Parsing chromosome '+this_chr)
                
            ## if we have a match to one of the weird chromosomes, just skip to the next entry
            if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                continue                    

            ## add all possibly overlapping entries for each type of element
            ## must have same chromosome and start before the end of the entry
            ## note that the bed files are 0-based, half-open so the start can be equal to the
            ## SNP position but the end must be strictly greater than
            ## positive strand:
            while (len(cur_pos_fp_exon) > 1
                   and cur_pos_fp_exon[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_pos_fp_exon[bed_idx['start']]) <= entry_pos):
                cur_pos_fp_exons.append(cur_pos_fp_exon)
                cur_pos_fp_exon = pos_fp_exons.readline().strip().split('\t')                
            while (len(cur_pos_fp_intron) > 1
                   and cur_pos_fp_intron[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_pos_fp_intron[bed_idx['start']]) <= entry_pos):
                cur_pos_fp_introns.append(cur_pos_fp_intron)
                cur_pos_fp_intron = pos_fp_introns.readline().strip().split('\t')
            while (len(cur_pos_tp_exon) > 1
                   and cur_pos_tp_exon[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_pos_tp_exon[bed_idx['start']]) <= entry_pos):
                cur_pos_tp_exons.append(cur_pos_tp_exon)
                cur_pos_tp_exon = pos_tp_exons.readline().strip().split('\t')                
            while (len(cur_pos_tp_intron) > 1
                   and cur_pos_tp_intron[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_pos_tp_intron[bed_idx['start']]) <= entry_pos):
                cur_pos_tp_introns.append(cur_pos_tp_intron)
                cur_pos_tp_intron = pos_tp_introns.readline().strip().split('\t')
            while (len(cur_pos_promoter) > 1
                   and cur_pos_promoter[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_pos_promoter[bed_idx['start']]) <= entry_pos):
                cur_pos_promoters.append(cur_pos_promoter)
                cur_pos_promoter = pos_promoters.readline().strip().split('\t')                
            while (len(cur_pos_exon) > 1
                   and cur_pos_exon[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_pos_exon[bed_idx['start']]) <= entry_pos):
                cur_pos_exons.append(cur_pos_exon)
                cur_pos_exon = pos_exons.readline().strip().split('\t')
            while (len(cur_pos_intron) > 1
                   and cur_pos_intron[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_pos_intron[bed_idx['start']]) <= entry_pos):
                cur_pos_introns.append(cur_pos_intron)
                cur_pos_intron = pos_introns.readline().strip().split('\t')            
            while (len(cur_pos_repeat) > 1
                   and cur_pos_repeat[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_pos_repeat[bed_idx['start']]) <= entry_pos):
                cur_pos_repeats.append(cur_pos_repeat)
                cur_pos_repeat = pos_repeats.readline().strip().split('\t')

            ## negative strand:
            while (len(cur_neg_fp_exon) > 1
                   and cur_neg_fp_exon[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_neg_fp_exon[bed_idx['start']]) <= entry_pos):
                cur_neg_fp_exons.append(cur_neg_fp_exon)
                cur_neg_fp_exon = neg_fp_exons.readline().strip().split('\t')                
            while (len(cur_neg_fp_intron) > 1
                   and cur_neg_fp_intron[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_neg_fp_intron[bed_idx['start']]) <= entry_pos):
                cur_neg_fp_introns.append(cur_neg_fp_intron)
                cur_neg_fp_intron = neg_fp_introns.readline().strip().split('\t')
            while (len(cur_neg_tp_exon) > 1
                   and cur_neg_tp_exon[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_neg_tp_exon[bed_idx['start']]) <= entry_pos):
                cur_neg_tp_exons.append(cur_neg_tp_exon)
                cur_neg_tp_exon = neg_tp_exons.readline().strip().split('\t')                
            while (len(cur_neg_tp_intron) > 1
                   and cur_neg_tp_intron[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_neg_tp_intron[bed_idx['start']]) <= entry_pos):
                cur_neg_tp_introns.append(cur_neg_tp_intron)
                cur_neg_tp_intron = neg_tp_introns.readline().strip().split('\t')
            while (len(cur_neg_promoter) > 1
                   and cur_neg_promoter[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_neg_promoter[bed_idx['start']]) <= entry_pos):
                cur_neg_promoters.append(cur_neg_promoter)
                cur_neg_promoter = neg_promoters.readline().strip().split('\t')                
            while (len(cur_neg_exon) > 1
                   and cur_neg_exon[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_neg_exon[bed_idx['start']]) <= entry_pos):
                cur_neg_exons.append(cur_neg_exon)
                cur_neg_exon = neg_exons.readline().strip().split('\t')
            while (len(cur_neg_intron) > 1
                   and cur_neg_intron[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_neg_intron[bed_idx['start']]) <= entry_pos):
                cur_neg_introns.append(cur_neg_intron)
                cur_neg_intron = neg_introns.readline().strip().split('\t')            
            while (len(cur_neg_repeat) > 1
                   and cur_neg_repeat[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_neg_repeat[bed_idx['start']]) <= entry_pos):
                cur_neg_repeats.append(cur_neg_repeat)
                cur_neg_repeat = neg_repeats.readline().strip().split('\t')
                
            ## remove entries that are before the current entry here:
            for i in range(len(cur_pos_fp_exons)-1, -1, -1):
                if cur_pos_fp_exons[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_pos_fp_exons[i][bed_idx['end']]) <= entry_pos:
                    cur_pos_fp_exons.pop(i)                    
            for i in range(len(cur_pos_fp_introns)-1, -1, -1):
                if cur_pos_fp_introns[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_pos_fp_introns[i][bed_idx['end']]) <= entry_pos:
                    cur_pos_fp_introns.pop(i)
            for i in range(len(cur_pos_tp_exons)-1, -1, -1):
                if cur_pos_tp_exons[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_pos_tp_exons[i][bed_idx['end']]) <= entry_pos:
                    cur_pos_tp_exons.pop(i)                    
            for i in range(len(cur_pos_tp_introns)-1, -1, -1):
                if cur_pos_tp_introns[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_pos_tp_introns[i][bed_idx['end']]) <= entry_pos:
                    cur_pos_tp_introns.pop(i)            
            for i in range(len(cur_pos_promoters)-1, -1, -1):
                if cur_pos_promoters[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_pos_promoters[i][bed_idx['end']]) <= entry_pos:
                    cur_pos_promoters.pop(i)
            for i in range(len(cur_pos_exons)-1, -1, -1):
                if cur_pos_exons[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_pos_exons[i][bed_idx['end']]) <= entry_pos:
                    cur_pos_exons.pop(i)
            for i in range(len(cur_pos_introns)-1, -1, -1):
                if cur_pos_introns[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_pos_introns[i][bed_idx['end']]) <= entry_pos:
                    cur_pos_introns.pop(i)
            for i in range(len(cur_pos_repeats)-1, -1, -1):
                if cur_pos_repeats[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_pos_repeats[i][bed_idx['end']]) <= entry_pos:
                    cur_pos_repeats.pop(i)                                        

            for i in range(len(cur_neg_fp_exons)-1, -1, -1):
                if cur_neg_fp_exons[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_neg_fp_exons[i][bed_idx['end']]) <= entry_pos:
                    cur_neg_fp_exons.pop(i)                    
            for i in range(len(cur_neg_fp_introns)-1, -1, -1):
                if cur_neg_fp_introns[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_neg_fp_introns[i][bed_idx['end']]) <= entry_pos:
                    cur_neg_fp_introns.pop(i)
            for i in range(len(cur_neg_tp_exons)-1, -1, -1):
                if cur_neg_tp_exons[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_neg_tp_exons[i][bed_idx['end']]) <= entry_pos:
                    cur_neg_tp_exons.pop(i)                    
            for i in range(len(cur_neg_tp_introns)-1, -1, -1):
                if cur_neg_tp_introns[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_neg_tp_introns[i][bed_idx['end']]) <= entry_pos:
                    cur_neg_tp_introns.pop(i)            
            for i in range(len(cur_neg_promoters)-1, -1, -1):
                if cur_neg_promoters[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_neg_promoters[i][bed_idx['end']]) <= entry_pos:
                    cur_neg_promoters.pop(i)
            for i in range(len(cur_neg_exons)-1, -1, -1):
                if cur_neg_exons[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_neg_exons[i][bed_idx['end']]) <= entry_pos:
                    cur_neg_exons.pop(i)
            for i in range(len(cur_neg_introns)-1, -1, -1):
                if cur_neg_introns[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_neg_introns[i][bed_idx['end']]) <= entry_pos:
                    cur_neg_introns.pop(i)
            for i in range(len(cur_neg_repeats)-1, -1, -1):
                if cur_neg_repeats[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_neg_repeats[i][bed_idx['end']]) <= entry_pos:
                    cur_neg_repeats.pop(i)                                        

            ## now compare the SNP (current entry) to each of the genomic loci
            ## track overlaps with all kinds of elements just in case they aren't mutually exclusive
            this_pos_fp_exon_bp = 0.0
            this_pos_fp_intron_bp = 0.0
            this_pos_tp_exon_bp = 0.0
            this_pos_tp_intron_bp = 0.0            
            this_pos_promoter_bp = 0.0
            this_pos_exon_bp = 0.0
            this_pos_intron_bp = 0.0
            this_pos_repeat_bp = 0.0

            this_neg_fp_exon_bp = 0.0
            this_neg_fp_intron_bp = 0.0
            this_neg_tp_exon_bp = 0.0
            this_neg_tp_intron_bp = 0.0
            this_neg_promoter_bp = 0.0
            this_neg_exon_bp = 0.0
            this_neg_intron_bp = 0.0
            this_neg_repeat_bp = 0.0

            ## calculate overlaps
            ## since we have positional SNPs, we can just ask if the SNP falls in the region
            for cfe in cur_pos_fp_exons:
                this_pos_fp_exon_bp += 1 if (int(cfe[bed_idx['start']]) <= entry_pos < int(cfe[bed_idx['end']])) else 0
            for cfi in cur_pos_fp_introns:
                this_pos_fp_intron_bp += 1 if (int(cfi[bed_idx['start']]) <= entry_pos < int(cfi[bed_idx['end']])) else 0
            for cte in cur_pos_tp_exons:
                this_pos_tp_exon_bp += 1 if (int(cte[bed_idx['start']]) <= entry_pos < int(cte[bed_idx['end']])) else 0
            for cti in cur_pos_tp_introns:
                this_pos_tp_intron_bp += 1 if (int(cti[bed_idx['start']]) <= entry_pos < int(cti[bed_idx['end']])) else 0
            for cp in cur_pos_promoters:
                this_pos_promoter_bp += 1 if (int(cp[bed_idx['start']]) <= entry_pos < int(cp[bed_idx['end']])) else 0
            for ce in cur_pos_exons:
                this_pos_exon_bp += 1 if (int(ce[bed_idx['start']]) <= entry_pos < int(ce[bed_idx['end']])) else 0
            for ci in cur_pos_introns:
                this_pos_intron_bp += 1 if (int(ci[bed_idx['start']]) <= entry_pos < int(ci[bed_idx['end']])) else 0
            for cr in cur_pos_repeats:
                this_pos_repeat_bp += 1 if (int(cr[bed_idx['start']]) <= entry_pos < int(cr[bed_idx['end']])) else 0

            for cfe in cur_neg_fp_exons:
                this_neg_fp_exon_bp += 1 if (int(cfe[bed_idx['start']]) <= entry_pos < int(cfe[bed_idx['end']])) else 0
            for cfi in cur_neg_fp_introns:
                this_neg_fp_intron_bp += 1 if (int(cfi[bed_idx['start']]) <= entry_pos < int(cfi[bed_idx['end']])) else 0
            for cte in cur_neg_tp_exons:
                this_neg_tp_exon_bp += 1 if (int(cte[bed_idx['start']]) <= entry_pos < int(cte[bed_idx['end']])) else 0
            for cti in cur_neg_tp_introns:
                this_neg_tp_intron_bp += 1 if (int(cti[bed_idx['start']]) <= entry_pos < int(cti[bed_idx['end']])) else 0
            for cp in cur_neg_promoters:
                this_neg_promoter_bp += 1 if (int(cp[bed_idx['start']]) <= entry_pos < int(cp[bed_idx['end']])) else 0
            for ce in cur_neg_exons:
                this_neg_exon_bp += 1 if (int(ce[bed_idx['start']]) <= entry_pos < int(ce[bed_idx['end']])) else 0
            for ci in cur_neg_introns:
                this_neg_intron_bp += 1 if (int(ci[bed_idx['start']]) <= entry_pos < int(ci[bed_idx['end']])) else 0
            for cr in cur_neg_repeats:
                this_neg_repeat_bp += 1 if (int(cr[bed_idx['start']]) <= entry_pos < int(cr[bed_idx['end']])) else 0
                
            ## go through and check for overlaps
            ## make list of overlapping elements
            pos_overlap_elements = []
            neg_overlap_elements = []
            
            if this_pos_fp_exon_bp > 0:
                pos_overlap_elements += ["fp_utr_exon"]
                element_counts['pos_fp_exon'] += 1.0
            if this_pos_fp_intron_bp > 0:
                pos_overlap_elements += ["fp_utr_intron"]
                element_counts['pos_fp_intron'] += 1.0
            if this_pos_tp_exon_bp > 0:
                pos_overlap_elements += ["tp_utr_exon"]
                element_counts['pos_tp_exon'] += 1.0
            if this_pos_tp_intron_bp > 0:
                pos_overlap_elements += ["tp_utr_intron"]
                element_counts['pos_tp_intron'] += 1.0
            if this_pos_promoter_bp > 0:
                pos_overlap_elements += ["promoter"]
                element_counts['pos_promoter'] += 1.0
            if this_pos_exon_bp > 0:
                pos_overlap_elements += ["exon"]
                element_counts['pos_exon'] += 1.0
            if this_pos_intron_bp > 0:
                pos_overlap_elements += ["intron"]
                element_counts['pos_intron'] += 1.0
            if this_pos_repeat_bp > 0:
                pos_overlap_elements += ["repeat"]
                element_counts['pos_repeat'] += 1.0
            if not pos_overlap_elements:
                pos_overlap_elements = ["intergenic"]
                element_counts['pos_intergenic'] += 1.0

            if this_neg_fp_exon_bp > 0:
                neg_overlap_elements += ["fp_utr_exon"]
                element_counts['neg_fp_exon'] += 1.0
            if this_neg_fp_intron_bp > 0:
                neg_overlap_elements += ["fp_utr_intron"]
                element_counts['neg_fp_intron'] += 1.0
            if this_neg_tp_exon_bp > 0:
                neg_overlap_elements += ["tp_utr_exon"]
                element_counts['neg_tp_exon'] += 1.0
            if this_neg_tp_intron_bp > 0:
                neg_overlap_elements += ["tp_utr_intron"]
                element_counts['neg_tp_intron'] += 1.0
            if this_neg_promoter_bp > 0:
                neg_overlap_elements += ["promoter"]
                element_counts['neg_promoter'] += 1.0
            if this_neg_exon_bp > 0:
                neg_overlap_elements += ["exon"]
                element_counts['neg_exon'] += 1.0
            if this_neg_intron_bp > 0:
                neg_overlap_elements += ["intron"]
                element_counts['neg_intron'] += 1.0
            if this_neg_repeat_bp > 0:
                neg_overlap_elements += ["repeat"]
                element_counts['neg_repeat'] += 1.0
            if not neg_overlap_elements:
                neg_overlap_elements = ["intergenic"]
                element_counts['neg_intergenic'] += 1.0
                                                 
            if len(pos_overlap_elements) > 1:
                logging_function("More than one overlapping positive strand element for SNP %s:%s!" % (entry_data[snp_idx['chr']], entry_data[snp_idx['pos']]))
                logging_function(pos_overlap_elements)
                logging_function("Writing out highest priority element: %s" % (pos_overlap_elements[0]))
            if len(neg_overlap_elements) > 1:
                logging_function("More than one overlapping negative strand element for SNP %s:%s!" % (entry_data[snp_idx['chr']], entry_data[snp_idx['pos']]))
                logging_function(neg_overlap_elements)
                logging_function("Writing out highest priority element: %s" % (neg_overlap_elements[0]))

            entry_out.write("\t".join(entry_data + [pos_overlap_elements[0], neg_overlap_elements[0]])+'\n')                                       
                                       # ",".join(pos_overlap_elements),
                                       # ",".join(neg_overlap_elements)])+"\n")

    with open(partition_outprefix+"_partition_summary.txt", "wb") as summary_out:
        summary_out.write("\t".join(['type', 'pos_number', 'pos_proportion', 'neg_number', 'neg_proportion'])+'\n')
        for x in element_list:
            pos_count = element_counts['pos_'+x]
            neg_count = element_counts['neg_'+x]
            ## TODO: fix roundoff of floats?
            summary_out.write("\t".join([x, str(pos_count), str(pos_count / entry_ctr),
                                         str(neg_count), str(neg_count / entry_ctr)])+"\n")
    
    return partition_outprefix+"_entrywise_partition.txt"

def calculate_unstranded_genomic_partition(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, partition_dir, logging_function):
    """
    determines what kind of genomic element each SNP falls in
    this uses annotation without regard to strand
    """
    try:
        os.makedirs(outdir+"/unstranded_genomic_partition/")
    except OSError:
      pass

    ## get all the file paths
    ## this assumes that the directory contains the fully parsed files
    fp_utr_exons_f = partition_dir+"/parsed_5utr_exons.merged.bed"
    fp_utr_introns_f = partition_dir+"/n5e_5utr_introns.bed"
    tp_utr_exons_f = partition_dir+"/n5e5i_3utr_exons.bed"
    tp_utr_introns_f = partition_dir+"/n5e5i3e_3utr_introns.bed"
    promoters_f = partition_dir+"/n5e5i3e3i_promoters.bed"
    exons_f = partition_dir+"/n5e5i3e3ip_exons.bed"
    introns_f = partition_dir+"/n5e5i3e3ipe_introns.bed"
    repeats_f = partition_dir+"/n5e5i3e3ipei_repeats.bed"

    partition_outprefix = outdir+"/unstranded_genomic_partition/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)+"_unstranded"

    # define a convenience dictionary to make subsetting clearer for bed files
    bed_idx = {} 
    bed_idx["chr"] = 0
    bed_idx["start"] = 1
    bed_idx["end"] = 2    
    bed_idx["name"] = 3
    
    ## write a single output file
    ## contains each entry with its positive and negative strand partition
    with open(ld_snp_file, 'rb') as snpin, open(fp_utr_exons_f, 'rb') as fp_exons, open(fp_utr_introns_f, 'rb') as fp_introns, open(tp_utr_exons_f, 'rb') as tp_exons, open(tp_utr_introns_f, 'rb') as tp_introns, open(promoters_f, 'rb') as promoters, open(exons_f, 'rb') as exons, open(introns_f, 'rb') as introns, open(repeats_f, 'rb') as repeats, open(partition_outprefix+"_entrywise_partition.txt", 'wb') as entry_out:
        # get the LD SNP header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}
        
        ## write the output header
        entry_out.write("\t".join(snp_header_data+['partition'])+"\n")

        ## store lists of each type of element
        cur_fp_exons = []
        cur_fp_introns = []
        cur_tp_exons = []
        cur_tp_introns = []        
        cur_promoters = []
        cur_exons = []
        cur_introns = []
        cur_repeats = []        
        
        ## start reading in these elements
        cur_fp_exon = fp_exons.readline().strip().split('\t')
        cur_fp_intron = fp_introns.readline().strip().split('\t')        
        cur_tp_exon = tp_exons.readline().strip().split('\t')
        cur_tp_intron = tp_introns.readline().strip().split('\t')        
        cur_promoter = promoters.readline().strip().split('\t')
        cur_exon = exons.readline().strip().split('\t')
        cur_intron = introns.readline().strip().split('\t')
        cur_repeat = repeats.readline().strip().split('\t')

        ## track which chromosome we're on
        this_chr = ''

        ## a regex to check that we are in a 'real' chromosome
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        

        ## track how many times each type of element was found
        element_list = ['fp_exon', 'fp_intron', 'tp_exon', 'tp_intron',
                        'promoter', 'exon', 'intron', 'repeat', 'intergenic']
        element_counts = {x:0.0 for x in element_list}

        entry_ctr = 0.0        
        for entry in snpin:
            entry_ctr += 1.0
            entry_data = entry.strip().split('\t')
            entry_pos = int(entry_data[snp_idx['pos']])
            this_rsID = entry_data[snp_idx['rsID']]
            this_tag_rsID = entry_data[snp_idx['tag_rsID']]
            
            ## while we're on a new chromosome, reset all the tracking variables:
            if entry_data[snp_idx['chr']] != this_chr:
                this_chr = entry_data[snp_idx['chr']]
                ## store lists of each type of element
                cur_fp_exons = []
                cur_fp_introns = []
                cur_tp_exons = []
                cur_tp_introns = []
                cur_promoters = []
                cur_exons = []
                cur_introns = []
                cur_repeats = []
                                
                ## read in data until we match chromosomes
                while (len(cur_fp_exon) > 1
                    and cur_fp_exon[bed_idx['chr']] != this_chr):
                    cur_fp_exon = fp_exons.readline().strip().split('\t')
                while (len(cur_fp_intron) > 1
                       and cur_fp_intron[bed_idx['chr']] != this_chr):
                    cur_fp_intron = fp_introns.readline().strip().split('\t')
                while (len(cur_tp_exon) > 1
                       and cur_tp_exon[bed_idx['chr']] != this_chr):
                    cur_tp_exon = tp_exons.readline().strip().split('\t')
                while (len(cur_tp_intron) > 1
                       and cur_tp_intron[bed_idx['chr']] != this_chr):
                    cur_tp_intron = tp_introns.readline().strip().split('\t')
                while (len(cur_promoter) > 1
                       and cur_promoter[bed_idx['chr']] != this_chr):
                    cur_promoter = promoters.readline().strip().split('\t')
                while (len(cur_exon) > 1
                       and cur_exon[bed_idx['chr']] != this_chr):
                    cur_exon = exons.readline().strip().split('\t')
                while (len(cur_intron) > 1
                       and cur_intron[bed_idx['chr']] != this_chr):
                    cur_intron = introns.readline().strip().split('\t')
                while (len(cur_repeat) > 1
                       and cur_repeat[bed_idx['chr']] != this_chr):
                    cur_repeat = repeats.readline().strip().split('\t')
                logging_function('Parsing chromosome '+this_chr)
                
            ## if we have a match to one of the weird chromosomes, just skip to the next entry
            if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                continue                    

            ## add all possibly overlapping entries for each type of element
            ## must have same chromosome and start before the SNP position
            ## note that the bed files are 0-based, half-open so the start can be equal to the
            ## SNP position but the end must be strictly greater than
            while (len(cur_fp_exon) > 1
                   and cur_fp_exon[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_fp_exon[bed_idx['start']]) <= entry_pos):
                if (int(cur_fp_exon[bed_idx['end']]) >= entry_pos):
                    cur_fp_exons.append(cur_fp_exon)
                cur_fp_exon = fp_exons.readline().strip().split('\t')                
            while (len(cur_fp_intron) > 1
                   and cur_fp_intron[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_fp_intron[bed_idx['start']]) <= entry_pos):
                if (int(cur_fp_intron[bed_idx['end']]) >= entry_pos):
                    cur_fp_introns.append(cur_fp_intron)
                cur_fp_intron = fp_introns.readline().strip().split('\t')
            while (len(cur_tp_exon) > 1
                   and cur_tp_exon[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_tp_exon[bed_idx['start']]) <= entry_pos):
                if (int(cur_tp_exon[bed_idx['end']]) >= entry_pos):
                    cur_tp_exons.append(cur_tp_exon)                   
                cur_tp_exon = tp_exons.readline().strip().split('\t')                
            while (len(cur_tp_intron) > 1
                   and cur_tp_intron[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_tp_intron[bed_idx['start']]) <= entry_pos):
                if (int(cur_tp_intron[bed_idx['end']]) >= entry_pos):
                    cur_tp_introns.append(cur_tp_intron)                   
                cur_tp_intron = tp_introns.readline().strip().split('\t')
            while (len(cur_promoter) > 1
                   and cur_promoter[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_promoter[bed_idx['start']]) <= entry_pos):
                if (int(cur_promoter[bed_idx['end']]) >= entry_pos):
                    cur_promoters.append(cur_promoter)                   
                cur_promoter = promoters.readline().strip().split('\t')                
            while (len(cur_exon) > 1
                   and cur_exon[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_exon[bed_idx['start']]) <= entry_pos):
                if (int(cur_exon[bed_idx['end']]) >= entry_pos):
                    cur_exons.append(cur_exon)
                cur_exon = exons.readline().strip().split('\t')
            while (len(cur_intron) > 1
                   and cur_intron[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_intron[bed_idx['start']]) <= entry_pos):
                if (int(cur_intron[bed_idx['end']]) >= entry_pos):
                    cur_introns.append(cur_intron)
                cur_intron = introns.readline().strip().split('\t')            
            while (len(cur_repeat) > 1
                   and cur_repeat[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_repeat[bed_idx['start']]) <= entry_pos):
                if (int(cur_repeat[bed_idx['end']]) >= entry_pos):
                    cur_repeats.append(cur_repeat)
                cur_repeat = repeats.readline().strip().split('\t')
                                
            ## remove entries that are before the current entry here:
            ## shouldn't be necessary..
            for i in range(len(cur_fp_exons)-1, -1, -1):
                if cur_fp_exons[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_fp_exons[i][bed_idx['end']]) <= entry_pos:
                    cur_fp_exons.pop(i)                    
            for i in range(len(cur_fp_introns)-1, -1, -1):
                if cur_fp_introns[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_fp_introns[i][bed_idx['end']]) <= entry_pos:
                    cur_fp_introns.pop(i)
            for i in range(len(cur_tp_exons)-1, -1, -1):
                if cur_tp_exons[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_tp_exons[i][bed_idx['end']]) <= entry_pos:
                    cur_tp_exons.pop(i)                    
            for i in range(len(cur_tp_introns)-1, -1, -1):
                if cur_tp_introns[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_tp_introns[i][bed_idx['end']]) <= entry_pos:
                    cur_tp_introns.pop(i)            
            for i in range(len(cur_promoters)-1, -1, -1):
                if cur_promoters[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_promoters[i][bed_idx['end']]) <= entry_pos:
                    cur_promoters.pop(i)
            for i in range(len(cur_exons)-1, -1, -1):
                if cur_exons[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_exons[i][bed_idx['end']]) <= entry_pos:
                    cur_exons.pop(i)
            for i in range(len(cur_introns)-1, -1, -1):
                if cur_introns[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_introns[i][bed_idx['end']]) <= entry_pos:
                    cur_introns.pop(i)
            for i in range(len(cur_repeats)-1, -1, -1):
                if cur_repeats[i][bed_idx['chr']] != entry_data[snp_idx['chr']] or int(cur_repeats[i][bed_idx['end']]) <= entry_pos:
                    cur_repeats.pop(i)                                        

            ## now compare the SNP (current entry) to each of the genomic loci
            ## track overlaps with all kinds of elements just in case they aren't mutually exclusive
            this_fp_exon_bp = 0.0
            this_fp_intron_bp = 0.0
            this_tp_exon_bp = 0.0
            this_tp_intron_bp = 0.0            
            this_promoter_bp = 0.0
            this_exon_bp = 0.0
            this_intron_bp = 0.0
            this_repeat_bp = 0.0

            ## calculate overlaps
            ## again, since we have positional data, we can just ask if it fits in the interval
            for cfe in cur_fp_exons:
                this_fp_exon_bp += 1 if (int(cfe[bed_idx['start']]) <= entry_pos < int(cfe[bed_idx['end']])) else 0
            for cfi in cur_fp_introns:
                this_fp_intron_bp += 1 if (int(cfi[bed_idx['start']]) <= entry_pos < int(cfi[bed_idx['end']])) else 0
            for cte in cur_tp_exons:
                this_tp_exon_bp += 1 if (int(cte[bed_idx['start']]) <= entry_pos < int(cte[bed_idx['end']])) else 0
            for cti in cur_tp_introns:
                this_tp_intron_bp += 1 if (int(cti[bed_idx['start']]) <= entry_pos < int(cti[bed_idx['end']])) else 0
            for cp in cur_promoters:
                this_promoter_bp += 1 if (int(cp[bed_idx['start']]) <= entry_pos < int(cp[bed_idx['end']])) else 0
            for ce in cur_exons:
                this_exon_bp += 1 if (int(ce[bed_idx['start']]) <= entry_pos < int(ce[bed_idx['end']])) else 0
            for ci in cur_introns:
                this_intron_bp += 1 if (int(ci[bed_idx['start']]) <= entry_pos < int(ci[bed_idx['end']])) else 0
            for cr in cur_repeats:
                this_repeat_bp += 1 if (int(cr[bed_idx['start']]) <= entry_pos < int(cr[bed_idx['end']])) else 0
                
            ## go through and check for overlaps
            ## make list of overlapping elements
            overlap_elements = []
            
            if this_fp_exon_bp > 0:
                overlap_elements += ["fp_utr_exon"]
                element_counts['fp_exon'] += 1.0
            if this_fp_intron_bp > 0:
                overlap_elements += ["fp_utr_intron"]
                element_counts['fp_intron'] += 1.0
            if this_tp_exon_bp > 0:
                overlap_elements += ["tp_utr_exon"]
                element_counts['tp_exon'] += 1.0
            if this_tp_intron_bp > 0:
                overlap_elements += ["tp_utr_intron"]
                element_counts['tp_intron'] += 1.0
            if this_promoter_bp > 0:
                overlap_elements += ["promoter"]
                element_counts['promoter'] += 1.0
            if this_exon_bp > 0:
                overlap_elements += ["exon"]
                element_counts['exon'] += 1.0
            if this_intron_bp > 0:
                overlap_elements += ["intron"]
                element_counts['intron'] += 1.0
            if this_repeat_bp > 0:
                overlap_elements += ["repeat"]
                element_counts['repeat'] += 1.0
            if not overlap_elements:
                overlap_elements = ["intergenic"]
                element_counts['intergenic'] += 1.0
                                                 
            if len(overlap_elements) > 1:
                logging_function("More than one overlapping element for %s:%s!" % (entry_data[snp_idx['chr']], entry_pos))
                logging_function(overlap_elements)
                logging_function("Writing out highest priority element: %s" % (overlap_elements[0]))

            entry_out.write("\t".join(entry_data + [overlap_elements[0]])+"\n")
                                       
                                       # ",".join(pos_overlap_elements),
                                       # ",".join(neg_overlap_elements)])+"\n")

    with open(partition_outprefix+"_partition_summary.txt", "wb") as summary_out:
        summary_out.write("\t".join(['type', 'number', 'proportion'])+'\n')
        for x in element_list:
            count = element_counts[x]
            ## TODO: fix roundoff of floats?
            summary_out.write("\t".join([x, str(count), str(count / entry_ctr)])+"\n")
    
    return partition_outprefix+"_entrywise_partition.txt"

# fantom5_dir = '/home/alexaml/data/FANTOM5/Enhancers/facet_expressed_enhancers/sorted/'
# enhancer_window = 1000
def compute_fantom5_midpoint_overlap(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, fantom5_dir, enhancer_window, skip_enh_summary, logging_function):
    """
    Calculates overlap of LD SNPs with all FANTOM5 data, using a window around the midpoint of
    transcription
    """
    try:
        os.makedirs(outdir+"/fantom5_overlap/")
    except OSError:
      pass

    ## prefix for the files generated
    overlap_outprefix = outdir+"/fantom5_overlap/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)+"_"+str(enhancer_window)+"bp_around_midpoint"

    # define a convenience dictionary to make subsetting clearer
    bed_idx = {} 
    bed_idx["chr"] = 0
    bed_idx["start"] = 1
    bed_idx["end"] = 2    
    bed_idx["name"] = 3

    # make a list to store all the cell / tissue types we have
    cell_tissue_types = []
    
    # make a dict for all the enhancer files
    enh_files = {}
    # also make a dict to store the lists of all current enhancer loci
    enh_loci = {}
    # finally, also have one to store the current enhancer for each tissue/cell type
    cur_enh_locus = {}
    for enh_f in glob.glob(fantom5_dir+"*.bed"):
        ## get just the ID and the name
        this_key = os.path.basename(enh_f).replace("_expressed_enhancers.bed", "")
        cell_tissue_types.append(this_key)
        
        ## open the file and associate it with the key
        enh_files[this_key] = open(enh_f, 'rb')
        ## initialize the list of enhancers
        enh_loci[this_key] = []
        ## initialize the current enhancer storage (this will get added to the list later)
        cur_enh_locus[this_key] = enh_files[this_key].readline().strip().split("\t")

    ## now sort the tissue types
    cell_tissue_types = sorted(cell_tissue_types)
        
    ## write two output files: one is a summary across the different cell types for each SNP, the other contains
    ## the specific enhancer loci that overlap
    with open(ld_snp_file, 'rb') as snpin,  open(overlap_outprefix+"_enh_overlaps.txt", 'wb') as enh_overlaps_out:
        if not skip_enh_summary:
            summary_out = open(overlap_outprefix+"_summary.txt", 'wb')
        # get the LD SNP header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

        if not skip_enh_summary:
            summary_out.write("\t".join(snp_header_data+[t+"_num_overlaps\t"+t+"_dists" for t in cell_tissue_types])+"\n")
        enh_overlaps_out.write("\t".join(snp_header_data + ['enh_source', 'enh_chr', 'enh_start', 'enh_end'])+"\n")

        ## track which chromosome we're on
        this_chr = ''

        ## a regex to check that we are in a 'real' chromosome
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        

        for entry in snpin:
            entry_data = entry.strip().split('\t')
            entry_pos = int(entry_data[snp_idx['pos']])
            this_rsID = entry_data[snp_idx['rsID']]
            this_tag_rsID = entry_data[snp_idx['tag_rsID']]
            
            ## if we're on a new chromosome, reset all the tracking variables:
            if entry_data[snp_idx['chr']] != this_chr:
                this_chr = entry_data[snp_idx['chr']]
                logging_function('Parsing chromosome '+this_chr)
                        
                for k in cell_tissue_types:
                    # reset the enhancer lists
                    enh_loci[k] = []

                if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                    ## for now, just skip this entry..
                    continue                    
                else:
                    # read in data until we match chromosomes
                    for k in cell_tissue_types:
                        while(len(cur_enh_locus[k]) > 1
                            and cur_enh_locus[k][bed_idx['chr']] != this_chr):
                            cur_enh_locus[k] = enh_files[k].readline().strip().split("\t")

            ## now we go through all the cell and tissue types
            ## store the number of overlaps and the distances to each
            num_enh_overlaps = {}
            enh_dists = {}
            for k in cell_tissue_types:                
                ## add all the new overlapping entries
                while(len(cur_enh_locus[k]) > 1
                      and cur_enh_locus[k][bed_idx['chr']]==entry_data[snp_idx['chr']]):
                    cur_enh_start = int(cur_enh_locus[k][bed_idx['start']])
                    cur_enh_end = int(cur_enh_locus[k][bed_idx['end']])
                    this_mid = (cur_enh_start + cur_enh_end)/2
                    
                    if (this_mid - enhancer_window) <= entry_pos:
                        ## only save it if it overlaps
                        if (this_mid + enhancer_window) >= entry_pos:
                            enh_loci[k].append(cur_enh_locus[k])
                        cur_enh_locus[k] = enh_files[k].readline().strip().split("\t")
                    # if we're not overlapping anymore, stop the loop
                    else:
                        break
                        
                ## remove all enhancers from the wrong chromosome or that end before this starts
                for i in range(len(enh_loci[k])-1, -1, -1):
                    if enh_loci[k][i][bed_idx['chr']] != this_chr:
                        enh_loci[k].pop(i)
                    else:                            
                        enh_start = int(enh_loci[k][i][bed_idx['start']])
                        enh_end = int(enh_loci[k][i][bed_idx['end']])
                        this_mid = (enh_start + enh_end)/2
                        ## <= because of the half-open interval
                        if (this_mid + enhancer_window) <= entry_pos:
                            enh_loci[k].pop(i)                            

                ## initialize the storage dicts
                num_enh_overlaps[k] = 0
                enh_dists[k] = []
                            
                ## now compare the SNP to the possible enhancers in this tissue type
                ## note that these are all guaranteed to overlap because of the above section
                ## i leave the redundant check in for now
                for enh in enh_loci[k]:
                    enh_start = int(enh[bed_idx['start']])
                    enh_end = int(enh[bed_idx['end']])
                    ## calculate the midpoint of transcription
                    ## this stays as an integer on purpose, since we can't have fractional position
                    this_mid = (enh_start + enh_end)/2

                    if (entry_pos >= this_mid - enhancer_window
                        and entry_pos <= this_mid + enhancer_window):
                        num_enh_overlaps[k] += 1
                        ## get the distance relative to the midpoint
                        enh_dists[k].append(str(this_mid - entry_pos))
                        ## write out the specific enhancer overlap
                        enh_overlaps_out.write("\t".join(entry_data + [k, enh[bed_idx['chr']],str(enh_start), str(enh_end)])+"\n")
            ## finally, write to the summary file
            if not skip_enh_summary:
                summary_out.write("\t".join(entry_data +
                                            [str(num_enh_overlaps[t])+"\t"+",".join(enh_dists[t])
                                            for t in cell_tissue_types])+"\n")
            
    # close all the enhancer files
    for k in enh_files.keys():
        enh_files[k].close()

    # close the summary, if it exists
    if not skip_enh_summary:
        summary_out.close()

    # return the file with one line per enhancer overlap
    return overlap_outprefix+"_enh_overlaps.txt"

def compute_fantom5_locus_overlap(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, fantom5_dir, enhancer_window, skip_enh_summary, logging_function):
    """
    Calculates overlap of LD SNPs with all FANTOM5 data, using a window around the actual
    originally transcribed locus
    """
    try:
        os.makedirs(outdir+"/fantom5_overlap/")
    except OSError:
      pass

    ## prefix for the files generated
    overlap_outprefix = outdir+"/fantom5_overlap/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)+"_"+str(enhancer_window)+"bp_around_orig_locus"

    # define a convenience dictionary to make subsetting clearer
    bed_idx = {} 
    bed_idx["chr"] = 0
    bed_idx["start"] = 1
    bed_idx["end"] = 2    
    bed_idx["name"] = 3

    # make a list to store all the cell / tissue types we have
    cell_tissue_types = []
    
    # make a dict for all the enhancer files
    enh_files = {}
    # also make a dict to store the lists of all current enhancer loci
    enh_loci = {}
    # finally, also have one to store the current enhancer for each tissue/cell type
    cur_enh_locus = {}
    for enh_f in glob.glob(fantom5_dir+"*.bed"):
        ## get just the ID and the name
        this_key = os.path.basename(enh_f).replace("_expressed_enhancers.bed", "")
        cell_tissue_types.append(this_key)

        ## open the file and associate it with the key        
        enh_files[this_key] = open(enh_f, 'rb')
        ## initialize the list of enhancers
        enh_loci[this_key] = []
        ## initialize the current enhancer storage (this will get added to the list later)
        cur_enh_locus[this_key] = enh_files[this_key].readline().strip().split("\t")

    ## now sort the tissue types
    cell_tissue_types = sorted(cell_tissue_types)
        
    ## write two output files: one is a summary across the different cell types, the other contains
    ## the specific enhancer loci that overlap

    with open(ld_snp_file, 'rb') as snpin, open(overlap_outprefix+"_enh_overlaps.txt", 'wb') as enh_overlaps_out:
        if not skip_enh_summary:
            summary_out = open(overlap_outprefix+"_summary.txt", 'wb')
        
        # get the LD SNP header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

        if not skip_enh_summary:
            summary_out.write("\t".join(snp_header_data+[t+"_num_overlaps\t"+t+"_dists" for t in cell_tissue_types])+"\n")
        enh_overlaps_out.write("\t".join(snp_header_data + ['enh_source', 'enh_chr', 'enh_start', 'enh_end'])+"\n")

        ## track which chromosome we're on
        this_chr = ''

        ## a regex to check that we are in a 'real' chromosome
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        

        for entry in snpin:
            entry_data = entry.strip().split('\t')
            entry_pos = int(entry_data[snp_idx['pos']])
            this_rsID = entry_data[snp_idx['rsID']]
            this_tag_rsID = entry_data[snp_idx['tag_rsID']]
            
            ## if we're on a new chromosome, reset all the tracking variables:
            if entry_data[snp_idx['chr']] != this_chr:
                this_chr = entry_data[snp_idx['chr']]
                logging_function('Parsing chromosome '+this_chr)
                        
                for k in cell_tissue_types:
                    # reset the enhancer lists
                    enh_loci[k] = []

                if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                    ## for now, just skip this entry..
                    continue                    
                else:
                    # read in data until we match chromosomes
                    for k in cell_tissue_types:
                        while(len(cur_enh_locus[k]) > 1
                            and cur_enh_locus[k][bed_idx['chr']] != this_chr):
                            cur_enh_locus[k] = enh_files[k].readline().strip().split("\t")

            ## now we go through all the cell and tissue types
            ## store the number of overlaps and the distances to each
            num_enh_overlaps = {}
            enh_dists = {}
            for k in cell_tissue_types:                
                ## add all the new overlapping entries
                while(len(cur_enh_locus[k]) > 1
                      and cur_enh_locus[k][bed_idx['chr']]==entry_data[snp_idx['chr']]):
                    cur_enh_start = int(cur_enh_locus[k][bed_idx['start']])
                    cur_enh_end = int(cur_enh_locus[k][bed_idx['end']])
                    
                    if (cur_enh_start - enhancer_window) <= entry_pos:
                        ## only save it if it might overlap
                        if (cur_enh_start + enhancer_window) >= entry_pos:
                            enh_loci[k].append(cur_enh_locus[k])
                        cur_enh_locus[k] = enh_files[k].readline().strip().split("\t")
                    # if we're not overlapping anymore, stop the loop
                    else:
                        break
                        
                ## remove all enhancers from the wrong chromosome or that end before this starts
                for i in range(len(enh_loci[k])-1, -1, -1):
                    if enh_loci[k][i][bed_idx['chr']] != this_chr:
                        enh_loci[k].pop(i)
                    else:                            
                        enh_start = int(enh_loci[k][i][bed_idx['start']])
                        enh_end = int(enh_loci[k][i][bed_idx['end']])
                        if (enh_end + enhancer_window) <= entry_pos:
                            enh_loci[k].pop(i)                            

                ## initialize the storage dicts
                num_enh_overlaps[k] = 0
                enh_dists[k] = []
                            
                ## now compare the SNP to the possible enhancers in this tissue type
                for enh in enh_loci[k]:
                    enh_start = int(enh[bed_idx['start']])
                    enh_end = int(enh[bed_idx['end']])
                    ## calculate the midpoint of transcription, to calculate distance
                    ## this stays as an integer on purpose, since we can't have fractional position
                    this_mid = (enh_start + enh_end)/2

                    if (entry_pos >= enh_start - enhancer_window
                        and entry_pos <= enh_end + enhancer_window):
                        num_enh_overlaps[k] += 1
                        ## get the distance relative to the midpoint
                        enh_dists[k].append(str(this_mid - entry_pos))
                        ## write out the specific enhancer overlap
                        enh_overlaps_out.write("\t".join(entry_data + [k, enh[bed_idx['chr']],str(enh_start), str(enh_end)])+"\n")                        
            ## finally, write to the summary file
            if not skip_enh_summary:
                summary_out.write("\t".join(entry_data +
                                            [str(num_enh_overlaps[t])+"\t"+",".join(enh_dists[t])
                                            for t in cell_tissue_types])+"\n")
            
    # close all the enhancer files
    for k in enh_files.keys():
        enh_files[k].close()

    # close the summary, if it exists
    if not skip_enh_summary:
        summary_out.close()
        
    # return the file with one line per enhancer overlap
    return overlap_outprefix+"_enh_overlaps.txt"

def compute_closest_enhancer(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, fantom5_dir, logging_function):
    """
    For each SNP, finds the closest enhancer across all FANTOM5 tissue types
    Distance is defined relative to the midpoint of enhancer transcription
    """
    try:
        os.makedirs(outdir+"/closest_fantom5_enh/")
    except OSError:
      pass

    ## prefix for the file(s) generated
    enh_outprefix = outdir+"/closest_fantom5_enh/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)

    # define a convenience dictionary to make subsetting clearer
    bed_idx = {} 
    bed_idx["chr"] = 0
    bed_idx["start"] = 1
    bed_idx["end"] = 2    
    bed_idx["name"] = 3

    # make a list to store all the cell / tissue types we have
    cell_tissue_types = []
    
    # make a dict for all the enhancer files
    enh_files = {}
    # also make a dict to store the lists of all current enhancer loci (by chromosome)
    enh_loci = {}
    # finally, also have one to store the current enhancer
    cur_enh_locus = {}
    for enh_f in glob.glob(fantom5_dir+"*.bed"):
        ## get just the ID and the name
        this_key = os.path.basename(enh_f).replace("_expressed_enhancers.bed", "")
        cell_tissue_types.append(this_key)
        
        ## open the file and associate it with the key
        enh_files[this_key] = open(enh_f, 'rb')
        ## initialize the list of enhancers
        enh_loci[this_key] = []
        ## initialize the current enhancer storage (this will get added to the list later)
        cur_enh_locus[this_key] = enh_files[this_key].readline().strip().split("\t")

    ## now sort the tissue types
    cell_tissue_types = sorted(cell_tissue_types)
        
    ## write one output file, containing the closest enhancer(s), the tissue(s) of origin,
    ## and the distances
    with open(ld_snp_file, 'rb') as snpin, open(enh_outprefix+"_closest_enhancers.txt", 'wb') as enh_out:
        # get the LD SNP header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}
        
        enh_out.write("\t".join(snp_header_data + ['enh_chr', 'enh_start', 'enh_end', 'enh_tissue', 'dist_to_enh_midpoint'])+'\n')

        ## track which chromosome we're on
        this_chr = ''

        ## a regex to check that we are in a 'real' chromosome
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        

        for entry in snpin:
            entry_data = entry.strip().split('\t')
            entry_pos = int(entry_data[snp_idx['pos']])
            this_rsID = entry_data[snp_idx['rsID']]
            this_tag_rsID = entry_data[snp_idx['tag_rsID']]
            
            ## if we're on a new chromosome, reset all the tracking variables
            ## since this is the closest enhancer analysis, we store all the ones on the
            ## current chromosome
            if entry_data[snp_idx['chr']] != this_chr:
                this_chr = entry_data[snp_idx['chr']]
                logging_function('Parsing chromosome '+this_chr)
                        
                for k in cell_tissue_types:
                    # reset the enhancer lists
                    enh_loci[k] = []

                if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                    ## for now, just skip this entry..
                    continue                    
                else:
                    # read through the enh files until we match chromosomes
                    for k in cell_tissue_types:
                        while(len(cur_enh_locus[k]) > 1
                            and cur_enh_locus[k][bed_idx['chr']] != this_chr):
                            cur_enh_locus[k] = enh_files[k].readline().strip().split("\t")
                        ## now add all enhancers on the same chr
                        while(len(cur_enh_locus[k]) > 1
                            and cur_enh_locus[k][bed_idx['chr']]==entry_data[snp_idx['chr']]):
                            enh_loci[k].append(cur_enh_locus[k])
                            cur_enh_locus[k] = enh_files[k].readline().strip().split("\t")

            ## now we go through all the cell and tissue types and look for the closest enh
            ## store the smallest current distance, and the enhancer(s) that gave that dist
            cur_min_enh_dist = "UNSET"
            closest_enhs = {'chr':[], 'start':[], 'end':[], 'tissue':[]}
            for k in cell_tissue_types:               
                ## compare the SNP to the enhancers in this tissue type
                for enh in enh_loci[k]:
                    enh_start = int(enh[bed_idx['start']])
                    enh_end = int(enh[bed_idx['end']])
                    ## calculate the midpoint of transcription
                    ## this stays as an integer on purpose, since we can't have fractional position
                    this_mid = (enh_start + enh_end)/2
                    ## calculate the relative distance to this enhancer midpoint
                    this_dist = this_mid - entry_pos
                    ## if this is the new closest one (by absolute value)
                    if(cur_min_enh_dist=="UNSET" or abs(this_dist) < abs(cur_min_enh_dist)):
                        cur_min_enh_dist = this_dist
                        ## store this enhancer in the list (overwriting any previous entries)
                        closest_enhs['chr'] = this_chr
                        closest_enhs['start'] = [str(enh_start)]
                        closest_enhs['end'] = [str(enh_end)]
                        closest_enhs['tissue'] = [k]
                    ## if we have a tie, add it to the storage list
                    elif(this_dist == cur_min_enh_dist):
                        ## only add it if we haven't seen it before
                        if(str(enh_start) not in closest_enhs['start'] and
                           str(enh_end) not in closest_enhs['end']):
                            closest_enhs['start'].append(str(enh_start))
                            closest_enhs['end'].append(str(enh_end))
                        ## always add the tissue
                        closest_enhs['tissue'].append(k)
                        ## might want to deal with possibility of equal distance on either side
                    ## if we're now farther away than the closest distance to the right, break
                    ## not using absolute value on this distance to encode the higher position
                    elif(this_dist > abs(cur_min_enh_dist)):
                        break
                        
            ## finally, write to the summary file
            enh_out.write("\t".join(entry_data + [closest_enhs['chr'],
                                                  ",".join(closest_enhs['start']),
                                                 ",".join(closest_enhs['end']),
                                                 ",".join(closest_enhs['tissue']),
                                                 str(cur_min_enh_dist)])+"\n")
            
    # close all the enhancer files
    for k in enh_files.keys():
        enh_files[k].close()

    return enh_outprefix+"_closest_enhancers.txt"

# correlation_file = '/home/alexaml/data/FANTOM5/Enhancers/enhancer_tss_associations.bed'
def correlation_enh_overlap_targets(outdir, outprefix, ld_threshold, ld_check_area, correlation_file, overlap_file, logging_function, overlap_type):
    """
    takes the enhancer overlap output file and annotates it with the
    correlation-based target genes of each enhancer
    """
    try:
        os.makedirs(outdir+"/correlation_enh_targets/")
    except OSError:
      pass

    ## prefix for the files generated
    target_outprefix = outdir+"/correlation_enh_targets/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)
  
    ## read in the correlation file, store in a dict for reference
    tss_associations = {}
    with open(correlation_file, 'rb') as enh_tss_associations:
        ## skip the header (it's just a bed format with a special name column)
        enh_tss_associations.next()
        ## now store all the enh-tss associations
        for line in enh_tss_associations:
            enh_tss_info = line.strip().split("\t")[3].split(";")
            ## use the enhancer info as the key, and store all the rest
            ## only store it if it has all the data we want
            ## (enhancer, refseq ID, symbol, R value, FDR value)
            if len(enh_tss_info) == 5:
                if enh_tss_info[0] in tss_associations:
                    tss_associations[enh_tss_info[0]].append(enh_tss_info[1:])
                else:
                    tss_associations[enh_tss_info[0]] = [enh_tss_info[1:]]

    ## analyze the overlap data:
    if overlap_type=="midpoint":
        overlaps = open(overlap_file, 'rb')
        target_out = open(target_outprefix+"_fantom5_midpoint_overlap_target_genes.txt", 'wb')
    elif overlap_type=="locus":
        overlaps = open(overlap_file, 'rb')
        target_out = open(target_outprefix+"_fantom5_locus_overlap_target_genes.txt", 'wb')
    else:
        logging_function("Overlap type not supported!")
        return

    ## store the header as a dict
    overlap_header = next(overlaps).strip().split("\t")
    overlap_dict = {overlap_header[x]:x for x in range(len(overlap_header))}
        
    ## write the output header
    target_out.write("\t".join([overlap_header[i] for i in range(len(overlap_header))] + 
                                ['num_assocs', 'refseqId', 'symbol', 'r', 'fdr'])+'\n')

    ## go through the enhancer overlaps and add target gene information
    ## count how many didn't match
    no_match_count = 0
    ## also track so we don't print enhancers with no matches
    no_match_enhs = []
    for line in overlaps:
        overlap_data = line.strip().split("\t")
        ## get the key to index the target gene dict
        this_enh_key = overlap_data[overlap_dict['enh_chr']]+":"+overlap_data[overlap_dict['enh_start']]+"-"+overlap_data[overlap_dict['enh_end']]

        if this_enh_key in tss_associations:
            this_assoc = tss_associations[this_enh_key]
            this_num_assoc = len(this_assoc)
            for assoc in this_assoc:
                target_out.write("\t".join(overlap_data+[str(this_num_assoc)]+assoc)+'\n')
        else:
            no_match_count += 1
            if this_enh_key not in no_match_enhs:
                logging_function("No match found for %s" % (this_enh_key))
                no_match_enhs.append(this_enh_key)
            target_out.write("\t".join(overlap_data+["0", "NA", "NA", "NA", "NA"])+'\n')

    ## close files
    overlaps.close()
    target_out.close()
            
def correlation_closest_enh_targets(outdir, outprefix, ld_threshold, ld_check_area, correlation_file, closest_enh_file, logging_function):
    """
    takes the closest enhancer output file and annotates it with the
    correlation-based target genes of each enhancer
    """
    try:
        os.makedirs(outdir+"/correlation_enh_targets/")
    except OSError:
      pass

    ## prefix for the files generated
    target_outprefix = outdir+"/correlation_enh_targets/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)
  
    ## read in the correlation file, store in a dict for reference
    tss_associations = {}
    with open(correlation_file, 'rb') as enh_tss_associations:
        ## skip the header
        enh_tss_associations.next()
        ## now store all the enh-tss associations
        for line in enh_tss_associations:
            enh_tss_info = line.strip().split("\t")[3].split(";")
            ## use the enhancer info as the key, and store all the rest
            ## only store it if it has all the data we want
            ## (enhancer, refseq ID, symbol, R value, FDR value)
            if len(enh_tss_info) == 5:
                if enh_tss_info[0] in tss_associations:
                    tss_associations[enh_tss_info[0]].append(enh_tss_info[1:])
                else:
                    tss_associations[enh_tss_info[0]] = [enh_tss_info[1:]]

    ## analyze the closest enhancer data
    with open(closest_enh_file, 'rb') as closest_enhs, open(target_outprefix+"_closest_enhancer_target_genes.txt", 'wb') as target_out:
        ## store the header as a dict
        closest_enh_header = next(closest_enhs).strip().split("\t")
        closest_enh_dict = {closest_enh_header[x]:x for x in range(len(closest_enh_header))}
        
        ## write the output header
        target_out.write("\t".join([closest_enh_header[i] for i in range(len(closest_enh_header))] +
                                    ['num_assocs', 'refseqId', 'symbol', 'r', 'fdr'])+'\n')

        ## go through the closest enhancers and add target gene information
        ## count how many didn't match
        no_match_count = 0
        ## also track it so we don't print too many
        no_match_enhs = []
        for line in closest_enhs:
            closest_enh_data = line.strip().split("\t")
            ## get the key to index the target gene dict
            this_enh_key = closest_enh_data[closest_enh_dict['enh_chr']]+":"+closest_enh_data[closest_enh_dict['enh_start']]+"-"+closest_enh_data[closest_enh_dict['enh_end']]
            
            if this_enh_key in tss_associations:
                this_assoc = tss_associations[this_enh_key]
                this_num_assoc = len(this_assoc)
                for assoc in this_assoc:
                    target_out.write("\t".join(closest_enh_data+[str(this_num_assoc)]+assoc)+'\n')
            else:
                no_match_count += 1
                if this_enh_key not in no_match_enhs:
                    logging_function("No match found for %s" % (this_enh_key))
                    no_match_enhs.append(this_enh_key)
                target_out.write("\t".join(closest_enh_data+["0", "NA", "NA", "NA", "NA"])+'\n')
                
# gtex_dir = "/home/alexaml/data/GTEx/single_cell_sig_eqtls_v6/"
def gtex_eqtl_overlap(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, gtex_dir, logging_function):
    """
    for each SNP, checks for overlap with GTEx eQTL signals across tissues
    """
    try:
        os.makedirs(outdir+"/gtex_eqtl_overlap/")
    except OSError:
      pass

    eqtl_outprefix = outdir+"/gtex_eqtl_overlap/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)
    
    with open(ld_snp_file, 'rb') as snpin, open(eqtl_outprefix+"_eqtl_overlaps_unsorted.txt", 'wb') as eqtl_overlap_out:
        ## write the header to the eQTL output file
        # start by getting the LD SNP header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}
        
        ## first define the columns we want to get from the eQTL file:
        eqtl_cols = ['gene', 'beta', 't_stat', 'se', 'p_value', 'nom_thresh', 'min_p',
                     'gene_emp_p', 'gene_q_value', 'beta_noNorm', 'minor_allele_samples',
                     'minor_allele_count', 'maf', 'has_best_p', 'is_chosen_snp',
                     'gene_name', 'gene_source', 'gene_chr', 'gene_start', 'gene_stop',
                     'orientation', 'gene_type', 'tss_distance']
        eqtl_overlap_out.write("\t".join(snp_header_data
                                         + ["tissue"] + eqtl_cols)+'\n')

        ## go through the eQTL files, open them, storing the file pointer in this dict
        eqtl_files = {}
        ## store the eQTL we're currently on
        cur_eqtl = {}
        ## note that we assume these files are all sorted correctly
        for eqtl_file in glob.glob(gtex_dir+"/*.snpgenes"):  
            this_tiss = eqtl_file.split("/")[-1].replace("_Analysis.snpgenes", "")
            eqtl_files[this_tiss] = open(eqtl_file, 'rb')
            ## read the header, and move on to the data
            ## we overwrite the header dict a bunch of times
            eqtl_header_data = eqtl_files[this_tiss].readline().strip().split("\t")            
            eqtl_idx = {eqtl_header_data[x]:x for x in range(len(eqtl_header_data))}
            ## also read in the first eQTL
            cur_eqtl[this_tiss] = eqtl_files[this_tiss].readline().strip().split("\t")

        ## read through the SNP list, store SNP information for each chromosome
        ## then read through the eQTL files and compare them
        this_chr = ""
        snp_ctr = 0
        for line in snpin:
            snp_data = line.strip().split('\t')
            ## if we're on a new chromosome, we will make a new dict to store these LD SNPs
            if snp_data[snp_idx['chr']] != this_chr:
                ## if we had an old chromosome, we can now go through and check eQTLs against it
                ## same as writing "if this_chr == '':"
                if this_chr:
                    logging_function("Analyzing eQTLs from chromosome %s" % (this_chr))
                    for this_tiss in eqtl_files.keys():
                        ## read through the eqtl file until we're on the right chromosome
                        while (len(cur_eqtl[this_tiss]) > 1
                               and "chr"+cur_eqtl[this_tiss][eqtl_idx["snp_chrom"]]!=this_chr):
                            cur_eqtl[this_tiss] = eqtl_files[this_tiss].readline().strip().split("\t")
                        ## now check the eqtl overlaps by position and alleles
                        ## check if this file is done being read
                        if len(cur_eqtl[this_tiss])==1:                   
                            continue
                        cur_eqtl_data = cur_eqtl[this_tiss]
                        eqtl_chr = "chr"+cur_eqtl_data[eqtl_idx["snp_chrom"]]
                        ## while we're on the right chromosome
                        while eqtl_chr == this_chr:
                            eqtl_key = ":".join([eqtl_chr,
                                                 cur_eqtl_data[eqtl_idx["snp_pos"]],
                                                 cur_eqtl_data[eqtl_idx["ref"]],
                                                 cur_eqtl_data[eqtl_idx["alt"]]])
                            if eqtl_key in this_chr_ld_snp_dict:
                                ## get all the LD SNPs
                                ld_snp_data = this_chr_ld_snp_dict[eqtl_key]
                                for ld_snp in ld_snp_data:
                                    ## we have an overlap, so write this SNP to the file
                                    outstring = "\t".join(ld_snp + [this_tiss] + 
                                                 [cur_eqtl_data[eqtl_idx[c]] for c in eqtl_cols])
                                    eqtl_overlap_out.write(outstring+"\n")
                            ## read the next eQTL from this tissue
                            cur_eqtl[this_tiss] = eqtl_files[this_tiss].readline().strip().split("\t")
                            if len(cur_eqtl[this_tiss]) > 1:
                                cur_eqtl_data = cur_eqtl[this_tiss]
                                eqtl_chr = "chr"+cur_eqtl[this_tiss][eqtl_idx["snp_chrom"]]
                            else:
                                break
                        
                ## now reset the dict and chromosome tracking variables
                this_chr_ld_snp_dict = {}
                this_chr = snp_data[snp_idx['chr']] 

            ## now read in this SNP and store it (by chr, position, and allele)
            this_key = ":".join([snp_data[snp_idx['chr']], snp_data[snp_idx['pos']],
                                 snp_data[snp_idx['ref']], snp_data[snp_idx['alt']]])
            if this_key not in this_chr_ld_snp_dict:
                this_chr_ld_snp_dict[this_key] = [snp_data]
            else:
                ## if it was already in there, just add it to the list of SNPs at the position
                this_chr_ld_snp_dict[this_key].append(snp_data)
                
        ## analyze eQTLs for the last chromosome
        logging_function("Analyzing eQTLs from chromosome %s" % (this_chr))
        for this_tiss in eqtl_files.keys():
            if len(cur_eqtl[this_tiss]) == 1:
                continue
            cur_eqtl_data = cur_eqtl[this_tiss]
            eqtl_chr = "chr"+cur_eqtl_data[eqtl_idx["snp_chrom"]]
            ## while we're on the right chromosome (the last one)
            while eqtl_chr == this_chr:
                eqtl_key = ":".join([eqtl_chr,
                     cur_eqtl_data[eqtl_idx["snp_pos"]],
                     cur_eqtl_data[eqtl_idx["ref"]],
                     cur_eqtl_data[eqtl_idx["alt"]]])
                                            
                if eqtl_key in this_chr_ld_snp_dict:
                    ## get all the LD SNPs
                    ld_snp_data = this_chr_ld_snp_dict[eqtl_key]
                    for ld_snp in ld_snp_data:
                        ## we have an overlap, so write this SNP to the file
                        outstring = "\t".join(ld_snp + [this_tiss] + 
                                     [cur_eqtl_data[eqtl_idx[c]] for c in eqtl_cols])
                        eqtl_overlap_out.write(outstring+"\n")
                ## read the next one
                cur_eqtl[this_tiss] = eqtl_files[this_tiss].readline().strip().split("\t")
                if len(cur_eqtl[this_tiss]) > 1:
                    cur_eqtl_data = cur_eqtl[this_tiss]
                    eqtl_chr = "chr"+cur_eqtl_data[eqtl_idx["snp_chrom"]]
                else:
                    break
            ## when we're done, close the file
            eqtl_files[this_tiss].close()
        
    ## now sort this file by chromosome, position, rsID, region name, tag rsID
    ## write the indexing strings here to make it neater
    chr_sort = "-k"+",".join([str(snp_idx['chr']+1)]*2)+"V"
    pos_sort = "-k"+",".join([str(snp_idx['pos']+1)]*2)+"n"
    rsid_sort = "-k"+",".join([str(snp_idx['rsID']+1)]*2)
    name_sort = "-k"+",".join([str(snp_idx['tag_name']+1)]*2)
    tag_rsid_sort = "-k"+",".join([str(snp_idx['tag_rsID']+1)]*2)

    logging_function("Sorting eQTL output")
    with open(eqtl_outprefix+"_eqtl_overlaps.txt", "wb") as eqtl_overlap_out:
        subprocess.call(['sort', chr_sort, pos_sort, rsid_sort, name_sort, tag_rsid_sort, eqtl_outprefix+"_eqtl_overlaps_unsorted.txt"], stdout=eqtl_overlap_out)
    subprocess.call(['rm', eqtl_outprefix+"_eqtl_overlaps_unsorted.txt"])

    return eqtl_outprefix+"_eqtl_overlaps.txt"

# factorbook_file = "/project/wang/alexaml/data/factorbook/wgEncodeRegTfbsClusteredWithCellsV3.sorted.bed"
def factorbook_overlap(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, tfbs_file, logging_function):
    """
    for each SNP, computes overlap with factorbook-defined TFBSs
    """
    try:
        os.makedirs(outdir+"/factorbook_overlap/")
    except OSError:
      pass
            
    tfbs_outf = outdir+"/factorbook_overlap/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)+"_tfbs_overlaps.txt"

    # define a convenience dictionary to make subsetting clearer
    bed_idx = {} 
    bed_idx["chr"] = 0
    bed_idx["start"] = 1
    bed_idx["end"] = 2    
    bed_idx["name"] = 3
    bed_idx["score"] = 4
    ## this is a special one for the factorbook data
    bed_idx["cell_types"] = 5
    
    ## we will write just the one file, containing all specific overlaps
    with open(ld_snp_file, 'rb') as snpin, open(tfbs_file, 'rb') as tf_sites, open(tfbs_outf, 'wb') as tf_overlap_out:
        # get the LD SNP header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

        # write a header
        tf_overlap_out.write("\t".join(snp_header_data +
                                       ['tfbs_chr', 'tfbs_start', 'tfbs_end', 'tf', 'score', 'cells'])+"\n")

        ## track which chromosome we're on
        this_chr = ''

        ## a regex to check that we are in a 'real' chromosome
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        

        ## initialize storage for the TFBSs
        cur_tf_sites = []
        cur_tfbs = tf_sites.readline().strip().split('\t')
        
        for entry in snpin:
            entry_data = entry.strip().split('\t')
            entry_pos = int(entry_data[snp_idx['pos']])
            this_rsID = entry_data[snp_idx['rsID']]
            this_tag_rsID = entry_data[snp_idx['tag_rsID']]
            
            ## if we're on a new chromosome, reset all the tracking variables:
            if entry_data[snp_idx['chr']] != this_chr:
                this_chr = entry_data[snp_idx['chr']]
                logging_function('Parsing chromosome '+this_chr)
                cur_tf_sites = []
                ## read in until we match
                while (len(cur_tfbs)>1 and cur_tfbs[bed_idx['chr']] != this_chr):
                    cur_tfbs = tf_sites.readline().strip().split('\t')
                    
            ## if we have a match to one of the weird chromosomes, just skip to the next entry
            if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                continue                    
                    
            ## now add all possibly overlapping entries
            while (len(cur_tfbs)>1 and cur_tfbs[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_tfbs[bed_idx['start']]) <= entry_pos):
                if (int(cur_tfbs[bed_idx['end']]) >= entry_pos):
                    cur_tf_sites.append(cur_tfbs)
                cur_tfbs = tf_sites.readline().strip().split('\t')

            ## remove all too-early entries here
            for i in range(len(cur_tf_sites)-1, -1, -1):
                if (cur_tf_sites[i][bed_idx['chr']] != entry_data[snp_idx['chr']]
                    or int(cur_tf_sites[i][bed_idx['end']]) <= entry_pos):
                    cur_tf_sites.pop(i)

            ## calculate overlaps, write them out directly
            tf_overlaps = 0.0
            for tfbs in cur_tf_sites:
                tf_overlaps += 1.0
                if (int(tfbs[bed_idx['start']]) <= entry_pos < int(tfbs[bed_idx['end']])):
                    ## just write out both entries
                    tf_overlap_out.write("\t".join(entry_data + tfbs)+"\n")

# roadmap_chromhmm_dir = "/home/alexaml/data/roadmap/chromHMM/sorted/"                    
def roadmap_chromhmm_overlap(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, roadmap_chromhmm_dir, logging_function):
    """
    for each SNP, finds the ChromHMM-defined chromatin state in all roadmap data sources
    """
    try:
        os.makedirs(outdir+"/roadmap_chromhmm_states/")
    except OSError:
      pass

    hmm_outprefix = outdir+"/roadmap_chromhmm_states/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)

    with open(ld_snp_file, 'rb') as snpin, open(hmm_outprefix+"_roadmap_chromHMM_states.txt", 'wb') as hmm_state_out:
        ## write the header to the chromHMM output file
        # start by getting the LD SNP header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

        ## now we need to go through and find all the Roadmap data sources that we have
        ## this dict contains the file pointers
        hmm_files = {}
        ## this dict points to the interval we are currently on
        cur_hmm = {}
        for hmm_file in glob.glob(roadmap_chromhmm_dir+"/*mnemonics.bed"):
            ## get the roadmap ID
            eid = hmm_file.split("/")[-1].split("_")[0]
            ## store the file pointer
            hmm_files[eid] = open(hmm_file, 'rb')
            ## read the first interval
            cur_hmm[eid] = hmm_files[eid].readline().strip().split("\t")

        ## make a list of the cell types we have
        roadmap_ids = sorted(hmm_files.keys())
            
        ## now write the header to the output file
        hmm_state_out.write("\t".join(snp_header_data + [k+"_state" for k in roadmap_ids])+'\n')

        ## read through the SNP list while reading through the HMM files
        this_chr = ""
        for line in snpin:
            snp_data = line.strip().split("\t")
            snp_pos = int(snp_data[snp_idx['pos']])
            snp_chr = snp_data[snp_idx['chr']]
            ## first write the SNP info to the file
            hmm_state_out.write("\t".join(snp_data))
            
            if snp_chr != this_chr:
                this_chr = snp_chr
                logging_function("Parsing chromosome "+this_chr)

                for eid in roadmap_ids:
                    ## keep reading through states until we're on the right chromosome
                    while(len(cur_hmm[eid]) > 1 and cur_hmm[eid][0] != this_chr):
                        cur_hmm[eid] = hmm_files[eid].readline().strip().split("\t")

            ## now go through the HMM states and find the one encompassing this position
            for eid in roadmap_ids:
                this_hmm = cur_hmm[eid]
                ## track whether we actually found the state for this SNP in this sample
                eid_found = False
                while(len(this_hmm) > 1 and this_hmm[0] == snp_chr and not eid_found):
                    ## check if the HMM window contains this SNP
                    if (int(this_hmm[1]) <= snp_pos < int(this_hmm[2])):
                        ## write the state out
                        hmm_state_out.write("\t"+this_hmm[3])
                        eid_found = True
                    else:
                        ## read the next file
                        cur_hmm[eid] = hmm_files[eid].readline().strip().split("\t")
                        this_hmm = cur_hmm[eid]
                    
                ## if we made it through without finding the SNP, print an error
                if not eid_found:
                    logging_function("Did not find state for snp %s in sample %s!" % (snp_data[snp_idx['rsID']], eid))
            ## write a newline to the output file
            hmm_state_out.write("\n")            

        ## finally, go through and close the files
        for eid in roadmap_ids:
            hmm_files[eid].close()

    return hmm_outprefix+"_roadmap_chromHMM_states.txt"

def homer_motif_overlap(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, motif_bed_file, motif_pwm_file, motif_seq_file, logging_function):
    """
    for each SNP, checks for overlap with HOMER-defined motifs, and if there is PWM information,
    calculate the effect on the PWM
    """
    try:
        os.makedirs(outdir+"/homer_motif_overlap/")
    except OSError:
      pass
  
    motif_outprefix = outdir+"/homer_motif_overlap/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)

    ## here, we only check that the PWM file exists, because to get to this point we assume
    ## that the sequence file also exists
    if motif_pwm_file:
        outf = motif_outprefix+"_motif_overlap_and_disruption.txt"
        ## read in the PWM information
        pwm_dict = {}
        with open(motif_pwm_file, 'rb') as motif_pwms:
            ## store the TF and PWM that we build up
            this_tf = ""
            this_pwm = []
            for line in motif_pwms:
                line_data = line.strip().split("\t")
                if line[0]==">":
                    ## if we have been building up a PWM, store it
                    if this_tf:
                        if this_tf not in pwm_dict:
                            ## store it in a list so that we can add others
                            pwm_dict[this_tf] = [this_pwm]
                        else:
                            logging_function("TF %s found multiple times in PWM data!" % (this_tf))
                            logging_function("Previous length: %d, current length: %d" % (len(pwm_dict[this_tf][0]), len(this_pwm)))
                            pwm_dict[this_tf].append(this_pwm)
                    ## reset the tracking variables
                    this_tf = line_data[1].split("/")[0]
                    this_pwm = []
                else:
                    ## add this PWM info to the list
                    this_pwm.append([float(x) for x in line_data])
            ## store the last PWM
            pwm_dict[this_tf] = this_pwm
        ## also open the sequence file for us to read concurrently with the bed file
        motif_seqs = open(motif_seq_file, 'rb')
        ## create a convenience dict for mapping bases to PWM indices
        pwm_map = {"A":0, "C":1, "G":2, "T":3}
        ## finally, create a dict for taking the complement of bases
        comp_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    else:
        outf = motif_outprefix+"_motif_overlap.txt"
        
    with open(ld_snp_file, 'rb') as snpin, open(motif_bed_file, 'rb') as motif_beds, open(outf, 'wb') as motif_out:
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

        ## manually define an index for the bed file
        motif_bed_idx = {"motif_chr":0, "motif_start":1, "motif_end":2, "tf_name":3, "log_odds_score":4, "strand":5}
        
        ## write the header to the output file, including motif information
        if motif_pwm_file:
            ## if we have PWM info, we output the reference, alternate, and change in PWM scores
            ## in addition, we report the scores for the specific nucleotides
            ## also output sequence information?
            motif_out.write("\t".join(snp_header_data+["motif_chr", "motif_start", "motif_end", "tf_name", "log_odds_score", "strand", "relative_snp_pos", "ref_pwm", "alt_pwm", "delta_pwm", "ref_prob", "alt_prob", "motif_seq"])+'\n')
            ## also store the current motif sequences
            all_cur_motif_seqs = []
            cur_seq = motif_seqs.readline().strip().split("\t")
        else:
            motif_out.write("\t".join(snp_header_data+["motif_chr", "motif_start", "motif_end", "tf_name", "log_odds_score", "strand", "relative_snp_pos"])+'\n')
        
        ## store all the current motifs (because they may overlap)
        all_cur_motifs = []
        cur_motif = motif_beds.readline().strip().split("\t")

        ## now read through the SNP list
        this_chr = ""
        for line in snpin:
            snp_data = line.strip().split("\t")
            snp_pos = int(snp_data[snp_idx['pos']])
            snp_chr = snp_data[snp_idx['chr']]
        
            ## if we're on a new chromosome, reset the tracking variables
            if snp_chr != this_chr:
                this_chr = snp_chr
                logging_function('Parsing chromosome '+this_chr)
                if motif_pwm_file:
                    all_cur_motif_seqs = []
                all_cur_motifs = []
                while(len(cur_motif) > 1 and cur_motif[motif_bed_idx['motif_chr']] != this_chr):
                    cur_motif = motif_beds.readline().strip().split("\t")
                    if motif_pwm_file:
                        cur_seq = motif_seqs.readline().strip().split("\t")
                        
            ## now go through all the possibly overlapping entries
            while(len(cur_motif) > 1 and cur_motif[motif_bed_idx['motif_chr']]==snp_chr and int(cur_motif[motif_bed_idx['motif_start']])<=snp_pos):
                ## only save the overlapping ones to save memory
                if(int(cur_motif[motif_bed_idx['motif_end']]) >= snp_pos):
                    all_cur_motifs.append(cur_motif)
                    if motif_pwm_file:
                        all_cur_motif_seqs.append(cur_seq)
                cur_motif = motif_beds.readline().strip().split("\t")
                if motif_pwm_file:
                    cur_seq = motif_seqs.readline().strip().split("\t")

            ## remove the early entries / ones from the wrong chromosome
            for i in range(len(all_cur_motifs)-1, -1, -1):
                if(all_cur_motifs[i][motif_bed_idx['motif_chr']] != this_chr or
                   int(all_cur_motifs[i][motif_bed_idx['motif_end']]) < snp_pos):
                   all_cur_motifs.pop(i)
                   if motif_pwm_file:
                       all_cur_motif_seqs.pop(i)

            ## finally, go through the current motifs and calculate overlaps
            for i in range(len(all_cur_motifs)):
                motif = all_cur_motifs[i]
                motif_start = int(motif[motif_bed_idx['motif_start']])
                motif_end = int(motif[motif_bed_idx['motif_end']])
                ## vcf files are 1-based and bed are 0-based half-open half-closed, but the
                ## HOMER motif "bed" files were found to be fully closed, so there's no
                ## coordinate correction
                if motif_start <= snp_pos <= motif_end:
                    motif_strand = motif[motif_bed_idx['strand']]                    
                    ## find the (1-based) position of the SNP relative to the motif
                    if motif_strand=="+":
                        relative_pos = snp_pos - motif_start + 1
                    else:
                        ## for negative strand hits, we need position relative to the end, but
                        ## we use the same correction
                        relative_pos = motif_end - snp_pos + 1

                    ## if we have PWM info, calculate the change
                    if motif_pwm_file:
                        this_tf = motif[motif_bed_idx['tf_name']]
                        ## make it uppercase to be consistent
                        this_seq = all_cur_motif_seqs[i][1].upper()

                        ## if we don't actually have alleles for this SNP, write dummy values
                        if snp_data[snp_idx['alt']].upper()=="NA" or snp_data[snp_idx['ref']].upper()=="NA":
                            logging_function("Missing alleles for SNP %s" % 
                                             (snp_data[snp_idx['rsID']]))
                            ## output dummy values
                            motif_out.write("\t".join(snp_data+motif+[str(relative_pos), "NA", "NA", "NA", "NA", "NA", this_seq])+"\n")
                            continue
                        try:
                            this_tf_pwms = pwm_dict[this_tf]
                        except KeyError:
                            logging_function("Error: TF %s not found in PWM dict!" % (this_tf))
                            ## output dummy values
                            motif_out.write("\t".join(snp_data+motif+[str(relative_pos), "NA", "NA", "NA", "NA", "NA", this_seq])+"\n")
                            continue

                        ## if this TF had more than one motif, we need to figure out which to use
                        if len(this_tf_pwms) > 1:
                            pwm_found = False
                            for pwm in this_tf_pwms:
                                if len(pwm)==len(this_seq):
                                    this_pwm = pwm
                                    pwm_found = True
                                    break
                            if not pwm_found:
                                logging_function("Error: no PWM matching length %d found for TF %s" % (len(this_seq), this_tf))
                        else:
                            this_pwm = this_tf_pwms[0]
                        
                        ## now calculate the PWM log-odds scores for reference and alternate
                        ref_score = 0
                        alt_score = 0
                        ## we also want to store the reference and allele probabilities
                        ref_prob = -1
                        alt_prob = -1
                        for i, c in enumerate(this_seq):
                            ## for now, just use an equiprobable distribution
                            try:
                                ref_score += math.log(this_pwm[i][pwm_map[c]]/0.25)
                            except IndexError:
                                logging_function("List index out of range, i=%s, c=%s for SNP %s and TF %s" % (str(i), c, snp_data[snp_idx['rsID']], this_tf))
                                continue
                            ## for debugging purposes, check whether our major allele matches
                            if i+1==relative_pos:
                                ## store the reference probability
                                ref_prob = this_pwm[i][pwm_map[c]]
                                ## grab the alternate allele
                                snp_alt = snp_data[snp_idx['alt']].upper()
                                if motif_strand=="+":
                                    snp_allele = snp_data[snp_idx['ref']].upper()
                                    if snp_allele != c:
                                        logging_function("Reference SNP allele %s for SNP %s doesn't match up with sequence base %s for TF %s on positive strand" % (snp_allele, snp_data[snp_idx['rsID']], c, this_tf))
                                else:
                                    snp_allele = comp_dict[snp_data[snp_idx['ref']].upper()]
                                    snp_alt = comp_dict[snp_alt]
                                    if snp_allele != c:
                                        logging_function("Reference SNP allele %s for SNP %s doesn't match up with sequence base %s for TF %s on negative strand" % (snp_allele, snp_data[snp_idx['rsID']], c, this_tf))
                                ## store the alternate probability
                                alt_prob = this_pwm[i][pwm_map[snp_alt]]
                                ## now update the alternate sequence score
                                alt_score += math.log(this_pwm[i][pwm_map[snp_alt]]/0.25)
                            else:
                                alt_score += math.log(this_pwm[i][pwm_map[c]]/0.25)

                        ## write out all the motif information with scoring (the score should
                        ## be relative to the reference, i.e. we want to measure the effect of
                        ## the alternate allele so if it increases binding, the delta score
                        ## should be positive)
                        motif_out.write("\t".join(snp_data+motif+[str(relative_pos), str(ref_score), str(alt_score), str(alt_score-ref_score), str(ref_prob), str(alt_prob), this_seq])+"\n")
                    else:
                        ## otherwise, just write out the relative position
                        motif_out.write("\t".join(snp_data+motif+[str(relative_pos)])+"\n")
                        
    if motif_pwm_file:
        motif_seqs.close()

def dashr_ncrna_overlap(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, dashr_locus_file, logging_function):
    """
    for each SNP, look for overlap with ncRNA loci defined in DASHR
    """
    try:
        os.makedirs(outdir+"/dashr_ncrna_loci_overlap/")
    except OSError:
      pass
  
    outf = outdir+"/dashr_ncrna_loci_overlap/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)+"_dashr_locus_overlap.txt"
    
    with open(ld_snp_file, 'rb') as snpin, open(dashr_locus_file, 'rb') as dashr_loci, open(outf, 'wb') as ncrna_out:
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}
        
        ## manually define the bed file index
        bed_idx = {"chr":0, "start":1, "end":2, "dashr_info":3, "rna_type":4, "strand":5}

        # write a header
        ncrna_out.write("\t".join(snp_header_data + ['ncrna_chr', 'ncrna_start', 'ncrna_end', 'ncrna_strand', 'ncrna_info', 'rna_type'])+"\n")

        ## track which chromosome we're on
        this_chr = ''

        ## a regex to check that we are in a 'real' chromosome
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        

        ## initialize storage for the ncRNA loci
        cur_ncrna_loci = []
        cur_ncrna = dashr_loci.readline().strip().split('\t')
        
        for entry in snpin:
            entry_data = entry.strip().split('\t')
            entry_pos = int(entry_data[snp_idx['pos']])
            this_rsID = entry_data[snp_idx['rsID']]
            this_tag_rsID = entry_data[snp_idx['tag_rsID']]
            
            ## if we're on a new chromosome, reset all the tracking variables:
            if entry_data[snp_idx['chr']] != this_chr:
                this_chr = entry_data[snp_idx['chr']]
                logging_function('Parsing chromosome '+this_chr)
                cur_ncrna_loci = []
                ## read in until we match
                while (len(cur_ncrna)>1 and cur_ncrna[bed_idx['chr']] != this_chr):
                    cur_ncrna = dashr_loci.readline().strip().split('\t')
                    
            ## if we have a match to one of the weird chromosomes, just skip to the next entry
            if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                continue                    
                    
            ## now add all possibly overlapping entries
            while (len(cur_ncrna)>1 and cur_ncrna[bed_idx['chr']]==entry_data[snp_idx['chr']]
                   and int(cur_ncrna[bed_idx['start']]) <= entry_pos):
                if (int(cur_ncrna[bed_idx['end']]) >= entry_pos):
                    cur_ncrna_loci.append(cur_ncrna)
                cur_ncrna = dashr_loci.readline().strip().split('\t')

            ## remove all too-early entries here
            for i in range(len(cur_ncrna_loci)-1, -1, -1):
                if (cur_ncrna_loci[i][bed_idx['chr']] != entry_data[snp_idx['chr']]
                    or int(cur_ncrna_loci[i][bed_idx['end']]) <= entry_pos):
                    cur_ncrna_loci.pop(i)

            ## calculate overlaps, write them out directly
            for locus in cur_ncrna_loci:
                if (int(locus[bed_idx['start']]) <= entry_pos < int(locus[bed_idx['end']])):
                    ## write out the info
                    ## need to reorder the locus to have chr,start,end,strand,info,type
                    ncrna_out.write("\t".join(entry_data +
                                              [str(locus[bed_idx[x]]) for x in ["chr", "start", "end", "strand", "dashr_info", "rna_type"]])+"\n")

    return outf

def targetscan_overlap(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, targetscan_dir, logging_function):
    """
    for each SNP, look for overlap with miRNA seed-based target sites from TargetScan
    """
    try:
        os.makedirs(outdir+"/targetscan_miRNA_overlap/")
    except OSError:
      pass

    outf = outdir+"/targetscan_miRNA_overlap/"+outprefix+"_"+str(ld_threshold)+"_ld_cutoff_snps_within_"+str(ld_check_area)+"_targetscan_overlap.txt"

    # define a convenience dictionary to make subsetting clearer
    # some of these are special fields for targetscan
    bed_idx = {} 
    bed_idx["chr"] = 0
    bed_idx["start"] = 1
    bed_idx["end"] = 2    
    bed_idx["target_info"] = 3
    bed_idx["context_score_percentile"] = 4
    bed_idx["strand"] = 5
    bed_idx["thickStart"] = 6
    bed_idx["thickEnd"] = 7
    bed_idx["itemRgb"] = 8
    bed_idx["blockCount"] = 9
    bed_idx["blockSizes"] = 10
    bed_idx["blockStarts"] = 11

    ## make a list to store the different targetScan conditions
    targetscan_types = []

    # make a dict for all the targetscan files
    targetscan_files = {}
    # also make a dict to store the lists of all current miRNA seeds
    targetscan_loci = {}
    # finally, also have one to store the current target in each condition
    cur_target = {}
    for targetscan_f in glob.glob(targetscan_dir+"*.bed"):
        ## get the condition (type of targetscan result) from this
        this_key = ".".join(os.path.basename(targetscan_f).split(".")[2:4])
        targetscan_types.append(this_key)

        ## open the file and associate it with the key        
        targetscan_files[this_key] = open(targetscan_f, 'rb')
        ## initialize the list of miRNA seed loci
        targetscan_loci[this_key] = []
        ## initialize the current miRNA locus storage (this will get added to the list later)
        cur_target[this_key] = targetscan_files[this_key].readline().strip().split("\t")

    ## for now, just write the one output file, containing the miRNA target overlaps
    with open(ld_snp_file, 'rb') as snpin, open(outf, 'wb') as targetscan_overlaps_out:        
        # get the LD SNP header
        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

        ## write the header
        targetscan_overlaps_out.write("\t".join(snp_header_data + ['targetscan_condition', 'target_chr', 'target_start', 'target_end', 'target_info', 'context_score_percentile', 'strand'])+"\n")
        
        ## track which chromosome we're on
        this_chr = ''

        ## a regex to check that we are in a 'real' chromosome
        nonref_chr_pattern = re.compile("^chr.*\_.*$")                        

        for entry in snpin:
            entry_data = entry.strip().split('\t')
            entry_pos = int(entry_data[snp_idx['pos']])
            this_rsID = entry_data[snp_idx['rsID']]
            this_tag_rsID = entry_data[snp_idx['tag_rsID']]
            
            ## if we're on a new chromosome, reset all the tracking variables:
            if entry_data[snp_idx['chr']] != this_chr:
                this_chr = entry_data[snp_idx['chr']]
                logging_function('Parsing chromosome '+this_chr)
                        
                for k in targetscan_types:
                    # reset the enhancer lists
                    targetscan_loci[k] = []

                if nonref_chr_pattern.match(this_chr) or this_chr=='chrM':
                    ## for now, just skip this entry..
                    continue                    
                else:
                    # read in data until we match chromosomes
                    for k in targetscan_types:
                        while(len(cur_target[k]) > 1
                            and cur_target[k][bed_idx['chr']] != this_chr):
                            cur_target[k] = targetscan_files[k].readline().strip().split("\t")

            ## now we go through all the targetscan conditions looking for overlaps
            ## store the number of overlaps in each condition (currently not using this)
            num_targetscan_overlaps = {}
                        
            for k in targetscan_types:                
                ## add all the new overlapping entries
                while(len(cur_target[k]) > 1
                      and cur_target[k][bed_idx['chr']]==entry_data[snp_idx['chr']]):
                    cur_target_start = int(cur_target[k][bed_idx['start']])
                    cur_target_end = int(cur_target[k][bed_idx['end']])
                    
                    if cur_target_start <= entry_pos:
                        if cur_target_end >= entry_pos:
                            targetscan_loci[k].append(cur_target[k])
                        cur_target[k] = targetscan_files[k].readline().strip().split("\t")
                    # if we're not overlapping anymore, stop the loop
                    else:
                        break
                        
                ## remove all loci from the wrong chromosome or that end before this starts
                for i in range(len(targetscan_loci[k])-1, -1, -1):
                    if targetscan_loci[k][i][bed_idx['chr']] != this_chr:
                        targetscan_loci[k].pop(i)
                    else:                            
                        targetscan_start = int(targetscan_loci[k][i][bed_idx['start']])
                        targetscan_end = int(targetscan_loci[k][i][bed_idx['end']])
                        if targetscan_end  <= entry_pos:
                            targetscan_loci[k].pop(i)                            

                ## initialize the storage dict for the number of overlaps in the condition
                num_targetscan_overlaps[k] = 0
                            
                ## now compare the SNP to the possible miRNA loci in this tissue type
                ## note that these are all guaranteed to overlap because of the above section
                ## i leave the redundant check in for now
                for target in targetscan_loci[k]:
                    target_start = int(target[bed_idx['start']])
                    target_end = int(target[bed_idx['end']])

                    if (entry_pos >= target_start and entry_pos <= target_end):
                        num_targetscan_overlaps[k] += 1
                        ## write out this overlap
                        targetscan_overlaps_out.write("\t".join(entry_data + [k, target[bed_idx['chr']],str(target_start), str(target_end)] + [target[bed_idx[x]] for x in ["target_info", "context_score_percentile", "strand"]])+"\n")
            
    # close all the targetscan files
    for k in targetscan_files.keys():
        targetscan_files[k].close()    
        
    # return the file with one line per targetscan overlap
    return outf
       
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Expand tag SNPs into LD block sets and annotate with regulatory information. Note that all annotation files MUST be sorted lexicographically i.e. using sort -k1,1V -k2,2n.")
    parser.add_argument("--kg_pop", default="EUR", help="The population from 1000 genomes to use")
    parser.add_argument("--skip_ld_expansion", action='store_true', help="If you don't want to do any LD expansion, give this flag to just do the functional overlaps of your input variants.")
    parser.add_argument("--loglevel", type=str.lower, choices=['full', 'print', 'save', 'quiet'], default="full", help="The amount of logging you want. Options: 'full' to print all messages and write to log file, 'print' to print only, 'save' to save output and only print high-level messages, and 'quiet' to skip all output")
    parser.add_argument("--ld_threshold", type=float, default=0.7, help="The threshold for LD with the tag SNPs")
    parser.add_argument("--ld_check_area", type=int, default=500000, help="How many base pairs in either direction to check for SNPs in LD with the tag SNPs")
    parser.add_argument("--ld_buffer_size", type=int, help="Number of SNPs to perform LD expansion on at a time. Useful if you have a long list of input SNPs.")
    parser.add_argument("--force_ld_recalc", action='store_true', help="Give this flag to force recalculation of the LD SNPs, even if you already ran this step.")
    parser.add_argument("--remove_ld_data", action='store_true', help="Give this flag if you want the pipeline to delete the extra data for LD expansion (the neighboring SNP file for each SNP and the full LD files). This saves a LOT of space.")
    parser.add_argument("--gene_bed_file", help="The bed file containing genic coordinates. Optional; if provided, closest genes and distance will be calculated for each SNP.")
    parser.add_argument("--kgxref_file", help="If computing closest genes, you can also give this optional file containing cross-references between IDs and gene names. NB: Only refseq cross-reference is currently supported!")
    parser.add_argument("--partition_dir", help="The directory containing the genomic partitioning files to be used, including 5' and 3' UTR exons and introns, mRNA exons and introns, promoters, repeats, and intergenic regions. If provided, the genomic localization of each SNP will be calculated")
    parser.add_argument("--unstranded_partition_dir", help="The directory containing the genomic partitioning files to be used without regard to strand, including 5' and 3' UTR exons and introns, mRNA exons and introns, promoters, repeats, and intergenic regions. If provided, the genomic localization of each SNP will be calculated. This can be provided in addition to the partition_dir argument")
    parser.add_argument("--bedtools_bin_dir", help="If you don't have bedtools in your path, you can provide the path to your own bin/ directory for bedtools with this argument.")
    parser.add_argument("--plink_path", help="If you don't have plink (>=v1.9) in your path, specify the direct path to the plink executable with this argument.")    
    parser.add_argument("--fantom5_dir", help="The path to the directory containing the sorted .bed files representing FANTOM5 enhancers.")
    parser.add_argument("--enhancer_midpoint_window", type=int, help="How far around the midopint of enhancer transcription you want to look (on each side of the midpoint).")
    parser.add_argument("--enhancer_locus_window", type=int, help="How far you want to extend the original enhancer loci for overlaps (set to 0 to use original loci).")
    parser.add_argument("--skip_closest_enh", action="store_true", help="If you don't want to compute the closest enhancers of each SNP.")
    parser.add_argument("--skip_enh_summary", action="store_true", help="If you want to skip generating the enhancer summary file (rows are SNPs, columns are # of enhancer overlaps and positions, for all FANTOM5 tissues, so it is an n x m table, where n is the number of SNPs and m is the number of FANTOM5 tissue/cell types (currently 116). This saves a lot of space, but means you lose some downstream analysis!")
    parser.add_argument("--fantom5_correlation_file", help="The path to the bed file containing the correlation-based target gene predictions for all FANTOM5 enhancers.")
    parser.add_argument("--gtex_dir", help="The path to the directory containing the sorted GTEx (v6) data of SNP-gene associations (files ending with .snpgenes).")
    parser.add_argument("--factorbook_file", help="The path to the file containing the FactorBook motifs (this should be the TFBS clusters with the input cell sources, i.e. wgEncodeRegTfbsClusteredWithCellsV3.bed).")
    parser.add_argument("--roadmap_chromhmm_dir", help="The path to the directory containing the ChromHMM Roadmap Mnemonic bed files (must all be sorted!)")
    parser.add_argument("--homer_motif_bed_file", help="The path to the sorted (with sort -k1,1V -k2,2n) bed file containing the results of HOMER motif scanning across the whole genome.")
    parser.add_argument("--homer_motif_pwm_file", help="The path to the file in HOMER motif format (http://homer.salk.edu/homer/motif/creatingCustomMotifs.html) containing the motifs and PWMs represented in the occurrence file. This is optional given the motif bed file but can be used to get more informative motif overlap results.")
    parser.add_argument("--homer_motif_seq_file", help="The file containing the sequences of the motifs from the HOMER motif file. Need both this file and the homer_motif_pwm_file to do PWM calculations! (ADD IN DESCRIPTION OF SCRIPT TO GENERATE FILE)")
    parser.add_argument("--dashr_locus_file", help="The bed file containing small ncRNA loci from DASHR. Must be sorted with -k1,1V -k2,2n")
    parser.add_argument("--targetscan_dir", help="The path to the directory containing the sorted (-k1,1V -k2,2n) .bed files containing targetScan predictions (http://www.targetscan.org/vert_71/vert_71_data_download/All_Target_Locations.hg19.bed.zip), named so the last two fields describe the miRNA family and site conservation.")
    parser.add_argument("kg_dir", help="The path to the directory containing 1000 genomes data")
    parser.add_argument("input_snp_list", help="The tab separated file containing the tag SNPs you want to analyze. Should be formatted with four columns: chromosome, rsID, region naming information, and position (in that order). IMPORTANT: Note that SNPs without dbSNP rsIDs should use 'chr-pos' naming format, not 'chr:pos', which is incompatible with this pipeline!")
    parser.add_argument("outdir", help="The directory to write the results to.")
    parser.add_argument("outprefix", help="The prefix of the filename to use to write results.")
    pargs = parser.parse_args()
            
    try:
        os.makedirs(pargs.outdir)
    except OSError:
        pass

    now = datetime.datetime.now()
        
    ## write the parameter file
    try:
        os.makedirs(pargs.outdir+"/parameters/")
    except OSError:
        pass

    date_string = "_".join([str("%02d" % now.month), str("%02d" % now.day), str(now.year)])+"_"+":".join(["%02d" % now.hour, "%02d" % now.minute, "%02d" % now.second])
    with open(pargs.outdir+"/parameters/"+date_string+"_parameters.txt", 'wb') as pargs_out:
        args_dict = vars(pargs)
        for k in args_dict.keys():
            pargs_out.write("\t".join([k, str(args_dict[k])])+'\n')

    with open(pargs.outdir+"/parameters/"+date_string+"_command_run.txt", "wb") as cmd_out:
        cmd_out.write(" ".join(["python"] + sys.argv)+"\n")
        
    ## also make the logging file    
    if pargs.loglevel != "print" and pargs.loglevel != "quiet":
        try:
            os.makedirs(pargs.outdir+"/logs/")
        except OSError:
            pass
        logfile = open(pargs.outdir+"/logs/"+date_string+"_log.txt", "wb")
    else:
        logfile = ""

    logging_function = create_logging_function(logfile, pargs.loglevel)
        
    ## first expand into LD SNPs:
    if pargs.skip_ld_expansion:
        if pargs.loglevel=="save":
            print "Skipping LD expansion! Setting LD parameters to R^2 = 1, distance = 0"
        logging_function("Skipping LD expansion! Setting LD parameters to R^2 = 1, distance = 0")
        pargs.ld_threshold = 1.0
        pargs.ld_check_area = 0
        ## definitely don't want to save the LD data, set this just in case
        pargs.remove_ld_data = True

        ## feed this to the normal LD expansion method, but it will skip some of the processing 
        ld_snp_file = expand_ld_snps(pargs.input_snp_list, pargs.kg_pop, pargs.kg_dir,
                                     pargs.ld_threshold, pargs.ld_check_area, pargs.outdir,
                                     pargs.outprefix, logging_function, pargs.ld_buffer_size,
                                     pargs.plink_path, pargs.force_ld_recalc, pargs.remove_ld_data)
    else:
        if pargs.loglevel=="save":
            print "Calculating LD for input SNPs, with threshold %s and distance %s" % (str(pargs.ld_threshold), str(pargs.ld_check_area))
        logging_function("Calculating LD for input SNPs, with threshold %s and distance %s" % (str(pargs.ld_threshold), str(pargs.ld_check_area)))
        ld_snp_file = expand_ld_snps(pargs.input_snp_list, pargs.kg_pop, pargs.kg_dir,
                                     pargs.ld_threshold, pargs.ld_check_area, pargs.outdir,
                                     pargs.outprefix, logging_function, pargs.ld_buffer_size,
                                     pargs.plink_path, pargs.force_ld_recalc, pargs.remove_ld_data)

    if pargs.gene_bed_file:
        if pargs.loglevel=="save":
            print "Computing closest genes"        
        logging_function("Computing closest genes.")
        ## assumes that the bedtools wrapper script is in the same directory (last arg)
        closest_gene_outf = compute_closest_genes(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.gene_bed_file, pargs.kgxref_file, pargs.bedtools_bin_dir, os.path.dirname(os.path.realpath(sys.argv[0])), logging_function)

    if pargs.partition_dir:
        if pargs.loglevel=="save":
            print "Computing genomic partition"        
        logging_function("Computing genomic partition")
        partition_outf = calculate_genomic_partition(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.partition_dir, logging_function)

    if pargs.unstranded_partition_dir:
        if pargs.loglevel=="save":
            print "Computing unstranded genomic partition"                
        logging_function("Computing unstranded genomic partition")
        unstranded_partition_outf = calculate_unstranded_genomic_partition(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.unstranded_partition_dir, logging_function)
        
    if pargs.fantom5_dir:
        if pargs.loglevel=="save":
            print "Performing FANTOM5 enhancer analysis"                
        logging_function("Performing FANTOM5 enhancer analysis")
        if pargs.enhancer_midpoint_window:
            if pargs.loglevel=="save":
                print "Computing FANTOM5 overlap using window around midpoint of transcription."
            logging_function("Computing FANTOM5 overlap using window around midpoint of transcription.")
            enh_midpoint_overlap_outf = compute_fantom5_midpoint_overlap(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.fantom5_dir, pargs.enhancer_midpoint_window, pargs.skip_enh_summary, logging_function)
            
        if pargs.enhancer_locus_window:
            if pargs.loglevel=="save":
                print "Computing FANTOM5 overlap using window around FANTOM5 transcription locus."
            logging_function("Computing FANTOM5 overlap using window around FANTOM5 transcription locus.")
            enh_locus_overlap_outf = compute_fantom5_locus_overlap(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.fantom5_dir, pargs.enhancer_locus_window, pargs.skip_enh_summary, logging_function)

        if not pargs.skip_closest_enh:
            if pargs.loglevel=="save":
                print "Computing closest enhancers"
            logging_function("Computing closest enhancers")
            closest_enh_outf = compute_closest_enhancer(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.fantom5_dir, logging_function)
            
        if pargs.fantom5_correlation_file:
            if not pargs.skip_closest_enh:
                if pargs.loglevel=="save":
                    print "Identifying correlation-based target genes for closest enhancers"
                logging_function("Identifying correlation-based target genes for closest enhancers")
                correlation_closest_enh_targets(pargs.outdir, pargs.outprefix, pargs.ld_threshold, pargs.ld_check_area, pargs.fantom5_correlation_file, closest_enh_outf, logging_function)
            ## now do analysis for whichever overlaps we had
            if pargs.enhancer_midpoint_window:
                if pargs.loglevel=="save":
                    print "Identifying correlation-based target genes for enhancer overlaps (by midpoint)"
                logging_function("Identifying correlation-based target genes for enhancer overlaps (by midpoint)")
                correlation_enh_overlap_targets(pargs.outdir, pargs.outprefix, pargs.ld_threshold, pargs.ld_check_area, pargs.fantom5_correlation_file, enh_midpoint_overlap_outf, logging_function, overlap_type="midpoint")
            if pargs.enhancer_locus_window:
                if pargs.loglevel=="save":
                    print "Identifying correlation-based target genes for enhancer overlaps (by window around locus of transcription)"                
                logging_function("Identifying correlation-based target genes for enhancer overlaps (by window around locus of transcription)")
                correlation_enh_overlap_targets(pargs.outdir, pargs.outprefix, pargs.ld_threshold, pargs.ld_check_area, pargs.fantom5_correlation_file, enh_locus_overlap_outf, logging_function, overlap_type="locus")
                
    if pargs.gtex_dir:
        if pargs.loglevel=="save":
            print "Computing overlap with GTEx eQTLs"
        logging_function("Computing overlap with GTEx eQTLs")
        eqtl_outf = gtex_eqtl_overlap(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.gtex_dir, logging_function)

    if pargs.factorbook_file:
        if pargs.loglevel=="save":
            print "Computing overlap with FactorBook TF Motifs"        
        logging_function("Computing overlap with FactorBook TF Motifs")
        factorbook_outf = factorbook_overlap(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.factorbook_file, logging_function)
                                       
    if pargs.roadmap_chromhmm_dir:
        if pargs.loglevel=="save":
            print "Computing overlap with ChromHMM states for Roadmap tissues"        
        logging_function("Computing overlap with ChromHMM states for Roadmap tissues")
        roadmap_chromhmm_outf = roadmap_chromhmm_overlap(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.roadmap_chromhmm_dir, logging_function)

    if pargs.homer_motif_bed_file:
        if pargs.loglevel=="save":
            print "Computing overlap with HOMER-predicted motifs"        
        logging_function("Computing overlap with HOMER-predicted motifs")
        if ((pargs.homer_motif_pwm_file and not pargs.homer_motif_seq_file) or
            (not pargs.homer_motif_pwm_file and pargs.homer_motif_seq_file)):
            logging_function("Error: Need both PWM and motif sequence file to do PWM calculations! Skipping PWM calculations")
            homer_motif_outf = homer_motif_overlap(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.homer_motif_bed_file, None, None, logging_function)
        ## otherwise, just use the arguments directly
        else:
            homer_motif_outf = homer_motif_overlap(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.homer_motif_bed_file, pargs.homer_motif_pwm_file, pargs.homer_motif_seq_file, logging_function)            

    if pargs.dashr_locus_file:
        if pargs.loglevel=="save":
            print "Computing overlap with DASHR ncRNA loci"        
        logging_function("Computing overlap with DASHR ncRNA loci")
        dashr_loci_outf = dashr_ncrna_overlap(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.dashr_locus_file, logging_function)

    if pargs.targetscan_dir:
        if pargs.loglevel=="save":
            print "Computing overlap with TargetScan miRNA targets"        
        logging_function("Computing overlap with TargetScan miRNA targets")
        targetscan_outf = targetscan_overlap(pargs.outdir, pargs.outprefix, ld_snp_file, pargs.ld_threshold, pargs.ld_check_area, pargs.targetscan_dir, logging_function)
                               
    ## close the log, if it exists
    if pargs.loglevel != "print" and pargs.loglevel != "quiet":
        logfile.close()
            
