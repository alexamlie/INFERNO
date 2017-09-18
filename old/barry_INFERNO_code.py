import pandas as pd, cPickle, time

def compute_fantom5_locus_overlap_hash(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, fantom5_dir, enhancer_window, skip_enh_summary, logging_function):
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
    bed_idx = dict()
    bed_idx["chr"] = 0
    bed_idx["start"] = 1
    bed_idx["end"] = 2    
    bed_idx["name"] = 3

    million = 1000000
    billion = 1000000000

    if os.path.exists("../pickle/encode_fantom5_enhancers.pickle"):

        with open('../pickle/encode_fantom5_enhancers.pickle', 'rb') as handle:
            saved = cPickle.load(handle)
            enh_dict = saved[0]
            enh_keys = saved[1]
            enh_keys_dict = saved[2]

    else:

        enh_dict = dict()
        enh_keys = []
        enh_keys_dict = dict()
        enh_keys_idx = 0

        for enh_f in glob.glob(fantom5_dir+"*.bed"):
            ## get just the ID and the name
            this_key = os.path.basename(enh_f).replace("_expressed_enhancers.bed", "")

            enh_keys.append(this_key)
            enh_keys_dict[this_key] = enh_keys_idx
            enh_keys_idx += 1

            # open the enhancer file
            enh_in = open(enh_f, 'rb')

            # read the first three cols into a pandas dataframe
            # enh_df = pd.read_csv(enh_in, sep='\t', header=None, names=["enh_chr", "enh_start", "enh_end"], usecols=[0,1,2])
            line = enh_in.readline()
            while( len(line) > 0 ):
                line = line.strip().split("\t")
                line[ bed_idx["name"] ] = this_key # overwrite

                chr_num = line[bed_idx["chr"]].replace("chr","")
                if chr_num == "X":
                    chr_num = 25
                elif chr_num == "Y":
                    chr_num = 26
                elif chr_num == "M":
                    chr_num = 27 # handle better?
                else:
                    chr_num = int(chr_num)

                low_hash = int(chr_num)*billion + int(math.floor( (float(line[bed_idx["start"]]) - enhancer_window) /million))*million
                if low_hash not in enh_dict:
                    enh_dict[low_hash] = [line[0:4]]
                else:
                    enh_dict[low_hash] = enh_dict[low_hash] + [line[0:4]]

                high_hash = int(chr_num)*billion + int(math.floor( (float(line[bed_idx["end"]]) + enhancer_window) /million))*million
                if high_hash != low_hash: #unusual but possible
                    if high_hash not in enh_dict:
                        enh_dict[high_hash] = [line]
                    else:
                        enh_dict[high_hash] = enh_dict[high_hash] + [line]

                line = enh_in.readline()

            enh_in.close()

        with open('../pickle/encode_fantom5_enhancers.pickle', 'wb') as handle:
            cPickle.dump([enh_dict, enh_keys, enh_keys_dict], handle, protocol=pickle.HIGHEST_PROTOCOL)

    #t0 = time.time()
    with open(ld_snp_file, 'rb') as snpin, open(overlap_outprefix+"_enh_overlaps.txt", 'wb') as enh_overlaps_out:

        if not skip_enh_summary:
             summary_out = open(overlap_outprefix+"_summary.txt", 'wb')

        snp_header_data = next(snpin).strip().split('\t')
        snp_idx = {snp_header_data[x]:x for x in range(len(snp_header_data))}

        if not skip_enh_summary:
            summary_out.write("\t".join(snp_header_data+[t+"_num_overlaps\t"+t+"_dists" for t in enh_keys])+"\n")
        enh_overlaps_out.write("\t".join(snp_header_data + ['enh_source', 'enh_chr', 'enh_start', 'enh_end'])+"\n")

        summ_line_size = 2*len(enh_keys_dict.keys())

        for entry in snpin:

            entry_data = entry.strip().split('\t')
            entry_chr = int(entry_data[snp_idx["chr"]].replace("chr",""))
            entry_pos = int(entry_data[snp_idx['pos']])
            this_rsID = entry_data[snp_idx['rsID']]
            this_tag_rsID = entry_data[snp_idx['tag_rsID']]

            entry_hash = entry_chr*billion + int(math.floor( entry_pos / million ) )*million

            enh_dists = []
            summary_line = ["0"]*summ_line_size

            if entry_hash in enh_dict:
                for enh in enh_dict[entry_hash]:

                    enh_start = int(enh[bed_idx['start']])
                    enh_end = int(enh[bed_idx['end']])

                    if (entry_pos >= enh_start - enhancer_window and entry_pos <= enh_end + enhancer_window):

                        midp = (enh_start + enh_end)/2
                        dist = midp - entry_pos

                        enh_source = enh[bed_idx["name"]]
                        enh_idx = 2*enh_keys_dict[ enh[bed_idx["name"]] ]

                        summary_line[enh_idx] = str(int(summary_line[enh_idx])+1) # update number of overlaps
                        if summary_line[enh_idx+1] != "0":
                            summary_line[enh_idx+1] = summary_line[enh_idx+1] + ","+str(dist)
                        else:
                            summary_line[enh_idx+1] = str(dist)

                        enh_overlaps_out.write("\t".join(entry_data + [enh_source, enh[bed_idx['chr']],str(enh_start), str(enh_end)])+"\n")                        

            if not skip_enh_summary:
                summary_out.write( "\t".join(entry_data + summary_line)+"\n" )

        snpin.close()
        if not skip_enh_summary:
            summary_out.close()

    #t1 = time.time()
    #print("\n" + str(t1-t0) + "\n")
    # return the file with one line per enhancer overlap
    return overlap_outprefix+"_enh_overlaps.txt"


def is_present(c):
    if pd.isnull(c): 
        return 0
    return 1

def compute_fantom5_locus_overlap_dframes(outdir, outprefix, ld_snp_file, ld_threshold, ld_check_area, fantom5_dir, enhancer_window, skip_enh_summary, logging_function):
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

    with open(ld_snp_file, 'rb') as snpin, open(overlap_outprefix+"_enh_overlaps.txt", 'wb') as enh_overlaps_out:

        if not skip_enh_summary:
             summary_out = open(overlap_outprefix+"_summary.txt", 'wb')

        # Note: to better manage memory we could have an outer loop where we read in 10K lines at a time from the snp file
        # For now we just read in the whole thing
        snps_df = pd.read_csv(snpin, sep='\t', header='infer')

        # Is "pos" a hard-coded column name for the snp offset?
        snps_df["pos"] = snps_df["pos"].astype(int)

        million = 1000000
        snps_df["pos_sub"] = [ round(c/million)*million for c in snps_df["pos"].astype(float) ]
        snps_df["idx"] = snps_df.index.copy()
        snps_df["idx2"] = snps_df.index.copy()

        # initialize variable
        all_enh_df = None
        enh_keys = []

        for enh_f in glob.glob(fantom5_dir+"*.bed"):
            ## get just the ID and the name
            this_key = os.path.basename(enh_f).replace("_expressed_enhancers.bed", "")

            # open the enhancer file
            enh_in = open(enh_f, 'rb')

            # read the first three cols into a pandas dataframe
            enh_df = pd.read_csv(enh_in, sep='\t', header=None, names=["enh_chr", "enh_start", "enh_end"], usecols=[0,1,2])

            # Add a new column named "source" with uniform value this_key for all entries
            enh_df.insert(0, "enh_source", [this_key]*len(enh_df["enh_chr"]))

            # concatenate this result to the accumulated record of results
            if all_enh_df is None: # first pass
                all_enh_df = enh_df.copy()
            else:
                all_enh_df = pd.concat([all_enh_df, enh_df])

            enh_keys.append(this_key)
            enh_in.close()

        # join the snps and enhancer dfs on the "chr" and "pos_sub" columns.
        all_enh_df["enh_start"] = all_enh_df["enh_start"].astype(int)
        all_enh_df["enh_end"] = all_enh_df["enh_end"].astype(int)

        all_enh_df["enh_start_sub"] = [ round(c/million)*million for c in all_enh_df["enh_start"].astype(float) - enhancer_window ]
        all_enh_df["enh_end_sub"] = [ round(c/million)*million for c in all_enh_df["enh_end"].astype(float) + enhancer_window ]

        result_df = snps_df.merge(all_enh_df, how='inner', left_on=["chr", "pos_sub"], right_on=["enh_chr", "enh_start_sub"])
        diff_df = all_enh_df[ all_enh_df["enh_start_sub"] != all_enh_df["enh_end_sub"] ]
        result_df = pd.concat([result_df, snps_df.merge(diff_df , how='inner', left_on=["chr", "pos_sub"], right_on=["enh_chr", "enh_end_sub"])] )

        # Limit the snp records in the joined table to the snps overlapping an enhancer window
        b1 = result_df["pos"] >= result_df["enh_start"] - enhancer_window
        b2 = result_df["pos"] <= result_df["enh_end"] + enhancer_window
        result_df = result_df[(b1 & b2)]

        result_df.insert(0, "chr_num", [c.replace("chr","") for c in result_df["chr"]])
        result_df["chr_num"] = result_df["chr_num"].astype(int)
        result_df = result_df.sort_values(["chr_num", "pos"])
        result_df = result_df.drop("chr_num", 1)

        summary_df = snps_df.copy()

        for this_key in enh_keys:

            merged_df = result_df[ result_df["enh_source"] == this_key ]

            #add to summary table
            #temp_df = snps_df[["rsID", "tag_rsID", "pos"]].merge(merged_df[["rsID","tag_rsID","enh_start","enh_end"]], how='left', on=["rsID", "tag_rsID"])
            temp_df = snps_df[["idx2", "rsID", "tag_rsID", "pos"]].join( merged_df[["idx", "enh_start","enh_end"]].set_index("idx") ) # merge on index

            dists_str = this_key+"_dists"

            temp_df[dists_str] = (temp_df["enh_start"] + temp_df["enh_end"])/2 - temp_df["pos"]
            temp_df[dists_str] = temp_df[dists_str].fillna(0).astype(int)

            overlaps_str = this_key+"_overlaps"

            temp_df[overlaps_str] = map(is_present, temp_df["enh_start"])
            temp_df[overlaps_str] = temp_df[overlaps_str].astype(int)

            temp_df_2 = temp_df[ temp_df[overlaps_str] > 0 ]
            temp_df_2 = temp_df_2.groupby(["rsID", "tag_rsID"]).count().reset_index()

            multiple = list(temp_df_2[ temp_df_2[overlaps_str] > 1 ]["rsID"].copy())

            if(len(multiple)==0):
                pass # no change to temp_df
            else:
                pt1_df = temp_df[ [i not in multiple for i in temp_df["rsID"]] ]
                
                pt2_df = temp_df[ [i in multiple for i in temp_df["rsID"]] ].groupby(["rsID", "tag_rsID"]).aggregate(lambda x: tuple(x)).reset_index()
                pt2_df[overlaps_str] = map(sum, pt2_df[overlaps_str])
                pt2_df[dists_str] = [ ",".join( map(str, tup) ) for tup in pt2_df[dists_str] ]
                pt2_df["idx2"] = [ tup[0] for tup in pt2_df["idx2"]]
                pt2_df["idx2"] = pt2_df["idx2"].astype(int)

                temp_df = pd.concat( [pt1_df, pt2_df] )

            temp_df = temp_df.set_index("idx2")

            # THIS IS THE SLOW PART remaining
            #summary_df = summary_df.merge(temp_df[["rsID", "tag_rsID", overlaps_str, dists_str]], how='inner', on=["rsID", "tag_rsID"])
            summary_df = summary_df.join( temp_df[[overlaps_str, dists_str]] )

        # drop added columns
        result_df = result_df.drop(["idx2", "idx", "enh_start_sub", "enh_end_sub", "pos_sub"], 1)
        summary_df = summary_df.drop(["idx2", "idx", "pos_sub"], 1)

        # write out results
        result_df.to_csv(enh_overlaps_out, sep='\t', index=False)
        enh_overlaps_out.close()

        if not skip_enh_summary:
            summary_df.to_csv(summary_out, sep='\t', index=False)
            summary_out.close()
        
    # return the file with one line per enhancer overlap
    return overlap_outprefix+"_enh_overlaps.txt"


