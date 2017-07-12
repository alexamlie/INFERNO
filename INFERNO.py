"""
INFERNO.py
Alex Amlie-Wolf, 07/10/17

This is the main driver method for running the INFERNO pipeline
"""

import argparse, subprocess, datetime, os, time, sys

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Driver script for the INFERNO pipeline")
    ## (optional) arguments for summary statistics analysis
    parser.add_argument("--rsid_column", type=int, help="The summary statistics column number containing the rsIDs")
    parser.add_argument("--pos_column", type=int, help="The summary statistics column number containing the positions")
    parser.add_argument("--pval_column", type=int, help="The summary statistics column number containing the p-values")
    parser.add_argument("--chr_column", type=int, help="The summary statistics column number containing the chromosomes")
    parser.add_argument("--allele1_column", type=int, help="The summary statistics column number containing allele 1 (TODO: CHECK THAT THIS IS MAJOR)")
    parser.add_argument("--allele2_column", type=int, help="The summary statistics column number containing allele 2 (TODO: CHECK THAT THIS IS MAJOR)")
    parser.add_argument("--maf_column", type=int, help="The summary statistics column number containing the minor allele frequency. Note that if this is ever >50% and a beta column is provided, the effect direction will be flipped to be defined relative to the minor allele.")
    parser.add_argument("--beta_column", type=int, help="The summary statistics column number containing the beta estimate (used for p-value expansion with consistent directions)")        
    ## for p-value expansion and pruning
    parser.add_argument("--sig_mult", type=float, help="The multiplier range for significance of the p-value expanded variant set (e.g. a value of 10 means one order of magnitude)")
    parser.add_argument("--skip_pval_expansion", action='store_true', help="If you don't want to do any expansion by p-values, give this flag to just expand your input variants when doing summary stats-based analyses.")
    parser.add_argument("--case_prop", type=float, help="The proportion of cases in the GWAS, for colocalization.")
    parser.add_argument("--sample_size", type=int, help="The total number of samples in the GWAS, for colocalization.")
    parser.add_argument("--summary_file", help="The path to the full summary statistics file.")
    ## required arguments 
    parser.add_argument("top_snpf", help="The tab separated file of the tag SNPs you want to analyze. Should be formatted with four columns: chromosome, rsID, region naming information, and position (in that order). IMPORTANT: Note that SNPs without dbSNP rsIDs should use 'chr-pos' naming format, not 'chr:pos', which is incompatible with this pipeline!")
    parser.add_argument("cfg_file", help="The configuration file containing paths to all the required functional annotation files. Should be formatted as a bash configuration file i.e. VARIABLE=DEFINTIION on each line.")
    parser.add_argument("outdir", help="The directory to write all the results to.")
    parser.add_argument("outprefix", help="The desired prefix for all the output file names")
    
    pargs = parser.parse_args()    

    ## source the variables from the config file
    config_vars = {}
    with open(pargs.cfg_file, 'rb') as cfg_file:
        for line in cfg_file:
            if line[0]!="#" and line.strip():
                this_vals = line.strip().split("=")
                ## add this variable to the dict, removing quotes around strings
                config_vars[this_vals[0]] = this_vals[1].strip('"\'')

    ## TODO: check that all the appropriate variables are defined
                
    try:
        os.makedirs(pargs.outdir+"/logs/")        
    except OSError:
        pass
    
    if len(sys.argv)==5:
        print "Running analysis without summary stats"
        ## first run the annotation script
        print "Running LD expansion and annotation"
        start_time = time.time()
        subprocess.call(["python", "./src/expand_and_annotate_snps.py", "--loglevel",
                        "full", "--kg_pop", config_vars["KG_POP"], "--ld_threshold",
                        config_vars["LD_THRESH"], "--ld_check_area", config_vars["LD_AREA"],
                        "--gene_bed_file", config_vars["GENE_BED_FILE"], "--kgxref_file",
                        config_vars["KGXREF_FILE"], "--unstranded_partition_dir",
                        config_vars["UNSTRANDED_PARTITION_DIR"], "--fantom5_dir",
                        config_vars["FANTOM5_DIR"], "--enhancer_locus_window",
                        config_vars["ENH_LOCUS_WINDOW"], "--fantom5_correlation_file",
                        config_vars["F5_CORR_FILE"], "--gtex_dir", config_vars["GTEX_DIR"],
                        "--factorbook_file", config_vars["FACTORBOOK_FILE"],
                        "--roadmap_chromhmm_dir", config_vars["ROADMAP_HMM_DIR"],
                        "--homer_motif_bed_file", config_vars["HOMER_BED_FILE"],
                        "--homer_motif_pwm_file", config_vars["HOMER_PWM_FILE"],
                        "--homer_motif_seq_file", config_vars["HOMER_SEQ_FILE"],
                        "--dashr_locus_file", config_vars["DASHR_LOCUS_F"],
                        "--targetscan_dir", config_vars["TARGETSCAN_DIR"], "--bedtools_bin_dir",
                        config_vars["BEDTOOLS_BIN_DIR"], config_vars["KG_DIR"],
                        pargs.top_snpf, pargs.outdir, pargs.outprefix])                    
        print "LD expansion and annotation took %.2f seconds" % (time.time()-start_time)
        
        print "Summarizing annotation results"
        subprocess.call(["./analysis_scripts/count_annotation_overlaps.sh",  "-l", config_vars["ENH_LOCUS_WINDOW"], "-f", config_vars["F5_CLASSES"], "-g", config_vars["GTEX_CLASSES"], "-r", config_vars["ROADMAP_CLASSES"], pargs.outdir, config_vars["LD_THRESH"], config_vars["LD_AREA"], pargs.outprefix, pargs.outdir+"/summaries/"])

        ## TODO: figure out how to make R loading generic
        ## find the parameter file from this most recent run
        param_file = subprocess.Popen('ls -t '+pargs.outdir+'/parameters/*parameters* | head -1', shell=True, stdout=subprocess.PIPE).stdout.read().strip()

        ## run the R analysis script
        print "Running R analysis script"
        start_time = time.time()
        subprocess.call("module load R/3.2.3; time Rscript ./analysis_scripts/Rscript_run_full_analysis.R ./analysis_scripts "+param_file+" "+pargs.outprefix+" TRUE", shell=True)
        print "R analysis took %.2f seconds" % (time.time()-start_time)
        
        ## finally submit a job to run bootstrapping
        print "Submitting job for enhancer bootstrapping analysis"
        subprocess.call(["bsub", "-M", "40000", "-J", pargs.outprefix+"_enh_bootstrapping", "-o",
                         pargs.outdir+"/logs/"+pargs.outprefix+"_enh_bootstrapping.o%J", "-e",
                         pargs.outdir+"/logs/"+pargs.outprefix+"_enh_bootstrapping.e%J",
                         "./bsub_wrappers/enhancer_bootstrap_bsub_wrapper.sh",
                         "./src/enhancer_only_sample_and_expand_matched_input_variants.R",
                         config_vars["NUM_SAMPLES"], config_vars["MAF_BIN_SIZE"],
                         config_vars["DIST_ROUND"], config_vars["DIST_THRESHOLD"],
                         config_vars["LD_PARTNER_THRESHOLD"], config_vars["BG_SNP_INFOF"],
                         config_vars["LD_SETS_DIR"], config_vars["REF_SUMMARY_DIR"], param_file])
    ## only run summary stats if there is a summary file
    elif pargs.summary_file:
        ## to do summary stats analysis, we need to check arguments at each step        
        ## first do p-value expansion, or not
        if pargs.skip_pval_expansion:
            print "Skipping p-value expansion"
            ## just run the annnotation script directly
            start_time = time.time()
            subprocess.call(["python", "./src/expand_and_annotate_snps.py", "--loglevel",
                            "full", "--kg_pop", config_vars["KG_POP"], "--ld_threshold",
                            config_vars["LD_THRESH"], "--ld_check_area", config_vars["LD_AREA"],
                            "--gene_bed_file", config_vars["GENE_BED_FILE"], "--kgxref_file",
                            config_vars["KGXREF_FILE"], "--unstranded_partition_dir",
                            config_vars["UNSTRANDED_PARTITION_DIR"], "--fantom5_dir",
                            config_vars["FANTOM5_DIR"], "--enhancer_locus_window",
                            config_vars["ENH_LOCUS_WINDOW"], "--fantom5_correlation_file",
                            config_vars["F5_CORR_FILE"], "--gtex_dir", config_vars["GTEX_DIR"],
                            "--factorbook_file", config_vars["FACTORBOOK_FILE"],
                            "--roadmap_chromhmm_dir", config_vars["ROADMAP_HMM_DIR"],
                            "--homer_motif_bed_file", config_vars["HOMER_BED_FILE"],
                            "--homer_motif_pwm_file", config_vars["HOMER_PWM_FILE"],
                            "--homer_motif_seq_file", config_vars["HOMER_SEQ_FILE"],
                            "--dashr_locus_file", config_vars["DASHR_LOCUS_F"],
                            "--targetscan_dir", config_vars["TARGETSCAN_DIR"],
                            "--bedtools_bin_dir", config_vars["BEDTOOLS_BIN_DIR"],
                            config_vars["KG_DIR"],pargs.top_snpf, pargs.outdir, pargs.outprefix])
            print "LD expansion and annotation took %.2f seconds" % (time.time()-start_time)
        else:
            ## do p-value expansion
            if pargs.sig_mult and pargs.rsid_column and pargs.pos_column and pargs.pval_column and pargs.chr_column:
                ## if we have beta and MAF columns, only allow consistent effects
                if pargs.beta_column and pargs.maf_column:
                    print "Performing p-value expansion with consistent effect directions"
                    start_time = time.time()
                    subprocess.call(["python", "./data_preprocessing/pval_expand_tagsnp_set.py", "--single_file", "--sig_multiplier", str(pargs.sig_mult), "--dist_threshold", config_vars["LD_AREA"], "--rsid_col", str(pargs.rsid_column), "--pos_col", str(pargs.pos_column), "--pval_col", str(pargs.pval_column), "--chr_col", str(pargs.chr_column), "--beta_col", str(pargs.beta_column), "--maf_col", str(pargs.maf_column), pargs.top_snpf, pargs.summary_file, pargs.outdir+"/"+pargs.outprefix+"_pval_expanded_snps.txt"])
                    print "P-value expansion took %.2f seconds" % (time.time()-start_time)
                else:
                    print "Performing p-value expansion with no consideration of effect direction"
                    start_time = time.time()
                    subprocess.call(["python", "./data_preprocessing/pval_expand_tagsnp_set.py", "--single_file", "--sig_multiplier", str(pargs.sig_mult), "--dist_threshold", config_vars["LD_AREA"], "--rsid_col", str(pargs.rsid_column), "--pos_col", str(pargs.pos_column), "--pval_col", str(pargs.pval_column), "--chr_col", str(pargs.chr_column), pargs.top_snpf, pargs.summary_file, pargs.outdir+"/"+pargs.outprefix+"_pval_expanded_snps.txt"])
                    print "P-value expansion took %.2f seconds" % (time.time()-start_time)
                
                ## do LD pruning
                print "Performing LD pruning"
                start_time = time.time()
                subprocess.call(["python", "./data_preprocessing/ld_prune_snp_set.py", config_vars["LD_THRESH"], pargs.outdir+"/"+pargs.outprefix+"_pval_expanded_snps.txt", config_vars["KG_DIR"]+"/"+config_vars["KG_POP"]+"/", pargs.outdir+"/"+pargs.outprefix+"_pruning/"])
                print "LD pruning took %.2f seconds" % (time.time()-start_time)
                
                ## now annotate these
                print "Running pruned variants through re-expansion and annotation"
                start_time = time.time()
                subprocess.call(["python", "./src/expand_and_annotate_snps.py", "--loglevel",
                                "full", "--kg_pop", config_vars["KG_POP"], "--ld_threshold",
                                config_vars["LD_THRESH"], "--ld_check_area", config_vars["LD_AREA"],
                                "--gene_bed_file", config_vars["GENE_BED_FILE"], "--kgxref_file",
                                config_vars["KGXREF_FILE"], "--unstranded_partition_dir",
                                config_vars["UNSTRANDED_PARTITION_DIR"], "--fantom5_dir",
                                config_vars["FANTOM5_DIR"], "--enhancer_locus_window",
                                config_vars["ENH_LOCUS_WINDOW"], "--fantom5_correlation_file",
                                config_vars["F5_CORR_FILE"], "--gtex_dir", config_vars["GTEX_DIR"],
                                "--factorbook_file", config_vars["FACTORBOOK_FILE"],
                                "--roadmap_chromhmm_dir", config_vars["ROADMAP_HMM_DIR"],
                                "--homer_motif_bed_file", config_vars["HOMER_BED_FILE"],
                                "--homer_motif_pwm_file", config_vars["HOMER_PWM_FILE"],
                                "--homer_motif_seq_file", config_vars["HOMER_SEQ_FILE"],
                                "--dashr_locus_file", config_vars["DASHR_LOCUS_F"],
                                "--targetscan_dir", config_vars["TARGETSCAN_DIR"],
                                "--bedtools_bin_dir", config_vars["BEDTOOLS_BIN_DIR"],
                                config_vars["KG_DIR"], pargs.outdir+"/"+pargs.outprefix+"_pruning/pruned_set_pipeline_input.txt", pargs.outdir, pargs.outprefix])
                print "LD expansion and annotation took %.2f seconds" % (time.time()-start_time)
            else:
                sys.exit("For p-value expansion, p-value multiplier (--sig_mult) and column numbers for rsID, position, pvalue, and chromosome must be provided")
                
        ## everything after this is independent of whether we did p-value expansion or not
        print "Summarizing annotation results"    
        subprocess.call(["./analysis_scripts/count_annotation_overlaps.sh",  "-l", config_vars["ENH_LOCUS_WINDOW"], "-f", config_vars["F5_CLASSES"], "-g", config_vars["GTEX_CLASSES"], "-r", config_vars["ROADMAP_CLASSES"], pargs.outdir, config_vars["LD_THRESH"], config_vars["LD_AREA"], pargs.outprefix, pargs.outdir+"/summaries/"])

        ## find the parameter file from this most recent run
        param_file = subprocess.Popen('ls -t '+pargs.outdir+'/parameters/*parameters* | head -1', shell=True, stdout=subprocess.PIPE).stdout.read().strip()

        ## run the R analysis script
        print "Running R analysis script"
        start_time = time.time()
        subprocess.call("module load R/3.2.3; time Rscript ./analysis_scripts/Rscript_run_full_analysis.R ./analysis_scripts "+param_file+" "+pargs.outprefix+" TRUE", shell=True)
        print "R analysis took %.2f seconds" % (time.time()-start_time)
        
        ## submit a job to run bootstrapping
        print "Submitting job for enhancer bootstrapping analysis"
        subprocess.call(["bsub", "-M", "40000", "-J", pargs.outprefix+"_enh_bootstrapping", "-o",
                         pargs.outdir+"/logs/"+pargs.outprefix+"_enh_bootstrapping.o%J", "-e",
                         pargs.outdir+"/logs/"+pargs.outprefix+"_enh_bootstrapping.e%J",
                         "./bsub_wrappers/enhancer_bootstrap_bsub_wrapper.sh",
                         "./src/enhancer_only_sample_and_expand_matched_input_variants.R",
                         config_vars["NUM_SAMPLES"], config_vars["MAF_BIN_SIZE"],
                         config_vars["DIST_ROUND"], config_vars["DIST_THRESHOLD"],
                         config_vars["LD_PARTNER_THRESHOLD"], config_vars["BG_SNP_INFOF"],
                         config_vars["LD_SETS_DIR"], config_vars["REF_SUMMARY_DIR"], param_file])
        
        ## submit a job for colocalization analysis
        ## first check that all the required arguments are present
        if pargs.rsid_column and pargs.pos_column and pargs.pval_column and pargs.chr_column and pargs.allele1_column and pargs.allele2_column and pargs.maf_column and pargs.case_prop and pargs.sample_size:
            print "Submitting co-localization analysis job"
            ## need to define the correct input file to use
            if pargs.skip_pval_expansion:
                coloc_input_f = pargs.top_snpf
            else:
                coloc_input_f = pargs.outdir+"/"+pargs.outprefix+"_pruning/pruned_set_pipeline_input.txt"
            subprocess.call(["bsub", "-M", "40000", "-J", pargs.outprefix+"_gtex_colocalization",
                             "-o", pargs.outdir+"/logs/"+pargs.outprefix+"_gtex_coloc.o%J",
                             "-e", pargs.outdir+"/logs/"+pargs.outprefix+"_gtex_coloc.e%J",
                             "./bsub_wrappers/gtex_coloc_bsub_wrapper.sh",
                             "./src/gtex_gwas_colocalization_analysis.R",
                             pargs.outdir+"/gtex_gwas_colocalization_analysis/",
                             param_file, config_vars["COLOC_H4_THRESH"],
                             config_vars["COLOC_ABF_THRESH"], coloc_input_f, 
                             pargs.summary_file, config_vars["COLOC_GTEX_DIR"],
                             config_vars["GTEX_SAMPLE_SIZEF"], config_vars["GTEX_CLASSES"],
                             config_vars["GTEX_RSID_MATCH"], config_vars["HG19_ENSEMBL_REF_FILE"],
                             config_vars["RELEVANT_CLASSES"], str(pargs.rsid_column), str(pargs.pos_column),
                             str(pargs.pval_column), str(pargs.chr_column), str(pargs.allele1_column),
                             str(pargs.allele2_column), str(pargs.maf_column),
                             str(pargs.case_prop), str(pargs.sample_size)])
            
            print "Submitting lncRNA correlation analysis job"
            subprocess.call(["bsub", "-M", "40000", "-J", pargs.outprefix+"_lncRNA_correlation",
                             "-o", pargs.outdir+"/logs/"+pargs.outprefix+"_gtex_lncRNA_corr.o%J",
                             "-e", pargs.outdir+"/logs/"+pargs.outprefix+"_gtex_lncRNA_corr.e%J",
                             "-w", "done("+pargs.outprefix+"_gtex_colocalization)",
                             "./bsub_wrappers/gtex_lncRNA_corr_bsub_wrapper.sh",
                             "./src/lncRNA_gtex_correlation.R",
                             pargs.outdir+"/gtex_lncRNA_correlation_analysis/",
                             pargs.outdir+"/gtex_gwas_colocalization_analysis/tables/"+pargs.outprefix+"_gtex_coloc_summaries.txt",
                             config_vars["GTEX_EXPR_DIR"], config_vars["SAMPLE_INFO_FILE"],
                             config_vars["GENCODE_LNCRNA_FILE"], config_vars["F5_CLASSES"],
                             config_vars["GTEX_CLASSES"], config_vars["ROADMAP_CLASSES"],
                             config_vars["COLOC_H4_THRESH"], config_vars["COR_THRESH"]])
                        
        else:
            print "Can't do colocalization without all the required columns"    
    else:
        print "Cannot perform p-value expansion or colocalization without summary statistics file."
