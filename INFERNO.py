"""
INFERNO.py
Alex Amlie-Wolf, 07/10/17

This is the main driver method for running the INFERNO pipeline. It can be run directly in this
python script or submitted through bsub/qsub using the INFERNO.sh wrapper script
"""

import argparse, subprocess, datetime, os, time, sys

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Driver script for the INFERNO pipeline")
    parser.add_argument("--skip_ld_expansion", action='store_true', help="Give this flag to skip any LD (and p-value-based) expansion and just run analysis directly on the provided list of input variants.")
    ## arguments for summary statistics analysis
    parser.add_argument("--rsid_column", type=int, help="The summary statistics column number containing the rsIDs")
    parser.add_argument("--pos_column", type=int, help="The summary statistics column number containing the positions")
    parser.add_argument("--pval_column", type=int, help="The summary statistics column number containing the p-values")
    parser.add_argument("--chr_column", type=int, help="The summary statistics column number containing the chromosomes")
    parser.add_argument("--allele1_column", type=int, help="The summary statistics column number containing allele 1, which should correspond to the major allele.")
    parser.add_argument("--allele2_column", type=int, help="The summary statistics column number containing allele 2, which should correspond to the minor allele.")
    parser.add_argument("--maf_column", type=int, help="The summary statistics column number containing the minor allele frequency. Note that if this is ever >50% and a beta column is provided, the effect direction will be flipped to be defined relative to the minor allele.")
    parser.add_argument("--beta_column", type=int, help="The summary statistics column number containing the beta estimate (used for p-value expansion with consistent directions)")        
    ## for p-value expansion and pruning
    parser.add_argument("--run_pval_expansion", action='store_true', help="If you want to do expansion by p-values when you have summary statistics, provide this flag. Otherwise, the top SNP file will be directly expanded.")
    parser.add_argument("--sig_mult", type=float, help="The multiplier range for significance of the p-value expanded variant set (e.g. a value of 10 means one order of magnitude)")
    parser.add_argument("--case_prop", type=float, help="The proportion of cases in the GWAS, for colocalization.")
    parser.add_argument("--sample_size", type=int, help="The total number of samples in the GWAS, for colocalization.")
    parser.add_argument("--summary_file", help="The path to the full summary statistics file, required for p-value expansion, colocalization analysis, and lncRNA target analysis.")
    ## flags for running the slow, bsub'd steps after expansion and annotation
    parser.add_argument("--run_enhancer_sampling", action='store_true', help="If you want to run the enhancer bootstrapping analysis, provide this flag.")
    parser.add_argument("--run_gtex_coloc", action="store_true", help="If you want to run COLOC analysis of your summary statistics against GTEx eQTLs from 44 tissues (requires summary statistics)")
    parser.add_argument("--run_lncrna_correlation", action="store_true", help="If you want to analyze expression correlation of any lncRNAs identified by COLOC analysis (--run_gtex_coloc flag) against all other RNAseq-based GTEx genes to find lncRNA targets.")
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
                config_vars[this_vals[0].upper()] = this_vals[1].strip('"\'')

    if not pargs.skip_ld_expansion:
        if "LD_THRESH" not in config_vars or "LD_AREA" not in config_vars:
            sys.exit("Must provide R^2 and distance thresholds for LD expansion")
            
    ## for expand_and_annotate_snps.py script, set up the flag arguments based on what's
    ## defined in the config file. to do this, make a list of the possible flags to the
    ## script. note that the config file variables must exactly match these (uppercase or not)
    annot_arg_options = ["gene_bed_file", "kgxref_file", "partition_dir", "unstranded_partition_dir", "fantom5_dir", "enhancer_midpoint_window", "enhancer_locus_window", "fantom5_correlation_file", "gtex_dir", "factorbook_file", "roadmap_chromhmm_dir", "homer_motif_bed_file", "homer_motif_pwm_file", "homer_motif_seq_file", "dashr_locus_file", "targetscan_dir", "bedtools_bin_dir", "plink_path"]
    
    ## TODO: deal with --skip_closest_enh, --skip_enh_summary
    annotation_arg_list = []
    for cfg_var in config_vars:
        if cfg_var.lower() in annot_arg_options:
            annotation_arg_list = annotation_arg_list + ["--"+cfg_var.lower(), config_vars[cfg_var]]
       
    try:
        os.makedirs(pargs.outdir+"/logs/")        
    except OSError:
        pass
    
    if not pargs.summary_file:
        print "Running analysis without summary stats"
        ## first run the annotation script
        start_time = time.time()
        if not pargs.skip_ld_expansion:
            print "Running LD expansion and annotation"        
            subprocess.call(["python", "./src/expand_and_annotate_snps.py", "--loglevel",
                            "full", "--kg_pop", config_vars["KG_POP"], "--ld_threshold",
                            config_vars["LD_THRESH"], "--ld_check_area", config_vars["LD_AREA"]] +
                            annotation_arg_list +
                            [config_vars["KG_DIR"], pargs.top_snpf, pargs.outdir, pargs.outprefix])
            print "LD expansion and annotation took %.2f seconds" % (time.time()-start_time)
        else:
            print "Running direct annotation (no LD expansion)"
            subprocess.call(["python", "./src/expand_and_annotate_snps.py", "--loglevel",
                            "full", "--skip_ld_expansion"] + 
                            annotation_arg_list +
                            [config_vars["KG_DIR"], pargs.top_snpf, pargs.outdir, pargs.outprefix])
            print "Direct annotation took %.2f seconds" % (time.time()-start_time)            
        
        print "Summarizing annotation results"
        ## TODO: check that the class files exist
        subprocess.call(["./analysis_scripts/count_annotation_overlaps.sh",  "-l", config_vars["ENHANCER_LOCUS_WINDOW"], "-f", config_vars["F5_CLASSES"], "-g", config_vars["GTEX_CLASSES"], "-r", config_vars["ROADMAP_CLASSES"], pargs.outdir, config_vars["LD_THRESH"], config_vars["LD_AREA"], pargs.outprefix, pargs.outdir+"/summaries/"])

        ## find the parameter file from this most recent run
        param_file = subprocess.Popen('ls -t '+pargs.outdir+'/parameters/*parameters* | head -1', shell=True, stdout=subprocess.PIPE).stdout.read().strip()

        ## run the R analysis script
        print "Running R analysis script"
        start_time = time.time()
        subprocess.call(["Rscript", "./analysis_scripts/Rscript_run_full_analysis.R", "./analysis_scripts", param_file, pargs.outprefix, "TRUE", config_vars["F5_CLASSES"], config_vars["GTEX_CLASSES"], config_vars["ROADMAP_CLASSES"]])
        print "R analysis took %.2f seconds" % (time.time()-start_time)
        
        ## finally submit a job to run bootstrapping, if we want to
        if pargs.run_enhancer_sampling:
            ## TODO: check that the config file has the right variables
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
    ## only run summary stats-based analyses if there is a summary file
    else:
        print "Running analysis using summary statistics file."
        ## to do summary stats analysis, we need to check arguments at each step
        ## first do p-value expansion, or not
        if not pargs.run_pval_expansion:
            if not pargs.skip_ld_expansion:
                print "Directly LD-expanding and annotating top variant file"
                ## just run the annnotation script directly
                start_time = time.time()
                subprocess.call(["python", "./src/expand_and_annotate_snps.py", "--loglevel",
                                "full", "--kg_pop", config_vars["KG_POP"], "--ld_threshold",
                                config_vars["LD_THRESH"], "--ld_check_area", config_vars["LD_AREA"]] +
                                annotation_arg_list + 
                                [config_vars["KG_DIR"],pargs.top_snpf, pargs.outdir, pargs.outprefix])
                print "LD expansion and annotation took %.2f seconds" % (time.time()-start_time)
            else:
                print "Directly LD-expanding and annotating top variant file"
                ## just run the annnotation script directly
                start_time = time.time()
                subprocess.call(["python", "./src/expand_and_annotate_snps.py", "--loglevel",
                                "full", "--skip_ld_expansion"] +
                                annotation_arg_list + 
                                [config_vars["KG_DIR"],pargs.top_snpf, pargs.outdir, pargs.outprefix])
                print "Direct annotation took %.2f seconds" % (time.time()-start_time)
        else:
            ## do p-value expansion.  note that i ignore the skip_ld_expansion flag; you can't
            ## skip the re-expansion for now..
            if pargs.sig_mult and pargs.rsid_column and pargs.pos_column and pargs.pval_column and pargs.chr_column:
                ## if we have beta and MAF columns, only allow consistent effects
                if pargs.beta_column and pargs.maf_column:
                    print "Performing p-value expansion with consistent effect directions"
                    start_time = time.time()
                    subprocess.call(["python", "./data_preprocessing/pval_expand_tagsnp_set.py", "--sig_multiplier", str(pargs.sig_mult), "--dist_threshold", config_vars["LD_AREA"], "--rsid_col", str(pargs.rsid_column), "--pos_col", str(pargs.pos_column), "--pval_col", str(pargs.pval_column), "--chr_col", str(pargs.chr_column), "--beta_col", str(pargs.beta_column), "--maf_col", str(pargs.maf_column), pargs.top_snpf, pargs.summary_file, pargs.outdir+"/"+pargs.outprefix+"_pval_expanded_snps.txt"])
                    print "P-value expansion took %.2f seconds" % (time.time()-start_time)
                else:
                    print "Performing p-value expansion with no consideration of effect direction"
                    start_time = time.time()
                    subprocess.call(["python", "./data_preprocessing/pval_expand_tagsnp_set.py", "--sig_multiplier", str(pargs.sig_mult), "--dist_threshold", config_vars["LD_AREA"], "--rsid_col", str(pargs.rsid_column), "--pos_col", str(pargs.pos_column), "--pval_col", str(pargs.pval_column), "--chr_col", str(pargs.chr_column), pargs.top_snpf, pargs.summary_file, pargs.outdir+"/"+pargs.outprefix+"_pval_expanded_snps.txt"])
                    print "P-value expansion took %.2f seconds" % (time.time()-start_time)
                
                ## do LD pruning
                print "Performing LD pruning of p-value expanded sets"
                start_time = time.time()
                subprocess.call(["python", "./data_preprocessing/ld_prune_snp_set.py", config_vars["LD_THRESH"], pargs.outdir+"/"+pargs.outprefix+"_pval_expanded_snps.txt", config_vars["KG_DIR"]+"/"+config_vars["KG_POP"]+"/", pargs.outdir+"/"+pargs.outprefix+"_pruning/"])
                print "LD pruning took %.2f seconds" % (time.time()-start_time)
                
                ## now annotate these
                print "Running pruned variants through re-expansion and annotation"
                start_time = time.time()
                subprocess.call(["python", "./src/expand_and_annotate_snps.py", "--loglevel",
                                "full", "--kg_pop", config_vars["KG_POP"], "--ld_threshold",
                                config_vars["LD_THRESH"], "--ld_check_area", config_vars["LD_AREA"]] +
                                annotation_arg_list + 
                                [config_vars["KG_DIR"], pargs.outdir+"/"+pargs.outprefix+"_pruning/pruned_set_pipeline_input.txt", pargs.outdir, pargs.outprefix])
                print "LD expansion and annotation took %.2f seconds" % (time.time()-start_time)
            else:
                sys.exit("For p-value expansion, p-value multiplier (--sig_mult) and column numbers for rsID, position, pvalue, and chromosome must be provided")
                
        ## everything after this is independent of whether we did p-value expansion or not
        print "Summarizing annotation results"
        ## TODO: check that the class files exist
        subprocess.call(["./analysis_scripts/count_annotation_overlaps.sh",  "-l", config_vars["ENHANCER_LOCUS_WINDOW"], "-f", config_vars["F5_CLASSES"], "-g", config_vars["GTEX_CLASSES"], "-r", config_vars["ROADMAP_CLASSES"], pargs.outdir, config_vars["LD_THRESH"], config_vars["LD_AREA"], pargs.outprefix, pargs.outdir+"/summaries/"])

        ## find the parameter file from this most recent run
        param_file = subprocess.Popen('ls -t '+pargs.outdir+'/parameters/*parameters* | head -1', shell=True, stdout=subprocess.PIPE).stdout.read().strip()

        ## run the R analysis script
        print "Running R analysis script"
        start_time = time.time()
        ## TODO: add class file arguments
        subprocess.call(["Rscript", "./analysis_scripts/Rscript_run_full_analysis.R", "./analysis_scripts", param_file, pargs.outprefix, "TRUE"])
        print "R analysis took %.2f seconds" % (time.time()-start_time)
        
        ## submit a job to run bootstrapping
        if pargs.run_enhancer_sampling:
            ## TODO: check config
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
        if pargs.run_gtex_coloc and pargs.rsid_column and pargs.pos_column and pargs.pval_column and pargs.chr_column and pargs.allele1_column and pargs.allele2_column and pargs.maf_column and pargs.case_prop and pargs.sample_size:
            ## TODO: check config
            print "Submitting co-localization analysis job"
            ## need to define the correct input file to use
            if pargs.run_pval_expansion:
                coloc_input_f = pargs.outdir+"/"+pargs.outprefix+"_pruning/pruned_set_pipeline_input.txt"
            else:
                coloc_input_f = pargs.top_snpf

            if "LOCUSZOOM_PATH" in config_vars:
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
                                 str(pargs.case_prop), str(pargs.sample_size),
                                 config_vars["LOCUSZOOM_PATH"]])
            else:
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

            if pargs.run_lncrna_correlation:
                ## TODO: check config
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
