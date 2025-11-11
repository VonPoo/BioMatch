# VCF_Geno_nonGeno.R - Backend processing script for Sample Comparison Analysis (DEBUG VERSION)
#
# Purpose
# - Process VCF inputs: convert to PLINK via PLINK2 with multiallelic filtering,
#   then perform downstream analysis via PLINK1 and deepKin.
# - Supports single or paired datasets (genomic/non-genomic) and robust validation.
#
# Inputs
# - `data_dir` containing `genomic.vcf.gz` and/or `non_genomic.vcf.gz` and `params.RData`.
# - `params.RData` must include a `res` list with keys `chr_set`, `file_type`, `dir`.
#
# Outputs
# - Converted PLINK files, deepKin summaries, relatedness outputs, plots, and CSV.
#
# Usage
#   Rscript VCF_Geno_nonGeno.R /path/to/data_dir
library(data.table)
suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(e1071)
  library(pls)
  library(deepKin)
})

# Enhanced logging function
debug_log <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
}

# Safe file reading function
safe_read_table <- function(file_path, ...) {
  debug_log(paste("Attempting to read file:", file_path))
  
  if (!file.exists(file_path)) {
    debug_log(paste("File does not exist:", file_path), "ERROR")
    return(NULL)
  }
  
  if (file.size(file_path) == 0) {
    debug_log(paste("File is empty:", file_path), "ERROR")
    return(NULL)
  }
  
  tryCatch({
    data <- read.table(file_path, ...)
    debug_log(paste("Successfully read file with", nrow(data), "rows and", ncol(data), "columns"))
    return(data)
  }, error = function(e) {
    debug_log(paste("Error reading file:", e$message), "ERROR")
    return(NULL)
  })
}

# Safe file counting function
safe_count_lines <- function(file_path) {
  if (!file.exists(file_path)) {
    debug_log(paste("File does not exist for counting:", file_path), "ERROR")
    return(0)
  }
  
  tryCatch({
    count <- length(readLines(file_path))
    debug_log(paste("File", basename(file_path), "has", count, "lines"))
    return(count)
  }, error = function(e) {
    debug_log(paste("Error counting lines in file:", e$message), "ERROR")
    return(0)
  })
}

# Override the calculate.me function to work with both PLINK1 and PLINK2
modified_calculate.me <- function(bfileprefix, method, plink1_path = NULL, plink2_path = NULL, gear_path = NULL, freq_path = NULL, pop_size = NULL, chr_set = NULL){
  debug_log("Starting modified_calculate.me function")
  debug_log(paste("Parameters: bfileprefix =", bfileprefix, ", method =", method))
  
  switch(method,
         "GRM" = {
           if(is.null(plink1_path) && is.null(plink2_path)) { 
             stop("Error in calculate.me(): No plink_path was specified!") 
           }
           
           # Check number of individuals in the .fam file
           fam_file <- paste0(bfileprefix, ".fam")
           debug_log(paste("Checking fam file:", fam_file))
           
           if(!file.exists(fam_file)) { 
             debug_log(paste("Fam file not found:", fam_file), "ERROR")
             stop("Error in calculate.me(): .fam file not found!") 
           }
           
           # Count individuals safely
           fam_data <- safe_read_table(fam_file)
           if (is.null(fam_data)) {
             stop("Error reading .fam file")
           }
           
           n_individuals <- nrow(fam_data)
           debug_log(paste("Number of individuals in merged dataset:", n_individuals))
           
           # For relatedness matrix, use PLINK2 if available, otherwise PLINK1
           if (!is.null(plink2_path) && file.exists(plink2_path)) {
             debug_log("Using PLINK2 for relatedness matrix calculation")
             
             # Conditionally add --bad-freqs if fewer than 50 individuals
             bad_freqs_flag <- ifelse(n_individuals < 50, "--bad-freqs", "")
             if(bad_freqs_flag != "") {
               debug_log("Using --bad-freqs flag due to small sample size (<50)")
             }
             
             # Add chr_set parameter to the command
             chr_set_flag <- ifelse(!is.null(chr_set), paste0("--chr-set ", chr_set), "")
             
             rel_file <- paste0(bfileprefix, ".rel.bin")
             if(!file.exists(rel_file)){
               debug_log("Creating relatedness matrix with PLINK2")
               cmd <- paste0(plink2_path, " --silent --bfile ", bfileprefix, 
                             " --make-rel triangle bin4 ", bad_freqs_flag, " ", chr_set_flag, " --out ", bfileprefix)
               debug_log(paste("Running command:", cmd))
               system(cmd)
               
               if (!file.exists(rel_file)) {
                 debug_log("Failed to create .rel.bin file", "ERROR")
                 stop("Failed to create relatedness matrix")
               }
             } else {
               debug_log("Relatedness matrix file already exists")
             }
             
             # If pop_size isn't provided, use the count from the fam file
             if(is.null(pop_size)) {
               pop_size <- n_individuals
             }
             
             debug_log(paste("Reading relatedness matrix for", pop_size, "individuals"))
             expected_size <- pop_size * (pop_size + 1) / 2
             debug_log(paste("Expected matrix size:", expected_size))
             
             grm = readBin(rel_file, what="numeric", n=expected_size, size=4)
             debug_log(paste("Successfully read", length(grm), "values from relatedness matrix"))
             
             if (length(grm) < expected_size) {
               debug_log(paste("Warning: Expected", expected_size, "values but got", length(grm)), "WARNING")
             }
             
             me = 1/var(grm[-cumsum(1:pop_size)])
             debug_log(paste("Calculated effective markers (PLINK2):", me))
             
           } else if (!is.null(plink1_path)) {
             debug_log("Using PLINK1 for relatedness matrix calculation")
             
             # Add chr_set parameter to the command
             chr_set_flag <- ifelse(!is.null(chr_set), paste0("--chr-set ", chr_set), "")
             
             genome_file <- paste0(bfileprefix, ".genome")
             # PLINK1 uses --genome instead of --make-rel
             if(!file.exists(genome_file)){
               debug_log("Creating genome file with PLINK1")
               cmd <- paste0(plink1_path, " --bfile ", bfileprefix, 
                             " --genome ", chr_set_flag, " --out ", bfileprefix)
               debug_log(paste("Running command:", cmd))
               system(cmd)
               
               if (!file.exists(genome_file)) {
                 debug_log("Failed to create .genome file", "ERROR")
                 stop("Failed to create genome file")
               }
             } else {
               debug_log("Genome file already exists")
             }
             
             # Read genome file and calculate me differently for PLINK1
             genome_data <- safe_read_table(genome_file, header = TRUE)
             if (is.null(genome_data)) {
               stop("Failed to read genome file")
             }
             
             debug_log(paste("Genome file has", nrow(genome_data), "rows"))
             
             if(nrow(genome_data) > 0 && "PI_HAT" %in% names(genome_data)) {
               me = 1/var(genome_data$PI_HAT)
               debug_log(paste("Calculated effective markers (PLINK1):", me))
             } else {
               debug_log("No PI_HAT data found, using default value", "WARNING")
               me = 1000  # Default value if no relationships detected
             }
           }
         },
         "RDM"  = {
           debug_log("Using RDM method")
           if(is.null(gear_path)) { stop("Error in calculate.me(): No gear_path was specified!") }
           if(!file.exists(paste0(bfileprefix,".it.me"))){
             system(paste0(gear_path, " --me --bfile ",bfileprefix," --iter 100 --out ",bfileprefix))
           }
           itme = read.table(paste0(bfileprefix,".it.me"), header = T)
           me = itme[nrow(itme),"Me"]
         }
  )
  
  debug_log(paste("Final effective markers value:", me))
  return(me)
}

# Override the original function
unlockBinding("calculate.me", getNamespace("deepKin"))
assignInNamespace("calculate.me", modified_calculate.me, "deepKin")


# Function to log errors
log_error <- function(data_dir, error_msg) {
  log_file <- file.path(data_dir, "error.log")
  debug_log(paste("Logging error to:", log_file))
  write(paste0(
    "Error occurred at ", Sys.time(), ":\n",
    error_msg, "\n",
    "--------------------\n"
  ), log_file, append = TRUE)
}

# Email functionality removed per BioMatch requirements

# Function to run system commands
run_command <- function(cmd) {
  debug_log(paste("Running system command:", cmd))
  result <- system(cmd, intern = TRUE)
  debug_log(paste("Command output:", paste(result, collapse = "; ")))
  return(result)
}

# Function to run PLINK commands with version awareness
run_plink <- function(cmd, plink_path, chr_set = NULL) {
  # Ensure --chr-set is present if provided
  if (!is.null(chr_set) && !grepl("--chr-set", cmd)) {
    cmd <- paste(cmd, "--chr-set", chr_set)
  }
  full_cmd <- paste(plink_path, cmd)
  debug_log(paste("Running PLINK command:", full_cmd))
  # 同样使用双引号包裹，避免内部参数中的单引号冲突
  wrapped <- sprintf('bash -lc "exec -a biomatch %s"', full_cmd)
  result <- system(wrapped, intern = TRUE)
  debug_log(paste("PLINK output:", paste(result, collapse = "; ")))
  return(result)
}

# Process VCF files if input type is VCF (混合策略：PLINK2转换，PLINK1分析)
process_vcf_files <- function(data_dir, chr_set, bcftools_path, plink2_path, plink1_path) {
  debug_log("Starting VCF file processing")
  
  # First, check if input files exist
  genomic_vcf <- file.path(data_dir, "genomic.vcf.gz")
  non_genomic_vcf <- file.path(data_dir, "non_genomic.vcf.gz")
  
  debug_log(paste("Checking for genomic VCF:", genomic_vcf))
  debug_log(paste("Checking for non-genomic VCF:", non_genomic_vcf))

  # 允许单数据集：若两者都不存在才报错
  genomic_exists <- file.exists(genomic_vcf)
  non_genomic_exists <- file.exists(non_genomic_vcf)
  if (!genomic_exists && !non_genomic_exists) {
    stop("Input VCF files not found in the data directory")
  }
  
  # Step 1: Convert to BGZF format
  debug_log("Converting to BGZF format...")
  
  # Function to convert to BGZF
  convert_to_bgzf <- function(input_vcf) {
    debug_log(paste("Converting to BGZF:", input_vcf))
    # Decompress first
    temp_vcf <- sub("\\.gz$", "", input_vcf)
    cmd <- paste("gunzip -c", input_vcf, ">", temp_vcf)
    run_command(cmd)
    
    # Compress with bgzip
    cmd <- paste("bgzip -c", temp_vcf, ">", input_vcf)
    run_command(cmd)
    
    # Clean up temporary file
    unlink(temp_vcf)
    debug_log(paste("BGZF conversion completed for:", input_vcf))
  }
  
  # Convert existing files
  if (genomic_exists) convert_to_bgzf(genomic_vcf)
  if (non_genomic_exists) convert_to_bgzf(non_genomic_vcf)
  
  # Step 2: Index BGZF files
  debug_log("Indexing VCF files...")
  
  # Index VCFs if present
  if (genomic_exists) {
    cmd <- paste(bcftools_path, "index", genomic_vcf)
    run_command(cmd)
  }
  if (non_genomic_exists) {
    cmd <- paste(bcftools_path, "index", non_genomic_vcf)
    run_command(cmd)
  }
  
  # Step 3: Add chromosome and position IDs
  debug_log("Adding chromosome and position IDs...")
  
  # Add IDs and index VCFs if present
  genomic_with_id <- file.path(data_dir, "genomic_with_id.vcf.gz")
  non_genomic_with_id <- file.path(data_dir, "non_genomic_with_id.vcf.gz")
  if (genomic_exists) {
    cmd <- paste(bcftools_path, "annotate --set-id '%CHROM\\_%POS'", genomic_vcf, "-Oz -o", genomic_with_id)
    run_command(cmd)
    cmd <- paste(bcftools_path, "index", genomic_with_id)
    run_command(cmd)
    # Validate: must include #CHROM header and at least one variant record
    hdr <- run_command(paste(bcftools_path, "view -h", genomic_with_id))
    has_chrom <- any(grepl("^#CHROM", hdr))
    first_var <- run_command(paste(bcftools_path, "view -H", genomic_with_id, "| head -n 1"))
    has_var <- length(first_var) > 0
    if (!has_chrom || !has_var) {
      debug_log("Annotated genomic VCF missing #CHROM header or variants", "ERROR")
      stop("Invalid annotated genomic VCF: missing header or variants")
    }
  }
  if (non_genomic_exists) {
    cmd <- paste(bcftools_path, "annotate --set-id '%CHROM\\_%POS'", non_genomic_vcf, "-Oz -o", non_genomic_with_id)
    run_command(cmd)
    cmd <- paste(bcftools_path, "index", non_genomic_with_id)
    run_command(cmd)
    hdr <- run_command(paste(bcftools_path, "view -h", non_genomic_with_id))
    has_chrom <- any(grepl("^#CHROM", hdr))
    first_var <- run_command(paste(bcftools_path, "view -H", non_genomic_with_id, "| head -n 1"))
    has_var <- length(first_var) > 0
    if (!has_chrom || !has_var) {
      debug_log("Annotated non-genomic VCF missing #CHROM header or variants", "ERROR")
      stop("Invalid annotated non-genomic VCF: missing header or variants")
    }
  }
  
  # Step 4: Convert VCF files using PLINK2 and filter multiallelic variants
  debug_log("Converting VCF files to PLINK format using PLINK2 (filtering multiallelic variants)...")
  
  # Convert genomic VCF to PLINK format using PLINK2 (if present)
  if (genomic_exists) {
    debug_log("Processing genomic VCF with PLINK2...")
    cmd <- sprintf("--vcf %s --max-alleles 2 --make-bed --out %s/genomic --chr-set %s --allow-no-sex",
                  genomic_with_id, data_dir, chr_set)
    run_plink(cmd, plink2_path)
  }
  
  # Convert non-genomic VCF to PLINK format using PLINK2 (if present)
  if (non_genomic_exists) {
    debug_log("Processing non-genomic VCF with PLINK2...")
    cmd <- sprintf("--vcf %s --max-alleles 2 --make-bed --out %s/non_genomic --chr-set %s --allow-no-sex",
                  non_genomic_with_id, data_dir, chr_set)
    run_plink(cmd, plink2_path)
  }
  
  # Verify conversion results
  genomic_files <- c(
    file.path(data_dir, "genomic.bed"),
    file.path(data_dir, "genomic.bim"), 
    file.path(data_dir, "genomic.fam")
  )
  
  non_genomic_files <- c(
    file.path(data_dir, "non_genomic.bed"),
    file.path(data_dir, "non_genomic.bim"),
    file.path(data_dir, "non_genomic.fam")
  )
  
  debug_log("Verifying conversion results...")
  if (genomic_exists && !all(file.exists(genomic_files))) {
    debug_log("Failed to create genomic PLINK files", "ERROR")
    debug_log(paste("Missing files:", paste(genomic_files[!file.exists(genomic_files)], collapse = ", ")))
    stop("Failed to create genomic PLINK files")
  }
  if (non_genomic_exists && !all(file.exists(non_genomic_files))) {
    debug_log("Failed to create non-genomic PLINK files", "WARNING")
    debug_log(paste("Missing files:", paste(non_genomic_files[!file.exists(non_genomic_files)], collapse = ", ")))
  }
  
  # Count number of variants after conversion
  genomic_variants <- if (genomic_exists) safe_count_lines(file.path(data_dir, "genomic.bim")) else 0
  non_genomic_variants <- if (non_genomic_exists) safe_count_lines(file.path(data_dir, "non_genomic.bim")) else 0
  
  debug_log("VCF to PLINK conversion completed successfully!")
  debug_log(paste("Genomic variants after filtering:", genomic_variants))
  debug_log(paste("Non-genomic variants after filtering:", non_genomic_variants))
  
  return(c(if (genomic_exists) genomic_with_id else NA,
           if (non_genomic_exists) non_genomic_with_id else NA))
}

# Main analysis function
main <- function() {
  tryCatch({
    debug_log("=== Starting VCF analysis ===")
    
    # Read command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    debug_log(paste("Command line arguments:", paste(args, collapse = ", ")))
    
    if (length(args) == 0) {
      stop("No data directory provided")
    }
    data_dir <- args[1]
    debug_log(paste("Data directory:", data_dir))
    
    # Ensure we're in the correct working directory (avoid hardcoded setwd)
    debug_log(paste("Working directory:", getwd()))
    
    # Resolve tool paths from environment PATH
    plink1_path <- Sys.which("plink")      # PLINK1 for downstream analysis
    plink2_path <- Sys.which("plink2")     # PLINK2 for VCF conversion
    bcftools_path <- Sys.which("bcftools")
    tabix_path <- Sys.which("tabix")
    
    debug_log("=== Tool verification ===")
    # Verify tools exist
    plink1_exists <- nzchar(plink1_path) && file.exists(plink1_path)
    plink2_exists <- nzchar(plink2_path) && file.exists(plink2_path)
    
    debug_log(paste("bcftools path:", bcftools_path, "- exists:", nzchar(bcftools_path) && file.exists(bcftools_path)))
    debug_log(paste("PLINK1 path:", plink1_path, "- exists:", plink1_exists))
    debug_log(paste("PLINK2 path:", plink2_path, "- exists:", plink2_exists))
    debug_log(paste("tabix path:", tabix_path, "- exists:", nzchar(tabix_path) && file.exists(tabix_path)))
    
    if (!nzchar(bcftools_path)) {
        stop("bcftools not found at: ", bcftools_path)
    }
    if (!plink2_exists) {
        stop("plink2 not found at: ", plink2_path, " (required for VCF processing)")
    }
    if (!plink1_exists) {
        stop("plink1 not found at: ", plink1_path)
    }
    if (!nzchar(tabix_path)) {
        stop("tabix not found at: ", tabix_path)
    }
    
    debug_log("All tools verified successfully")
    
    # Load and verify parameters
    debug_log("=== Loading parameters ===")
    params_file <- file.path(data_dir, "params.RData")
    debug_log(paste("Params file:", params_file))
    
    if (!file.exists(params_file)) {
      stop("Cannot find params.RData in ", data_dir)
    }
    load(params_file)
    
    # Verify required parameters
    if (!exists("res")) {
      stop("params.RData does not contain 'res' object")
    }
    
    debug_log(paste("Loaded parameters. res object class:", class(res)))
    debug_log(paste("res object names:", paste(names(res), collapse = ", ")))
    
    required_params <- c("chr_set", "file_type", "dir")
    missing_params <- required_params[!required_params %in% names(res)]
    if (length(missing_params) > 0) {
      stop("Missing required parameters: ", paste(missing_params, collapse = ", "))
    }
    
    chr_set <- res$chr_set
    debug_log(paste("Using chromosome set:", chr_set))
    
    # Process VCF files if input type is VCF (use PLINK2 for conversion)
    if (res$file_type == "vcf") {
      debug_log("=== VCF Processing Phase ===")
      vcf_files <- process_vcf_files(data_dir, chr_set, bcftools_path, plink2_path, plink1_path)
    } else {
      stop("This script is designed for VCF input. Please use appropriate script for other formats.")
    }
    
    # === Use PLINK1 for downstream analysis from here ===
    debug_log("=== Starting PLINK1 Analysis Phase ===")
    
    # Step 1: Verify input files exist
    genomic_fam <- file.path(data_dir, "genomic.fam")
    non_genomic_fam <- file.path(data_dir, "non_genomic.fam")
    
    debug_log(paste("Checking genomic fam file:", genomic_fam))
    debug_log(paste("Checking non-genomic fam file:", non_genomic_fam))
    
    # Relaxed check: allow single dataset (continue if either exists)
    if (!file.exists(genomic_fam) && !file.exists(non_genomic_fam)) {
      stop("PLINK files not found after conversion. Neither genomic nor non_genomic dataset exists.")
    }
    
    # Step 2: Merge or skip merge (single dataset)
    debug_log("=== Combining datasets (skip if single dataset) ===")

    genomic_exists <- file.exists(file.path(data_dir, "genomic.fam")) &&
                      file.exists(file.path(data_dir, "genomic.bed")) &&
                      file.exists(file.path(data_dir, "genomic.bim"))
    non_genomic_exists <- file.exists(file.path(data_dir, "non_genomic.fam")) &&
                          file.exists(file.path(data_dir, "non_genomic.bed")) &&
                          file.exists(file.path(data_dir, "non_genomic.bim"))

    if (genomic_exists && non_genomic_exists) {
      # Both datasets exist: normalize missing variant IDs first to avoid '.' causing duplicate IDs and multiallelic warnings
      debug_log("Normalizing missing variant IDs for both datasets")
      genomic_norm_cmd <- sprintf(
        "--bfile %s/genomic --set-missing-var-ids @:# --make-bed --out %s/genomic_norm --chr-set %s --allow-no-sex",
        data_dir, data_dir, chr_set
      )
      run_plink(genomic_norm_cmd, plink1_path)
      non_genomic_norm_cmd <- sprintf(
        "--bfile %s/non_genomic --set-missing-var-ids @:# --make-bed --out %s/non_genomic_norm --chr-set %s --allow-no-sex",
        data_dir, data_dir, chr_set
      )
      run_plink(non_genomic_norm_cmd, plink1_path)

      # 规范后执行合并
      merge_cmd <- sprintf(
        "--bfile %s/genomic_norm --bmerge %s/non_genomic_norm --make-bed --out %s/combined_data --chr-set %s --allow-no-sex",
        data_dir, data_dir, data_dir, chr_set
      )
      merge_success <- tryCatch({
        run_plink(merge_cmd, plink1_path)
        TRUE
      }, error = function(e) {
        debug_log(paste("Direct merge failed:", e$message), "WARNING")
        FALSE
      })
      if (!merge_success || !file.exists(file.path(data_dir, "combined_data.fam"))) {
        debug_log("Performing SNP consistency check with PLINK1...")
        merge_check_cmd <- sprintf(
          "--bfile %s/genomic_norm --bmerge %s/non_genomic_norm --merge-mode 6 --out %s/merge_check --chr-set %s --allow-no-sex",
          data_dir, data_dir, data_dir, chr_set
        )
        tryCatch({
          run_plink(merge_check_cmd, plink1_path)
        }, error = function(e) {
          debug_log("Merge check completed with conflicts")
        })
        missnp_file <- file.path(data_dir, "merge_check.missnp")
        if (file.exists(missnp_file)) {
          debug_log("Found inconsistent SNPs, cleaning data with PLINK1...")
          genomic_clean_cmd <- sprintf(
            "--bfile %s/genomic_norm --exclude %s/merge_check.missnp --make-bed --out %s/genomic_clean --chr-set %s --allow-no-sex",
            data_dir, data_dir, data_dir, chr_set
          )
          run_plink(genomic_clean_cmd, plink1_path)
          non_genomic_clean_cmd <- sprintf(
            "--bfile %s/non_genomic_norm --exclude %s/merge_check.missnp --make-bed --out %s/non_genomic_clean --chr-set %s --allow-no-sex",
            data_dir, data_dir, data_dir, chr_set
          )
          run_plink(non_genomic_clean_cmd, plink1_path)
          final_merge_cmd <- sprintf(
            "--bfile %s/genomic_clean --bmerge %s/non_genomic_clean --make-bed --out %s/combined_data --chr-set %s --allow-no-sex",
            data_dir, data_dir, data_dir, chr_set
          )
          run_plink(final_merge_cmd, plink1_path)
          excluded_snps <- safe_count_lines(missnp_file)
          debug_log(paste("Number of inconsistent SNPs removed:", excluded_snps))
        } else {
          stop("Merge failed but no inconsistent SNPs found. Check input files.")
        }
      }
      combined_fam <- file.path(data_dir, "combined_data.fam")
      if (!file.exists(combined_fam)) {
        stop("Dataset merge failed. Please check input files.")
      }
      debug_log("Dataset merge completed successfully!")
      # 基于 keep-base 进行位点过滤（在合并之后）
      if (!is.null(res$keep_base) && nchar(res$keep_base) > 0) {
        # BIM ID 标准化已在 Python 过滤脚本中完成，无需在此处 awk 处理
        debug_log("=== Base-set filtering on combined_data BIM ===")
        snp_list <- file.path(data_dir, "combined_data_filtered_snps.txt")
        filter_py <- res$script_filter_py
        cmd <- sprintf(
          "python '%s' '%s/combined_data.bim' --output '%s' --keep-base '%s'",
          filter_py, data_dir, snp_list, res$keep_base
        )
        debug_log(paste("Running base-set filter:", cmd))
        system(cmd)
        if (file.exists(snp_list) && file.info(snp_list)$size > 0) {
          debug_log("Applying PLINK extract to create filtered combined dataset")
          filt_cmd <- sprintf(
            "--bfile %s/combined_data --extract %s --make-bed --out %s/combined_data_filtered --chr-set %s --allow-no-sex",
            data_dir, snp_list, data_dir, chr_set
          )
          run_plink(filt_cmd, plink1_path)
          bfileprefix <- file.path(data_dir, "combined_data_filtered")
        } else {
          debug_log("No SNPs remained after base-set filtering; proceeding with unfiltered combined_data", "WARNING")
          bfileprefix <- file.path(data_dir, "combined_data")
        }
      } else {
        bfileprefix <- file.path(data_dir, "combined_data")
      }
      # 计算频率（基于最终使用的 bfileprefix）
      freq_cmd <- sprintf(
        "--bfile %s --freq --out %s --chr-set %s",
        bfileprefix, bfileprefix, chr_set
      )
      run_plink(freq_cmd, plink1_path)
      debug_log(paste("bfileprefix set to:", bfileprefix))
      n_individuals <- safe_count_lines(paste0(bfileprefix, ".fam"))
      n_variants <- safe_count_lines(paste0(bfileprefix, ".bim"))
      debug_log(paste("Combined dataset: ", n_individuals, " individuals, ", n_variants, " variants"))
    } else if (genomic_exists || non_genomic_exists) {
      debug_log("Single dataset detected. Skipping merge and using it directly.")
      single_prefix <- if (genomic_exists) "genomic" else "non_genomic"
      file.copy(file.path(data_dir, paste0(single_prefix, ".bed")), file.path(data_dir, "combined_data.bed"), overwrite = TRUE)
      file.copy(file.path(data_dir, paste0(single_prefix, ".bim")), file.path(data_dir, "combined_data.bim"), overwrite = TRUE)
      file.copy(file.path(data_dir, paste0(single_prefix, ".fam")), file.path(data_dir, "combined_data.fam"), overwrite = TRUE)
      # 基于 keep-base 进行位点过滤（单数据集路径）
      if (!is.null(res$keep_base) && nchar(res$keep_base) > 0) {
        # BIM ID 标准化已在 Python 过滤脚本中完成，无需在此处 awk 处理
        debug_log("=== Base-set filtering on single dataset BIM ===")
        snp_list <- file.path(data_dir, "combined_data_filtered_snps.txt")
        filter_py <- res$script_filter_py
        cmd <- sprintf(
          "python '%s' '%s/combined_data.bim' --output '%s' --keep-base '%s'",
          filter_py, data_dir, snp_list, res$keep_base
        )
        debug_log(paste("Running base-set filter:", cmd))
        system(cmd)
        if (file.exists(snp_list) && file.info(snp_list)$size > 0) {
          debug_log("Applying PLINK extract to create filtered dataset")
          filt_cmd <- sprintf(
            "--bfile %s/combined_data --extract %s --make-bed --out %s/combined_data_filtered --chr-set %s --allow-no-sex",
            data_dir, snp_list, data_dir, chr_set
          )
          run_plink(filt_cmd, plink1_path)
          bfileprefix <- file.path(data_dir, "combined_data_filtered")
        } else {
          debug_log("No SNPs remained after base-set filtering; proceeding with unfiltered combined_data", "WARNING")
          bfileprefix <- file.path(data_dir, "combined_data")
        }
      } else {
        bfileprefix <- file.path(data_dir, "combined_data")
      }
      # 计算频率（基于最终使用的 bfileprefix）
      freq_cmd <- sprintf(
        "--bfile %s --freq --out %s --chr-set %s",
        bfileprefix, bfileprefix, chr_set
      )
      run_plink(freq_cmd, plink1_path)
      debug_log(paste("bfileprefix set to:", bfileprefix))
      n_individuals <- safe_count_lines(paste0(bfileprefix, ".fam"))
      n_variants <- safe_count_lines(paste0(bfileprefix, ".bim"))
      debug_log(paste("Single dataset: ", n_individuals, " individuals, ", n_variants, " variants"))
    } else {
      stop("PLINK files not found after VCF conversion. Expected genomic/non_genomic datasets or a single dataset in ", data_dir)
    }
    
    # Step 4: 计算有效标记数
    debug_log("=== Starting kinship analysis ===")
    n <- n_individuals
    npairs <- n * (n - 1) / 2
    debug_log(paste("Number of pairs to analyze:", npairs))
    
    debug_log("Calculating minimum effective markers...")
    me.min <- me.min(theta = (1 / 2)^(0:2), alpha = 0.05, beta = 0.1)
    debug_log("Minimum effective markers calculated")
    print(me.min)
    
    # Step 5: Calculate effective markers (使用修改后的函数)
    debug_log("=== Calculating effective markers ===")
    
    me <- modified_calculate.me(
      bfileprefix = bfileprefix,
      method = "GRM",
      plink1_path = plink1_path,
      plink2_path = if(plink2_exists) plink2_path else NULL,
      pop_size = n,
      chr_set = chr_set
    )
    debug_log(paste("Final effective markers value:", me))
    
    # Step 6: Generate summary report
    debug_log("=== Generating deepKin summary ===")
    dK <- deepKin(n = n, me = me, alpha = 0.05, beta = 0.1, max.degree = 2)
    summary_text <- deepKin.summary(dK)
    debug_log("deepKin summary generated")
    
    # Save summary to file
    writeLines(summary_text, file.path(data_dir, "deepkin_summary.txt"))
    
    # Save plots
    # debug_log("Saving plots...")
    # pdf(file.path(data_dir, "Minimum_Effective_Markers_by_Degree.pdf"))
    # plot(dK$me.min$degree, dK$me.min$Me.min,
    #      ylab = "me", xlab = "Degree",
    #      main = "Minimum Effective Markers by Degree")
    # dev.off()
    
    # pdf(file.path(data_dir, "Maximum_Power_by_Degree.pdf"))
    # plot(dK$power.max$degree, dK$power.max$Power.max,
    #      ylab = "Power", xlab = "Degree",
    #      main = "Maximum Power by Degree")
    # abline(h = 0.9, col = "blue")
    # abline(v = dK$delta, col = "red")
    # dev.off()
    
    # Get critical values
    theta.min <- dK$theta.min
    degree.deep <- dK$delta
    debug_log(paste("Critical values - Theta minimum:", theta.min, ", Degree deep:", degree.deep))
    
    # Step 7: Perform deepKin estimation
    debug_log("=== Performing deepKin estimation ===")
    
    # 添加额外的参数检查
    if (is.null(me) || !is.numeric(me) || is.na(me)) {
      stop("Invalid effective markers value: ", me)
    }
    
    deepkin <- deepKin.estimation(bfileprefix, plink1_path, me, xcohort = FALSE, pop_size1 = n)
    debug_log(paste("deepKin estimation completed. Results have", nrow(deepkin), "rows"))
    
    write.table(
      deepkin,
      file = file.path(data_dir, "results.deepkin"),
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE
    )
    
    # Step 8: Perform classification
    debug_log("=== Performing classification ===")
    
    if (nrow(deepkin) == 0) {
      debug_log("No deepkin results to classify", "WARNING")
      deepkin.qc.rst <- data.frame(ID1 = character(0), ID2 = character(0), theta = numeric(0), class = numeric(0))
    } else {
      deepkin.qc <- deepkin[which(deepkin$theta > theta.min), ]
      debug_log(paste("Found", nrow(deepkin.qc), "pairs above theta minimum"))
      
      if (nrow(deepkin.qc) > 0) {
        id <- safe_read_table(paste0(bfileprefix, '.fam'), header = FALSE)
        if (is.null(id)) {
          stop("Failed to read .fam file for ID extraction")
        }
        
        index <- as.numeric(rownames(deepkin.qc))
        deepkin.qc.id <- extract.indi.id(id = id[,2], index = index, xcohort = FALSE)
        
        deepkin.qc.class <- deepKin.classification(deepkin.qc$theta, me = me, alpha = 0.05)
        deepkin.qc.rst <- cbind(deepkin.qc.id, deepkin.qc.class)
        debug_log(paste("Classification completed. Final results have", nrow(deepkin.qc.rst), "pairs"))
      } else {
        debug_log("No pairs above theta minimum", "WARNING")
        deepkin.qc.rst <- data.frame(ID1 = character(0), ID2 = character(0), theta = numeric(0), class = numeric(0))
      }
    }
    
    # Step 9: 准备相关性结果（不绘制heatmap）
    debug_log("=== Preparing relatedness outputs ===")

    # 修正：使用 <= 1 过滤显著配对
    if (nrow(deepkin.qc.rst) > 0) {
      significant_pairs <- deepkin.qc.rst[deepkin.qc.rst$class <= 1, ]
      debug_log(paste("Found", nrow(significant_pairs), "significant pairs (class <= 1)"))
    } else {
      significant_pairs <- deepkin.qc.rst
      debug_log("No pairs to filter for significance")
    }

    # Write out both full and significant pair datasets (for reference/debugging)
    write.table(
      deepkin.qc.rst, 
      file = paste0(bfileprefix, '.related'), 
      quote = FALSE, 
      col.names = TRUE, 
      row.names = FALSE
    )
    write.table(
      significant_pairs, 
      file = paste0(bfileprefix, '_significant_pairs.related'), 
      quote = FALSE, 
      col.names = TRUE, 
      row.names = FALSE
    )


    
    # Step 10: Generate analysis summary
    debug_log("=== Generating analysis summary ===")
    fam <- safe_read_table(file.path(data_dir, "combined_data.fam"))
    bim <- safe_read_table(file.path(data_dir, "combined_data.bim"))
    
    if (is.null(fam) || is.null(bim)) {
      stop("Failed to read combined data files for summary")
    }
    
    analysis_summary <- paste0(
      "Analysis Summary:\n",
      "----------------\n",
      "Input type: VCF files\n",
      "VCF processing: PLINK2 (with --max-alleles 2)\n",
      "Merging and analysis: PLINK1\n",
      "Total samples: ", nrow(fam), "\n",
      "Total SNPs: ", nrow(bim), "\n",
      "Species: ", res$species, "\n",
      "Chromosome set: ", chr_set, "\n",
      "Effective markers: ", me, "\n",
      "Minimum theta: ", theta.min, "\n",
      "Deep degree: ", degree.deep, "\n",
      "Significant pairs found: ", nrow(significant_pairs), "\n"
    )
    
    writeLines(analysis_summary, file.path(data_dir, "analysis_summary.txt"))
    
    # Also save as CSV if needed
    write.table(
      significant_pairs,
      file = file.path(data_dir, "final_results.related.csv"),
      sep = ",",
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE
    )

    # 移除打包与清理步骤（命令行版本不需要打包）
    debug_log("=== Analysis completed successfully! ===")
    debug_log("=== Final Summary ===")
    debug_log("VCF processing: PLINK2 (with --max-alleles 2)")
    debug_log("SNP merging and analysis: PLINK1")
    debug_log(paste("Relatedness matrix:", ifelse(plink2_exists, "PLINK2", "PLINK1")))
    debug_log(paste("Final dataset:", n_individuals, "individuals,", n_variants, "biallelic variants"))
    debug_log(paste("Significant pairs found:", nrow(significant_pairs)))
    debug_log(paste("Output files created in:", data_dir))

    
  }, error = function(e) {
    error_msg <- paste("Error in analysis:", e$message)
    debug_log(error_msg, "ERROR")
    debug_log(paste("Traceback:", paste(traceback(), collapse = "\n")), "ERROR")
    
    # Log the error
    if (exists("data_dir")) {
      log_error(data_dir, error_msg)
    }
    
    # Email sending removed
    
    # Re-throw the error to ensure the script exits with an error status
    stop(error_msg)
  })
}

# Run the analysis
main()
