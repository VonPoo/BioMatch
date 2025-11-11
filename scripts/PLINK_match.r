# PLINK_Geno_nonGeno.R - Backend processing script for Sample Comparison Analysis (DEBUG VERSION)
#
# Purpose
# - Perform sample comparison analysis using PLINK1/PLINK2 and deepKin.
# - Support merging two datasets (genomic and non-genomic) or processing a single dataset.
# - Optionally filter SNPs by a base-set (keep-base) prior to deepKin analysis.
#
# Inputs
# - A `data_dir` containing PLINK files (`genomic.*` and/or `non_genomic.*`) or converted outputs.
# - `params.RData` with a `res` list including at least: `chr_set`, `file_type`, `dir`.
# - Optional: `res$keep_base` and `res$script_filter_py` for base-set SNP filtering.
#
# Outputs
# - deepKin summary (`deepkin_summary.txt`), relatedness results (`*.related`),
#   significant pairs (`*_significant_pairs.related`), analysis summary (`analysis_summary.txt`),
#   plots, and CSV output.
#
# Usage
#   Rscript PLINK_Geno_nonGeno.R /path/to/data_dir
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

# Function to run PLINK commands with biomatch argv0
run_plink <- function(cmd, plink_path, chr_set = NULL) {
  # Ensure we always include --chr-set if it's provided and not already in the command
  if (!is.null(chr_set) && !grepl("--chr-set", cmd)) {
    cmd <- paste(cmd, "--chr-set", chr_set)
  }
  full_cmd <- paste(plink_path, cmd)
  debug_log(paste("Running PLINK command:", full_cmd))
  # Use double quotes to avoid conflicts with single quotes in arguments
  wrapped <- sprintf('bash -lc "exec -a biomatch %s"', full_cmd)
  result <- system(wrapped, intern = TRUE)
  debug_log(paste("PLINK output:", paste(result, collapse = "; ")))
  return(result)
}

# Main analysis function
main <- function() {
  tryCatch({
    debug_log("=== Starting PLINK analysis ===")
    
    # Read command line arguments for the data directory
    args <- commandArgs(trailingOnly = TRUE)
    debug_log(paste("Command line arguments:", paste(args, collapse = ", ")))
    
    if (length(args) == 0) {
      stop("No data directory provided")
    }
    data_dir <- args[1]
    debug_log(paste("Data directory:", data_dir))
    
    # Keep current working directory; avoid hardcoded setwd
    debug_log(paste("Working directory:", getwd()))
    
    # Resolve tool paths from PATH
    plink1_path <- Sys.which("plink")      # PLINK1用于合并和大部分分析
    plink2_path <- Sys.which("plink2")     # PLINK2用于特定功能
    
    debug_log("=== Tool verification ===")
    # Verify tools exist
    plink1_exists <- nzchar(plink1_path) && file.exists(plink1_path)
    plink2_exists <- nzchar(plink2_path) && file.exists(plink2_path)
    
    debug_log(paste("PLINK1 path:", plink1_path, "- exists:", plink1_exists))
    debug_log(paste("PLINK2 path:", plink2_path, "- exists:", plink2_exists))
    
    if (!plink1_exists) {
      stop("PLINK1 not found at: ", plink1_path)
    }
    
    debug_log("Using tools:")
    debug_log(paste("PLINK1:", plink1_path, "(", ifelse(plink1_exists, "Found", "NOT FOUND"), ")"))
    debug_log(paste("PLINK2:", plink2_path, "(", ifelse(plink2_exists, "Found", "NOT FOUND"), ")"))
    
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
    
    # === Start analysis: PLINK1 for most operations ===
    debug_log("=== Starting analysis with PLINK1/PLINK2 hybrid approach ===")
    
    # Step 1: Verify input files exist
    genomic_fam <- file.path(data_dir, "genomic.fam")
    non_genomic_fam <- file.path(data_dir, "non_genomic.fam")
    
    debug_log(paste("Checking genomic fam file:", genomic_fam))
    debug_log(paste("Checking non-genomic fam file:", non_genomic_fam))
    
    # Relaxed check: exit only if both datasets are missing; allow single dataset
    if (!file.exists(genomic_fam) && !file.exists(non_genomic_fam)) {
      stop("PLINK files not found. Expected at least one of genomic.fam or non_genomic.fam in ", data_dir)
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

      # Perform merge after normalization
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

      # If direct merge fails, perform conflict check and cleanup
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

      # Validate merge result
      combined_fam <- file.path(data_dir, "combined_data.fam")
      if (!file.exists(combined_fam)) {
        stop("Dataset merge failed. Please check input files.")
      }

      debug_log("Dataset merge completed successfully!")
      # Filter SNPs by keep-base (after merge)
      if (!is.null(res$keep_base) && nchar(res$keep_base) > 0) {
        # BIM ID normalization is handled in the Python filter script; no awk needed here
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

      # Calculate frequencies (based on the final `bfileprefix`)
      debug_log("=== Calculating frequencies with PLINK1 ===")
      freq_cmd <- sprintf(
        "--bfile %s --freq --out %s --chr-set %s",
        bfileprefix, bfileprefix, chr_set
      )
      run_plink(freq_cmd, plink1_path)

      debug_log(paste("bfileprefix set to:", bfileprefix))
      n_individuals <- safe_count_lines(paste0(bfileprefix, ".fam"))
      n_variants <- safe_count_lines(paste0(bfileprefix, ".bim"))
      debug_log(paste("Combined dataset:", n_individuals, "individuals,", n_variants, "variants"))
    } else if (genomic_exists || non_genomic_exists) {
      # Single dataset detected; skip merge and use directly as combined_data
      debug_log("Single dataset detected. Skipping merge and using it directly.")
      single_prefix <- if (genomic_exists) "genomic" else "non_genomic"
      # 拷贝为 combined_data 前缀
      file.copy(file.path(data_dir, paste0(single_prefix, ".bed")), file.path(data_dir, "combined_data.bed"), overwrite = TRUE)
      file.copy(file.path(data_dir, paste0(single_prefix, ".bim")), file.path(data_dir, "combined_data.bim"), overwrite = TRUE)
      file.copy(file.path(data_dir, paste0(single_prefix, ".fam")), file.path(data_dir, "combined_data.fam"), overwrite = TRUE)

      # Filter SNPs by keep-base (single dataset path)
      if (!is.null(res$keep_base) && nchar(res$keep_base) > 0) {
        # BIM ID normalization is handled in the Python filter script; no awk needed here
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

      # Calculate frequencies (based on the final `bfileprefix`)
      freq_cmd <- sprintf(
        "--bfile %s --freq --out %s --chr-set %s",
        bfileprefix, bfileprefix, chr_set
      )
      run_plink(freq_cmd, plink1_path)

      debug_log(paste("bfileprefix set to:", bfileprefix))
      n_individuals <- safe_count_lines(paste0(bfileprefix, ".fam"))
      n_variants <- safe_count_lines(paste0(bfileprefix, ".bim"))
      debug_log(paste("Single dataset:", n_individuals, "individuals,", n_variants, "variants"))
    } else {
      stop("PLINK files not found. Expected genomic/non_genomic datasets or a single dataset in ", data_dir)
    }

    # Step 4: Compute effective markers (prefer PLINK2, fallback to PLINK1)
    debug_log("=== Starting kinship analysis ===")
    n <- n_individuals
    npairs <- n * (n - 1) / 2
    debug_log(paste("Number of pairs to analyze:", npairs))
    
    debug_log("Calculating minimum effective markers...")
    me.min <- me.min(theta = (1 / 2)^(0:2), alpha = 0.05, beta = 0.1)
    debug_log("Minimum effective markers calculated")
    print(me.min)
    
    # Step 5: Calculate effective markers
    debug_log("=== Calculating effective markers ===")
    
    # Use modified calculate.me supporting both PLINK1 and PLINK2
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
    debug_log("Saving plots...")
    pdf(file.path(data_dir, "Minimum_Effective_Markers_by_Degree.pdf"))
    plot(dK$me.min$degree, dK$me.min$Me.min,
         ylab = "me", xlab = "Degree",
         main = "Minimum Effective Markers by Degree")
    dev.off()
    
    pdf(file.path(data_dir, "Maximum_Power_by_Degree.pdf"))
    plot(dK$power.max$degree, dK$power.max$Power.max,
         ylab = "Power", xlab = "Degree",
         main = "Maximum Power by Degree")
    abline(h = 0.9, col = "blue")
    abline(v = dK$delta, col = "red")
    dev.off()
    
    # Get critical values
    theta.min <- dK$theta.min
    degree.deep <- dK$delta
    debug_log(paste("Critical values - Theta minimum:", theta.min, ", Degree deep:", degree.deep))
    
    # Step 7: Perform deepKin estimation
    debug_log("=== Performing deepKin estimation ===")
    
    # Additional parameter checks
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

    # Step 9: Prepare relatedness outputs (no heatmap)
    debug_log("=== Preparing relatedness outputs ===")

    # Filter significant pairs using class <= 1
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
      "Input type: PLINK files\n",
      "Analysis approach: Hybrid PLINK1/PLINK2\n",
      "Merging and most analysis: PLINK1\n",
      "Relatedness matrix: ", ifelse(plink2_exists, "PLINK2", "PLINK1"), "\n",
      "Total samples: ", nrow(fam), "\n",
      "Total SNPs: ", nrow(bim), "\n",
      "Species: ", res$species, "\n",
      "Chromosome set: ", chr_set, "\n",
      "Effective markers: ", me, "\n",
      "Minimum theta: ", theta.min, "\n",
      "Deep degree: ", degree.deep, "\n",
      "Significant pairs found: ", nrow(significant_pairs), "\n",
      "PLINK1 available: ", plink1_exists, "\n",
      "PLINK2 available: ", plink2_exists, "\n"
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

    # Remove packaging/cleanup steps (CLI version does not require packaging)
    debug_log("=== Analysis completed successfully! ===")
    debug_log("=== Final Summary ===")
    debug_log("Analysis approach: Hybrid PLINK1/PLINK2")
    debug_log("Merging and most analysis: PLINK1")
    debug_log(paste("Relatedness matrix:", ifelse(plink2_exists, "PLINK2", "PLINK1")))
    debug_log(paste("Final dataset:", n_individuals, "individuals,", n_variants, "variants"))
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
