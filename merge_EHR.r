#!/usr/bin/env Rscript

warning("This script is a WIP. It is designed for use with updates of the electronic health record (EHR) data from UKB. It takes the field IDs from the merged phenotype file that already exist in the update of the EHR data, and merges these fields together in a new file. Unlike regular phenotype updates, the EHR data has MANY columns (>24K), which makes working with the data difficulty if it's all in one file. What is the best way to work with this data?")

##############################
# VARIABLES AND CHECKS
##############################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Incorrect number of arguments, please specify: <date for version of old phenotype file as YYYY_MM_DD> </path/to/new/phenotype/file/>")
}
old_date <- args[1] # date for version of old files that you want to merge into (e.g., "2020_07_24")
new_filePath <- args[2] # file path for data that you want to merge into the old files (e.g., "/hightide/UKB_PHENOS/Data_Portal_Downloads/Death_Register/deathRegistry_wide_2020_07_10.txt")
if (!file.exists(paste0("/hightide/UKB_PHENOS/UKB_MERGED_PHENO_BRITISH_", old_date, "_FINAL_INCLUDE_PCs.txt")) && !file.exists(paste0("/hightide/UKB_PHENOS/UKB_MERGED_PHENO_BRITISH_", old_date, "_FINAL_INCLUDE_PCs.txt.gz"))) {
    stop("The specified date does not have any files associated with it.")
}
if (!file.exists(new_filePath)) {
    stop("The specified file to merge in does not exist.")
}
withdrawn_consent_file <- "/hightide/UKB_PHENOS/ids_withdrawnConsent_2020_08_20.txt"
if (!file.exists(withdrawn_consent_file)) {
    stop(paste0("ERROR: Input file does not exist! ", withdrawn_consent_file))
}

##############################
# LOAD PACKAGES
##############################
require(data.table)
require(dplyr)

##############################
# LOAD DATA
##############################
start_time <- Sys.time()
message("************************************************\nStarted merging script at ", start_time, "\n************************************************")

setwd("/hightide/UKB_PHENOS/")

updated_df <- fread(new_filePath, header = TRUE, colClasses = "character", na.strings = c("", "NA", NA))
regex<-paste0("^", unique(gsub("-.*", "", names(updated_df)[-1])), collapse = "|")

for (input_file in list.files(pattern = paste0("UKB_MERGED_PHENO_.*_", old_date, "_FINAL_INCLUDE_PCs.txt.*"))) {
    message("Updating phenotypes for... ", input_file)
    current_date <- gsub("-", "_", Sys.Date())
    output_file <- paste0("UKB_EHR_Fields_", sub("_[0-9][0-9][0-9][0-9].*", "", sub("UKB_MERGED_PHENO_", "", basename(input_file))), "_", current_date)

    original_df <- fread(input_file, header = TRUE, colClasses = "character", na.strings = c("", "NA", NA))
    exclude_df <- fread(withdrawn_consent_file, header = FALSE, colClasses = "character", na.strings = c("", "NA", NA))
    # extract columns belonging to the fields you want to merge from existing data
    original_df <- original_df[, .SD, .SDcols = names(original_df) %like% paste("eid", regex, sep = "|")]
    # merge with old phenotype file
    output_df <- merge(original_df, updated_df, by = "eid", all.x = TRUE)
    # check number of SNPs
     if (nrow(output_df) != nrow(original_df)) {
        stop("ERROR: New dataframe doesn't have the same number of rows as the original dataframe! There may be duplicate IDs, troubleshooting is needed...")
    }
    exclude_rows <- which(output_df$eid %in% exclude_df[[1]])
    if (length(exclude_rows) > 0) {
        output_df <- output_df[-exclude_rows, ] # completely excludes individuals that withdrew consent; this may results in incompatibility for downstream analyses
    }
    old_columns <- grep(".x$", colnames(output_df))
    if (length(old_columns) > 0) {
        warning("There were ", length(old_columns), " old columns overwritten with new data.")
        output_df[, grep(".x$", colnames(output_df)) := NULL]
        names(output_df) <- sub("\\.y", "", names(output_df))
    } else {
        message("There were 0 old columns overwritten with new data.")
    }
    if (!dir.exists("temp_updatedPhenotypes/")) {
        dir.create("temp_updatedPhenotypes/")
    }
    fwrite(output_df, paste0("temp_updatedPhenotypes/", output_file, ".txt.gz"), sep = "\t", quote = FALSE, na = "NA")
    print(Sys.time() - start_time)
    start_time <- Sys.time() # reset timer for next subset
}

message("************************************************\nFinished merging script at ", Sys.time(), "\n************************************************")