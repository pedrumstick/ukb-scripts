#!/usr/bin/env Rscript

##############################
# VARIABLES AND CHECKS
##############################
warning("This script will overwrite any old columns. Therefore, you shouldn't use it to update if your new data is a subset of the broader phenotype file with the same column names. For instance, if the data was downloaded from the Data Portal (e.g., death registry). You will need to either assign new column names, or use a different script to avoid losing columns with the same name in the previous phenotype file.")
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {stop("Incorrect number of arguments, please specify: <date for version of old phenotype file as YYYY_MM_DD> </path/to/new/phenotype/file/>")}
old_date <- "2020_07_24" # date for version of old files that you want to merge into
if (!file.exists(paste0("/hightide/UKB_PHENOS/UKB_MERGED_PHENO_BRITISH_",old_date,"_FINAL_INCLUDE_PCs.txt")) && !file.exists(paste0("/hightide/UKB_PHENOS/UKB_MERGED_PHENO_BRITISH_",old_date,"_FINAL_INCLUDE_PCs.txt.gz"))) {stop("The specified date does not have any files associated with it.")}
new_filePath <- "/hightide/UKB_PHENOS/Data_Portal_Downloads/Death_Register/deathRegistry_wide_2020_07_10.txt" # file path for data that you want to merge into the old files
if (!file.exists(new_filePath)) {stop("The specified file to merge in does not exist.")}
withdrawn_consent_file <- "/hightide/UKB_PHENOS/ids_withdrawnConsent_2020_02_05.txt"
if (!file.exists(withdrawn_consent_file)) {stop(paste0("ERROR: Input file does not exist! ", withdrawn_consent_file))}

##############################
# LOAD PACKAGES
##############################
require(data.table)
require(dplyr)

##############################
# LOAD DATA FRAMES
##############################
start_time <- Sys.time()
message("************************************************\nStarted merging script at ", start_time,"\n************************************************")

setwd("/hightide/UKB_PHENOS/")
for (input_file in list.files(pattern=paste0("UKB_MERGED_PHENO_.*_",old_date,"_FINAL_INCLUDE_PCs.txt.*"))) {
    message("Updating phenotypes for... ", input_file)
    current_date <- gsub("-", "_", Sys.Date())
    output_file <- sub("\\..*", "", sub("[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9]", current_date, basename(input_file))) 

    updated_df <- fread(new_filePath, header=TRUE, stringsAsFactors=FALSE, na.strings=c("", "NA", NA))
    message("There are ", ncol(updated_df)-1, " new columns of data in this update.")
    original_df <- fread(input_file, header=TRUE, stringsAsFactors=FALSE, na.strings=c("", "NA", NA))
    exclude_df <- fread(withdrawn_consent_file, header=FALSE, stringsAsFactors=FALSE, na.strings=c("", "NA", NA))

    temp <- left_join(data.frame(original_df, check.names=FALSE), updated_df, by="eid") #need to keep the order the same as the original file
    # temp<-merge(data.frame(original_df, check.names=FALSE), updated_df, by="eid", all.x=TRUE) #need to keep the order the same as the original file
    if (nrow(temp)!=nrow(original_df)) {stop("ERROR: New dataframe doesn't have the same number of rows as the original dataframe! There may be duplicate IDs, troubleshooting is needed...")}

    exclude_rows <- which(temp$eid %in% exclude_df[[1]])
    if (length(exclude_rows)>0) {
        temp <- temp[-exclude_rows, ] # completely excludes individuals that withdrew consent; this may results in incompatibility for downstream analyses
    }
    # temp[exclude_rows, (names(temp)[2:ncol(temp)]) := .SD[NA], .SDcols=names(temp)[2:ncol(temp)]] # alternative where we set individuals that withdrew consent to NA for all variables, but keep their IDs to maintain compatibility with previous files or scripts

    old_columns <- grep(".x$", colnames(temp))
    if (length(old_columns)>0) {
        warning("There were ", length(old_columns), " old columns overwritten with new data.")
        new_df <- temp[,-old_columns] # new_df[, grep(".x$", colnames(new_df)):=NULL]
        names(new_df) <- sub("\\.y","", names(new_df))
    } else {
        message("There were 0 old columns overwritten with new data.")
        new_df <- temp
    }

    if (!dir.exists("temp_updatedPhenotypes/")) {dir.create("temp_updatedPhenotypes/")}

    fwrite(new_df, paste0("temp_updatedPhenotypes/",output_file,".txt.gz"), sep="\t", quote=FALSE, na="NA")
    print(Sys.time()-start_time)
    start_time <- Sys.time() # reset timer for next subset
}

message("************************************************\nFinished merging script at ", Sys.time(),"\n************************************************")