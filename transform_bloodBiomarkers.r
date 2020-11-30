#!/usr/bin/env Rscript

##############################
# VARIABLES AND CHECKS
##############################
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Incorrect number of arguments, please specify: </path/to/phenotype/directory/> <date for version of phenotype file as YYYY_MM_DD> </path/to/output/directory/>")
}
pheno_dir <- args[1]
pheno_date <- args[2]
output_dir <- args[3]
if (!file.exists(paste0(pheno_dir,"/UKB_MERGED_PHENO_BRITISH_", pheno_date, "_FINAL_INCLUDE_PCs.txt.gz"))) {
    stop("The specified phenotype directory and date does not have any files associated with it: ", paste0(pheno_dir,"/UKB_MERGED_PHENO_BRITISH_", pheno_date, "_FINAL_INCLUDE_PCs.txt"))
}

if (!dir.exists(output_dir)){warning("Specified output directory does not exist. Creating now..."); dir.create(output_dir)}

#############################################
# Functions
#############################################
winsorize <- function(x, sd_threshold=5) {
	x<-as.numeric(x)
    mean <- mean(x, na.rm=TRUE)
    upper_thresh <- mean + sd_threshold*sd(x, na.rm=TRUE)
    lower_thresh <- mean - sd_threshold*sd(x, na.rm=TRUE)
    x[which(x>upper_thresh)]<-upper_thresh
	x[which(x<lower_thresh)]<-lower_thresh
	return(x)
}
quantile_normalize <- function(x) {
	return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
}

#############################################
# Extract and clean blood biomarker columns
# Extract relevant columns, calculate change in biomarkers between visits, winsorize, and quantile normalize biomarkers for all ethnicities
# To get difference between columns programmatically (https://stackoverflow.com/questions/35785270/data-table-difference-set-of-columns-from-another-column)
#############################################
setwd(output_dir)

start_time <- Sys.time()
message("************************************************\nStarted script at ", start_time, "\n************************************************")
message("Files are being output to ",output_dir," using phenotype files from date ",)

final_output <- data.table() # initialize empty data.table
for (ethnicity in c("ALL_AFRICAN", "BRITISH", "CARRIBEAN", "IRISH", "NONBRIT_CAUCASIAN", "OTHER_WHITE", "SOUTH_ASIAN")) {
    message("Running for ",ethnicity,"...")
    phenotype_file<-paste0(pheno_dir,"/UKB_MERGED_PHENO_",ethnicity,"_",pheno_date,"_FINAL_INCLUDE_PCs.txt.gz")
    raw_df <- fread(phenotype_file, sep="\t", header=TRUE, na.strings = c("","NA"))
    output <- data.table(SAMPLE=raw_df$eid)
    for (biomarker in legend_df$Field_ID) {
        message("Extracting ",biomarker,"...")
        regex <- paste0("^",biomarker)
        biomarker_name <- legend_df[Field_ID==biomarker, "Description"][[1]]
        extract_df <- raw_df[, .SD, .SDcols = names(raw_df) %like% paste("eid", regex, sep="|")]
        colnames(extract_df)<-c("SAMPLE", "baseline", "followup1")
        extract_df[, diff:=followup1-baseline]
        # If testosterone, estradiol, or SHBG, then perform transformations separately for each sex to account for sex differences in biomarker levels
        if (biomarker=="30800" | biomarker=="30830" | biomarker=="30850") {
            message("Stratifying transformations by sex...")
            sex_df <- raw_df[, .SD, .SDcols = names(raw_df) %like% paste("eid|^22001")]
            colnames(sex_df)<-c("SAMPLE", "sex")
            extract_df <- merge(extract_df, sex_df, by="SAMPLE")
            # Winsorize (5SD)
            extract_df[, baseline_winsorize:=winsorize(baseline), by=sex]
            extract_df[, followup1_winsorize:=winsorize(followup1), by=sex]
            # Normalize
            extract_df[, baseline_norm:=quantile_normalize(baseline), by=sex]
            extract_df[, followup1_norm:=quantile_normalize(followup1), by=sex]
            extract_df[, sex:=NULL]
        } else {
            # Winsorize (5SD)
            extract_df[, baseline_winsorize:=winsorize(baseline)]
            extract_df[, followup1_winsorize:=winsorize(followup1)]
            # Normalize
            extract_df[, baseline_norm:=quantile_normalize(baseline)]
            extract_df[, followup1_norm:=quantile_normalize(followup1)]
        }
        colnames(extract_df)[2:ncol(extract_df)]<-paste(biomarker_name,colnames(extract_df)[2:ncol(extract_df)], sep="_")
        output <- merge(output, extract_df, by="SAMPLE", all.x=TRUE)
    }
    fwrite(output, paste0(output_dir,"/UKB_",ethnicity,"_biomarkersNormalized.txt.gz"), sep="\t", row.names=FALSE, quote=FALSE, na="NA")
    output$ethnicity <- ethnicity
    final_output <- rbind(final_output, output)
}
fwrite(final_output, paste0(output_dir,"/UKB_biomarkersNormalized.txt.gz"), sep="\t", row.names=FALSE, quote=FALSE, na="NA")

print(Sys.time() - start_time)
message("************************************************\nFinished script at ", Sys.time(), "\n************************************************")