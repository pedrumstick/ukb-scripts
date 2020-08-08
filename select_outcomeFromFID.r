#!/usr/bin/env Rscript

require(data.table)

# USER-DEFINED VARIABLES
args = commandArgs(trailingOnly=TRUE)
phenotype_name<-as.character(args[1]) # disease name
fid<-as.character(args[2]) # regex corresponding to column IDs in UK Biobank
icd_code<-as.character(args[3]) # ICD codes for disease definition (if applicable)
algorithmic_status<-as.character(args[4]) # is it an algorithmically-defined outcome?
if (any(c("Yes","YES","Y","True","TRUE","T")==algorithmic_status)) {
    algorithmic_status<-TRUE
} else {
    algorithmic_status<-FALSE
}
transform<-as.character(args[5]) # NA, log, or normalize

# FUNCTIONS
extract_cols <- function(df, regex) {
    temp <- df[, .SD, .SDcols = names(df) %like% paste("eid",regex, sep="|")]
    # temp <- temp[, lapply(.SD, function(x) gsub("^$|^ $", NA, x))] # set any empty cells as NA
    colnames(temp)<-paste0("X", colnames(temp)) # R doesn't like numeric column names
    return(temp)
}
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

# LOAD DATA
raw_df <- fread("/hightide/UKB_PHENOS/UKB_MERGED_PHENO_BRITISH_2020_01_28_FINAL_INCLUDE_PCs.txt", sep="\t", header=TRUE, na.strings = c("","NA"))

# EXTRACT RELEVANT COLUMNS BASED ON FID
message("Extracting ",phenotype_name," using FID columns matching ", fid)
extract_df<-extract_cols(raw_df, fid)

if (ncol(extract_df)<2) {
    stop("ERROR: The specified column ",fid," coding for ", phenotype_name," does not exist.")
}

############################################
# Format Phenotype
############################################
if (algorithmic_status) {
    # If Algorithmic Outcome
    output<-extract_df
    output[!is.na(output[[2]]),2]<-1
    output[is.na(output[[2]]),2]<-0
    phenotype_name<-paste(phenotype_name,"algo",sep="_")
} else if (ncol(extract_df)>2 & (!is.na(icd_code) | icd_code=="NA")) {
    # If ICD code
    output<-data.table(extract_df[, "Xeid"], case_status=Reduce(`|`, lapply(extract_df, `%like%`, icd_code))) # https://stackoverflow.com/questions/45827337/in-r-check-if-string-appears-in-row-of-dataframe-in-any-column
    output[, case_status:=as.numeric(case_status)]
} else {
    output<-extract_df
}

# Rename column names and check if there is unexpected output
if (ncol(output)==2) {
    names(output)[1]<-"SAMPLE"
    names(output)[2]<-phenotype_name
} else {
    stop("There are NOT 2 columns in the phenotype file. Please check the headers of the output file.")
}

############################################
# Data Transformation (if excluding any samples, please consider whether this should be done before transforming data)
############################################
if (length(unique(output[[2]])) > 3 ) {
    output[, 2] <- winsorize(output[[2]]) # do I winsorize, then transform? or vice-versa?
    if (transform == "log") {
        output[,2] <- log(output[[2]]-(min(output[[2]], na.rm=TRUE) - 1)) # add constant to avoid problems log-transforming 0 or negative numbers
        phenotype_name<-paste(phenotype_name,"log",sep="_")
    } else if (transform == "normalize") {
        output[,2] <- quantile_normalize(output[[2]])
        phenotype_name<-paste(phenotype_name,"norm",sep="_")
    }
    names(output)[2]<-phenotype_name
}

############################################
# Output
############################################
output_dir<-paste0("/widetide/mohamp1/analysis/BOLT_GWAS/OUTCOMES/",phenotype_name,"/")
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}
fwrite(output, paste0(output_dir,"/UKB_pheno.txt"), sep=" ", row.names=FALSE, quote=FALSE, na="NA")