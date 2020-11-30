################################################################################
# VARIABLES
################################################################################
ukb_pheno_file <- "/storage/UKBIOBANK/UKB_PHENOS/UKB_MERGED_PHENO_BRITISH_2020_08_28_FINAL_INCLUDE_PCs.txt.gz"
ehr_dir <- "/storage/UKBIOBANK/UKB_PHENOS/common_references/Showcase_EHR/"
output_dir <- "/storage/UKBIOBANK/mohamp1/for_Sukrit/angina"
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<2) {stop("Please provide:\n1: name of phenotype\n2: path to file containing list of ICD-10 codes (one per line; no header)\n3: path to file containing list of OPCS4 codes (one per line; no header) (optional)")}
phenotype_name <- args[1] # "angina"
if (file.exists(args[3])) {
    icd_file <- args[2] # paste0("/storage/UKBIOBANK/mohamp1/for_Sukrit/angina/",phenotype_name,"_codes.txt")
    temp <- data.table::fread(icd_file, header = FALSE)
    icd_list <- as.list(temp[[1]])
} else {
    stop("No valid ICD code file path provided")
}
if (length(args)>2 && file.exists(args[3])) {
    opcs_file <- args[3]
    temp2 <- data.table::fread(opcs_file, header = FALSE)
    opcs_list <- temp2[[1]]
} else {
    warning("No valid OPCS code file path provided.")
    opcs_list <- "NA"
}
current_date <- gsub("-", "_", Sys.Date())

################################################################################
# LOAD PACKAGES
################################################################################
library(parallel)
library(data.table)
threads <- floor(detectCores()*0.125)
setDTthreads(threads=threads)

################################################################################
# INITIALIZE FUNCTIONS
################################################################################
extract_cols <- function(df, regex) {
    temp <- df[, .SD, .SDcols = names(df) %like% paste("eid", regex, sep = "|")]
    colnames(temp) <- paste0("X", colnames(temp))
    return(temp)
}
generate_df <- function(file1, file2, CAD_LIST, pheno_dir = ehr_dir) {
    # Determine index of pertinent phenotypes to extract corresponding date of event
    message("STEP: Determining index of specified phenotypes and matching to corresponding dates of each event...")
    ICD10 <- as.matrix(fread(
        file = paste0(pheno_dir, file1, ".txt.gz"),
        header = TRUE
    ))
    ICD10_DATES <- as.matrix(fread(
        file = paste0(pheno_dir, file2, ".txt.gz"),
        header = TRUE
    ))
    rownames(ICD10) <- paste0("X", ICD10[, 1])
    ICD10 <- ICD10[, -1]
    rownames(ICD10_DATES) <- paste0("X", ICD10_DATES[, 1])
    ICD10_DATES <- ICD10_DATES[, -1]
    CAD_INDEX <- lapply(
        CAD_LIST,
        function(x) which(ICD10 == x, arr.ind = TRUE))
    CAD_DFS <- lapply(
        CAD_INDEX,
        function(x) data.frame(eid = rownames(x), AFF = rep(1, nrow(x))))
    CAD_DATE_DFS <- lapply(
        CAD_INDEX,
        function(x) data.frame(eid = rownames(x), DOE = ICD10_DATES[x]))
    # Generate a dataframe of just UKB IDs
    message("STEP: Generating dataframe of unaffected IDs...")
    UKB_eid <- rownames(ICD10)
    outersect <- function(x, y) {
        sort(c(
            x[!x %in% y],
            y[!y %in% x]
        ))
    }
    NOCAD_DFS <- lapply(
        CAD_DFS,
        function(x) {
            data.frame(
                eid = outersect(x[, 1], UKB_eid),
                AFF = rep(0, length(outersect(x[, 1], UKB_eid)))
            )
        }) # create dataframe of unaffected IDs
    NOCAD_DATE_DFS <- lapply(
        NOCAD_DFS,
        function(x) {
            data.frame(eid = x[, 1], DOE = rep(NA, nrow(x)))
        })
    f1 <- function(x, y) {
        z <- rbind(x, y)
        return(z)
    }
    message("STEP: Merging unaffected and affected dataframes...")
    CAD_NOCAD_LIST <- mapply(f1, CAD_DFS, NOCAD_DFS, SIMPLIFY = FALSE)
    CAD_NOCAD_DATE_LIST <- mapply(f1,
        CAD_DATE_DFS,
        NOCAD_DATE_DFS,
        SIMPLIFY = FALSE)
    f2 <- function(x, y) {
        out <- setNames(x, c("eid", y))
        return(out)
    }
    f3 <- function(x, y) {
        out <- setNames(x, c("eid", paste0("DOE_", y)))
        return(out)
    }
    CAD_NOCAD_LIST_named <- mapply(f2,
        CAD_NOCAD_LIST,
        CAD_LIST,
        SIMPLIFY = FALSE)
    CAD_NOCAD_DATE_LIST_named <- mapply(f3,
        CAD_NOCAD_DATE_LIST,
        CAD_LIST,
        SIMPLIFY = FALSE)
    f4 <- function(x, y) {
        out <- merge(x, y, by = "eid")
        return(out)
    }
    ICD10_fin_merge <- mapply(f4,
        CAD_NOCAD_LIST_named,
        CAD_NOCAD_DATE_LIST_named,
        SIMPLIFY = FALSE)
    output <- lapply(ICD10_fin_merge, function(x) {
        x[, 1] <- as.numeric(gsub("X", "", x[, 1]))
        return(x)
    })

    return(output)
}
get_prevalent_cases <- function(dt, recruitment_df, CAD_LIST, phenotype_name) {
    message("STEP: Getting year of event...")
    ICD10_fin_merge_with_recruitment <- lapply(
        dt,
        function(x) {
            merge(recruitment_df, transform(x, DOE_year = substring(as.character(x[[3]]), 1, 4)), by = "eid")
        })
    message("STEP: Setting incident cases as date, prevalent cases as 1, and others as 0...")
    ICD10_fin_merge_with_disease <- lapply(
        ICD10_fin_merge_with_recruitment,
        function(x) {
            transform(x,
                prevalent = ifelse(as.numeric(as.character(x[[3]])) == 0, 0,
                    ifelse(as.numeric(as.character(x[[5]])) <= as.numeric(as.character(x[[2]])), 1, as.character(x[[4]]))
                )
            )
        }) # set incident cases as date of event, set cases with date before recruitment as 1, and set all others as 0

    ICD10_fin_merge_disease_sub <- lapply(
        ICD10_fin_merge_with_disease,
        function(x) x[, c(1, 6)])
    f4 <- function(x, y) {
        out <- setNames(x, c("eid", paste0("prevalent_", y)))
        return(out)
    }
    ICD10_fin_merge_disease_sub_named <- mapply(f4,
        ICD10_fin_merge_disease_sub,
        CAD_LIST,
        SIMPLIFY = FALSE)

    f5 <- function(x) {
        # If there is a patient with multiple instances of an ICD code, we need to pick the earliest one in order to cleanly merge downstream
        col_name <- names(x)[2]
        names(x)[2] <- "pheno"
        # Get all cases with prevalent disease
        out <- unique(x[pheno=="1", ])
        # Get all cases with incident diseases that don't have prevalent disease; pick the earliest date for each ID (https://stackoverflow.com/questions/24070714/extract-row-corresponding-to-minimum-value-of-a-variable-by-group)
        temp <- x[!(eid %in% out$eid) & pheno!="0", .SD[which.min(lubridate::ymd(pheno))], by = eid]
        out <- rbind(out, temp)
        # Get all remaining controls that don't have prevalent (1) or incident disease (date)
        temp <- x[!(eid %in% out$eid) & pheno=="0", ]
        out <- rbind(out, temp)
        names(out)[2] <- col_name
        return(out)
    }
    message("STEP: Getting earliest instance of event in the case of multiple instances...")
    ICD10_fin_merge_disease_sub_named_clean <- lapply(ICD10_fin_merge_disease_sub_named,
        f5
    )
    message("STEP: Merging all ICD codes...")
    prevalent_df <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "eid"), ICD10_fin_merge_disease_sub_named_clean)
    message("STEP: Combining across all ICD codes to get any individuals with 1 or more disease...")
    prevalent_df[, eval(phenotype_name) := rowSums(prevalent_df[, 2:ncol(prevalent_df)] == "1")] # if there is a case present in any column for any of the diseases, set it as 1
    prevalent_df[[ncol(prevalent_df)]][prevalent_df[[ncol(prevalent_df)]] > 1 & !is.na(prevalent_df[[ncol(prevalent_df)]])] <- 1 # set any with individuals with more than 1 prevalent case of a disease to 1

    # set all 0 values (i.e., controls) to NA
    for (i in 2:(ncol(prevalent_df) - 1)) {
        set(prevalent_df, i = which(prevalent_df[[i]] == 0), j = i, value = NA)
    }
    message("STEP: Set individuals with incident disease as date of earliest event...")
    # set individuals with incident disease as date of earliest event after recruitment
    x <- names(prevalent_df)[2:(ncol(prevalent_df) - 1)]
    incident_df <- prevalent_df[rowSums(!is.na(prevalent_df[, 2:(ncol(prevalent_df) - 1)])) > 0 & rowSums(prevalent_df[, 2:(ncol(prevalent_df) - 1)] == "1", na.rm = TRUE) == 0, ]
    incident_df[, eval(phenotype_name) := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = x]

    # select relevant columns for merging
    prevalent_df <- prevalent_df[, .SD, .SDcols = c(1, ncol(prevalent_df))]
    incident_df <- incident_df[, .SD, .SDcols = c(1, ncol(incident_df))]
    output <- merge(prevalent_df, incident_df, by = "eid", all.x = TRUE, suffixes = c("_prevalent", "_incident"))
    output[[2]][!is.na(output[[3]])] <- NA # set incident cases as NA in prevalent column
    output[[3]][output[[2]] == 0] <- 0 # set controls as 0 in incident column, but leave prevalent cases as NA

    return(output)
}
################################################################################
# FIND INSTANCES OF CODE FROM EHR DATA based on prevalent or incident event
################################################################################
setwd(ehr_dir)
# Load UKB assessment date file
start <- Sys.time()
message("Starting script at ", start)
assessment_date_df <- fread("date_assessment.txt.gz",
    sep = "\t",
    colClasses = "character",
    na.strings = c("", "NA")
)
assessment_date_df[, eid:=as.numeric(eid)]

# Get instances of cases based on your specified regular expression
# Prevalent case is coded as 0 (controls), or 1 (case)
# Incident case is coded as 0 (controls), date of event (case), or NA (if prevalent)
code_file <- "event_deathRegister"
date_file <- "date_deathRegister"
message("Extracting from death registry...")
date_df <- generate_df(code_file, date_file, icd_list)
message("Identifying prevalent cases...")
output2 <- get_prevalent_cases(date_df, assessment_date_df, icd_list, phenotype_name)
names(output2)[2:3] <- paste0(names(output2)[2:3], 2)
fwrite(output2, paste0("temp_death_",phenotype_name,".txt.gz"), na="NA", sep="\t", quote=FALSE)
print(Sys.time() - start)

code_file <- "event_cancerRegister"
date_file <- "date_cancerRegister"
message("Extracting from cancer registry...")
date_df <- generate_df(code_file, date_file, icd_list)
message("Identifying prevalent cases...")
output3 <- get_prevalent_cases(date_df, assessment_date_df, icd_list, phenotype_name)
names(output3)[2:3] <- paste0(names(output3)[2:3], 3)
fwrite(output3, paste0("temp_cancer_",phenotype_name,".txt.gz"), na="NA", sep="\t", quote=FALSE)
print(Sys.time() - start)

code_file <- "event_ICD10"
date_file <- "date_ICD10"
message("Extracting from ICD10...")
date_df <- generate_df(code_file, date_file, icd_list)
message("Identifying prevalent cases...")
output1 <- get_prevalent_cases(date_df, assessment_date_df, icd_list, phenotype_name)
names(output1)[2:3] <- paste0(names(output1)[2:3], 1)
fwrite(output1, paste0("temp_icd_",phenotype_name,".txt.gz"), na="NA", sep="\t", quote=FALSE)
print(Sys.time() - start)

if (opcs_list != "NA") {
    code_file <- "event_OPCS4"
    date_file <- "date_OPCS4"
    message("Extracting from OPCS...")
    date_df <- generate_df(code_file, date_file, opcs_list)
    message("Identifying prevalent cases...")
    output4 <- get_prevalent_cases(date_df, assessment_date_df, opcs_list, phenotype_name)
} else {
    output4 <- data.table(eid = output1$eid)
    output4[, paste0(eval(phenotype_name),"_prevalent"):=0]
    output4[, paste0(eval(phenotype_name),"_incident"):=0]
}
names(output4)[2:3] <- paste0(names(output4)[2:3], 4)
fwrite(output4, "temp_opcs.txt.gz", na="NA", sep="\t", quote=FALSE)
print(Sys.time() - start)

output1 <- fread(paste0("temp_icd_",phenotype_name,".txt.gz"), sep="\t", header=TRUE)
output2 <- fread(paste0("temp_death_",phenotype_name,".txt.gz"), sep="\t", header=TRUE)
output3 <- fread(paste0("temp_cancer_",phenotype_name,".txt.gz"), sep="\t", header=TRUE)

message("Merging prevalent cases...")
merged_output_prevalent <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "eid"), 
list(output1[, c(1, 2)], output2[, c(1, 2)], output3[, c(1, 2)], output4[, c(1, 2)])
)
merged_output_prevalent[, eval(phenotype_name) := rowSums(merged_output_prevalent[, 2:ncol(merged_output_prevalent)], na.rm=TRUE)] # exclude any NA (i.e., incident cases)
merged_output_prevalent[[ncol(merged_output_prevalent)]][merged_output_prevalent[[ncol(merged_output_prevalent)]] > 1 & !is.na(merged_output_prevalent[[ncol(merged_output_prevalent)]])] <- 1 # set any with individuals with more than 1 prevalent case of a disease to 1
merged_output_prevalent <- merged_output_prevalent[, .SD, .SDcols = c(1, ncol(merged_output_prevalent))]
fwrite(merged_output_prevalent,
    file = paste0(output_dir, "/prevalent_", phenotype_name, "_", current_date, ".txt.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)
print(Sys.time() - start)

message("Merging incident cases...")
merged_output_incident <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "eid"), list(output1[, c(1, 3)], output2[, c(1, 3)], output3[, c(1, 3)], output4[, c(1, 3)]))
# Convert incident dates from character to date format
col_name <- names(merged_output_incident)[2:ncol(merged_output_incident)]
merged_output_incident[, (col_name) := lapply(.SD, lubridate::ymd), .SDcols = col_name]
# Find minimum/earliest date in each row; remember controls and prevalent cases are NA right now after converting the columns to date format
merged_output_incident[, eval(phenotype_name) := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = col_name]
# Convert date column back to character
merged_output_incident[, (phenotype_name) := lapply(.SD, as.character), .SDcols = phenotype_name]
# Assign any remaining individuals as controls (except recall our prevalent cases are NA as well)
merged_output_incident[is.na(get(phenotype_name)), eval(phenotype_name) := 0]
# Identify individuals with prevalent cases and set as NA
prevalent_ids <- merged_output_prevalent[get(phenotype_name)==1, ][[1]]
merged_output_incident[eid %in% prevalent_ids, eval(phenotype_name) := NA]
merged_output_incident <- merged_output_incident[, .SD, .SDcols = c(1, ncol(merged_output_incident))]
fwrite(merged_output_incident,
    file = paste0(output_dir, "/incident_", phenotype_name, "_", current_date, ".txt.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

message("FINISHED!")
print(Sys.time() - start)

# Cleanup
file.remove(c(paste0("temp_death_",phenotype_name,".txt.gz"), paste0("temp_cancer_",phenotype_name,".txt.gz"), "temp_opcs.txt.gz", paste0("temp_icd_",phenotype_name,".txt.gz")))