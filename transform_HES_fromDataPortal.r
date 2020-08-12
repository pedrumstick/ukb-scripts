#!/usr/bin/env Rscript

##############################
# VARIABLES AND CHECKS
##############################
input_file_diag <- "/hightide/UKB_PHENOS/Data_Portal_Downloads/Inpatient_Hospital/hesin_diag.txt" # diagnosis
if (!file.exists(input_file_diag)) {
  stop(paste0("ERROR: There is no specified input file ", input_file_diag))
}
input_file_dates <- "/hightide/UKB_PHENOS/Data_Portal_Downloads/Inpatient_Hospital/hesin.txt" # dates
if (!file.exists(input_file_dates)) {
  stop(paste0("ERROR: There is no specified input file with dates ", input_file_dates))
}
input_file_oper <- "/hightide/UKB_PHENOS/Data_Portal_Downloads/Inpatient_Hospital/hesin_oper.txt" # operation
if (!file.exists(input_file_oper)) {
  stop(paste0("ERROR: There is no specified input file ", input_file_oper))
}
output_dir <- "/hightide/UKB_PHENOS/Data_Portal_Downloads/Inpatient_Hospital/"
if (!dir.exists(output_dir)) {
  stop(paste0("ERROR: There is no specified output directory ", output_dir))
}

##############################
# LOAD PACKAGES
##############################
require(data.table)
require(tidyr)
require(dplyr)
require(lubridate)

##############################
# LOAD DATA FRAMES
##############################
# http://biobank.ndph.ox.ac.uk/showcase/showcase/docs/HospitalEpisodeStatistics.pdf
setwd(output_dir)
start <- Sys.time()
# DATES
# Based on the UKB data dictionary, the date field is taken as the value of the epistart field for the corresponding record, or where this is unknown the admidate field
# eid and ins_index uniquely identifies this record, i.e. eid & ins_index together form a primary key for this table. But the ins_index will not generally index records in chronological order
df <- fread(input_file_dates, header = TRUE, na.strings = c("", "NA", NA)) # load dates
# Based on the UKB data dictionary, the date field is taken as the value of the epistart field for the corresponding record, or where this is unknown the admidate field
df_dates <- df[, .SD, .SDcols = c("eid", "ins_index", "epistart", "admidate")]
df_dates[is.na(epistart), epistart := admidate] # for instances without an episode start date, we replace with the admission date
df_dates[, admidate := NULL] # remove admission date column
df_dates[, epistart := dmy(epistart)]
df_dates[, epistart := as.character(epistart)]

# DIAGNOSTIC CODES
df2 <- fread(input_file_diag, header = TRUE, na.strings = c("", "NA", NA)) # diagnosis codes
# 12,627,730 total rows (i.e., some type of diagnosis/case)

# OPERATIONS CODES
# Based on the UKB data dictionary, the opdate field in the HESIN_OPER table has proved unreliable in a number of cases, and it should not be preferred over the epistart or admidate fields in the DATES table.
df3 <- fread(input_file_oper, header = TRUE, na.strings = c("", "NA", NA)) # operation codes
# 7,169,074 total rows (i.e., some type of diagnosis/case)

##############################
# ICD-9 (41271)
##############################
icd9_df <- df2[!is.na(diag_icd9), .SD, .SDcols = c("eid", "ins_index", "arr_index", "diag_icd9")]
setorderv(icd9_df, c("ins_index", "arr_index"))
# Merge with Dates (41281)
icd9_df <- setDT(left_join(icd9_df, df_dates, by = c("eid", "ins_index")))
# Transform to wide data
icd9_diagnoses <- dcast(icd9_df, eid ~ ins_index + arr_index, value.var = "diag_icd9", sep = ".")
icd9_dates <- dcast(icd9_df, eid ~ ins_index + arr_index, value.var = "epistart", sep = ".")
# Shift non-NA cells to the left of the dataframe
# https://stackoverflow.com/questions/49079789/using-r-to-shift-values-to-the-left-of-data-frame
icd9_diagnoses_clean <- setDT(as.data.frame(t(apply(icd9_diagnoses, 1, function(x) {
  c(x[!is.na(x)], x[is.na(x)])
}))))
icd9_dates_clean <- setDT(as.data.frame(t(apply(icd9_dates, 1, function(x) {
  c(x[!is.na(x)], x[is.na(x)])
}))))
# Delete columns with all NA values
icd9_diagnoses_clean <- icd9_diagnoses_clean[, which(unlist(lapply(icd9_diagnoses_clean, function(x) !all(is.na(x))))), with = F]
icd9_dates_clean <- icd9_dates_clean[, which(unlist(lapply(icd9_dates_clean, function(x) !all(is.na(x))))), with = F]
# Check
if (!identical(dim(icd9_diagnoses_clean), dim(icd9_dates_clean))) {
  stop("ERROR: Dimensions of date and diagnosis codes are different!")
}
# Output
names(icd9_diagnoses_clean) <- c("eid", paste0("41271-1.", as.character(seq(2, ncol(icd9_diagnoses_clean)) - 2)))
names(icd9_dates_clean) <- c("eid", paste0("41281-1.", as.character(seq(2, ncol(icd9_dates_clean)) - 2)))
write.table(icd9_diagnoses_clean, "temp_icd9_diag.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(icd9_dates_clean, "temp_icd9_dates.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##############################
# OPCS-3 (41273)
##############################
opcs3_df <- df3[!is.na(oper3), .SD, .SDcols = c("eid", "ins_index", "arr_index", "oper3")]
setorderv(opcs3_df, c("ins_index", "arr_index"))
# Merge with Dates (41283)
opcs3_df <- setDT(left_join(opcs3_df, df_dates, by = c("eid", "ins_index")))
# Transform to wide data
opcs3_diagnoses <- dcast(opcs3_df, eid ~ ins_index + arr_index, value.var = "oper3", sep = ".")
opcs3_dates <- dcast(opcs3_df, eid ~ ins_index + arr_index, value.var = "epistart", sep = ".")
names(opcs3_diagnoses)[2:ncol(opcs3_diagnoses)] <- paste0("V", as.character(names(opcs3_diagnoses)[2:ncol(opcs3_diagnoses)]))
names(opcs3_dates)[2:ncol(opcs3_dates)] <- paste0("V", as.character(names(opcs3_dates)[2:ncol(opcs3_dates)]))
# Shift non-NA cells to the left of the dataframe
# https://stackoverflow.com/questions/49079789/using-r-to-shift-values-to-the-left-of-data-frame
opcs3_diagnoses_clean <- setDT(as.data.frame(t(apply(opcs3_diagnoses, 1, function(x) {
  c(x[!is.na(x)], x[is.na(x)])
}))))
opcs3_dates_clean <- setDT(as.data.frame(t(apply(opcs3_dates, 1, function(x) {
  c(x[!is.na(x)], x[is.na(x)])
}))))
# Delete columns with all NA values
opcs3_diagnoses_clean <- opcs3_diagnoses_clean[, which(unlist(lapply(opcs3_diagnoses_clean, function(x) !all(is.na(x))))), with = F]
opcs3_dates_clean <- opcs3_dates_clean[, which(unlist(lapply(opcs3_dates_clean, function(x) !all(is.na(x))))), with = F]
# Check
if (!identical(dim(opcs3_diagnoses_clean), dim(opcs3_dates_clean))) {
  stop("ERROR: Dimensions of date and diagnosis codes are different!")
}
# Output
names(opcs3_diagnoses_clean) <- c("eid", paste0("41273-1.", as.character(seq(2, ncol(opcs3_diagnoses_clean)) - 2)))
names(opcs3_dates_clean) <- c("eid", paste0("41283-1.", as.character(seq(2, ncol(opcs3_dates_clean)) - 2)))
write.table(opcs3_diagnoses_clean, "temp_opcs3_diag.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(opcs3_dates_clean, "temp_opcs3_dates.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##############################
# OPCS-4 (41272)
##############################
opcs4_df <- df3[!is.na(oper4), .SD, .SDcols = c("eid", "ins_index", "arr_index", "oper4")]
setorderv(opcs4_df, c("ins_index", "arr_index"))
# Merge with Dates (41282)
opcs4_df <- setDT(left_join(opcs4_df, df_dates, by = c("eid", "ins_index")))
# Transform to wide data
# test <- opcs4_df[sample(x=1:nrow(opcs4_df), size=10000, replace=FALSE), ]
temp_split <- split(opcs4_df, opcs4_df$eid) # dataframe is too big; need to breakup into chunks based on ID
opcs4_diagnoses <- lapply(temp_split, function(x) {
  message("Transforming codes from sample ID: ", x$eid, " at ", Sys.time())
  dcast.data.table(x, formula = eid ~ ins_index + arr_index, value.var = "oper4", sep = ".")
})
opcs4_dates <- lapply(temp_split, function(x) {
  message("Transforming dates from sample ID: ", x$eid, " at ", Sys.time())
  dcast.data.table(x, formula = eid ~ ins_index + arr_index, value.var = "epistart", sep = ".")
})
# First Output
max_cols <- as.data.frame(t(paste0("V", seq_len(max(sapply(opcs4_dates, ncol)))))) # get maximum number of columns from list of data.tables
write.table(max_cols, "temp_opcs4_diag.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(max_cols, "temp_opcs4_dates.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
lapply(opcs4_diagnoses, function(x) {
  write.table(x, "temp_opcs4_diag.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
})
lapply(opcs4_dates, function(x) {
  write.table(x, "temp_opcs4_dates.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
})
# Read data back in to assign column names, and fill in any columns with "NA"
opcs4_diagnoses <- fread("temp_opcs4_diag.txt", fill = TRUE, na.strings = c("", "NA"))
opcs4_dates <- fread("temp_opcs4_dates.txt", fill = TRUE, na.strings = c("", "NA"))
# Output
names(opcs4_diagnoses) <- c("eid", paste0("41272-1.", as.character(seq(2, ncol(opcs4_diagnoses)) - 2)))
names(opcs4_dates) <- c("eid", paste0("41282-1.", as.character(seq(2, ncol(opcs4_dates)) - 2)))
fwrite(opcs4_diagnoses, "temp_opcs4_diag.txt", sep = "\t", quote = FALSE, na = "NA")
fwrite(opcs4_dates, "temp_opcs4_dates.txt", sep = "\t", quote = FALSE, na = "NA")

##############################
# ICD-10 (41270)
##############################
icd10_df <- df2[!is.na(diag_icd10), .SD, .SDcols = c("eid", "ins_index", "arr_index", "diag_icd10")]
setorderv(icd10_df, c("ins_index", "arr_index"))
# Merge with Dates (41280)
icd10_df <- setDT(left_join(icd10_df, df_dates, by = c("eid", "ins_index")))
# Transform to wide data
temp_split <- split(icd10_df, icd10_df$eid) # dataframe is too big; need to breakup into chunks based on ID
icd10_diagnoses <- lapply(temp_split, function(x) {
  message("Transforming codes from sample ID: ", x$eid, " at ", Sys.time())
  dcast.data.table(x, formula = eid ~ ins_index + arr_index, value.var = "diag_icd10", sep = ".")
})
icd10_dates <- lapply(temp_split, function(x) {
  message("Transforming dates from sample ID: ", x$eid, " at ", Sys.time())
  dcast.data.table(x, formula = eid ~ ins_index + arr_index, value.var = "epistart", sep = ".")
})
# First Output
max_cols <- as.data.frame(t(paste0("V", seq_len(max(sapply(icd10_dates, ncol)))))) # get maximum number of columns from list of data.tables
write.table(max_cols, "temp_icd10_diag.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(max_cols, "temp_icd10_dates.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
lapply(icd10_diagnoses, function(x) {
  write.table(x, "temp_icd10_diag.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
})
lapply(icd10_dates, function(x) {
  write.table(x, "temp_icd10_dates.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
})
# Read data back in to assign column names, and fill in any columns with "NA"
icd10_diagnoses <- fread("temp_icd10_diag.txt", fill = TRUE, na.strings = c("", "NA"))
icd10_dates <- fread("temp_icd10_dates.txt", fill = TRUE, na.strings = c("", "NA"))
# Output
names(icd10_diagnoses) <- c("eid", paste0("41270-1.", as.character(seq(2, ncol(icd10_diagnoses)) - 2)))
names(icd10_dates) <- c("eid", paste0("41280-1.", as.character(seq(2, ncol(icd10_dates)) - 2)))
fwrite(icd10_diagnoses, "temp_icd10_diag.txt", sep = "\t", quote = FALSE, na = "NA")
fwrite(icd10_dates, "temp_icd10_dates.txt", sep = "\t", quote = FALSE, na = "NA")

##############################
# MERGE ALL
##############################
icd10_diagnoses <- fread("temp_icd10_diag.txt", fill = TRUE, na.strings = c("", "NA"))
icd9_diagnoses <- fread("temp_icd9_diag.txt", fill = TRUE, na.strings = c("", "NA"))
opcs4_diagnoses <- fread("temp_opcs4_diag.txt", fill = TRUE, na.strings = c("", "NA"))
opcs3_diagnoses <- fread("temp_opcs3_diag.txt", fill = TRUE, na.strings = c("", "NA"))
icd10_dates <- fread("temp_icd10_dates.txt", fill = TRUE, na.strings = c("", "NA"))
icd9_dates <- fread("temp_icd9_dates.txt", fill = TRUE, na.strings = c("", "NA"))
opcs4_dates <- fread("temp_opcs4_dates.txt", fill = TRUE, na.strings = c("", "NA"))
opcs3_dates <- fread("temp_opcs3_dates.txt", fill = TRUE, na.strings = c("", "NA"))
output <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "eid"), list(icd10_diagnoses, icd9_diagnoses, opcs4_diagnoses, opcs3_diagnoses, icd10_dates, icd9_dates, opcs4_dates, opcs3_dates))
fwrite(output, paste0("hesin_wide_", gsub("-", "_", Sys.Date()), ".txt.gz"), sep = "\t", quote = FALSE, na = "NA")

message("FINISHED!")
print(Sys.time() - start)