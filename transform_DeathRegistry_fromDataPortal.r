#!/usr/bin/env Rscript

##############################
# VARIABLES AND CHECKS
##############################
input_file<-"/hightide/UKB_PHENOS/Data_Portal_Downloads/Death_Register/death_cause.txt" # cause of death
if (!file.exists(input_file)) {stop(paste0("ERROR: There is no specified input file ", input_file))}
input_file_dates<-"/hightide/UKB_PHENOS/Data_Portal_Downloads/Death_Register/death.txt" # dates of death
if (!file.exists(input_file_dates)) {stop(paste0("ERROR: There is no specified input file with dates ", input_file_dates))}

##############################
# LOAD PACKAGES
##############################
require(reshape2)
require(data.table)

##############################
# LOAD DATA FRAMES AND CONVERT TO "WIDE" FORMAT
##############################
# Primary Causes of Death
df<-fread(input_file, header=TRUE)
primary_df<-df[level==1, ] # extract only primary causes of death
primary_df[order(primary_df$eid), ]
primary_df <- dcast(primary_df, eid ~ ins_index, value.var="cause_icd10")
names(primary_df)[2:ncol(primary_df)] <- paste0("40001-",as.numeric(names(primary_df)[2:ncol(primary_df)]),".0")

# Dates of Death
df2<-fread(input_file_dates, header=TRUE)
df_dates <- dcast(df2, eid ~ ins_index, value.var="date_of_death")
names(df_dates)[2:ncol(df_dates)] <- paste0("40000-",as.numeric(names(df_dates)[2:ncol(df_dates)]),".0")

# Secondary Causes of Death
secondary_df<-df[level==2, ]
secondary_df[, arr_index:=arr_index-1]
secondary_df<-dcast(secondary_df, eid ~ ins_index + arr_index, value.var="cause_icd10", sep=".")
names(secondary_df)[2:ncol(secondary_df)] <- paste0("40002-", as.character(names(secondary_df)[2:ncol(secondary_df)]))

# Merge everything into one dataset with dates
merged_df <- Reduce(function(x, y) merge(x, y, by="eid", all=TRUE), list(df_dates, primary_df, secondary_df))

##############################
# OUTPUT
##############################
if (length(unique(c(df$eid, df2$eid)))!=nrow(merged_df)) {stop("Incorrect number of rows after transforming to wide format.")}
write.table(merged_df, paste0("/hightide/UKB_PHENOS/Data_Portal_Downloads/Death_Register/deathRegistry_wide_",gsub("-", "_", Sys.Date()),".txt"), sep="\t", row.names=FALSE, quote=FALSE)