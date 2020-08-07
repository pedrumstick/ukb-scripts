#!/usr/bin/env Rscript

##############################
# VARIABLES AND CHECKS
##############################
input_file<-"/hightide/UKB_PHENOS/Data_Portal_Downloads/COVID19/covid19_result.txt"
if (!file.exists(input_file)) {stop(paste0("ERROR: There is no specified input file ", input_file))}

##############################
# LOAD PACKAGES
##############################
require(reshape2)
require(lubridate)

##############################
# LOAD DATA FRAMES
##############################
df<-read.table(input_file, header=TRUE)
df$specdate <- dmy(df$specdate)
df<-df[order(df$specdate),]
df$specdate<-as.character(df$specdate) # convert back to character due to downstream issues with dates in dcast (https://stackoverflow.com/questions/12289731/posixct-values-become-numeric-in-reshape2-dcast)
df$count<-with(df, ave(eid, eid, FUN = seq_along)) # https://stackoverflow.com/questions/38335099/counting-duplicates-in-r

# Convert to wide form data for merging with broader UK Biobank phenotype files
wide_df<-unique(df[,"eid", drop=FALSE])
for (i in names(df)[2:(ncol(df)-1)]) {
        print(i)
        temp_wide<-dcast(df, eid ~ count, value.var = i)
        names(temp_wide)[2:ncol(temp_wide)]<-paste0(i,"-",names(temp_wide)[2:ncol(temp_wide)])
        wide_df<-merge(wide_df, temp_wide, by="eid", all.x=TRUE)
}

if (length(unique(df$eid))!=nrow(wide_df)) {stop("Incorrect number of rows after transforming to wide format.")}
if (max(table(df$eid))!=length(grep("specdate", names(wide_df)))) {stop("Incorrect number of columns after transforming to wide format.")}
write.table(wide_df, paste0("/hightide/UKB_PHENOS/Data_Portal_Downloads/COVID19/covid19_result_wide_",gsub("-", "_", Sys.Date()),".txt"), sep="\t", row.names=FALSE, quote=FALSE)

df<-df[,-ncol(df)]
df<-df[order(df$eid, df$specdate),]
write.table(df, paste0("/hightide/UKB_PHENOS/Data_Portal_Downloads/COVID19/covid19_result_long_",gsub("-", "_", Sys.Date()),".txt"), sep="\t", row.names=FALSE, quote=FALSE)