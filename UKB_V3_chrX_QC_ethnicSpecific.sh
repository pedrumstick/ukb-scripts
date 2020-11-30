#!/bin/bash

##################### QUALITY CONTROL FOR X CHROMOSOME IN UK BIOBANK BRITISH INDIVIDUALS #####################
chr=X # change to XY for pseudoautosomal regions that can be QC'd similar to autosomes
ethnicity=$1
IMP_DIR="/londonbridge/UKB_TRANSFER_2019_12_13/xfer.crg.eu"
bgenix="/hightide/bgen_tools/bgen/gavinband-bgen-798eca81c0fa/final_bgenie_program/bin/bgenix"
cat_bgen="/hightide/bgen_tools/bgen/gavinband-bgen-798eca81c0fa/final_bgenie_program/bin/cat-bgen"
# qctool="/hightide/chongm/programs/qctool_v2.0-rc2-linux-x86_64/qctool"
qctool="/hightide/software/qctool_v2.0.8-CentOS_Linux7.6.1810-x86_64/qctool"
analysis_dir="/londonbridge/UKB_FULL_v3_${ethnicity}_QC"
plink2="/hightide/software/plink_2/plink2"
OUTPUT_BGEN="/widetide/UKB_FULL_v3_${ethnicity}_X/" # Decide where to output data
FINAL_OUTPUT="/widetide/UKB_FULL_v3_${ethnicity}_8BIT_forBOLT/"
OUTPUT_PLINK="/widetide/UKB_FULL_v3_${ethnicity}_QC_PLINK/"

mkdir ${OUTPUT_BGEN}
mkdir ${FINAL_OUTPUT}

echo "************************************************************"
echo "Started QC for ${ethnicity} at: "`date +"%m/%d/%Y %H:%M:%S"`
echo "************************************************************"

# CREATE PROPERLY SYNTAXED BGI FILES
# A) INFO > 0.30 (Using UK biobank field)
cd $analysis_dir
if [[ ! -f ${analysis_dir}/chr.${chr}.UKB_v3_INFO_0.30.bgen.bgi ]]
then
awk 'BEGIN{FS="\t";OFS="\t"} $8 > 0.30 && ($6>0.0001 && $6 < 0.999){print $2}' ${IMP_DIR}/ukb_mfi_chr${chr}_v3.txt > chr${chr}.v3_INFO_0.30_VARIANT_IDS.txt
$bgenix -g ${IMP_DIR}/ukb_imp_chr${chr}_v3.bgen \
-incl-rsids chr${chr}.v3_INFO_0.30_VARIANT_IDS.txt \
> chr.${chr}.UKB_v3_INFO_0.30.bgen
$bgenix -g chr.${chr}.UKB_v3_INFO_0.30.bgen -index -clobber
fi

####################################################################################
# 0 - Create related and unrelated BGEN files
####################################################################################
# https://stackoverflow.com/questions/40377623/bash-wait-command-waiting-for-more-than-1-pid-to-finish-execution
cd $analysis_dir
SECONDS=0

${qctool} -g chr.${chr}.UKB_v3_INFO_0.30.bgen \
-s ${IMP_DIR}/ukb15255_imp_chr${chr}_v3_s486663.sample \
-incl-samples ${analysis_dir}/UKB_FULL_${ethnicity}_POTENTIALLY_RELATED_QC_PASS_IDS.txt \
-og chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.bgen \
-os chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.sample &
pids+=($!)
${qctool} -g chr.${chr}.UKB_v3_INFO_0.30.bgen \
-s ${IMP_DIR}/ukb15255_imp_chr${chr}_v3_s486663.sample \
-incl-samples ${analysis_dir}/UKB_FULL_${ethnicity}_UNRELATED_QC_PASS_IDS.txt \
-og chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.bgen \
-os chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.sample &
pids+=($!)

wait "${pids[@]}"

$bgenix -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.bgen -index -clobber &
pids=($!)
$bgenix -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.bgen -index -clobber &
pids+=($!)

wait "${pids[@]}"

duration=$SECONDS
echo "Creating related and unrelated BGEN files took: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"

####################################################################################
# 1 - Separate into males and females
####################################################################################
cd $analysis_dir

SECONDS=0
for sex in Male Female; do
echo "Creating "$sex"-specific files..."
${qctool} -g chr.${chr}.UKB_v3_INFO_0.30.bgen \
-s ${IMP_DIR}/ukb15255_imp_chr${chr}_v3_s486663.sample \
-incl-samples ${analysis_dir}/UKB_FULL_${ethnicity}_POTENTIALLY_RELATED_QC_PASS_IDS.txt \
-excl-samples /hightide/mohamp1/UKB_IDs/ids_exclude_for${sex}.txt \
-og chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED_${sex}.bgen \
-os chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED_${sex}.sample &
${qctool} -g chr.${chr}.UKB_v3_INFO_0.30.bgen \
-s ${IMP_DIR}/ukb15255_imp_chr${chr}_v3_s486663.sample \
-incl-samples ${analysis_dir}/UKB_FULL_${ethnicity}_UNRELATED_QC_PASS_IDS.txt \
-excl-samples /hightide/mohamp1/UKB_IDs/ids_exclude_for${sex}.txt \
-og chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.bgen \
-os chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.sample &
done
wait

echo "$(wc -l chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED_Male.sample)"
echo "$(wc -l chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED_Female.sample)"
echo "$(wc -l chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_Male.sample)"
echo "$(wc -l chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_Female.sample)"

for sex in Male Female; do
for subset in UNRELATED RELATED; do
if [ -f chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}_${sex}.bgen ]; then
$bgenix -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}_${sex}.bgen -index -clobber
echo "Total elapsed time for BGEN indexing: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
else
echo "ERROR: The specified file does not exist!"
fi
done
done
wait

duration=$SECONDS
echo "Separating sexes took: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"

####################################################################################
# 2 - Exclude variants based on MAF  > 0.001 in males and MAF > 0.001 + HWE > 1x10-10 in females within the unrelated dataset
# Principles of HWE don't apply to males for the X chromosome
####################################################################################
cd $analysis_dir

SECONDS=0
# In males
sex=Male
${qctool} -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.bgen \
-s chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.sample \
-snp-stats -osnp chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.SNP_STATS.txt
# Find MAF passing SNPs
grep -v "#" chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.SNP_STATS.txt | awk 'NR>1{if ($14>0.001) print $2}' > chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt
# In females
sex=Female
${qctool} -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.bgen \
-s chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.sample \
-snp-stats -osnp chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.SNP_STATS.txt
# Find MAF and HWE passing SNPs
grep -v "#" chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.SNP_STATS.txt | awk 'NR>1{if ($8>1e-10 && $14>0.001) print $2}' > chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_${sex}.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt
# Merge male and female variants to include only variants passing QC in both
awk 'NR==FNR {a[$1];next;} $1 in a {print}' chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_Male.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED_Female.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt > chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt

# IN THE MERGED SEX FILE, INCLUDE ONLY SNPs THAT ARE PASSING QC IN BOTH SEXES
$bgenix -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.bgen -incl-rsids chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt > ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.MAF_0.001.HWE.FINAL.bgen &
pids+=($!)
$bgenix -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.bgen -incl-rsids chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt > ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.MAF_0.001.HWE.FINAL.bgen &
pids+=($!)

wait "${pids[@]}"

cp chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.sample ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.MAF_0.001.HWE.FINAL.sample 
cp chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.sample ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.MAF_0.001.HWE.FINAL.sample 

duration=$SECONDS
echo "FINISHED! Variant exclusion took: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"

####################################################################################
# 3 - Recode males as diploid (0/2)
####################################################################################
# If you wish to analyze X chromosome data, you will need to create a BGEN v1.2 file in which all genotypes are coded as diploid. By default, plink2 will code males as haploid, but you can force it to create diploid X chromosome data by setting the sex of all individuals to female before converting.
# Starting with v2.3.2, BOLT-LMM accepts X chromosome genotypes for both model-fitting (via --bfile or --bed/bim/fam PLINK-format input) and association testing on imputed variants (e.g., in BGEN files). Males should be coded as diploid (as PLINK does for chromosome code 23 = X non-PAR), such that male genotypes are coded as 0/2 and female genotypes are coded as 0/1/2 (corresponding to a model of random X inactivation) 

cd ${OUTPUT_BGEN}

SECONDS=0
for subset in UNRELATED RELATED; do
awk 'NR<3 {$1=$1; print} NR>=3 {$4=2; print}' ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.sample > ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL_recodedMales.sample # code all samples as female to ensure PLINK codes males as diploid (0/2)
$plink2 --bgen ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.bgen \
--sample ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL_recodedMales.sample \
--export bgen-1.2 bits=8 \
--out ${OUTPUT_BGEN}/temp

mv ${OUTPUT_BGEN}/temp.bgen ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.bgen
mv ${OUTPUT_BGEN}/temp.sample ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.sample

# CONVERT FINAL FILES TO PLINK FORMAT 
$plink2 --bgen ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.bgen \
--sample ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.sample \
--make-bed \
--out ${OUTPUT_PLINK}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL

# Add haploid BGEN files to final output folder
mv ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.bgen ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.HaploidMales.bgen
mv ${OUTPUT_BGEN}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.sample ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.HaploidMales.sample

duration=$SECONDS
echo "FINISHED! Recoding males took: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"

####################################################################################
# 4 - Final BGEN indexing
####################################################################################
cd ${FINAL_OUTPUT}
SECONDS=0
for subset in UNRELATED RELATED; do
$bgenix -g ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.bgen -index -clobber &
pids+=($!)
done
wait "${pids[@]}"
duration=$SECONDS
echo "FINISHED! Elapsed time for final BGEN indexing: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"

####################################################################################
# PSEUDOAUTOSOMAL REGION?
####################################################################################
chr=XY # change to XY for pseudoautosomal regions that can be QC'd similar to autosomes
################### UK BIOBANK FULL 500K (v3-IMPUTED SNPS ONLY) QC ######################
############# VARIANT QC #############
# Create properly syntaxed bgi files 
# for chr in {1..22} ; do
# cp ${IMP_DIR}/ukb_bgi_chr${chr}_v2.bgi ${IMP_DIR}/ukb_imp_chr${chr}_v2.bgen.bgi
# done
# A) INFO > 0.30 (Using UK biobank field)
cd $analysis_dir
# Also exclude ultra rare variants - eliminates about 1/3 of variants
SECONDS=0
awk 'BEGIN{FS="\t";OFS="\t"} $8 > 0.30 && ($6>0.0001 && $6 < 0.999){print $2}' ${IMP_DIR}/ukb_mfi_chr${chr}_v3.txt > chr${chr}.v3_INFO_0.30_VARIANT_IDS.txt
$bgenix -g ${IMP_DIR}/ukb_imp_chr${chr}_v3.bgen \
-incl-rsids chr${chr}.v3_INFO_0.30_VARIANT_IDS.txt \
> chr.${chr}.UKB_v3_INFO_0.30.bgen
$bgenix -g chr.${chr}.UKB_v3_INFO_0.30.bgen -index -clobber

duration=$SECONDS
echo "FINISHED! Elapsed time for creating proper BGEN chromosome file: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"

# C) Sample Exclusions
# British w/ potentially related
SECONDS=0
${qctool} -g chr.${chr}.UKB_v3_INFO_0.30.bgen \
-s ${IMP_DIR}/ukb15255_imp_chr${chr}_v3_s486349.sample \
-incl-samples ${analysis_dir}/UKB_FULL_${ethnicity}_POTENTIALLY_RELATED_QC_PASS_IDS.txt \
-og chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.bgen \
-os chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.sample &

# British unrelated 
${qctool} -g chr.${chr}.UKB_v3_INFO_0.30.bgen \
-s ${IMP_DIR}/ukb15255_imp_chr${chr}_v3_s486349.sample \
-incl-samples ${analysis_dir}/UKB_FULL_${ethnicity}_UNRELATED_QC_PASS_IDS.txt \
-og chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.bgen \
-os chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.sample &
wait

duration=$SECONDS
echo "FINISHED! Elapsed time for creating related and unrelated subsets: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"

# D) MAF > 0.001 (Re-calculate within unrelated dataset)
# E) HWE > 1 x 10-10 (Recalculate within unrelated dataset only)
SECONDS=0
${qctool} -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.bgen \
-s chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.sample \
-snp-stats -osnp chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.SNP_STATS.txt
# Find MAF and HWE passing SNPs
grep -v "#" chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.SNP_STATS.txt | awk 'NR>1{if ($8>1e-10 && $14>0.001) print $2}' > chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt

# Variant Subset for both Related and Unrelated Subsets
# Unrelated
$bgenix -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.bgen -index -clobber &
$bgenix -g chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.bgen -index -clobber &
wait
duration=$SECONDS
echo "FINISHED! Elapsed time for finding variants that failed QC: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"

SECONDS=0
$plink2 --bgen chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.bgen \
--sample chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.sample \
--extract chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt \
--export bgen-1.2 bits=8 \
--out ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.MAF_0.001.HWE.FINAL

$plink2 --bgen chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.bgen \
--sample chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.sample \
--extract chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.MAF_0.001.HWE.SNPS_TO_INCLUDE.txt \
--export bgen-1.2 bits=8 \
--out ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.MAF_0.001.HWE.FINAL

$bgenix -g ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.UNRELATED.bgen -index -clobber &
$bgenix -g ${FINAL_OUTPUT}/chr.${chr}.UKB_v3_INFO_0.30.${ethnicity}.RELATED.bgen -index -clobber &
wait
duration=$SECONDS
echo "FINISHED! Elapsed time for variant-level QC: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"

################################################################
# MERGE PAR (PSEUDOAUTOSOMAL REGION) INTO X CHROMOSOME
################################################################
SECONDS=0
# Imputed X chromosome SNPs can also be included in BOLT-LMM association tests; again, males should be coded as diploid in one of the currently-supported formats (e.g., BGEN v1.1 or 8-bit BGEN v1.2). (BGEN v1.2 includes a data format that natively encodes a mixture of haploid and diploid SNPs, but BOLT-LMM currently does not support this format.) Chromosomes named 23, X, XY, PAR1, and PAR2 are all acceptable.
# There is no need to separate chrX into PAR and non-PAR; for PLINK input, you should simply merge PAR and non-PAR SNPs into a single “chromosome 23” using PLINK --merge-x
cd $FINAL_OUTPUT

for subset in UNRELATED RELATED; do
# Find samples common to both chromosome X and XY (PAR) for merging
grep -F -f <(awk '{print $1,$2}' chr.XY.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.sample) <(awk '{print $1, $2}' chr.X.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.sample) | awk 'NR>2 {print}' > temp_common_chrX_samples.txt

# Change chromosome XY code to chromosome X
$plink2 --bgen chr.XY.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.bgen \
--sample chr.XY.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.sample \
--merge-x \
--export bgen-1.2 bits=8 \
--out temp_chrXY_forMerge_${subset}

# Exclude excess individuals from chromosome X that aren't in chromosome XY
$plink2 --bgen chr.X.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.bgen \
--sample chr.X.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.sample \
--keep temp_common_chrX_samples.txt \
--export bgen-1.2 bits=8 \
--out temp_chrX_forMerge_${subset}

# Concatenate bgen files together (https://bitbucket.org/gavinband/bgen/wiki/cat-bgen)
# assumes the bgen files all have the same header information (i.e. they are encoded with the same version of the BGEN format, they have the same number of samples, etc.)
$cat_bgen -g temp_chrX_forMerge_${subset}.bgen temp_chrXY_forMerge_${subset}.bgen -og temp_chr23_${subset}.bgen
$plink2 --bgen temp_chr23_${subset}.bgen \
--sample temp_chrX_forMerge_${subset}.sample  \
--export bgen-1.2 bits=8 \
--out chr.23.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL

$bgenix -g ${FINAL_OUTPUT}/chr.23.UKB_v3_INFO_0.30.${ethnicity}.${subset}.MAF_0.001.HWE.FINAL.bgen -index -clobber

# CLEAN UP
rm temp*
done

duration=$SECONDS
echo "FINISHED EVERYTHING! Elapsed time for merging PAR into X chromosome: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"