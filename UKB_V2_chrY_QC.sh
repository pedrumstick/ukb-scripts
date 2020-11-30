# Y chromosome QC for Rob
# https://biobank.ctsu.ox.ac.uk/crystal/crystal/docs/ukbgene_instruct.html
# https://www.ahajournals.org/doi/10.1161/ATVBAHA.119.312405 (813 and 1041 SNPs on Y chromosome; 162 after QC; 168,686 males of European ancestry with available data; Male-specific region of the Y chromosome (MSY) constitutes 95% of the entire length of the human Y chromosome)

# DEFINE VARIABLES
plink2="/hightide/software/plink_2/plink2"
chr=Y
output_dir="/londonbridge/mohamp1/chrY_for_Rob/"
geno_dir="/londonbridge/UKB_TRANSFER_2017_06_09/xfer.crg.eu/EGAD00010001226/001/"
fam_file="/londonbridge/UKB_TRANSFER_2017_06_09/xfer.crg.eu/OTHER_UKB_FILES/ukb1525_cal_v2_s488374.fam" # confirm whether this is the right file????!!!!!

######################################################
# SPLIT INTO RELATED AND UNRELATED
######################################################
cd ${output_dir}
SECONDS=0
# RELATED
${plink2} --bed ${geno_dir}/ukb_cal_chr${chr}_v2.bed \
    --bim ${geno_dir}/ukb_snp_chr${chr}_v2.bim \
    --fam ${fam_file} \
    --keep-fam /londonbridge/UKB_FULL_HRC_BRIT_QC/UKB_FULL_BRIT_POTENTIALLY_RELATED_QC_PASS_IDS.2018_09_15.txt \
    --make-bed \
    --out ${output_dir}/ukb_chr${chr}_v2_BRITISH_RELATED

# UNRELATED
${plink2} --bed ${geno_dir}/ukb_cal_chr${chr}_v2.bed \
    --bim ${geno_dir}/ukb_snp_chr${chr}_v2.bim \
    --fam ${fam_file} \
    --keep-fam /londonbridge/UKB_FULL_HRC_BRIT_QC/UKB_FULL_BRIT_UNRELATED_QC_PASS_IDS.2018_09_15.txt \
    --make-bed \
    --out ${output_dir}/ukb_chr${chr}_v2_BRITISH_UNRELATED

echo "There are "$(cat ${output_dir}/ukb_chr${chr}_v2_BRITISH_RELATED.fam | wc -l)" White British related individuals."
echo "There are "$(cat ${output_dir}/ukb_chr${chr}_v2_BRITISH_UNRELATED.fam | wc -l)" White British unrelated individuals."

######################################################
# INCLUDE MALES ONLY
######################################################
sex=Male
for subset in RELATED UNRELATED; do
    ${plink2} --bfile ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset} \
        --remove-fam /hightide/mohamp1/UKB_IDs/ids_exclude_for${sex}.txt \
        --make-bed \
        --out ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}
done

echo "There are "$(cat ${output_dir}/ukb_chr${chr}_v2_BRITISH_RELATED_Male.fam | wc -l)" White British related males."
echo "There are "$(cat ${output_dir}/ukb_chr${chr}_v2_BRITISH_UNRELATED_Male.fam | wc -l)" White British unrelated males."

######################################################
# EXCLUDE RARE SNPs (MAF < 0.1%)
######################################################
# Principles of HWE don't apply to Y chromosome since it's haploid
for subset in RELATED UNRELATED; do
    ${plink2} --bfile ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex} \
        --maf 0.001 \
        --make-bed \
        --out ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}_MAF_0dot001
done

echo "There are "$(cat ${output_dir}/ukb_chr${chr}_v2_BRITISH_RELATED_Male_MAF_0dot001.bim | wc -l)" genotype SNPs on Y chromosome with MAF > 0.01% in the related subset."
echo "There are "$(cat ${output_dir}/ukb_chr${chr}_v2_BRITISH_UNRELATED_Male_MAF_0dot001.bim | wc -l)" genotype SNPs on Y chromosome with MAF > 0.01% in the unrelated subset."

######################################################
# CHECK FOR SNPs IN PSEUDOAUTOSOMAL REGION (https://en.wikipedia.org/wiki/Pseudoautosomal_region)
######################################################
# par1_start=10001
# par1_end=2649520
# par2_start=59034050
# par2_end=59363566

# ${plink2} --bfile ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex} \
# --chr Y \
# --write-snplist \
# --from-bp ${par1_start} \
# --to-bp ${par1_end}

# ${plink2} --bfile ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex} \
# --chr Y \
# --write-snplist \
# --from-bp ${par2_start} \
# --to-bp ${par2_end}

######################################################
# RECODE MALES AS DIPLOID (0/2)
######################################################
# When analyzing X chromosome data, some tools (e.g., BOLT) will want you to explicitly code male genotypes as diploid (0/2). By default, plink2 will code males as haploid, but you can force it to call diploid genotypes by setting the sex of all individuals to female before converting.
# I don't know if plink2 will do the same for Y chromosome, but I've created a diploid version of the files if there are issues.

for subset in RELATED UNRELATED; do
    awk 'BEGIN {FS=OFS="\t"} {$5=2; print}' ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}_MAF_0dot001.fam >${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}_MAF_0dot001_recodeDiploid.fam # code all samples as female to ensure PLINK codes males as diploid (0/2)
    ${plink2} --bed ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}_MAF_0dot001.bed \
        --bim ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}_MAF_0dot001.bim \
        --fam ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}_MAF_0dot001_recodeDiploid.fam \
        --make-bed \
        --out ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}_MAF_0dot001_recodeDiploid

    mv ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}_MAF_0dot001.fam ${output_dir}/ukb_chr${chr}_v2_BRITISH_${subset}_${sex}_MAF_0dot001_recodeDiploid.fam # recode fam file as males
done

duration=$SECONDS
echo "FINISHED! Elapsed time: $(($duration / 3600)) hours $((($duration - $duration / 3600 * 3600) / 60)) min $((duration % 60)) sec"
