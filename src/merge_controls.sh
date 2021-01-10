# MERGING BLOOD DONOR DATA TO HSP, IBD, AND HSCT DATA 

# variant names have been unified with ./src/unify_snp_names.R and files moved to appropriate locations beforehand

# unifying file names for handling the blood donor files all at once
for chr in {11..20}; do
    mv /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr"${chr}"_BLOOD_SERVICE_extracted_20200326_rid.bed /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr"${chr}"_BLOOD_SERVICE_extracted_20200325_rid.bed;
    mv /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr"${chr}"_BLOOD_SERVICE_extracted_20200326_rid.bim /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr"${chr}"_BLOOD_SERVICE_extracted_20200325_rid.bim;
    mv /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr"${chr}"_BLOOD_SERVICE_extracted_20200326_rid.fam /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr"${chr}"_BLOOD_SERVICE_extracted_20200325_rid.fam;
done

mv /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr22_BLOOD_SERVICE_extracted_20200327_rid.bed /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr22_BLOOD_SERVICE_extracted_20200325_rid.bed;
mv /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr22_BLOOD_SERVICE_extracted_20200327_rid.bim /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr22_BLOOD_SERVICE_extracted_20200325_rid.bim;
mv /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr22_BLOOD_SERVICE_extracted_20200327_rid.fam /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr22_BLOOD_SERVICE_extracted_20200325_rid.fam;

mv /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chrX_BLOOD_SERVICE_extracted_20200327_rid.bed /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr23_BLOOD_SERVICE_extracted_20200325_rid.bed;
mv /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chrX_BLOOD_SERVICE_extracted_20200327_rid.bim /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr23_BLOOD_SERVICE_extracted_20200325_rid.bim;
mv /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chrX_BLOOD_SERVICE_extracted_20200327_rid.fam /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr23_BLOOD_SERVICE_extracted_20200325_rid.fam;


# extracting common snips (between blood donors and HSP, IBD, and HSCT data) from blood donor files to reduce the number of variants in the files and for easier handling
for chr in {1..23}; do
    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/yhdessa/finngen_R4_bb_chr"${chr}"_BLOOD_SERVICE_extracted_20200325_rid \
      --extract ./data/merge_controls/common_snps_chr"${chr}".txt \
      --make-bed \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/filtered/filtered_chr"${chr}";
done

# blood donor files for merging listed to /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/files_for_merging.txt

# merging blood donor data -> all chromosomes into one file:
plink \
  --merge-list /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/files_for_merging.txt \
  --make-bed \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/merged_1_23

# -> Error: 548 multiallelic snips (listed in .missnp file) -> excluding them to enable merging files
for chr in {1..23}; do
    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/filtered/filtered_chr"${chr}" \
      --exclude /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/merged_1_23-merge.missnp \
      --make-bed \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/filtered/multiallelic_filtered_chr"${chr}";
done

# blood donor files for merging listed to /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/files_for_merging_2.txt

# merging blood donor data -> all chromosomes into one file:
plink \
  --merge-list /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/files_for_merging_2.txt \
  --make-bed \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/merged_1_23

# checking sex information:
# for HSP, IBD, and HSCT samples:
plink \
  --bfile ./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped_ambiguous_flipped \
  --check-sex
  --out ./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped_ambiguous_flipped
# sex ok for HSP, IBD, and HSCT samples: information in .fam compatible with genetic information, no ambiguous sexes

# checking sex information:
# imputing sex for blood donor data:
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/merged_1_23 \
  --impute-sex \
  --make-bed \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/imputed_sex_merged_1_23
# successful for 19535/19543 samples -> remove individuals with ambiguous sex: 
# samples with no sex listed in /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/no_sex.txt

# removing individuals with ambiguous sex information from blood donor data:
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/imputed_sex_merged_1_23 \
  --remove /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/no_sex.txt \
  --make-bed \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/sex_ok_merged_1_23

# cleaning merged blood donor data:
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/sex_ok_merged_1_23 \
  --geno 0.05 \
  --mind 0.1 \
  --hwe 1E-6 \
  --maf 0.01 \
  --make-bed \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/cleaned_merged_1_23

# finding common snips for blood donor data and HSP, IBD, and HSCT data
# blood donor data was cleaned to include only SNPs in HSP, IBD, and HSCT data but some snips have been removed from blood donor data while cleaning the data and merging files. Listing common snips again ensures all snips in the final data set are not missing from any of the data sets
sort /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/cleaned_merged_1_23.bim ./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped_ambiguous_flipped.bim | uniq -d | cut -f2 > ./data/merge_controls/common_snps_filtered_controls.txt

# merging blood donor data with HSP, IBD, and HSCT data
plink \
  --bfile ./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped_ambiguous_flipped \
  --bmerge /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/cleaned_merged_1_23 \
  --extract ./data/merge_controls/common_snps_filtered_controls.txt \
  --make-bed \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/cleaned_merged_1_23






