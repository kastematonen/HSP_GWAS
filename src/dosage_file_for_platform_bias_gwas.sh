# dosage file for platform bias GWAS 

# extract snips for each chromosome & create dosage file for each chromosome
for chr in {8..23}; do

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_cleaned_merged_1_23 \
      --extract /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snips_chr_"${chr}".txt \
      --make-bed \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/filtered_chr_"${chr}";

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/filtered_chr_"${chr}" \
      --recode A \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/dosage_chr_"${chr}";
done

# extract snips for each chromosome & create dosage file for each chromosome (half a chromosome at a time -> smaller files)
for chr in {1..7}; do

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_cleaned_merged_1_23 \
      --extract /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snips_chr_"${chr}"_beginning.txt \
      --make-bed \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/filtered_chr_"${chr}"_beginning;

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_cleaned_merged_1_23 \
      --extract /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snips_chr_"${chr}"_end.txt \
      --make-bed \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/filtered_chr_"${chr}"_end;

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/filtered_chr_"${chr}"_beginning \
      --recode A \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/dosage_chr_"${chr}"_beginning;

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/filtered_chr_"${chr}"_end \
      --recode A \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/dosage_chr_"${chr}"_end;
done

#-------------------------------------------------------------------------------------------------------------------------------------

# platform bias GWAS with ./src/platform_bias_gwas.R

#--------------------------------------------------------------------------------------------------------------------------------------

# removing SNPs associated with platform (both from the entire data set as well as from individual chromosome files)

plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_cleaned_merged_1_23 \
      --exclude /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snps_for_removal.txt \
      --make-bed \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final


for chr in {1..7}; do

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/filtered_chr_"${chr}"_beginning \
      --exclude /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snps_for_removal.txt \
      --make-bed \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_chr_"${chr}"_beginning;

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/filtered_chr_"${chr}"_end \
      --exclude /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snps_for_removal.txt \
      --make-bed \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_chr_"${chr}"_end;


done

for chr in {8..23}; do

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/filtered_chr_"${chr}" \
      --exclude /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snps_for_removal.txt \
      --make-bed \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_chr_"${chr}";
done

