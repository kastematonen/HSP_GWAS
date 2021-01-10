
# dosage file for GWAS 

# create dosage file for each chromosome
for chr in {8..23}; do

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_chr_"${chr}" \
      --recode A \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/dosage_chr_"${chr}";
done

# create dosage file for each chromosome (half a chromosome at a time -> smaller files)
for chr in {1..7}; do

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_chr_"${chr}"_beginning \
      --recode A \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/dosage_chr_"${chr}"_beginning;

    plink \
      --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_chr_"${chr}"_end \
      --recode A \
      --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/dosage_chr_"${chr}"_end;
done





