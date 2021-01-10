#merging files
sort ./data/hsp_genotypes/VPU_ILLUMINA_APR_2020.bim ./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.bim | uniq -d | cut -f2 > ./data/common_snps.txt
plink --bfile ./data/hsp_genotypes/VPU_ILLUMINA_APR_2020 --bmerge ./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.bed ./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.bim ./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.fam --extract ./data/common_snps.txt --make-bed --out ./data/cleaned_genotypes/merged_data

