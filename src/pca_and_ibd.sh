# LD pruning the data, PCA, and IBD

# finding LD-based variants
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/cleaned_merged_1_23 \
  --indep-pairwise 50 5 0.5 \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/snps_for_pruning

# excracting LD variants
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/cleaned_merged_1_23 \
  --extract /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/snps_for_pruning.prune.in \
  --make-bed \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/ld_pruned_cleaned_merged_1_23

# PCA
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/ld_pruned_cleaned_merged_1_23 \
  --pca header var-wts \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/pca_ld_pruned_cleaned_merged_1_23

# IBD 
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/ld_pruned_cleaned_merged_1_23 \
  --genome \
  --min 0.2 \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/ibd_ld_pruned_cleaned_merged_1_23

# identifying PCA outliers and listing relatives based on IBD for removal: pca_plots_and_ibd_list.R

# IBD relatives: removing individual with more missingness from each relative pair
# missingness:

plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/cleaned_merged_1_23 \
  --keep /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/relatives.txt \
  --missing \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/relatives_missingness

# removing related individuals & pca outliers:
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/cleaned_merged_1_23 \
  --remove /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/remove_relatives.txt \
  --make-bed \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_cleaned_merged_1_23

# PCA for final data set (before platform bias analysis):
# finding LD-based variants
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_cleaned_merged_1_23 \
  --indep-pairwise 50 5 0.5 \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/final_snps_for_pruning

# excracting LD variants
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_cleaned_merged_1_23 \
  --extract /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/final_snps_for_pruning.prune.in \
  --make-bed \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/final_ld_pruned_cleaned_merged_1_23

# PCA
plink \
  --bfile /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/final_ld_pruned_cleaned_merged_1_23 \
  --pca header var-wts \
  --out /media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/final_pca_ld_pruned_cleaned_merged_1_23





