# cleaning HSP, IBD and HSCT data, removing samples based on genotyping reports:

plink \
  --bfile ./data/cleaned_genotypes/merged_data \
  --geno 0.05 \
  --mind 0.1 \
  --hwe 1E-6 \
  --maf 0.01 \
  --remove ./data/remove_samples.txt \
  --make-bed \
  --out ./data/cleaned_genotypes/merged_cleaned_data

# removing duplicate variants:
# uses a duplicate list made with ./src/find__duplicates.R

plink \
  --bfile ./data/cleaned_genotypes/merged_cleaned_data \
  --exclude ./data/cleaned_genotypes/duplicates_cleaned_merged.txt \
  --make-bed \
  --out ./data/cleaned_genotypes/merged_cleaned_no_duplicates

