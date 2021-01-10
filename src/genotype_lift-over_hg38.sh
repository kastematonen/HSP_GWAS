# Genotype lift-over from hg19 to hg38 for merged HSP, IBD and HSCT data

# scripts based on protocol Genotyping chip data lift-over to reference genome build GRCh38/hg38 V.2

#############################################################################################################

# lift-over : step 2

# Getting chip information and script for lift-over:
#wget http://www.well.ox.ac.uk/~wrayner/strand/GSA-24v2-0_A1-b38.strand.zip
#wget http://www.well.ox.ac.uk/~wrayner/strand/update_build.sh

# Uncompress the files
#unzip /media/tk-tutk/TOSHIBA_3T_BU01/HSP/referenssigenomi/GSA-24v2-0_A1-b38-strand.zip

# -> unifying SNP names in .strand file with ./src/lift-over_snp_names.R
# -> new .strand file ./GSA-24v2-0_A1-b38_mod.strand

DATASET=./data/cleaned_genotypes/merged_cleaned_no_duplicates
OUTPUT=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod


#Run the script
./src/update_build.sh \
    ${DATASET} \
    GSA-24v2-0_A1-b38_mod.strand \
    ${OUTPUT}

###############################################################################################################

# lift-over : step 3.1 lift-over verification

FASTA=/media/tk-tutk/TOSHIBA_3T_BU01/HSP/reference-data-full/reference-data/fasta/Homo_sapiens_assembly38.fasta
BCFTOOLS_PLUGINS=/usr/lib/x86_64-linux-gnu/bcftools/

DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod # from step 2
#DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped # from step 5 
#DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped_ambiguous_flipped # from step 6

# Convert to VCF utilizing fasta file
plink2 \
    --bfile ${DATASET} \
    --recode vcf id-paste=iid bgz id-delim="*"\
    --ref-from-fa \
    --fa ${FASTA} \
    --output-chr chrM \
    --out ${DATASET}

###################################################################################################

# lift-over : step 3.2

BCFTOOLS_PLUGINS=/usr/lib/x86_64-linux-gnu/bcftools/

DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod # from step 2
VCF=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod.vcf.gz # from step 3.1

#DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped # from step 5
#VCF=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped.vcf.gz # from step 3.1

#DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped_ambiguous_flipped # from step 6
#VCF=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped_ambiguous_flipped.vcf.gz # from step 3.1

# Extract frequencies with BCFtools
# Requires environment variable BCFTOOLS_PLUGINS
export BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS}
bcftools +fill-tags ${VCF}\
     -Oz -o ${DATASET}_AF.vcf.gz -- -t AF
bcftools query \
    -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\n' \
    ${DATASET}_AF.vcf.gz | \
sed '1iCHR\tSNP\tREF\tALT\tAF' > ${DATASET}.frq

##################################################################################################

# lift-over : step 4 - was lift-over successful or not

DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod # from step 2 
#DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped # from step 5
#DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped_ambiguous_flipped # from step 6

Rscript --no-save ./src/compare_AF.R \
    ${DATASET} \
    ${DATASET}.frq \
    /media/tk-tutk/TOSHIBA_3T_BU01/HSP/reference-data-full/reference-data/frq/1000GP_EUR_GRCh38.frq \
    1000GP_EUR_GRCh38 \
    0.1

# genome build lift-over was successful the first time -> continuing to step 5 to fix more variants
# genome build lift-over was successful after step 5 -> continue to step 6  to fix more variants
# genome build lift-over was successful after step 6

########################################################################################################

# lift-over : step 5

FASTA=/media/tk-tutk/TOSHIBA_3T_BU01/HSP/reference-data-full/reference-data/fasta/Homo_sapiens_assembly38.fasta
DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod # from step 2 

# Align to reference fasta file and 
# keep only RSIDs from the multiallelic sites
bcftools norm -f ${FASTA} -c ws ${DATASET}.vcf.gz -Ou | \
bcftools view -m 3 -Ou | \
bcftools query -f '%ID\n' \
> ${DATASET}_nonreference_alleles.rsid

plink \
    --bfile ${DATASET} \
    --flip ${DATASET}_nonreference_alleles.rsid \
    --make-bed \
    --out ${DATASET}_nonreference_flipped

#-> continue to step 3 and 4 for lift-over vericication

############################################################################################################

# lift-over : step 6 

DATASET=./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped

Rscript --no-save ./src/find_ambiguous.R \
    ${DATASET} \
    ${DATASET}.frq \
    /media/tk-tutk/TOSHIBA_3T_BU01/HSP/reference-data-full/reference-data/frq/1000GP_EUR_GRCh38.frq

bcftools view -T ${DATASET}_flippable_AF.txt \
    ${DATASET}.vcf.gz | \
bcftools query -f '%ID\n' \
    > ${DATASET}_ambiguous_for_flipping.rsid

plink \
    --bfile ${DATASET} \
    --flip ${DATASET}_ambiguous_for_flipping.rsid \
    --make-bed \
    --out ${DATASET}_ambiguous_flipped

#-> continue to step 3 and 4 for lift-over vericication



