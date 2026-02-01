# vcftools Filter
vcftools --vcf chicken.ref_population.vcf --max-missing 0.9 --maf 0.05 --min-alleles 2 --max-alleles 2 --recode --stdout > chicken.ref_population.filter.vcf

# beagle filling
java -jar -Xmx200000m /gpfs/users/zhangz/biodata/2023/wiseed_GS/software/beagle.jar gt=chicken.ref_population.filter.vcf out=chicken.ref_population.filter_impute

# Plink Processing
plink --vcf chicken.ref_population.filter_impute.vcf.gz --double-id --keep-allele-order --allow-extra-chr --chr-set 90 --set-missing-var-ids @:# --make-bed --out chicken.ref_population.filter_impute_012 
plink --bfile chicken.ref_population.filter_impute_012 --keep-allele-order --chr-set 90 --allow-extra-chr --recodeA --out chicken.ref_population.filter_impute_012
cat chicken.ref_population.filter_impute_012.raw|cut -d\" \" -f1,7- |sed '1s/_[A-Z,0]//g' >chicken.ref_population.filter_impute_012.txt