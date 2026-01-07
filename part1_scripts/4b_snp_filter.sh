#trial 1: strict - ideal scenario according to Brenden
vcftools --vcf final_combined.vcf --remove-indels --maf 0.05 --max-missing 0.50 --minQ 30 --min-meanDP 5 --max-meanDP 40 --minDP 10 --maxDP 40 --recode --recode-INFO-all --out final_strict
#kept 593/593 individuals
#kept 542/4461904

#trial 2: Brenden's parameters
vcftools --vcf final_combined.vcf --remove-indels --maf 0.05 --max-missing 0.50 --minQ 30 --min-meanDP 5 --max-meanDP 80 --minDP 5 --maxDP 80 --recode --recode-INFO-all --out final_b
#kept 593/593 individuals
#kept 1326/4461904

#based on the vcf results
vcftools --vcf final_combined.vcf --remove-indels --maf 0.1 --max-missing 0.80 --minQ 30 --min-meanDP 5 --max-meanDP 40 --minDP 5 --maxDP 40 --recode --recode-INFO-all --out final_4
#this only kept 299 sites
