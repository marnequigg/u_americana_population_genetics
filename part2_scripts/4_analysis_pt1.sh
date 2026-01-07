############################################################################################################################################
##### VCF FILTERING #####
#trial 1: strict - ideal scenario according to Brenden
vcftools --vcf final_2n.vcf --remove-indels --maf 0.05 --max-missing 0.50 --minQ 30 --min-meanDP 5 --max-meanDP 40 --minDP 10 --maxDP 40 --recode --recode-INFO-all --out final_2n_strict
#kept 78/78 individuals
#kept 406/561497

#trial 2: Brenden's parameters
vcftools --vcf final_2n.vcf --remove-indels --maf 0.05 --max-missing 0.50 --minQ 30 --min-meanDP 5 --max-meanDP 80 --minDP 5 --maxDP 80 --recode --recode-INFO-all --out final_2n_lazy
#kept 78/78 individuals
#kept 998/561497

#based on the vcf results
vcftools --vcf final_2n.vcf --remove-indels --maf 0.1 --max-missing 0.80 --minQ 30 --min-meanDP 5 --max-meanDP 40 --minDP 5 --maxDP 40 --recode --recode-INFO-all --out final_2n_ideal
#kept 78/78 individuals
#kept 304/561497
###using this one moving forward

############################################################################################################################################
##### LINKAGE PRUNING #####
cd ~/gatk_again/diploids/
mkdir 05.linkage_prune
cd /home/imh4101/gatk_again/diploids/05.linkage_prune

###I want to prune any possible sites with linkage from all 3 vcfs: strict, laxy, ideal

##start with strict
#set the variables
VCF=/home/imh4101//home/imh4101/gatk_again/diploids/03.vcfs/final_2n_strict.recode.vcf

#linkage pruning
plink --vcf /home/imh4101/gatk_again/diploids/03.vcfs/final_2n_strict.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out 2n_strict
#pruned 295 out of 406
#111 left

plink --vcf /home/imh4101/gatk_again/diploids/03.vcfs/final_2n_lazy.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out 2n_lazy
#pruned 734 out of 998
#264 left

plink --vcf /home/imh4101/gatk_again/diploids/03.vcfs/final_2n_ideal.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out 2n_ideal
#pruned 191 of 304
#113 left

############################################################################################################################################
##### RUN PCA #####
cd /home/imh4101/gatk_again/diploids
mkdir 06.PCA
cd /home/imh4101/gatk_again/diploids/06.PCA/

#start with the lazy one (most SNPs)
plink --vcf /home/imh4101/gatk_again/diploids/03.vcfs/final_2n_lazy.recode.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# --mind 0.5 \
--extract /home/imh4101/gatk_again/diploids/05.linkage_prune/2n_lazy.prune.in \
--make-bed --pca --out 2n_lazy
#removed 9 trees due to lack of data

#and onwards to the strict one
plink --vcf /home/imh4101/gatk_again/diploids/03.vcfs/final_2n_strict.recode.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# --mind 0.5 \
--extract /home/imh4101/gatk_again/diploids/05.linkage_prune/2n_strict.prune.in \
--make-bed --pca --out 2n_strict
#removed 14 trees

#now finally the ideal one
plink --vcf /home/imh4101/gatk_again/diploids/03.vcfs/final_2n_ideal.recode.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# --mind 0.5 \
--extract /home/imh4101/gatk_again/diploids/05.linkage_prune/2n_ideal.prune.in \
--make-bed --pca --out 2n_ideal
#removed 0 trees

############################################################################################################################################
##### RUN ADMIXTURE #####
conda activate admixture

cd /home/imh4101/gatk_again/diploids/07.admixture

#generate the files for admixture
plink --bfile /home/imh4101/gatk_again/diploids/06.PCA/2n_ideal/2n_ideal --recode12 --out ideal_2n --allow-extra-chr --double-id

#then we use this line of code to call the script from github and run it
wget https://raw.githubusercontent.com/dportik/admixture-wrapper/master/admixture-wrapper.py

#time to run it
screen -L python3 /home/imh4101/gatk_again/diploids/07.admixture/admixture-wrapper.py -i /home/imh4101/gatk_again/diploids/07.admixture --kmin 2 --kmax 25 --reps 10 -t 40 --cv 10
cd /home/imh4101/gatk_again/diploids/07.admixture/Outputs-ideal_2n
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > ideal_2n.cv.error
