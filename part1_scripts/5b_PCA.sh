#ok first things first lets do some linkage pruning
cd ~/
mkdir plink
cd ~/plink

#set the variables
VCF=/home/imh4101/cohort_gvcfs/final_4.recode.vcf

#linkage pruning
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out elm
#wow that remove A LOT of SNPs
#now were down to 128

#now to run the PCA
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --mind 0.5 \
--extract elm.prune.in \
--make-bed --pca --out elm
#--mind gets rid of individuals who are missing more then 80% of data
#that removed 112 individuals
#lets slide that back to 0.5
#ok now only 29 went missing
#kept 128 variants
