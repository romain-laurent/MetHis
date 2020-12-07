
# this README explains how to produce MetHis source population genetic reservoirs from a VCF file
# the files needed are:
#  - s1.txt -> list of samples for source population 1
#  - s2.txt -> list of samples for source population 2
#  - original.vcf -> a VCF file containing genotypes for samples in s1.txt and s2.txt

# extract samples from source population
# this will produce file "source_pops_s1_s2.vcf"
vcftools --vcf original.vcf --keep s1.txt --keep s2.txt --recode --out source_pops_s1_s2

# find independant SNPs
# this requires plink (https://www.cog-genomics.org/plink2) to be installed
# this will produce files "indep_SNPs.prune.in" and "indep_SNPs.prune.out"
plink --vcf source_pops_s1_s2.vcf --const-fid --indep-pairwise 100 10 0.1 --out indep_SNPs

# subsample SNPs
# you can set the number of SNPs wanted with variable "nb_wanted" in the Python script
# this will produce the file "wanted_SNPs.txt"
./subsample_SNPs.py > wanted_SNPs.txt

# extract SNPs from VCF
# this will produce file "wanted_SNPs_s1_s2.vcf"
vcftools --vcf source_pops_s1_s2.vcf --snps wanted_SNPs.txt --recode --out wanted_SNPs_s1_s2

# compute allele frequencies in source populations
# this will produce files "sfs_s1.frq" and "sfs_s2.frq"
vcftools --vcf wanted_SNPs_s1_s2.vcf --keep s1.txt --freq --out sfs_s1
vcftools --vcf wanted_SNPs_s1_s2.vcf --keep s2.txt --freq --out sfs_s2

# generate data
# this will produce file "MetHis_input.arp" containing 20000 gametes for each source population
# filenames and number of gametes can be modified in the Python script in the "main" function at the end of the file
./generate_data_from_sfs.py
