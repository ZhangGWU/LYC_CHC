# LYC_CHC by Linyi Zhang 2023 #

### 1. File preparation ###
Path to the data 
```
cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/lycaeides_chc_experiment/fastq/alignment
```

Convert library barcodes from mac to Unix format 
```
dos2unix lyc_barcodeKey_L1.csv
mac2unix lyc_barcodeKey_L1.csv 

dos2unix  lyc_barcodeKey_L2.csv
mac2unix lyc_barcodeKey_L2.csv 
```

### 2. Split fastq files ####
This is because the parsing time for the single big fastq file is too long, well beyond wall time 
```
split -l 90000000 ../Gomp032_S1_L001_R1_001.fastq
split -l 90000000 ../Gomp033_S2_L002_R1_001.fastq
```

### 3. Parse ### 
Removing barcodes and cut sites from raw reads and attaching ind. IDs

This is the code file: parse.sh
```
#!/bin/sh 
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=L1lycparse

cd library1

module load perl

perl ../RunParseFork.pl x*
```
RunParseForkl.pl is attached in the depository 

#### Get individual ID file 
```
cut -d',' -f3- lyc_barcodeKey_L1.csv >lycID_L1.csv ### remove the first two columns
sed 1d lycID_L1.csv> lycID_L1new.csv

cut -d',' -f3- lyc_barcodeKey_L2.csv >lycID_L2.csv
sed 1d lycID_L2.csv> lycID_L2new.csv
```
#### combine all fastq files #######
```
cat parsed_x* > parsed_comb_L1.fastq
cat parsed_x* > parsed_comb_L2.fastq
```

### 4. Alignment and variant calling ###
##### aligning individual reads to the reference genome and get SNP calls ##### 
conda activate pipeline-structural-variation
######## create index files from reference genome #######
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=Index

bwa index -p lycCHC_ref -a is /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta

########## generate sam files for variant calling ###########
########## change ref.ID in runbwa.pl ###########
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=samfile

module load perl
perl runbwa.pl *fastq

nano sam2bam.sh 
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=samtobam

module load perl
perl sam2bam.pl *.sam

########## call SNPS #############
nano variantcall.sh
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
bcftools mpileup -d 8000 -o lyc_CHC.bcf -O b -I -f /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta aln*sorted.bam 

nano variantcall2.sh
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
bcftools call -c -V indels -v -p 0.05 -P 0.001 -o CHCvariants.vcf lyc_CHC.bcf

### 5. Filtering ######
sbatch filter1.sh
sbatch filter2.sh

########### r script , get the depth filter stats #
####### first round filtering: remove alleles that are fixed (afreqMin= 0.001, afreqMax=0.999); set p value for bqrs, mqrs, rprs: 0.00001 ####
############ depth filtering: maxCoverage= mean + 3sd : 49181.89###

perl vcfFilter_CCN_1.9 more.pl  filtered_firstRound_variants.vcf  ### make sure everything is the same as vcfFilter_CCN_1.9.pl, expect for depth coverage ###
mv filtered_secondRound_filtered_firstRound_variants.vcf doubleFiltered_variants.vcf

### 6. quantify coverage ###
sbatch coverage.sh  #### output file: total_coverage.txt
######### remove individuals with <2x depth ########

### 7. redo variant calling and filtering ####

### 8. prepare files for Entropy ###
###### genotype likelihoods from the variants ########
perl vcf2mpgl_CCN_1.9.pl doubleFiltered_variants.vcf  ### need to change the expression to match the header of vcf ### this generate doubleFiltered_variantsnew.mpgl

cat header_ids.txt doubleFiltered_variantsnew.mpgl >lyc_variantsnew.gl

############ remove one individual which is lack of individual id ###################
cut -d' ' -f1-1762,1766-2294 doubleFiltered_variants.mpgl >doubleFiltered_variantsnew.mpgl

###### transform the genotype likelihood files into a genotype matrix (point estimate of genotype) ####
sbatch mpgl2peg.sh ######## perl gl2genest.pl doubleFiltered_variantsnew.mpgl ###########
######## this will generate gl_doubleFiltered_variantsnew.mpgl, this genotype matrix is used for generating ldak files ######

sbatch run_entropy_lmel.sh



