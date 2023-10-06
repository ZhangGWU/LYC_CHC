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
`RunParseForkl.pl` is attached in the depository 

#### Get individual ID file 
```
cut -d',' -f3- lyc_barcodeKey_L1.csv >lycID_L1.csv ### remove the first two columns
sed 1d lycID_L1.csv> lycID_L1new.csv

cut -d',' -f3- lyc_barcodeKey_L2.csv >lycID_L2.csv
sed 1d lycID_L2.csv> lycID_L2new.csv
```
#### Combine all fastq files #######
```
cat parsed_x* > parsed_comb_L1.fastq
cat parsed_x* > parsed_comb_L2.fastq
```

### 4. Alignment and variant calling
Aligning individual reads to the reference genome and get SNP calls 
```
conda activate pipeline-structural-variation ### this is to activate bwa function 
```
#### Create index files from reference genome
```
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=Index

bwa index -p lycCHC_ref -a is /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta
```

Generate sam files for variant calling 
change ref.ID in the file `runbwa.pl`
```
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=samfile

module load perl
perl runbwa.pl *fastq
```

Convert sam file into bam file
code file is sam2bam.sh 
```
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=samtobam

module load perl
perl sam2bam.pl *.sam
```
#### Call SNPS 
Code file is variantcall.sh
```
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
bcftools mpileup -d 8000 -o lyc_CHC.bcf -O b -I -f /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta aln*sorted.bam
```
Code file is variantcall2.sh
```
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
bcftools call -c -V indels -v -p 0.05 -P 0.001 -o CHCvariants.vcf lyc_CHC.bcf
```
### 5. Filtering ######
First round of filtering include a series of parameters including 
+ minCoverage = 2* number of individuals
+ allele frequency min = 0.001
+ allele frequency max = 0.999
+ missing number of individuals with no data = 153 ### 20% of individuals ###

Filtering round 1 code file: filter1.sh
```
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=filter1CHC

perl vcfFilter_CCN_1.9.pl variants.vcf
```

Get the coverage depth for each loci, filter2.sh
```
#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=filter2CHC

perl depthCollector_1.9.pl filtered_firstRound_variants.vcf
```

Use R script to get the depth filter stats # The strategy is to  lter those loci which have extremely high sequencing depth 
because they are likely to be assemblies of paralogous loci.
```
dat<-read.table("depth_filtered_firstRound_variants.vcf",header=F)
dim(dat)
dat[1:5,]
summary(dat$V3)
sd(dat$V3)
max<-mean(dat$V3)+2*sd(dat$V3)
length(which(dat$V3 > max))

```

Filtering round 2 code file: filter3.sh ### depth filtering: maxCoverage= mean + 2sd : 37760.82, make sure everything is the same as vcfFilter_CCN_1.9.pl, expect for depth coverage
```
#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=filter3CHC

perl vcfFilter_CCN_1.9_more.pl filtered_firstRound_variants.vcf
```

```
mv filtered_secondRound_filtered_firstRound_variants.vcf doubleFiltered_variants.vcf
```
### 6. Quantify coverage ###
Code file for coverage counting: sbatch coverage.sh
```
#!/bin/sh
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert
#SBATCH --partition=kingspeak
#SBATCH --job-name=EScov2


# Directory containing BAM files
bam_dir="/uufs/chpc.utah.edu/common/home/gompert-group3/data/lycaeides_chc_experiment/fastq/alignment/bamfile/"

# Output file name
output_file="total_coverage.txt"

# Get a list of all BAM files in the directory
bam_files=("$bam_dir"/*.sorted.bam)

# Iterate over BAM files
for bam_file in "${bam_files[@]}"; do
    # Extract individual name from the BAM file name
    individual=$(basename "$bam_file" .bam)

    # Use samtools to calculate coverage and extract relevant information
    coverage=$(samtools depth -a "$bam_file" | awk '{sum += $3} END {print sum}')

    # Print individual coverage to the output file
    echo -e "$individual\t$coverage" >> "$output_file"
done
```
Output file: total_coverage.txt
Remove individuals with <2x depth

### 7. Redo variant calling and filtering ####

### 8. Prepare files for Entropy ###
##### directory ###
```
/uufs/chpc.utah.edu/common/home/gompert-group3/data/lycaeides_chc_experiment/fastq/alignment/entropynew/

```
#### Genotype likelihoods from the variants
```
perl vcf2mpgl_CCN_1.9.pl doubleFiltered_variants.vcf  ### need to change the expression in vcf2mpgl_CCN_1.9.pl to match the header of vcf

```
This generate doubleFiltered_variants.mpgl, remove one individual which is lack of individual id 
```
cut -d' ' -f1-1762,1766-2294 doubleFiltered_variants.mpgl >doubleFiltered_variantsnew.mpgl
```
Add the header, this is input file for entropy run

```
cat header_ids.txt doubleFiltered_variantsnew.mpgl >lyc_variantsnew.gl
```

#### Transform the genotype likelihood files into a genotype matrix (point estimate of genotype) 
This is the code file mpgl2peg.sh 
```
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=peg

perl gl2genest.pl doubleFiltered_variantsnew.mpgl
```
This will generate gl_doubleFiltered_variantsnew.mpgl, this genotype matrix is used for generating ldak files. R script for ldak files is in the depository

### 9. Run entropy: run_entropy_lmel.sh
```
#!/bin/sh 
#SBATCH --time=172:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=entropy
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=linyizhangecnu@gmail.com

module purge
module load gcc/8.5.0 hdf5/1.10.7

cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/lycaeides_chc_experiment/fastq/alignment/entropynew/

perl forkEntropy_lmel.pl

```

### 10. Get the genotype estimates from entropy: run CHCgprob.sh ####
```
#!/bin/sh 
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=CHCgprob


module purge
module load gcc/8.5.0 hdf5/1.10.7

/uufs/chpc.utah.edu/common/home/u6000989/bin/estpost.entropy out_lmel_k2_ch1.hdf5 out_lmel_k2_ch2.hdf5 out_lmel_k2_ch3.hdf5 out_lmel_k2_ch4.hdf5 out_lmel_k2_ch5.hdf5 out_lmel_k2_ch6.hdf5 out_lmel_k2_ch7.hdf5 out_lmel_k2_ch8.hdf5 out_lmel_k2_ch9.hdf5 out_lmel_k2_ch10.hdf5  out_lmel_k3_ch1.hdf5 out_lmel_k3_ch2.hdf5 out_lmel_k3_ch3.hdf5 out_lmel_k3_ch4.hdf5 out_lmel_k3_ch5.hdf5 out_lmel_k3_ch6.hdf5 out_lmel_k3_ch7.hdf5 out_lmel_k3_ch8.hdf5 out_lmel_k3_ch9.hdf5 out_lmel_k3_ch10.hdf5 out_lmel_k4_ch1.hdf5 out_lmel_k4_ch2.hdf5 out_lmel_k4_ch3.hdf5 out_lmel_k4_ch4.hdf5 out_lmel_k4_ch5.hdf5 out_lmel_k4_ch6.hdf5 out_lmel_k4_ch7.hdf5 out_lmel_k4_ch8.hdf5 out_lmel_k4_ch9.hdf5 out_lmel_k4_ch10.hdf5 -o CHC_gprobnew.txt -p gprob -s 0

```
