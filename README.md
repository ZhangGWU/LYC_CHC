# LYC_CHC by Linyi Zhang 2023 #
## Get the genotype estimates from Entropy from ddRAD seq data 
### 1. File preparation ###
Path to the data 
```
cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/lycaeides_chc_experiment/fastq/
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
### get individual fastq files ###
```
#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=L1splitfastq

## This script is to run the parse_barcodes.pl script on the cluster


cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/lycaeides_chc_experiment/fastq/parsed/library1/

module load perl
perl ../splitFastq.pl ../lycID_L1new.csv parsed_comb_L1.fastq
```
Note: you can find [splitFastq.pl](splitFastq.pl) here. The script here is written to match the individual id like this: "VIC-L-20-17". If your individual ID doesn't look like this, you might need to adjust the individual ID row to match your individual id with the splitFastq.pl script.
This is the code line in splitFastq.pl that you need to adjust to match you individual ID: if (/^\@([A-Za-z0-9\-_]+)/).Here A-Za-z0-9\-_ means any uppercase letter (A-Z), any lowercase letter (a-z), any digit (0-9), the hyphen (-), or the underscore (_). You can ask chatGPT how to revise the code to match your individual id in perl language. 

### You need a reference genome to do the part 4 analysis. If you don't have a assemble genome ready for analysis, you want to do a de Novo assembly based on your sequences (please follow steps in section 3 in the [NextGenNotes_v4_2019.pdf](NextGenNotes_v4_2019.pdf). 
Here is the list of the scripts you need to make a de novo genome assembly: [ddocent.sh](ddocent.sh), [ddocentPart2.sh](ddocentPart2.sh), [tapewormRemover.pl](tapewormRemover.pl). Note: Parallel comuputing are required in the ddocent script, you might want to check how to do the parallel computing in your cluster setting. 

### 4. Alignment and variant calling 
Aligning individual reads to the reference genome and get SNP calls 
```
conda activate YOURENVIRONMENTNAME ### this is to activate your environment where you loaded bwa, bcftools function.
```
Instructions on how to create your conda environment, install and load the programs you need, please ask chatgpt or any AI assistant. 
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

perl vcfFilter_CCN_1.9.pl variants.vcf #### need to chech the perl script to make sure the line search for loci matches the scaffold heads in your genome. 
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
bam_dir="/uufs/chpc.utah.edu/common/home/gompert-group3/data/lycaeides_chc_experiment/fastq/vcall"

# Output file name
output_file="total_coverage.txt"

# Get a list of all BAM files in the directory
bam_files=("$bam_dir"/*.sorted.bam)

# Iterate over BAM files
calculate_coverage () {

    bam_file="$1"
    # Extract individual name from the BAM file name
    individual=$(basename "$bam_file" .bam)

    # Use samtools to calculate coverage and extract relevant information
    coverage=$(samtools depth -a "$bam_file" | awk '{sum += $3} END {print sum}')

    # Print individual coverage to the output file
    echo -e "$individual\t$coverage"

}

# Export the function so it's accessible to parallel
export -f calculate_coverage

parallel -j 24 calculate_coverage ::: "${bam_files[@]}" > "$output_file"
```

Remove individuals with <2x depth (2*211376,211376 is the number of loci), remove 15 more individuals, leave 750 individuals

### 7. Redo variant calling and filtering (if you use individual sequences to do de novo genome assembly, you also redo the de novo genome assembly after excluding these low coverage individuals)####
in the directory vrecall
```
sbatch variantcall.sh
sbatch filter1.sh
sbatch filter2.sh
sbatch filter3.sh ### depth filtering: total loci： 217166， maxCoverage= mean + 2sd : 37413.41, make sure everything is the same as vcfFilter_CCN_1.9.pl, expect for depth coverage

mv filtered_secondRound_filtered_firstRound_variants.vcf doubleFiltered_variants.vcf ### rename the file ##
```
### 8. Prepare files for Entropy ###
#### Directory ###
```
/uufs/chpc.utah.edu/common/home/gompert-group3/data/lycaeides_chc_experiment/fastq/entropy/
```
#### Genotype likelihoods from the variants
```
perl vcf2mpgl_CCN_1.9.pl doubleFiltered_variants.vcf  ### need to change the expression in vcf2mpgl_CCN_1.9.pl to match the header of vcf

```
This generate doubleFiltered_variants.mpgl

###### You want to do another round of filtering if your genome assembly is a de novo assembly, in this filtering step, you are reducing the degree of linkage disequilibrium by only selection one locus per scaffold. Following the steps in [NextGenNotes_v4_2019.pdf](NextGenNotes_v4_2019.pdf).

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

perl gl2genest.pl doubleFiltered_variants.mpgl
```
This will generate gl_doubleFiltered_variants.mpgl, this genotype matrix is used for generating ldak files. R script for ldak files is `generate ldak.R`

#### Add the header, this is input file for entropy run ###

```
awk -F, '{print $1}' indiv_ids.txt |tr '\n' ' '>ids_col1.txt
awk -F- '{print $1}' indiv_ids.txt |tr '\n' ' '>ids_col2.txt
cat ids_col1.txt ids_col2.txt >header_ids.txt
cp header_ids.txt header_T.txt
## add #ind.s #loci and \1 to header_T.txt 
cat header_T.txt doubleFiltered_variants.mpgl >lyc_entropy.mpgl

```

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

cd /uufs/chpc.utah.edu/common/home/gompert-group3/data/lycaeides_chc_experiment/fastq/entropy/

perl forkEntropy_lmel.pl

```
### 10.1 Check if the MCMC chains converge ###

```
/uufs/chpc.utah.edu/common/home/u6000989/bin/estpost.entropy ES_out__k2_ch1.hdf5 ES_out__k2_ch2.hdf5 ES_out__k2_ch3.hdf5 ES_out__k2_ch4.hdf5 ES_out__k2_ch5.hdf5 ES_out__k2_ch6.hdf5 ES_out__k2_ch7.hdf5 ES_out__k2_ch8.hdf5 -p q -s 4 -o ESk2dig.txt
```
###### check if in the output file ESk2dig.txt, the values for each loci are close to 1, the closer it is, the more converge the chains are ###

### 10.2 Get the genotype estimates from entropy: run CHCgprob.sh ####
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

/uufs/chpc.utah.edu/common/home/u6000989/bin/estpost.entropy out_lmel_k2_ch1.hdf5 out_lmel_k2_ch2.hdf5 out_lmel_k2_ch3.hdf5 out_lmel_k2_ch4.hdf5 out_lmel_k2_ch5.hdf5 out_lmel_k2_ch6.hdf5 out_lmel_k2_ch7.hdf5 out_lmel_k2_ch8.hdf5 out_lmel_k2_ch9.hdf5 out_lmel_k2_ch10.hdf5  out_lmel_k3_ch1.hdf5 out_lmel_k3_ch2.hdf5 out_lmel_k3_ch3.hdf5 out_lmel_k3_ch4.hdf5 out_lmel_k3_ch5.hdf5 out_lmel_k3_ch6.hdf5 out_lmel_k3_ch7.hdf5 out_lmel_k3_ch8.hdf5 out_lmel_k3_ch9.hdf5 out_lmel_k3_ch10.hdf5 out_lmel_k4_ch1.hdf5 out_lmel_k4_ch2.hdf5 out_lmel_k4_ch3.hdf5 out_lmel_k4_ch4.hdf5 out_lmel_k4_ch5.hdf5 out_lmel_k4_ch6.hdf5 out_lmel_k4_ch7.hdf5 out_lmel_k4_ch8.hdf5 out_lmel_k4_ch9.hdf5 out_lmel_k4_ch10.hdf5 -o CHC_gprob.txt -p gprob -s 0

```
