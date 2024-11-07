# 1. Prepare the file for GWAS run.
> ### Here we have individuals reared in high (36/30) and low temperatures (26/20), we want to run GWAS seperately for these individuals. Therefore, we want to generate two genotype files BSTh_G.txt and BSTl_G.txt. The genoetype file looks like this, each row is an individual locus with their genotype for each individual. 

<img width="500" alt="image" src="https://github.com/user-attachments/assets/cc58d896-b997-4b6e-a143-c0ac25e53793">

The R code for preparing the Phenotype and Genotype file for GWAS is here [GWAS_file_prep.R](GWAS_file_prep.R) . 

The output files BSTh_G.txt and BSTl_G.txt from R has quotes for the first three columns. We need to run this command to remove them:
```
awk -F ' ' '{gsub(/"/, "", $1); gsub(/"/, "", $2); gsub(/"/, "", $3); print}' OFS=' ' BSTl_G.txt >BSTl_Gnew.txt
```

Phenotype file has only one column, see example [BSTh_growth.txt](BSTh_growth.txt). 

# 2. Estimating PVE (heritability) and polygenemic scores using BSLMMs with gemma

### perl script for running gemma forkRunGemmaH_growth.pl  

```
#!/usr/bin/perl
#
# fit gemma BSLMM for L. melissa, BST population, high temp treatment
#

use Parallel::ForkManager;
my $max = 50;
my $pm = Parallel::ForkManager->new($max);

$g = "BSTh_G.txt";
#$dir = "output_ranchem_$ph_ub";
$dir = "output_BSThigh";
 
		foreach $ch (0..9){
			sleep 2;
			$pm->start and next;
			$o = "BSTh_$base"."growth"."_ch$ch";
    			system "gemma -g $g -p BSTh_growth.txt -bslmm 1 -o $o -maf 0 -w 200000 -s 1000000 -outdir $dir\n";
			$pm->finish;
		}
	
$pm->wait_all_children;
```

### unix script for running gemma 
```
#!/bin/sh 
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=gemma


module load gemma
# 0.95a

perl forkRunGemmaH_growth.pl BSTh_growth.txt 
```
### MCMC convergence check 
[calcDiags.R](calcDiags.R)

### Summarize posterior results and obtain estimates of model-averaged effects
Growth_sum_post.sh
```
#!/bin/sh 
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=sumgemma


for dir in output_BST*
#for dir in output_chem_*
do 
cd $dir
echo "perl ../calpost.pl $dir *ch0*hyp.txt"
perl ../calpost.pl $dir BST*_growth_ch0.hyp.txt 
cd ..
done
```

Growth_sum_effects.sh
```
#!/bin/sh 
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=sumgemma


cd output_BSThigh

perl ../grabMavEffects.pl BSTh_growth_ch*.param.txt
```

You can find perl script [calpost.pl](calpost.pl) and [grabMavEffects.pl](grabMavEffects.pl) here. 

### summarize the GWAS result including polygenic score 

[GWASanalysis.R](GWASanalysis.R)

