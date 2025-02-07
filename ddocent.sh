#!/bin/bash
#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --account=gompert
#SBATCH --partition=kingspeak
#SBATCH --job-name=Aqlddocent
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=linyi.zhang@usu.edu
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(\%z)) {print "$v\t$k\n";}'

cat namelist | parallel --no-notice -j 8 "zcat {}.fastq | mawk '$AWK1' | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"

mkdir zippedFastqs
mv *.fastq.gz zippedFastqs

cat *.uniq.seqs > uniq.seqs

parallel --no-notice -j 8 mawk -v x=4 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp: $z{$_}++;} while (($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv

for ((i = 2; i <= 10; i++));
do
echo $i >> ufile
done

cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data

rm ufile

