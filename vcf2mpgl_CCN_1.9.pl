#!/usr/bin/perl
#
# Time-stamp: <Tuesday, 27 May 2014, 11:23 MDT -- cbuerkle>

# 22ix19 CCN modifications (after other modifications) to use with vcf files from samtools/bcftools 1.9
# - the way I am running bcftools, the vcf does not include the genotype quality score that was reported at the end of the likelihoods in older versions

# This scripts converts a vcf file to a simpler format for downstream
# analysis. I am calling this format multiple population genetoype
# likelihood (mpgl). The first line lists: the number of individuals,
# loci.  There is a following line, which we do not use.  Instead,
# look in the indiv_ids.txt for the order and identifiers for inds.
# with one entry per individual. This is followed by one line per SNP
# that gives the SNP id (scaffold, position) and the phred-scaled
# genotype likelihoods, three per individual.

#
# USAGE: vcf2mpgl.pl in.vcf
# Follow with cat header_mpglvariants_0.9.txt variants_0.9.mpgl > entropy_in_variants_0.9_all.mpgl

use warnings;

my @line = ();
my $word;
my %keepind;
my $nind = 0;
my $nloc = 0;
my $nline = 0;


# my $idin = shift @ARGV;
# open (ID, "$idin") or die "Could not read the id file\n";
# ## individual names we want to keep
# chomp(<ID>);
# @line = split / /;
# foreach $word (@line){
#     $keepind{$word} = 1;
# }
# push (@ids,$_);
# while (<ID>){
#     chomp;
#     push (@ids,$_);
#     $nline++;
# }
# close (ID);


my $in = shift (@ARGV);
open (IN, $in) or die "Could not read the vcf file\n";
my $out = $in;
if ($out =~ s/vcf/mpgl/){
    open (OUT, "> $out") or die "Could not write $out\n";
}
else {
    die "File extension was not vcf\n";
}

while (<IN>){
    chomp;

    ## read genetic data lines, write gl
    ## match both scaffold111 and C113 and exclude loci with multiple alternative alleles
    if (m/^(\w+\d+)\s+(\d+)\s+[.\w]+\s+[AGCT]\s+[AGCT]\s+/){ 
	print OUT "$1:$2 ";   
	@line = split /\s+/;
	$i = 0;
	foreach $word (@line){
	    if ($word =~ s/^\d\/\d\://){ # && ($word =~ s/\:\d+$//)){
		$word =~ s/,/ /g;
		#if (exists $keepind{$inds[$i]}){		    
		print OUT " $word";
		#}
		$i++; #This was originally in the wrong place
	    }
	}
	if($i > 0){
	    print OUT "\n";
	    $nloc++;
	}
    }
    ## get individual ids
    elsif (m/^#CHROM/){
	   @line = split /\s+/;
	   foreach $word (@line){
	       if ($word =~ s/\.sorted\.bam$//){
		   # $word =~ s/.*aln(.*)/$1/; ### <<<< addition to strip leading path and "aln_" CCN modified here delete _
		 $word =~ s/aln//;   # simplified for 1.9
		 push (@inds, $word);
	       }
	   }
	   $nind = scalar @inds;
       }
}
close (OUT);

open (INDS, "> indiv_ids.txt") or die "Could not write the id file\n";
print INDS "ind,pop,taxon\n";
foreach $word (@inds){
    $word =~ m/([A-Za-z\d-]+).([A-Za-z\d-]+)/;
    print INDS "$word,$1,$2\n";
}
close (INDS);

$out =~ s/mpgl$/txt/;
open (INDS, "> header_mpgl$out") or die "Could not write the id file\n";
print INDS "$nind $nloc\n";
print INDS "discarded ind line\n";
close (INDS);

print "Number of loci: $nloc; number of individuals $nind\n";
