#!/usr/bin/perl
# fork script for entropy 

use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);


foreach $k (2..4){
	foreach $ch (1..10){ 
		sleep 2;
		$pm->start and next;
		$out = "out_lmel_k$k"."_ch$ch".".hdf5";
		system "/uufs/chpc.utah.edu/common/home/u6000989/bin/entropy_mp -i lyc_entropy.mpgl -n 2 -m 1 -l 2000 -b 1000 -t 5 -k $k -o $out -q ldak$k"."_lmel.txt -s 20\n";
		$pm->finish;
	}
}
$pm->wait_all_children;


