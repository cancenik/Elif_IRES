#!/usr/bin/perl -w

use strict; 
use Statistics::Basic;
use constant REN_KNOCK => 0.5; 

## This script will remove renilla and firefly reps if
## Renilla firefly less than REN_KNOCK; 

my $sirna = "072814_heksiRNAreplicates.csv"; 
open (SIRNA, $sirna); 
while (<SIRNA>) {
    chomp;
    my @F = split(/,/); 
    my $gene = shift (@F);
    print "$gene"; 
    my $ren_sum = 0; 
    my $firefly_sum = 0; 
    my $num = 0; 
    my @REN;
    my @FIRE;
    for (my $i =1; $i < scalar(@F) ; $i = $i+2) { 
	if ($F[$i] <  REN_KNOCK ) {
	    push (@REN, $F[$i] );
	    push (@FIRE, $F[$i - 1 ] );
	    $num++; 
	    $ren_sum += $F[$i]; 
	    $firefly_sum += $F[$i - 1]; 
	    print "\t$F[$i - 1 ]\t$F[$i]"; 
	}
	else { 

	}
    }
    unless ( $num == 0 ) { 
	my $ren_mean = $ren_sum /$num; 
	my $firefly_mean = $firefly_sum / $num; 
	print "\t$firefly_mean\t$ren_mean\n";
	my $t1= mean(@REN); 
	print "$t1\n"; 
    }
    else { 
	print "\n"; 
    }

} 
