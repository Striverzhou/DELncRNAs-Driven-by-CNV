#!/usr/bin/perl

use strict;
use warnings;

my %hash=();

open(RF,"diffCNV_lncRNAs.txt") or die $!;
while(my $line=<RF>){
	next if($.==1);
	chomp($line);
	my @arr=split(/\t/,$line);
	$hash{$arr[0]}=1;
}
close(RF);

my $normalCount=0;
open(RF,"cnvMatrix.txt") or die $!;
while(my $line=<RF>){
	next if($.==1);
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1){
	  for(my $i=1;$i<=$#arr;$i++){
			my @sampleArr=split(/\-/,$arr[$i]);
			unless($sampleArr[3]=~/^0/){
				$normalCount++;
			}
		}
		next;
	}
	if(exists $hash{$arr[0]}){
	  my %cnvHash=();
	  for(my $i=($normalCount+1);$i<=$#arr;$i++){
      if($arr[$i] != 0){
    	  $cnvHash{$arr[$i]}++;
      }
	  }
	  MARK:foreach my $key (sort{$cnvHash{$b}<=>$cnvHash{$a}}(keys %cnvHash)){
	  	$hash{$arr[0]}=$key;
	  	#print $key . "\n";
	  	last MARK;
	  }
  }
}
close(RF);

open(RF,"lncRNAPos.txt") or die $!;
open(GENE,">lncRNA_Label.txt") or die $!;
open(SCATTER,">lncRNA_scatter.txt") or die $!;
print SCATTER "chromosome\tstart\tstop\tseg.mean\n";
print GENE "Chromosome\tchromStart\tchromEnd\tGene\n";
while(my $line=<RF>){
	next if($.==1);
	chomp($line);
	my @arr=split(/\t/,$line);
	if(exists $hash{$arr[0]}){
		print GENE "chr$arr[1]\t$arr[2]\t$arr[3]\t$arr[0]\n";
		print SCATTER "chr$arr[1]\t$arr[2]\t$arr[3]\t$hash{$arr[0]}\n";
		delete($hash{$arr[0]});
	}
}
close(SCATTER);
close(GENE);
close(RF);
