#!/usr/bin/perl

use strict;
use warnings;

my $normalCount=0;
my $tumorCount=0;
open(RF,"cnvMatrix.txt") or die $!;
open(WF,">chiInput.txt") or die $!;
print WF "Gene\tnormalCNV\tnormalWild\ttumorCNV\ttumorWild\n";
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1){
		for(my $i=1;$i<=$#arr;$i++){
			my @sampleArr=split(/\-/,$arr[$i]);
			if($sampleArr[3]=~/^0/){
				$tumorCount++;
			}
			else{
				$normalCount++;
			}
		}
		next;
	}
	my $normalWild=0;
	my $tumorWild=0;
	for(my $i=1;$i<=$normalCount;$i++){
		if($arr[$i]==0){
		  $normalWild++;
		}
	}
	for(my $i=($normalCount+1);$i<=$#arr;$i++){
		if($arr[$i]==0){
		  $tumorWild++;
	  }
	}
	my $normalCNV=$normalCount-$normalWild;
	my $tumorCNV=$tumorCount-$tumorWild;
	unless(($normalCNV==0) && ($tumorCNV==0)){
	  print WF "$arr[0]\t$normalCNV\t$normalWild\t$tumorCNV\t$tumorWild\n";
	}
}
close(WF);
close(RF);
