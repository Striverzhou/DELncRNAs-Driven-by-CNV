#!/usr/bin/perl
use strict;
use warnings;

my $file="sample.tsv"; # sample.tsv is a file that contains only two columns of data, the id of segmentation file and the entity submitter id derived from TCGA database.

my %hash=();
my @normalSamples=();
my @tumorSamples=();

open(TF,"$file") or die $!;
while (my $Line=<TF>)
{
    next if($.==1);
    next if($Line=~/^\n/);
    chomp($Line);
    my @Arr=split(/\t/,$Line);
    my $file_name=$Arr[0];
	my $entity_submitter_id=$Arr[1];
	$file_name=$file_name . ".txt";
	#print "$file_name\n$entity_submitter_id\n";
	if((defined $file_name) && (-f $file_name))
    {
		my @idArr=split(/\-/,$entity_submitter_id);
		if($idArr[3]=~/^0/)
        {
            push(@tumorSamples,$entity_submitter_id);
        }
        else
        {
            push(@normalSamples,$entity_submitter_id);
        }        	
        open(RF,"$file_name") or die $!;
        while(my $line=<RF>)
        {
            next if($.==1);
            next if($line=~/^\n/);
            next if($line=~/^\_/);
            chomp($line);
            my @arr=split(/\t/,$line);
            ${$hash{$arr[6]}}{$entity_submitter_id}=$arr[5];
        }
        close(RF);
    }
}
close(TF);

open(WF,">cnvMatrix.txt") or die $!;
my $normalCount=$#normalSamples+1;
my $tumorCount=$#tumorSamples+1;

if($normalCount==0)
{
				print WF "id";
}
else
{
				print WF "id\t" . join("\t",@normalSamples);
}
print WF "\t" . join("\t",@tumorSamples) . "\n";

foreach my $key(keys %hash)
{
				print WF $key;
				foreach my $normal(@normalSamples)
				{
					      if(exists ${$hash{$key}}{$normal}){
					      	    my $copyNum=${$hash{$key}}{$normal};
					      	    $copyNum=2**(1+$copyNum);
					      	    my $copyOut=0;
					      	    if($copyNum<0.5){
					      	    	$copyOut=-2;
					      	    }
					      	    elsif($copyNum<1.5){
					      	    	$copyOut=-1;
					      	    }
					      	    elsif($copyNum>3.5){
					      	    	$copyOut=2;
					      	    }
					      	    elsif($copyNum>2.5){
					      	    	$copyOut=1;
					      	    }
								      print WF "\t" . $copyOut;
							  }
							  else{
							  	    print WF "\t0";
							  }
				}
				foreach my $tumor(@tumorSamples)
				{
					      if(exists ${$hash{$key}}{$tumor}){
					      	    my $copyNum=${$hash{$key}}{$tumor};
					      	    $copyNum=2**(1+$copyNum);
					      	    my $copyOut=0;
					      	    if($copyNum<0.5){
					      	    	$copyOut=-2;
					      	    }
					      	    elsif($copyNum<1.5){
					      	    	$copyOut=-1;
					      	    }
					      	    elsif($copyNum>3.5){
					      	    	$copyOut=2;
					      	    }
					      	    elsif($copyNum>2.5){
					      	    	$copyOut=1;
					      	    }
								      print WF "\t" . $copyOut;
							  }
							  else{
							  	    print WF "\t0";
							  }
				}
				print WF "\n";
}
close(WF);

print "normal count: $normalCount\n";
print "tumor count: $tumorCount\n";
