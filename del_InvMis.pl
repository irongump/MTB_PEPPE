#!/usr/bin/perl
use strict; 
use warnings;

#exclude loci at which N(>threshold) strains have ambiguous base
#by using fasta file which composed by all of the strain.

die "usage:perl $0 <cfa_file> <threshold>\n" if @ARGV==0;

open (FA, "<$ARGV[0]");

my $hpt = $ARGV[0];
$hpt =~ s/\.fa|.fasta|.fna//;

my $thr = 0.05;
if (scalar @ARGV == 2){
	$thr = $ARGV[1]/100;
}
my $per = $thr * 100;

open (HPT, ">$hpt"."_del_InvMisF$per.fa");
open (LOC, ">$hpt"."_del_InvMisF$per.loc");

my %seq;
my $str;
while(<FA>){
	chomp;
	$_ =~ s/\s+//g;
	if ($_ =~ /^>/){
		$str = $_;
	}elsif($_){
		$seq{$str} = $_; 
	}
}
close FA;

my %loc;
my %ref;
my $n = 0;

my $sn = keys %seq;
my @loci; 
my %len;
for my $loc (0..4411531){
	my %gt=();
	$gt{N}=0;
	foreach my $k(keys %seq){
		my $nuc = substr ($seq{$k}, $loc, 1);
        if($nuc =~ m/A|T|G|C|a|t|g|c|n|N/){
			$nuc= uc($nuc);    #capitalise "a t g c n"
       	}
        if($nuc =~ m/N|-|\?/){     
			$gt{N}++;
		}
		elsif(!exists $gt{$nuc}){
			$gt{$nuc}=1; 
		}else{
			$gt{$nuc}++;
		}
	}
	my $frqmis = $gt{N}/$sn;
	my $gtn = keys %gt; 
	my $valid_gtn = $gtn - (exists $gt{N} ? 1 : 0); 
	if($frqmis < $thr && $valid_gtn >= 2){    #nucletide type should be at least 2, if only 1, no useful information provided.
        push @loci, $loc;
	}
}

foreach my $k(keys %seq){
        my $j;
        $j=$k;
        $j=~s/\.cns//;
        print HPT $j, "\n";
        foreach (@loci){
                my $nuc = substr ($seq{$k},$_,1);
                print HPT $nuc;
        }
        print HPT "\n";
}

print scalar(@loci)." check loci length\n";

foreach(@loci){
    print LOC $_+1,"\n";
}

close HPT;
close LOC;
