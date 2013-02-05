#!/usr/bin/perl

use strict;
#use warnings;
######## This script creates a binary version of the fasta file, coded in 0's and 1's for bases occuring in coding and noncoding regions. You must run this script before genomic_agg.
#This script uses the glob utility, which attempts to access folders or files in the specific file path. After accessing genomic features files (.gb), this script will create a file
#ending in _crs (mentioned above)


my@files = glob("/home/user/documents/*"); ##You will need to specify the file path
my$name;

foreach (@files){

if (/\/home\/user\/documents\/(\S+)/){  ##You will need to specify the file path
  $name = $1;
}

print "\n$name\n";
my$name_crs = $name."_crs";
my$namegb = $name.".gb";
my@start = ();
my@end = ();

my$L;

open (FILE, "/home/user/documents/$name/$namegb") or die;   ##You will need to specify the file path
open (CONTIGS, ">/home/user/documents/$name/$name_crs");    ##You will need to specify the file path

while (<FILE>){
    
    # loop grabs the length
    if (/LOCUS\s+\S+\s+(\d+)/){
	$L = $1;
    }
    # get protein coding gene coordinates from CDS tags in the data
    elsif (/CDS\s+(\d+)\S\S(\d+)/){
	push (@start,$1);
	push (@end,$2);
	}
    elsif (/CDS\s+complement\S(\d+)\S\S(\d+)\S/){
	push(@start,$1);
	push(@end,$2);
	}
    elsif (/CDS\s+join\W(\d+)\S\S(\d+)\W(\d+)\S\S(\d+)\S/){
	push(@start,$1);
	push(@end,$2);
	push(@start,$3);
	push(@end,$4);
	}
}

my$E = @end;
my$S = @start;
my$C = 0;
my$N = 0;

print "Length: $L : $E,$S\n";

while($end[$N]){
	$C += $end[$N] - $start[$N] + 1;
	$N++;
	}

print "Coding: $C, $N\n";
	my@start = sort{$a <=> $b} @start;
	my@end = sort{$a <=> $b} @end;

my$A = 0;
my$C = 1;
my$TC = 0; my$TN = 0;
my@code = ();
print "start: $start[$A], end: $end[$A]\n";
while ($C <= $L) {	
	if ($C < $start[$A]) {
	push(@code,0);
	$TN++;
	}
	elsif (($C >= $start[$A]) && ($C < $end[$A])) {
	push(@code,1);
	$TC++;
	}	
	elsif (($C == $end[$A]) && ($end[$A] != $end[$S-1])) {
	push(@code,1);
	$TC++;
	$A++;
	}
	elsif ($C == $end[$S-1]) {
	push(@code,1);
	$TC++;
	}
	elsif (($C > $end[$S-1]) && ($C <= $L)) {
	push(@code,0);
	$TN++;
	}
$C++;
}

print "Coding: $TC; Noncod: $TN; $A\n";
print CONTIGS @code;	
}
