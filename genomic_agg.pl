#!/usr/local/bin/perl 

#use warnings;
#use strict; this has been blocked to allow the script to create files 'on the fly'

#The following script does several tasks:

#1. It calculates the average measure of Morisita's Index of aggregation for Pyrimidines within coding and noncoding regions,
#and for Purines within coding and noncoding regions.
#2. It creates 100 randomized versions of each genome and performs the above operations. 
#3. It creates files to hold every estimate of aggregation from real genomes, and 100 measurements of average aggregation from randomized genomes.
#4. It creates a plain text document of the genome in binary code representing coding and noncoding regions.

#The glob utility uses a specific path to access a folder containing subfolders of two genome files. One file is a genomic features file (.gb) conforming to the 
#NCBI microbial genome database format, the other file is the .fasta file corresponding to the feature file.
 
#The script reads the fasta file and recodes it in binary form for purines (A,G) and pyrimidines (C,T). The script then reads the genome features file (.gb) and
#creates a plain text document (ending in _crs) to hold a binary representation of coding and noncoding regions. It then reads both binary files and separates 
#the binary sequence file (purines and pyrimidines) into two arrays. One represents purines or pyrimidines within coding regions. The other represents 
#purines or pyrimidines in noncoding regions. 

#A window (containing number of quadrats and number of bases within quadrats) is then moved along the arrays using while loops. Morisita's I is calculated for each window.
#Since each array (coding and noncoding) will contain thousands of windows, thousands of measurements of Morisita's I are averaged.
#The script then creates 100 randomizations per genome and performs the above operations (lines 22-23) on each.

#Upon completion, each species folder will contain 11 files: a features file (.gb), a .fasta file, a binary coding/noncoding sequence file (_crs),
#and 8 files holding measurements of Morosita's I.
# The following abbreviations denote files the script will create:
#C_pur: purines within coding regions
#C_pyr: pyrimidines within coding regions
#NC_pur: purines within noncoding regions
#NC_pyr: pyrimidines within noncoding regions
#RC_pur: purines within randomized coding regions
#RC_pyr: pyrimidines within randomized coding regions
#RN_pur: purines within randomized noncoding regions
#RN_pyr: pyrimidines within randomized noncoding regions

##                FEATURE FILES, FASTA FILES, AND FILE FOLDERS MUST HAVE IDENTICAL NAMES!

##                FILE PATHS WILL LIKELY NEED TO BE SPECIFIED FOR YOUR FILE SYSTEM ON THESE LINES:
#                  55, 58, 65, 88, 94, 95, 216, 217, 399, 400, 403, 537, 538
# 
#
#
#
#
use List::Util 'shuffle';  #a utility to randomize coding and noncoding regions
my$name;
my$quads = 10;  #Determines the size of the window (e.g. $quads = 10, $bases = 10 represents a window of 10 quadrats and 10 bases per quadrat).	
my$bases = 10;  #The maximum occupancy of a quadrat is $bases.

my($Cavg,$NCavg,$countC,$countN,$simcount,$codreg,$nonreg,$countR,$Ni,$end,$Ci);
my@miniN = ();my@miniC = ();my@mIC = ();my@mIN = ();my@aQS = 10;

my@files = glob("/home/kenlocey/bioinfo/genomes/NOT_DONE/*"); #This file path should be rescripted to the user's system.

foreach (@files){
if (/\/home\/kenlocey\/bioinfo\/genomes\/NOT_DONE\/(\S+)/){  #This file path should be rescripted to the user's system
	$name = $1;
}
my$name_crs = $name . "_crs";
print " 	           $name\n";


open (FILE, "NOT_DONE/$name/$name.fasta") or die;  #This file path should be rescripted to the user's system
my$file = <FILE>;
open (AGG, ">agg") or die;

while ($file = <FILE>) { 
	chomp($file);	
	$file =~ s/(<br>)//g;  # binary coding assignment
	$file =~ s/c/1/ig;
	$file =~ s/t/1/ig;
	$file =~ s/a/0/ig;
	$file =~ s/g/0/ig; 
	print AGG $file; 
}
close(FILE);
close(AGG);
#print "A "; ##prints to screen to indicate progress (tends to 'choke' here on genomes exceeding 13mbp)

open(IN, "agg");
my$seq = <IN>;
my@seq2 = split(//,$seq); #split is the utility called to 'split' the binary sequence file into coding and noncoding arrays
my$seqsize = @seq2;
#print "seqsize: $seqsize\n";

open(IN2, "NOT_DONE/$name/$name_crs") or die;
my$cod = <IN2>;
my@code = split(//,$cod); #splits the binary file of coding and noncoding seqs into coding and noncoding arrays  
my$codesize = @code;

#print "B\n"; ##prints to screen to indicate progress
open(OUTN, ">NOT_DONE/$name/N_pyr");
open(OUTC, ">NOT_DONE/$name/C_pyr");


my$count=0;my$counts=0;my$C=0;my$SS=0;my$avg=0;my$mtotalC=0;my$N=0;my$nW=0;my$qtotal=0;my$wtotal=0;my$pos=0;my$count1=0;my$count2=0;my$mI=0;my$mtotalN=0;my$mavgN=0;my$mavgC=0;my$mTN=0;my$mTC=0;my$mavgRC=0;my$mavgRN=0;my$mtotalRN=0;my$mtotalRC=0;

#######################################################################
####################################################################### FOR PYRIMIDINES  #########################################################################

$end = $bases-1;
@miniN = ();
@miniC = ();

while($C < $codesize) {
if ($code[$C] == 0) {
push(@miniN,$seq2[$C]);
}
else {
push(@miniC,$seq2[$C]);	
}
$C++;
}
$C = 0;
$Ni = @miniN;
$N = @miniN;
	#print "Noncoding: $N ; ";
	$nW = ($N-($N%100))/100;
	while ($count1 < $nW) { 
		while ($count2 < 10){ 
			#print "$mini[$pos]\n";
	  		while ($pos <= $end) { 
				 if ($miniN[$pos] == 1) {    
				$qtotal++;                   
				$aQS[$count2] = $qtotal;  
				$wtotal++; 
				#print "$wtotal\n";   
				}
			$pos++;  
			}
		$end = $end + 10;  
		$count2++;         
		$qtotal = 0;       
		}
	
	#print "total 1's: $wtotal\n"; 	
	$avg = $wtotal/10;
	#print "average: $avg\n";
		$count = 0;
		while ($count < 10) {  
		$SS = $SS + (($aQS[($count2 - $count -1)] - $avg)**2);
		$count++;
		}
		$count = 0;


	$mI = ($wtotal/($wtotal-1))*(1/$avg)*((($SS/($quads-1))/$avg) + $avg - 1);	
	print OUTN "$mI\n";
	$mtotalN = $mtotalN + $mI;
	$count1++;
	$wtotal=0;$count2=0;$SS=0;
		 
}##          END OF WINDOW LOOP

	@miniN = ();	
	$mavgN = $mtotalN/$nW;
	
	$count1=0;$pos=0;$mI=0;$mtotalN=0;$nW=0;$N=0;$count=0;$avg=0;
	$end = 9; 

$Ci = @miniC;	
$N = @miniC;

	$nW = ($N-($N%100))/100;
	
	while ($count1 < $nW) { 
		while ($count2 < 10){ 
	  		while ($pos <= $end) {  
	    			if ($miniC[$pos] == 1) {    
				$qtotal++;
				                   
				$aQS[$count2] = $qtotal;  
				$wtotal++;    
				}
			$pos++;
			  
			}
		$end = $end + 10;  
		$count2++;         
		$qtotal = 0;       
		}
	#print "total 11's: $wtotal\n"; 	
	$avg = $wtotal/10;
	$count = 0;
		while ($count < 10) {  
		$SS = $SS + (($aQS[($count2 - $count -1)] - $avg)**2);
		$count++;
		}
		$count = 0;
	if($wtotal == 1) {		
	print "@mini\n$avg\n";
	}
	$mI = ($wtotal/($wtotal-1))*(1/$avg)*((($SS/($quads-1))/$avg) + $avg - 1);	
	print OUTC "$mI\n";	
	$mtotalC = $mtotalC + $mI;
	
	$count1++;
	$wtotal=0;$count2=0;$SS=0;$mI=0;
	
}##          END OF WINDOW LOOP
	@miniC= ();
	$mavgC = $mtotalC/$nW;
	$end = 9; 

$mtotalN=0;$mtotalC=0;
close(OUTN);
close(OUTC);

################################################################################### For Random Pyrimidines
#########################################################################################################################################################################################
                                                                                
#print "Random Pyrimidines\n"; ##prints to screen to indicate progress

open(OUT4C, ">NOT_DONE/$name/RC_pyr");
open(OUT4N, ">NOT_DONE/$name/RN_pyr");

$count=0;$C=0;$avg=0;$N=0;$nW=0;$qtotal=0;$pos=0;$count1=0;$mI=0;$simcount=0;


while($simcount < 100){
print "$name PYR $simcount\n";
#print "count1: $count1\n";
$codreg=0;$nonreg=0;$countR=0;
my@AoA = ();
my@mini2=();

while($countR < $codesize) {
	if($code[$countR] == 1){
	while(($code[$countR] == 1) && ($countR < $codesize)){
	push(@mini2,$seq2[$countR]);
	$countR++;	
	}
	@shuf = shuffle(@mini2);
	push (@AoA, @shuf);
	@mini2 = ();
	$codreg++;	
	}
	else{
	while(($code[$countR] == 0) && ($countR < $codesize)){
	push(@mini2,$seq2[$countR]);
	$countR++;
	}	
	@shuf = shuffle(@mini2);
	push (@AoA, @shuf);
	@mini2 = ();
	$nonreg++;	
	}
}
$AoAsize = @AoA;
@mini2 =();
$countR = 0;
#print "AoAsize: $AoAsize\n";
#$N = @miniN;	
#print "miniN: $N\n";

$C = 0;
while($C < $codesize) {
if ($code[$C] == 0) {
push(@miniN,$AoA[$C]);
}
else {
push(@miniC,$AoA[$C]);	
}
$C++;#print "N: $N\n";
}
$C = 0;
@AoA = ();
#print "$simcount\n";
$N = @miniN;
#print "miniN: $N\n";
	$nW = ($N-($N%100))/100;
#	print "number of windows: $nW\n";
	while ($count1 < $nW) { 
		while ($count2 < 10){ 
			#print "$mini[$pos]\n";
	  		while ($pos <= $end) { 
				if ($miniN[$pos] == 1) {    
				$qtotal++;                   
				$aQS[$count2] = $qtotal;  
				$wtotal++; 
				#print "$wtotal\n";   
				}
			$pos++;  
			}
		$end = $end + 10;  
		$count2++;         
		$qtotal = 0;       
		}
	
	#print "total 1's: $wtotal\n"; 	
	$avg = $wtotal/10;
	#print "average: $avg\n";
		$count = 0;
		while ($count < 10) {  
		$SS = $SS + (($aQS[($count2 - $count -1)] - $avg)**2);
		$count++;
		}
		$count = 0;

	$mI = ($wtotal/($wtotal-1))*(1/$avg)*((($SS/($quads-1))/$avg) + $avg - 1);	
	$mtotalRN = $mtotalRN + $mI;
	$count1++;	
	$mI=0;$wtotal=0;$count2=0;$SS=0;
		
}##          END OF WINDOW LOOP
	#print "windows: $nW\n";	
	$mavgRN = $mtotalRN/$nW;
	@miniN= ();	
	push (@mIN,$mavgRN);
	print OUT4N "$mavgRN\n";	
	$mTN = $mTN + $mavgRN;
	$mavgRN=0;$mtotalRN=0;$count1=0;$pos=0;$nW=0;$N=0;$avg=0;
	$end = 9; 
		
$N = @miniC;
	#print "$N\n";
	$nW = ($N-($N%100))/100;
	#print "number of windows: $nW\n";
	while ($count1 < $nW) { 
		while ($count2 < 10){ 
	  		while ($pos <= $end) {  
	    			if ($miniC[$pos] == 1) {    
				$qtotal++;                   
				$aQS[$count2] = $qtotal;  
				$wtotal++;    
				}
			$pos++;  
			}
		$end = $end + 10;  
		$count2++;         
		$qtotal = 0;       
		}
	#print "total 11's: $wtotal\n"; 	
	$avg = $wtotal/10;
	$count = 0;
		while ($count < 10) {  
		$SS = $SS + (($aQS[($count2 - $count -1)] - $avg)**2);
		$count++;
		}
		$count = 0;
	#if($wtotal == 1) {		
	#print "@miniC\n$avg\n";
	#}
	$mI = ($wtotal/($wtotal-1))*(1/$avg)*((($SS/($quads-1))/$avg) + $avg - 1);	
	$mtotalRC = $mtotalRC + $mI;
	$count1++;
	$mI=0;$wtotal=0;$count2=0;$SS=0;
	
}##          END OF WINDOW LOOP
	@miniC= ();
	$mavgRC = $mtotalRC/$nW;
	print OUT4C "$mavgRC\n";
	push (@mIC,$mavgRC);
	$mTC = $mTC + $mavgRC;
	$mavgRC=0;$mtotalRC=0;$C=0;$avg=0;$N=0;$nW=0;$qtotal=0;$pos=0;$count1=0;
	$end = 9; 

$simcount++;

}
$Cavg = $mTC/$simcount;
$NCavg = $mTN/$simcount;

$Ccount=0;$Ccounts=0;
foreach(@mIC){
	if ($mavgC > $_){
	$Ccount++;
	}
	if ($mavgC < $_){
	$Ccounts++;
	}
}
$Ncount=0;$Ncounts=0;
foreach(@mIN){
	if ($mavgN > $_){
	$Ncount++;
	}
	if ($mavgN < $_){
	$Ncounts++;
	}
}

open (REPORT, ">>report2_2-OCT-2010") or die;
print REPORT "   	           $name\nQ$quads,B$bases  Seq:$seqsize Cod:$codesize C: $Ci N: $Ni\n\nPYR: mI Coding: $mavgC; mI Noncoding: $mavgN\nRand PYR: C: $Cavg, B$Ccount,S$Ccounts  N: $NCavg, B$Ncount,S$Ncounts\n";
close(REPORT);
close(OUT4C);
close(OUT4N);
close(IN);
@mIC = ();
@mIN = ();
$mavgN=0;$mavgC=0;


##########################################################################  For Purines
###################################################################################################################################

open(OUT2C, ">NOT_DONE/$name/C_pur");
open(OUT2N, ">NOT_DONE/$name/N_pur");


open (FILEE, "NOT_DONE/$name/$name.fasta") or die;
my$filee = <FILEE>;
open (AGG, ">agg") or die; 
while ($filee = <FILEE>) { 
	chomp($filee);	
	$filee =~ s/(<br>)//g;  # substituting 0's for pyrimidines and 1's for purines
	$filee =~ s/c/0/ig;
	$filee =~ s/t/0/ig;
	$filee =~ s/a/1/ig;
	$filee =~ s/g/1/ig; 
	print AGG $filee; 
}
close(FILEE);
close(AGG);

open(IN3, "agg");
my$seqq = <IN3>;
my@seq3 = split(//,$seqq);
my$seqsize = @seq3;
#print "$seqsize\n";

$wtotal=0;
my$end = $bases-1;
my@mini = ();

while($C < $codesize) {
if ($code[$C] == 0) {
push(@miniN,$seq3[$C]);
}
else{
push(@miniC,$seq3[$C]);	
}
$C++;
}

$N = @miniN;
	#print "$N\n";
	$nW = ($N-($N%100))/100;
	#print "windows: $nW\n";
	#print "number of windows: $nW\n";
	while ($count1 < $nW) { 
		while ($count2 < 10){ 
			#print "$mini[$pos]\n";
	  		while ($pos <= $end) { 
				 
	    			if ($miniN[$pos] == 1) {    
				$qtotal++;                   
				$aQS[$count2] = $qtotal;  
				$wtotal++; 
				#print "$wtotal\n";   
				}
			$pos++;  
			}
		$end = $end + 10;  
		$count2++;         
		$qtotal = 0;       
		}
	
	#print "total 1's: $wtotal\n"; 	
	$avg = $wtotal/10;
	#print "average: $avg\n";
		
		
		$count = 0;
		while ($count < 10) {  
		$SS = $SS + (($aQS[($count2 - $count -1)] - $avg)**2);
		$count++;
		}
		$count = 0;


	$mI = ($wtotal/($wtotal-1))*(1/$avg)*((($SS/($quads-1))/$avg) + $avg - 1);	
	print OUT2N "$mI\n";
	$mtotalN = $mtotalN + $mI;
	
	$count1++;
	$wtotal=0;$count2=0;$SS=0;
		 
}##          END OF WINDOW LOOP
	$mavgN = $mtotalN/$nW;
	$end = 9; 
	$count1=0;$pos=0;$mI=0;$nW=0;$N=0;$count=0;$avg=0;$C=0;
	@miniN= ();


$N = @miniC;
	#print "$N\n";
	$nW = ($N-($N%100))/100;
	#print "number of windows: $nW\n";
	while ($count1 < $nW) { 
		while ($count2 < 10){ 
	  		while ($pos <= $end) {  
	    			if ($miniC[$pos] == 1) {    
				$qtotal++;                   
				$aQS[$count2] = $qtotal;  
				$wtotal++;    
				}
			$pos++;  
			}
		$end = $end + 10;  
		$count2++;         
		$qtotal = 0;       
		}
	#print "total 11's: $wtotal\n"; 	
	$avg = $wtotal/10;
	$count = 0;
		while ($count < 10) {  
		$SS = $SS + (($aQS[($count2 - $count -1)] - $avg)**2);
		$count++;
		}
		$count = 0;
	if($wtotal == 1) {		
	print "@mini\n$avg\n";
	}
	$mI = ($wtotal/($wtotal-1))*(1/$avg)*((($SS/($quads-1))/$avg) + $avg - 1);	
	print OUT2C "$mI\n";
	$mtotalC = $mtotalC + $mI;
	
	$count1++;
	$wtotal=0;$count2=0;$SS=0;
	

}##          END OF WINDOW LOOP
	$mavgC = $mtotalC/$nW;
	$end = 9;
	@miniC = ();
	$count1=0;$pos=0;$mI=0;$nW=0;$N=0;$count=0;$avg=0;$C=0;

close(OUT2C);
close(OUT2N);

##################################################################################### FOR RANDOM PURINES
################################################################################################################################################################
#print "\n$seqsize\n";
open(OUT3C, ">NOT_DONE/$name/RC_pur");
open(OUT3N, ">NOT_DONE/$name/RN_pur");

$simcount=0;$mTN=0;$mTC=0;$qtotal=0;

while($simcount < 100){
print "$name PUR $simcount\n";

$codreg=0;$nonreg=0;$countR=0;
@AoA = ();
@mini2=();
#$AoAsize = @AoA;
#print "AoAsize: $AoAsize\n";
#print "codesize: $codesize\n";	
while($countR < $codesize) {
	if($code[$countR] == 1){
	while(($code[$countR] == 1) && ($countR < $codesize)){
	push(@mini2,$seq3[$countR]);
	$countR++;	
	}
	@shuf = shuffle(@mini2);
	push (@AoA, @shuf);
	@mini2 = ();
	$codreg++;	
	}
	else{
	while(($code[$countR] == 0) && ($countR < $codesize)){
	push(@mini2,$seq3[$countR]);
	$countR++;
	}	
	@shuf = shuffle(@mini2);
	push (@AoA, @shuf);
	@mini2 = ();
	$nonreg++;	
	}
}
$AoAsize = @AoA;
#print "AoAsize: $AoAsize\n";

while($C < $codesize) {
if ($code[$C] == 0) {
push(@miniN,$AoA[$C]);
}
elsif($code[$C] == 1) {
push(@miniC,$AoA[$C]);	
}
$C++;
}

my$countR=0;
my@AoA = ();
my@mini2=();

print "$simcount\n";

$N = @miniN;
#print "@miniN\n";
	#print "N: $N\n";
	$nW = ($N-($N%100))/100;
	while ($count1 < $nW) { 
		while ($count2 < 10){ 
	  		while ($pos <= $end) { 
				if ($miniN[$pos] == 1) {    
				$qtotal++;                   
				$aQS[$count2] = $qtotal;  
				$wtotal++; 
				}
			$pos++;  
			}
		$end = $end + 10;  
		$count2++;         
		$qtotal = 0;       
		}
	
	#print "total 1's: $wtotal\n"; 	
	$avg = $wtotal/10;
	#print "average: $avg\n";
				
		$count = 0;
		while ($count < 10) {  
		$SS = $SS + (($aQS[($count2 - $count -1)] - $avg)**2);
		$count++;
		}
		$count = 0;

	$mI = ($wtotal/($wtotal-1))*(1/$avg)*((($SS/($quads-1))/$avg) + $avg - 1);	
	$mtotalRN = $mtotalRN + $mI;
	$count1++;
	$wtotal=0;$count2=0;$SS=0;
		 
}##          END OF WINDOW LOOP
	#print "windows: $nW\n";	
	$mavgRN = $mtotalRN/$nW;
	push (@mIN,$mavgRN);
	print OUT3N "$mavgRN\n";	
	
	$mTN = $mTN + $mavgRN;
	$mtotalRN=0;$mavgRN=0;$count1=0;$pos=0;$mI=0;$nW=0;$N=0;$avg=0;
	$end = 9; 
	@miniN= ();
	
$N = @miniC;
	#print "$N\n";
	$nW = ($N-($N%100))/100;
	#print "number of windows: $nW\n";
	while ($count1 < $nW) { 
		while ($count2 < 10){ 
	  		while ($pos <= $end) {  
	    			if ($miniC[$pos] == 1) {    
				$qtotal++;                   
				$aQS[$count2] = $qtotal;  
				$wtotal++;    
				}
			$pos++;  
			}
		$end = $end + 10;  
		$count2++;         
		$qtotal = 0;       
		}
	#print "total 11's: $wtotal\n"; 	
	$avg = $wtotal/10;
	$count = 0;
		while ($count < 10) {  
		$SS = $SS + (($aQS[($count2 - $count -1)] - $avg)**2);
		$count++;
		}
		$count = 0;
	if($wtotal == 1) {		
	print "@mini\n$avg\n";
	}
	$mI = ($wtotal/($wtotal-1))*(1/$avg)*((($SS/($quads-1))/$avg) + $avg - 1);	
	$mtotalRC = $mtotalRC + $mI;
	$count1++;
	$wtotal=0;$count2=0;$SS=0;
	
}##          END OF WINDOW LOOP
	@miniC= ();
	$mavgRC = $mtotalRC/$nW;
	push (@mIC,$mavgRC);
	print OUT3C "$mavgRC\n";
	$mTC = $mTC + $mavgRC;
	$end = 9; 

$mavgRC=0;$mtotalRC=0;$C=0;$avg=0;$N=0;$nW=0;$pos=0;$count1=0;$mI=0;$qtotal=0;
$simcount++;
}
close(IN3);
close(IN2);
close(OUT3C);
close(OUT3N);
my$Cavg = $mTC/$simcount;
my$NCavg = $mTN/$simcount;

$Ccount=0;$Ccounts=0;
foreach(@mIC){
	if ($mavgC > $_){
	$Ccount++;
	}
	if ($mavgC < $_){
	$Ccounts++;
	}
}
$Ncount=0;$Ncounts=0;
foreach(@mIN){
	if ($mavgN > $_){
	$Ncount++;
	}
	if ($mavgN < $_){
	$Ncounts++;
	}
}
@mIN = ();
@mIC = ();

open (REPORT, ">>report2_2-OCT-2010") or die;
print REPORT "\nPUR: mI Coding: $mavgC; mI Noncoding: $mavgN\n";
print REPORT "R PUR: mI C: $Cavg, B$Ccount,S$Ccounts N: $NCavg, B$Ncount,S$Ncounts\n\n\n";
close(REPORT);
}
################################################################# The End ###############################


