#!/usr/bin/perl
use strict; 
use subs;
use Compress::Zlib;
use List::Util qw(sum);
use List::MoreUtils qw(first_index);
use warnings;

##### Author: Josh Schmidt Jan 2015 @ MPI EVA #############
##### Version: 0.1:27.01.15.
###### From ARGV: 0: vcf file; 1: Pop_1.txt; 2: Pop_2.txt; chr; outdirectory 
###### Lets do some pair-wise FST calculations !!!! yeah, you know you always wanted to....
###### open GZ VCF file, ignore all meta-info lines. capture the header-line, create an array from this line
###### which will include the literal names of the all the individuals in the VCF. Open two population files,
###### that contain the pop name, then each individual, presented as a tab separated line. i.e. ONE LINE per file. 
###### We'll then use the first_index function from List::MoreUtils. This gives us the index for each of the individuals in the VCF.
###### We can now go line by line, counting ref alleles etc for both populations.
###### Where p = ref allele frequency, n = number of individuals with genotype data. IF the VCF is unfiltered (missingness etc) n can/will differ.
###### Calculate Pattersons FST as derived in Bhatia et al Genome Res. 2013. 23: 1514-1521.
###### Output: 	Chr	Site	p_Pop1	p_Pop2	n_Pop1	n_Pop2	nHat dHat FST
###### Example: 1	10642	0.99	0.98	61	61	-6.606647e-05	1.216071e-02	-5.432781e-03


my @indPop1=(); # will put the names of the individuals of pop1 in here
my @indPop2=(); # will put the names of the individuals of pop2 in here
my @pop1=(); #will put indexes in here
my @pop2=(); # same as line above, but pop2
my $pop1Name = (); # the first column in our Pop file.
my $pop2Name = ();
my $maxPossiblePop1 =(); # how many IDs are in the popfile.
my $maxPossiblePop2 = ();
my $chr = $ARGV[3];
my $path =  $ARGV[4];
my @pop1FileName = split (/[\/.]/, $ARGV[1]);
my @pop2FileName = split (/[\/.]/, $ARGV[2]);

# open pop file 1, get Popname (should be column 1), columns2..n contain the names.
open (POPFILE1, "<$ARGV[1]") or die "Couldn't open the population1 file";
while (<POPFILE1>){
	chomp $_;
	if ($.==1){
		@indPop1 = split ("\t", $_);
		$pop1Name = shift (@indPop1);
		}
}
# open pop file 2, get Popname (should be column 1), columns2..n contain the names.
open (POPFILE2, "<$ARGV[2]") or die "Couldn't open the population2 file";
while (<POPFILE2>){
	chomp $_;
	if ($.==1){
		@indPop2 = split ("\t", $_);
		$pop2Name = shift (@indPop2);
		}
}
close POPFILE1;
close POPFILE2;
my $outFH = join ('_', $pop1FileName[-2], "VS",$pop2FileName[-2],"chr",$chr,"bhatia_FST.txt");
# lets open the vcf. can be gz.! sweet....
my $gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[0]: $gzerrno\n" ;
my $gzLine;
open (OUTFILE, ">$path/$outFH")  or die "Couldn't open the output file";
print "\n$outFH";
print OUTFILE "Chr\tPos\tp_"."$pop1Name"."\tp_"."$pop2Name"."\tn_"."$pop1Name"."\tn_"."$pop2Name"."\tNhat\tDhat\tFST_bahtia\n";
while($gz->gzreadline($gzLine) > 0 ){
	chomp $gzLine;
	my @outputLine = (); # this guy will hold our FST related info. We print this to output.
	if ($gzLine =~ m/##/){;} # skip meta-info lines which are marked with double hash ##
	elsif ($gzLine =~ m/#CHROM/){ # single hash: this line contains that actual IDS of the individuals. We would like to find the indexes
		my @lineArray = split ("\t", $gzLine);
		for my $ind (0 .. $#indPop1){
      			my $index = first_index { $_ eq $indPop1[$ind] } @lineArray;
      			push (@pop1, $index);
		}
		for my $ind2 (0 .. $#indPop2){
      			my $index2 = first_index { $_ eq $indPop2[$ind2] } @lineArray;
      			push (@pop2, $index2);
		}
	# for your delectation. Lets print some info to screen, about the populations: Name, Indvs, and their index in the VCF lines. Who knows, maybe something weird has happened.
	print "$pop1Name";
	print "@indPop1";
	print "@pop1\n";
	print "$pop2Name";
	print "@indPop2";
	print "@pop2\n";
	$maxPossiblePop1 = scalar(@pop1); 	# need to define these two variables now. They will help with error catching in UNFILTERED files, when either no genotype data exists for a site, 
									#or a single individual in a pop. Both are big no-noes generating 0 divisors in calculations.
	$maxPossiblePop2 = scalar(@pop2);
	print "$maxPossiblePop1\n$maxPossiblePop2";
	} 
	else {
	
		# define a whole lot of empties. Is this how the professionals do it? I doubt it.
		my $refAllelePop1 = 0; # count of ref (0) alleles
		my $altAllelePop1 = 0; # count of alt (1) alleles
		my $noGenotypePop1 = 0; # count of missing (.) 'alleles'
		my $refAllelePop2 = 0;
		my $altAllelePop2 = 0;
		my $noGenotypePop2 = 0;


		#split the VCF input line on tab.
		my @lineArray = split ("\t", $gzLine);
		
		#now, in inelegant fashion are two FOR loops
		#iterating through each index array for each Population
		# we look at the genotype calls for each relevant individual, count homo 0 as 2 ref
		# 0/1 or 1/0 as 1 ref and 1 alt, 1/1 as 2 alt and any matching . (dot) as 2 missing.
		# 0/1 or 1/0 are also tallied once as an observed het.
		# what an ind is.... 0/0:3:3:119:0:0:0,-0.664216,-3.5708 (in the chimp case).
		# as genotype is / divided remember to escape \/
		# also for phased vcf, | and not / is the genotype divider. updated now to search for either case.
		
		for my $i (0 .. $#pop1){
				my $slot = $pop1[$i];
				my @indData = split(":", $lineArray[$slot]);
				# ref homozygous
				if ($indData[0] =~ m/0\/0|0\|0/){
					$refAllelePop1 = $refAllelePop1+2;
					}
				# heterozygous
				elsif($indData[0] =~ m/0\/1|1\/0|0\|1|1\|0/){
					$refAllelePop1 = $refAllelePop1+1; ## $refAllelePop1++; could work to add one, but this makes it more explicit about how allele freqs are being derived.
					$altAllelePop1 = $altAllelePop1+1;
					}
				# alt homozygous
				elsif ($indData[0] =~ m/1\/1|1\|1/){
					$altAllelePop1 = $altAllelePop1+2;
					}
				# naughty - no genotype data
				elsif ($indData[0]!~ m/1|0/g){
					$noGenotypePop1 = $noGenotypePop1+1;
					}			
				else {;}
				}
		# same as the sequence of loops immediately above, just for pop2 this time.
		for my $j (0 .. $#pop2){
				my $slot2 = $pop2[$j];
				my @indData2 = split(":", $lineArray[$slot2]);
				if ($indData2[0] =~ m/0\/0|0\|0/){
					$refAllelePop2 = $refAllelePop2+2;
					}
				elsif($indData2[0] =~ m/0\/1|1\/0|0\|1|1\|0/){
					$refAllelePop2 = $refAllelePop2+1;
					$altAllelePop2 = $altAllelePop2+1;
					}
				elsif ($indData2[0] =~ m/1\/1|1\|1/){
					$altAllelePop2 = $altAllelePop2+2;
					}
				elsif ($indData2[0] !~ m/1|0/g){ # use global string match here, just to make sure we catch the non 0/1 (.)
					$noGenotypePop2 = $noGenotypePop2+1;
					}			
				else {;}
				}

#### some of these thoughts have been implemented....
# think about a suitable flow control before too many class are done: i.e. pop is totally missing
# monomorphic across pops. we dont want undef errors etc. beware bugs! 
# oh and a check does total alleles and missing alleles equal 2 times the array length? It should. 
# Otherwise there is an error######

# note, if we are running on an unfiltered vcf, there could be sites with all genotypes missing, or only 1 individual from each pop that has genotype data. Both cases are undef, and 
# perl will abort the first time such a case is seen. so lets skip them, and output a line saying why it was skipped: "All_missing_or_single"	
		if (($maxPossiblePop1 - $noGenotypePop1)<$maxPossiblePop1 or ($maxPossiblePop2 - $noGenotypePop2)<$maxPossiblePop2){;}
		# note WC 1984 suggest that for invariant alleles, FST should be undefined - that is what we have here.
		# however this is the change for PBS4. where we need to come out with a vector of all pairwise distances. FST=0 is therefore informative in this case!!! so put out a dummy line!
		elsif ($altAllelePop1 == 0 and $altAllelePop2 == 0){
		  push (@outputLine, $lineArray[0]);
		  push (@outputLine, $lineArray[1]);
		  push (@outputLine,sprintf '%.2f', 1.00);
		  push (@outputLine,sprintf '%.2f', 1.00);
		  my $pop1N=($refAllelePop1+$altAllelePop1)/2; #pop1 sample size. equals total number of alleles divided by 2, because diploid
		  push (@outputLine, $pop1N);
		  my $pop2N=($refAllelePop2+$altAllelePop2)/2; #pop2 sample size
		  push (@outputLine, $pop2N);
		  push (@outputLine, sprintf '%.6e', 0);
		  push (@outputLine, sprintf '%.6e',0);
		  push (@outputLine, sprintf '%.6e',0);
		  #print OUTFILE join("\t",@outputLine);
		  my $output = join("\t",@outputLine);
                  print OUTFILE "$output\n";  
		}
		# more rarely both pops are fixed for the alternate allele (though it is all relative)
		elsif ($refAllelePop1 == 0 and $refAllelePop2 == 0){
		  push (@outputLine, $lineArray[0]);
		  push (@outputLine, $lineArray[1]);
		  push (@outputLine,sprintf '%.2f', 0.00);
		  push (@outputLine,sprintf '%.2f', 0.00);
		  my $pop1N=($refAllelePop1+$altAllelePop1)/2; #pop1 sample size. equals total number of alleles divided by 2, because diploid
		  push (@outputLine, $pop1N);
		  my $pop2N=($refAllelePop2+$altAllelePop2)/2; #pop2 sample size
		  push (@outputLine, $pop2N);
		  push (@outputLine, sprintf '%.6e', 0);
		  push (@outputLine, sprintf '%.6e', 0);
		  push (@outputLine, sprintf '%.6e', 0);
		  #print OUTFILE join("\t",@outputLine);
		  my $output = join("\t",@outputLine);
                  print OUTFILE "$output\n";
		}
		
	# you get to here, you should be in the clear....
		else {
			my $alleleFreqPop1 = $refAllelePop1/($refAllelePop1+$altAllelePop1); # alleleFreq - number of ref alleles/(total number of alleles)
			my $alleleFreqPop2 = $refAllelePop2/($refAllelePop2+$altAllelePop2);
		# put this data in our output array	
					#push onto the outputLine array Chr and site. these are the first two columns.
			push (@outputLine, $lineArray[0]);
			push (@outputLine, $lineArray[1]);
			push (@outputLine,sprintf '%.2f', $alleleFreqPop1);
			push (@outputLine,sprintf '%.2f', $alleleFreqPop2);
			 my $pop1N=($refAllelePop1+$altAllelePop1)/2; #pop1 sample size. equals total number of alleles divided by 2, because diploid
			push (@outputLine, $pop1N);
			my $pop2N=($refAllelePop2+$altAllelePop2)/2; #pop2 sample size
			push (@outputLine, $pop2N);
		 
			my $pop1Acount=$altAllelePop1+$refAllelePop1; ## count of the alleles, fr each population
			my $pop2Acount=$altAllelePop2+$refAllelePop2;
			##
			my $h1=($altAllelePop1*($pop1Acount-$altAllelePop1))/($pop1Acount*($pop1Acount-1));
			my $h2=($altAllelePop2*($pop2Acount-$altAllelePop2))/($pop2Acount*($pop2Acount-1));

			my $nHat=( (($altAllelePop1/$pop1Acount )-($altAllelePop2/$pop2Acount))*(($altAllelePop1/$pop1Acount )-($altAllelePop2/$pop2Acount))  ) -($h1/$pop1Acount)-($h2/$pop2Acount);
			my $dHat=$nHat+$h1+$h2;
			my $fst=$nHat/$dHat;
			push (@outputLine, sprintf '%.6e', $nHat);
			push (@outputLine, sprintf '%.6e', $dHat);
			push (@outputLine, sprintf '%.6e', $fst);
			#print OUTFILE join("\t",@outputLine);
			my $output = join("\t",@outputLine);
			print OUTFILE "$output\n";
		      }	
	      }
      }
close OUTFILE;

