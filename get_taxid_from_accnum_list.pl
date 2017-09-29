#By Rodrigo G. 
#Created: Sept 1st, 2016
#Update 1: Jan 19th, 2017
#The script was modified to avoid sorting (large DBs have problems with large acc files)
#It now uses a single column table containing a list of acc_nums
#Last update Sept 18, 2017
#Some variables names were changed
#!/usr/bin/perl
###########################IMPORTANT###########################
#As of September 2016, GIs will be phased out so we are changing our procedure to use Acc_Num ids instead. However, some variable names remain as "gi"
#This program reads a list of gis, a series of taxonomy dmps and gets the corresponding taids for each entry
#It was designed to parse dmp files: nodes.dmp, and names.dmp
#The input files are:
#	1.- List of accession numbers (1 column only)
#	2.- Accession2taxid.gz file containing the Accession and taxids
#In order to avoid using too much memory, the program stores the input accnums into memory and passes only once through the whole acc2taxid file
#The output file is a 2 column tab file containing an Acc_Num and taxid (terminal nodes)
#run as follows:
#	perl get_taxid_from_accnum_list.pl <acc_num_list> <Accession2taxid.gz>
use strict;
use warnings;
my $file = $ARGV[0];
my $acc2taxid = $ARGV[1];
my @line=();
my %acc=();
my ($i,$j,$items)="";

load_accnum();#first we read the gi list, keeping only derep indexed list
parse_acc2taxid();#now we parse the gi taxid dmp file (this is the slowest part
print_missing();#Print the ones that failed at the end.

sub load_accnum{#read the gi list, keeping only derep indexed list
open(ACC,"$file")|| die "Couldnt read acc list file\n";
	while (<ACC>){
#  		print $_;
		chomp ($_);
		$acc{$_}=1 if ! exists $acc{$_};#this also avoids repetition
	}
	close(ACC);
	$items=scalar(keys(%acc));
	print "Total different acc in input file: $items\n";
}
sub parse_acc2taxid{#parse the gi taxid dmp file (this is the slowest part
	$i=0;
	open(OUT,">$file.acc_taxid.txt")|| die "Couldn't create acc2taxid file\n";
	open (FH, "gunzip -c $acc2taxid | ")|| die "Couldn't find dmp file $acc2taxid\n";
	while(<FH>){  ## read single line from the file
		$j++;
		@line = split(/\t/,$_,2);
		if(exists $acc{$line[0]}){#if the key exists (it is a target acc_num)
			print OUT $_;
			$i++;
			print "Found Item $i in line $j\n";
			delete $acc{$line[0]};#the ones that are found are removed from the hash
		}else{
			next;
		}
		if ($i==$items){print "Total accs with taxid: $i\nAll taxids found\n";exit};#exit if all accnums are found
	}
	close (FH);
	close (OUT);
}
sub print_missing{#Print the ones that failed at the end.
	my $bad=0;
	open(FAIL,">$file.accs_with_no_taxid.txt")|| die "Couldn't create notaxid file\n";
	foreach $j (keys(%acc)){
		print FAIL "$j\n";
		$bad++;
	}
	close(FAIL);
	print "Total accs with taxid: $i\n$bad numbers had no taxid\n";
}