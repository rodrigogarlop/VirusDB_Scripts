#By Rodrigo G. 
#Last Update: Sept 09, 2017
#Some minor code depurations
#Last Update: Nov 28th, 2016
#The script no longer assumes the acc_num has a version (it may not have the "."). Also, I changed the separator character of the sptaxid table to "\t". Added alternative terminal node taxid for splitting instead of the sptaxid (optional; as a comment on the same line where it %taxid is filled).
#Original: 2016-09-1
#This script separates a multi fasta DB file into multiple parts, one per species so they can be processed separately downstream.
#It was originally created to split a bacterial DB and separate plasmids but may be used regardless of the source.
#Contrary to the no_Ns version, created for removing Ns in the sequences this does not modify the sequences in any sense.
#It loads an acc_num to species_taxid table in a hash and then parses a fasta with headers matching those acc_nums so that the sequences are splitted accordingly.
#Those sequences with an Acc_num that cannot be traced to a taxid are forwarded to a different file ("<file>_missing.fasta") and are not processed at all, just copied.
#Plasmids are separated into another file but only if their acc_num can be associated to a taxid.
#It may be modified to cut by using the terminal taxid rather than the actual species taxid
#It is executed as follows:
#	perl split_fastaDB_by_taxid.pl <acc_num2sptaxid_table> <input_fasta>
#the input table has 3 columns: Acc_num taxid_entry taxid_species (we use the 1st and 3rd)
use strict;
# use warnings;
my ($total,$i,$single,$header,$sequence,$acc,$missing,$cur_tax,$part,$temp,$no_id_flag)="";
my %taxid=();
my (@subseqs,@line)=();
my $acc2tax=$ARGV[0];#gets the desired output name
my $fasta=$ARGV[1];#gets the input fasta file
my $folder = "03_BuildDB";
mkdir($folder, 0777) unless(-d $folder );
open(IN,"$acc2tax")|| die "Couldn't acc2taxid table $acc2tax\n";

load_acc2taxid();#store the expected taxids in a hash
split_fasta();#use the headers and the hash to separate sequences by the provided taxids

sub load_acc2taxid{
	while(<IN>){    #read input table
		chomp $_;
		@subseqs=split(/\t/,$_);
		$taxid{$subseqs[0]}=$subseqs[1];#This may be modified to cut by using the terminal taxid rather than the actual species taxid. Change to $taxid{$subseqs[0]}=$subseqs[2];
		#print "$subseqs[0]\t$taxid{$subseqs[0]}\n";
	}
	close(IN);
}
$total=0;
sub split_fasta{
	open(IN,"gunzip -c $fasta |")|| die "Couldn't find input fasta file\n";
	while(<IN>){    #read sequencially
		chomp $_;#remove trailing "\n"
		if ($_ =~ /^>/){#every time we get a new sequence
			$no_id_flag=0;
			$total++;
			$acc="";
			$header=$_;#now we process the new header
			$cur_tax="$fasta\_missing";#we asume no taxid is found by default
			$sequence=$missing="";
			@line = split(/\s/,$_,3);
			$acc=$line[1];#extract the Acc_num
			$acc=~s/\..*//;#ignor the version
# 			print "$acc\n";
# 			$_ =~/(\w*)\.?\w*\s/;#We extract the acc_num
			if (! exists $taxid{$acc}){#If no taxid is available we print the header as it is
				open(SUB,">>$folder/$cur_tax.fasta")|| die "Couldn't open output fasta file $cur_tax.fasta\n";#the file is created/opened
				print SUB "$header\n";
				next;
			}
			else{ #Alternatively we store the corresponding taxid
				$cur_tax=$taxid{$acc};
				$_=~s/\t/ /g;#remove tabs if present
				$_=~s/[^a-zA-Z\d\s\(\)\[\]\.\-\_\/\:]//g; #delete all rare signs
				$_=~tr/ //s; #squash spaces
				#$_=~s/ /\*__\*/g; #[Optional]replace spaces with a fixed pattern to prevent downstream clustering from truncating the original names
				$header=">$cur_tax|$_";
				#if ($header=~/lasmid/){$cur_tax="Plasmids"}#plasmids are sent to a different file
				open(SUB,">>$folder/$cur_tax.fasta")|| die "Couldn't open output fasta file $cur_tax.fasta\n";#the file is created/opened
				print SUB "$header\n";
			}
		}
		else{
			print SUB "$_\n";
		}
	}
	close(IN);
	print "$total\n";
}
