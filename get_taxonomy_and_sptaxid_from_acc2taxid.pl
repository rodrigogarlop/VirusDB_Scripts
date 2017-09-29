#Update Sept 19th, 2017
#Fixed path problems. The nodes and names files must be in the same folder as the input
#Restored the whole taxonomy print, this time by species taxid
#Update November 25th, 2016
#The script now only prints one species-level taxid for each input taxid
#Update Sept 15th, 2016
#Apart of getting the whole 22 categories table, this now prints a taxid to species_lvl_taxid to the stdout to use it for splitting a fasta by sptaxid (originally for the whole 200GB BacteriaDB_2016_08_08) These are also printed in the last 2 columns of the 24 categories table (columns 23 and 24). Column 23 is the initial taxid corresponding to the Acc_num and column 24 has the actual species taxid (most times they are the same but not always).
#Created Aug 31th, 2016
#Just changed this version of the script so that not fasta is requested. Thus, this prints a list of non repeated acc_nums to taxonomies.
#This program reads a list of accs, a series of dmps and a database in fasta format and gets the corresponding taxonomy for each entry in the DB
#It also requires the names.dmp and nodes.dmp files included in the taxdump.tar.gz pack to be in the same folder
#This can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
#The dmp is read only once and no multiple comparisons are made
#The input is a tab separated file containing acc numbers, version, taxid and gi in each row (only the first and third column are used)
#e.g.
#	KY921644 KY921644.1  1384672 1179780879
#the result would be (Acc_num, species_node, terminal_node):
#	KY921644 1384672  1979165
#run as follows:
#	perl get_taxonomy_and_sptaxid_from_acc2taxid.pl <input_taxid_file(acc_taxid)>
#!/usr/bin/perl
use File::Basename;
use strict;
use warnings;
my $file = $ARGV[0];
my $dir=dirname($file);
my $flag=0;
my @line=();
my @gilist=();
my %nodes=();
my %taxid=();
my %gis=();
my %names=();
my %rank=();
my %actual_ranks=();
my %sptaxid2tax=();
my $species=();
# my @rank_names=("superkingdom","phylum","class","order","family","genus","species");
my @rank_names=("superkingdom","kingdom","subkingdom","superphylum","phylum","subphylum","superclass","class","subclass","superorder","order","suborder","superfamily","family","subfamily","tribe","subtribe","genus","subgenus","species","subspecies");
my @keys=();
my ($temp,$current,$next,$norank)="";
my ($i,$j)="";

load_nodes_and_ranks(); #First load everything into hashes (relatively small files)
load_names();
load_taxids();
build_taxonomy(); #Build the 22 level taxonomy for each sptaxid
create_acc2sptaxid_table(); #Get a table with acc_nums to species taxids

sub load_nodes_and_ranks{ #First, the Nodes.dmp file is indexed using two ref tables (hash)
	open(NODES,"$dir/nodes.dmp")|| die "Couldn't read nodes.dmp. Expected in directory $dir\n";
	while (<NODES>){
		chomp ($_);
		@line = split(/\t\|\t/,$_,4);
		$nodes{$line[0]}=$line[1]; #This stores the nodes (keys are taxids)
		$rank{$line[0]}=$line[2]; #This stores the ranks
	}
	close(NODES);
}
sub load_names{ #Then, the scientific names are loaded as well
	open(NAMES,"$dir/names.dmp")|| die "Couldn't read names.dmp. Expected in directory $dir\n";
	while (<NAMES>){
		if ($_ =~ /scientific name/){ #only use scientific names
			chomp ($_);
			@line = split(/\t\|\t/,$_,3);
			$names{$line[0]}=$line[1]; #store taxids and sci_names
		}
	}
	close(NAMES);
}
sub load_taxids{ #now we get the actual taxids in out dataset
	open(ACTAX,"$file")|| die "Couldn't read $file\n";
	while (<ACTAX>){
		chomp($_);
		@line = split(/\t/,$_);
		$taxid{$line[2]}=1;
	}
	close(ACTAX);
}
sub build_taxonomy{
	open(TAX,">$file.pretaxonomy")|| die "Couldn't create output taxonomy file\n";
	foreach $i (sort(keys(%taxid))){ #for each taxid
		#print "\n";
		$taxid{$i}="root\t"; #start filling the taxonomy string while recycling the hash
		%actual_ranks=(); foreach $j (@rank_names){$actual_ranks{$j}="n";} #The taxonomy is reset for each taxid (default to n)
		$current=$i; $next=$nodes{$current}; #we set the first nodes for comparison
		#The acc2taxid file has some errors (not our fault) so we separate these exceptions
 		if((! exists $nodes{$current})||($current eq "0")){#when the next node is empty or we have a 0 as taxid
			open(BAD,">>$file.bad")|| die "Couldn't create/append file\n";
			print BAD "No information found for taxid: $current\n";
			close(BAD);
			delete $taxid{$i};
			next;
		}
		if($rank{$current} eq "no rank"){$norank="$names{$current}";}#when rank is "no rank" save it for extra info
		else{$actual_ranks{$rank{$current}}="$names{$current}";}#in any other case, save the name associated to the taxid into a hash matching each rank and print the node to the node map
		while ($current != $next){#this will advance node by node and will stop when taxid reaches root taxonomic level (1)
			$temp=$norank;#just for avoiding overwriting this value
			if ($rank{$current} eq "no rank"){$norank="$temp|$names{$current}";}#when rank is "no rank" append it to "extra information"
			else{$actual_ranks{$rank{$current}}="$names{$current}";}# any no rank information is saved in a cummulative string
			if ($rank{$current} eq "species"){$species=$current;}
			$current=$next;#we advance one node at a time
			$next=$nodes{$current};#and, again, get the next one for comparison
		}
		foreach($j=0;$j<scalar(@rank_names);$j++){$taxid{$i}= "$taxid{$i}$actual_ranks{$rank_names[$j]}\t";}#finally we print the complete taxonomy (including 'n's where no corresponding rank was found
		$norank=~s/^\|//;
 	 	if (! exists $sptaxid2tax{$species}){
			print TAX "$species\t$taxid{$i}$norank\n";
			$sptaxid2tax{$species}=undef;
			
		}
		$taxid{$i}="$species\t$i";#Recycle the hash to print the acc2sptaxid table
		$norank="";
	}
	close(TAX);
	close(SP);
	(%nodes,%names,%rank,%actual_ranks)=();#free some memory
}
sub create_acc2sptaxid_table {#we go through the acc_taxid file again, this time to construct a 3-column table containing the Acc_num, the species taxid and the terminal taxid node
	open(SP,">$file.sptaxid")|| die "Couldn't create output tax file\n";
	open(ACTAX,"$file")|| die "Couldnt read acc_taxid file\n";
	while (<ACTAX>){
		chomp($_);
		@line = split(/\t/,$_);
		next if ! exists $taxid{$line[2]};
		print SP "$line[0]\t$taxid{$line[2]}\n";
		print "No taxonomy for $_\n" if ! exists $taxid{$line[2]};
	}
	close(SP);
	close(ACTAX);
	%taxid=();
}