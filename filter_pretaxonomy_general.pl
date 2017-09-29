#Author: Rodrigo G.
#Last Update Sept 20th, 2017
#Updated Nov 29th,2016
#Adapted the script for different type of viruses
#Updated Oct 10th, 2016
#Modified the output for stdout
#Updated Sept 09,2016
#Modified for viral taxonomy
#First: July 4th, 2016
#This program was created to filter a 22-column bacterial/archaeal or viral taxonomy and match it to corresponding taxids
#Two files are expected to be inputed as arguments:
# 	-A file containing the taxonomy, the pretaxonomy output of get_taxonomy_and_sptaxid_from_acc2taxid.pl or equivalent. This has 22 taxonomic ranks.
# 	 The 2nd column has the acc_num, followed by 22 taxonomic labels, most of which are commonly just "n"s. The last 2 columns are the taxid and a common name that matches that of  	genbank
# 	-A file containing the headers in the first column (this should match headers found in the DB and an acc_num in the 2nd.
# Missing taxonomy tags can occur in most levels. In order to avoid counting missing categories together, we add the tag of its nearest neighbour, in a case controlled scenario where each permutation is accounted for separately, depending which levels are missing. 
# This may seem overkill as some cases will never be present but they are addressed here nonetheless.
# Execute as follows:
#	perl filter_pretaxonomy_general.pl <taxonomy_file> <headers_file>
#Modifications were made to avoid using a second file. This should only be used when taxonomies are not obtained using header files and thus are unrepeated. Instead, this version just prints the corresponding filtered taxonomy for each line.
#New modifications include a few tweaks into viral classification: the first one to recover Baltimore's classification (which will be included in the place of phyla) and a 2nd one to get subfamily. This will be inserted in its place (after families, before genera).
#Fixed for tables not containing the first "header column" by checking if root is in the second column or not.

use strict;
use Switch;
my $taxonomy = $ARGV[0];
my $header2accnum = $ARGV[1];
my $i="";
my @line=();
my %filtered_tax=();
my %baltimore= ('dsDNA viruses, no RNA stage'=>'I:dsDNA','ssDNA viruses'=>'II:ssDNA','dsRNA viruses'=>'III:dsRNA','ssRNA positive-strand viruses, no DNA stage'=>'IV:(+)ssRNA','ssRNA negative-strand viruses'=>'V:(-)ssRNA','Retro-transcribing viruses'=>'VI-VII:RT','Satellites'=>'Satellites');#This is for the Baltimore classification of viruses
#Important: Viruses having a viral-encoded retrotranscriptase may belong to groups VI and VII but there is no way to automatize this separation other than using the families, which must be added manually as they grow in number
#Satellites were added as a separate category in the baltimor classification group
my ($flags,$is_virus,$has_baltimore)="";
open(TAX,"$taxonomy")|| die "Couldn't read nodes file\n";
while (<TAX>){
	chomp ($_);
	@line = split(/\t/,$_);
	unshift(@line, "") if $line[1] eq "root";
#  	print "$line[3]\n";
	$flags="";
	next if $line[3] eq "n";
# 	next if $line[3] eq "Eukaryota";
	if ($line[3] eq "Viruses"){
		$is_virus=1;#Make the program know the current registry is a virus
		$line[10]=$line[13];#shift the taxonomy to include subfamily
		$line[13]=$line[16];
		$line[16]=$line[17];
		if ($_ =~ /[pP]hage/){$line[3]="Bacteriophages"}#We will now identify phages by prefixes (this one is perhaps a poor criterium)
		elsif ($_ =~ /Caudovirales/){$line[3]="Bacteriophages"}#This is for the main phage category
		elsif ($_ =~ /Myovir/){$line[3]="Bacteriophages"}#and these next three are for the families within Caudovirales (if any remained unclassified)
		elsif ($_ =~ /Siphovir/){$line[3]="Bacteriophages"}
		elsif ($_ =~ /Podovir/){$line[3]="Bacteriophages"}
		elsif ($_ =~ /Corticovir/){$line[3]="Bacteriophages"}#The next are for also bacteriophages 
		elsif ($_ =~ /Cystovir/){$line[3]="Bacteriophages"}
		elsif ($_ =~ /Inovir/){$line[3]="Bacteriophages"}
		elsif ($_ =~ /Levivir/){$line[3]="Bacteriophages"}
		elsif ($_ =~ /Microvir/){$line[3]="Bacteriophages"}
		elsif ($_ =~ /Tectivir/){$line[3]="Bacteriophages"}
		elsif ($_ =~ /Plasmavir/){$line[3]="Bacteriophages"}
		if ($_ =~ /Ligamenvir/){$line[3]="Archaeophages"}#These next are given the arhchaeophage label
		elsif ($_ =~ /Ampullavir/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Bicaudavir/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Clavavir/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Fusellovir/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Globulovir/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Guttavir/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Lipothrixvir/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Rudivir/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Haloarcula/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Halorubrum/){$line[3]="Archaeophages"}
		elsif ($_ =~ /Nitrososphaera/){$line[3]="Archaeophages"}
		if ($_ =~ /virophage/){$line[3]="Virophages";$line[13]="Lavidaviridae"}#Also, these next are to give virophages their own category (Lavidaviridae is just a proposed family but I'm not including the sputnikvirus, the proposed genus for these.
		elsif ($_ =~ /Mavirus/){$line[3]="Virophages";$line[13]="Lavidaviridae"}
		if ($_ =~ /[pP]rophage/){$line[3]="Prophages"}#some special consideration in the case that actual prophages are found
		if ($_ =~ /Aureococcus anophagefferens virus/){$line[3]="Viruses"}#Finally, revert some unwanted changes for undetermined or badly classified phages. A generic Viruses label is applied.
		elsif ($_ =~ /Phage NCTB/){$line[3]="Viruses"}
		elsif ($_ =~ /Phage NST1/){$line[3]="Viruses"}
		elsif ($_ =~ /ssRNA phage DC/){$line[3]="Viruses"}
		foreach $i (keys(%baltimore)){
			if ($line[24] =~ /$i/){
				$line[7]=$baltimore{$i};
# 				print "$line[1]\t$line[24]\t$line[7]\n";
				last;
			}
		}
	}
	########################################## End of Modification ###########################################
	if($line[7] eq "n"){$flags=$flags."0"}else{$flags=$flags."1"}#flag 0 if no tag is found for Phylum in Bacteria, Baltimore Class in Viruses
	if($line[10] eq "n"){$flags=$flags."0"}else{$flags=$flags."1"}#flag 0 if no tag is found for Class in Bacteria, Order Class in Viruses
	if($line[13] eq "n"){$flags=$flags."0"}else{$flags=$flags."1"}#flag 0 if no tag is found for Order in Bacteria, Family Class in Viruses
	if($line[16] eq "n"){$flags=$flags."0"}else{$flags=$flags."1"}#flag 0 if no tag is found for Family in Bacteria, Subfamily Class in Viruses
	if($line[20] eq "n"){$flags=$flags."0"}else{$flags=$flags."1"}#flag 0 if no tag is found for Genus in Bacteria and Viruses
 	if($line[22] eq "n"){$flags=$flags."0"}else{$flags=$flags."1"}#flag 0 if no tag is found for Species in Bacteria and Viruses
	switch ($flags){#Depending on the permutation, perform accordingly (total permutations: 64)
		#####################################################################################################################################################
		# Depending on which levels that are missing, an especific procedure is performed to fill gaps (to make them distinct "undefined tags".		    #
		# "1" Represents those present.															    #
		# "0" are those that are missing. 														    #
		# Basically, if one is missing it takes the name of the higher level and an "n" (remains undefined but it has information of its nearest neighbour. #
		# Only in specific cases will the higher taxonomy be inherited. 										    #
		#####################################################################################################################################################
		case "000000" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[13]\t$line[16]\t$line[20]\t$line[22]\n"; next}
		case "000001" {print "$line[1]\t$line[3]\tn_n_n_n_n_$line[22]\tn_n_n_n_$line[22]\tn_n_n_$line[22]\tn_n_$line[22]\tn_$line[22]\t$line[22]\n"; next}
		case "000010" {print "$line[1]\t$line[3]\tn_n_n_n_$line[20]\tn_n_n_$line[20]\tn_n_$line[20]\tn_$line[20]\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "000011" {print "$line[1]\t$line[3]\tn_n_n_n_$line[20]\tn_n_n_$line[20]\tn_n_$line[20]\tn_$line[20]\t$line[20]\t$line[22]\n"; next}
		case "000100" {print "$line[1]\t$line[3]\tn_n_n_$line[16]\tn_n_$line[16]\tn_$line[16]\t$line[16]\t$line[16]_n\t$line[16]_n_n\n"; next} #not found in our current DB
		case "000101" {print "$line[1]\t$line[3]\tn_n_n_$line[16]\tn_n_$line[16]\tn_$line[16]\t$line[16]\t$line[16]_n\t$line[22]\n"; next} #not found in our current DB
		case "000110" {print "$line[1]\t$line[3]\tn_n_n_$line[16]\tn_n_$line[16]\tn_$line[16]\t$line[16]\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "000111" {print "$line[1]\t$line[3]\tn_n_n_$line[16]\tn_n_$line[16]\tn_$line[16]\t$line[16]\t$line[20]\t$line[22]\n"; next} #not found in our current DB
		case "001000" {print "$line[1]\t$line[3]\tn_n_$line[13]\tn_$line[13]\t$line[13]\t$line[13]_n\t$line[13]_n_n\t$line[13]_n_n_n\n"; next} #not found in our current DB
		case "001001" {print "$line[1]\t$line[3]\tn_n_$line[13]\tn_$line[13]\t$line[13]\t$line[13]_n\t$line[13]_n_n\t$line[22]\n"; next} #not found in our current DB
		case "001010" {print "$line[1]\t$line[3]\tn_n_$line[13]\tn_$line[13]\t$line[13]\t$line[13]_n\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "001011" {print "$line[1]\t$line[3]\tn_n_$line[13]\tn_$line[13]\t$line[13]\t$line[13]_n\t$line[20]\t$line[22]\n"; next} #not found in our current DB
		case "001100" {print "$line[1]\t$line[3]\tn_n_$line[13]\tn_$line[13]\t$line[13]\t$line[16]\t$line[16]_n\t$line[16]_n_n\n"; next} #not found in our current DB
		case "001101" {print "$line[1]\t$line[3]\tn_n_$line[13]\tn_$line[13]\t$line[13]\t$line[16]\t$line[16]_n\t$line[22]\n"; next} #not found in our current DB
		case "001110" {print "$line[1]\t$line[3]\tn_n_$line[13]\tn_$line[13]\t$line[13]\t$line[16]\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "001111" {print "$line[1]\t$line[3]\tn_n_$line[13]\tn_$line[13]\t$line[13]\t$line[16]\t$line[20]\t$line[22]\n"; next} #not found in our current DB
		case "010000" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[10]_n\t$line[10]_n_n\t$line[10]_n_n_n\t$line[10]_n_n_n_n\n"; next} #not found in our current DB
		case "010001" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[10]_n\t$line[10]_n_n\t$line[10]_n_n_n\t$line[22]\n"; next} #not found in our current DB
		case "010010" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[10]_n\t$line[10]_n_n\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "010011" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[10]_n\t$line[10]_n_n\t$line[20]\t$line[22]\n"; next} #not found in our current DB
		case "010100" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[10]_n\t$line[16]\t$line[16]_n\t$line[16]_n_n\n"; next} #not found in our current DB
		case "010101" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[10]_n\t$line[16]\t$line[16]_n\t$line[22]\n"; next} #not found in our current DB
		case "010110" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[10]_n\t$line[16]\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "010111" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[10]_n\t$line[16]\t$line[20]\t$line[22]\n"; next} #not found in our current DB
		case "011000" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[13]\t$line[13]_n\t$line[13]_n_n\t$line[13]_n_n_n\n"; next} #not found in our current DB
		case "011001" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[13]\t$line[13]_n\t$line[13]_n_n\t$line[22]\n"; next} #not found in our current DB
		case "011010" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[13]\t$line[13]_n\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "011011" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[13]\t$line[13]_n\t$line[20]\t$line[22]\n"; next} #not found in our current DB
		case "011100" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[13]\t$line[16]\t$line[16]_n\t$line[16]_n_n\n"; next} #not found in our current DB
		case "011101" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[13]\t$line[16]\t$line[16]_n\t$line[22]\n"; next} #not found in our current DB
		case "011110" {print "$line[1]\t$line[3]\tn_$line[10]\t$line[10]\t$line[13]\t$line[16]\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "011111" {print "$line[1]\t$line[3]\tn_$line[13]\t$line[10]\t$line[13]\t$line[16]\t$line[20]\t$line[22]\n"; next} #not found in our current DB
		case "100000" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[7]_n_n\t$line[7]_n_n_n\t$line[7]_n_n_n_n\t$line[7]_n_n_n_n_n\n"; next} #not found in our current DB
		case "100001" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[7]_n_n\t$line[7]_n_n_n\t$line[7]_n_n_n_n\t$line[22]\n"; next}
		case "100010" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[7]_n_n\t$line[7]_n_n_n\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "100011" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[7]_n_n\t$line[7]_n_n_n\t$line[20]\t$line[22]\n"; next}
		case "100100" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[7]_n_n\t$line[16]\t$line[16]_n\t$line[16]_n_n\n"; next} #not found in our current DB
		case "100101" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[7]_n_n\t$line[16]\t$line[16]_n\t$line[22]\n"; next} #not found in our current DB
		case "100110" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[7]_n_n\t$line[16]\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "100111" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[7]_n_n\t$line[16]\t$line[20]\t$line[22]\n"; next} #not found in our current DB
		case "101000" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[13]\t$line[13]_n\t$line[13]_n_n\t$line[13]_n_n_n\n"; next} #not found in our current DB
		case "101001" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[13]\t$line[13]_n\t$line[13]_n_n\t$line[22]\n"; next}
		case "101010" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[13]\t$line[13]_n\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "101011" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[13]\t$line[13]_n\t$line[20]\t$line[22]\n"; next}
		case "101100" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[13]\t$line[16]\t$line[16]_n\t$line[16]_n_n\n"; next} #not found in our current DB
		case "101101" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[13]\t$line[16]\t$line[16]_n\t$line[22]\n"; next} #not found in our current DB
		case "101110" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[13]\t$line[16]\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "101111" {print "$line[1]\t$line[3]\t$line[7]\t$line[7]_n\t$line[13]\t$line[16]\t$line[20]\t$line[22]\n"; next}
		case "110000" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[10]_n\t$line[10]_n_n\t$line[10]_n_n_n\t$line[10]_n_n_n_n\n"; next} #not found in our current DB
		case "110001" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[10]_n\t$line[10]_n_n\t$line[10]_n_n_n\t$line[22]\n"; next}
		case "110010" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[10]_n\t$line[10]_n_n\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "110011" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[10]_n\t$line[10]_n_n\t$line[20]\t$line[22]\n"; next}
		case "110100" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[10]_n\t$line[16]\t$line[16]_n\t$line[16]_n_n\n"; next} #not found in our current DB
		case "110101" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[10]_n\t$line[16]\t$line[16]_n\t$line[22]\n"; next} #not found in our current DB
		case "110110" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[10]_n\t$line[16]\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "110111" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[10]_n\t$line[16]\t$line[20]\t$line[22]\n"; next}
		case "111000" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[13]\t$line[13]_n\t$line[13]_n_n\t$line[13]_n_n_n\n"; next} #not found in our current DB
		case "111001" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[13]\t$line[13]_n\t$line[13]_n_n\t$line[22]\n"; next}
		case "111010" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[13]\t$line[13]_n\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "111011" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[13]\t$line[13]_n\t$line[20]\t$line[22]\n"; next}
		case "111100" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[13]\t$line[16]\t$line[16]_n\t$line[16]_n_n\n"; next} #not found in our current DB
		case "111101" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[13]\t$line[16]\t$line[16]_n\t$line[22]\n"; next}
		case "111110" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[13]\t$line[16]\t$line[20]\t$line[20]_n\n"; next} #not found in our current DB
		case "111111" {print "$line[1]\t$line[3]\t$line[7]\t$line[10]\t$line[13]\t$line[16]\t$line[20]\t$line[22]\n"; next}
	}
}
close(TAX);