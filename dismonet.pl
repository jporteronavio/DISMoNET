#########################################
#	DISMoNET: DISEASE MODULE NETWORK	#
#########################################

# GOAL: analysis of diseasome based on modules 
#########################################################################################

# PARAMETERS
#########################################################################################
# $ARGV[0]--> Three options: 
# 1) -s single disease: introduce the name of the disease
# 2) -m multiple diseases: introduce the disease names separated by spaces (\s)
# and between double quotes (" ") 
# 3) -f given file: a file with your collection of gene-disease data  
# $ARGV[1]--> disease name(s)
# $ARGV[2]--> q-value (e.g. 1.00e-5)
# $ARGV[3]--> name to identify result files 

# e.g. perl -w dismonet.pl -s lung_neoplasms 1.00e-3 lung_neoplasms_trial
# e.g. perl -w dismonet.pl -m "malaria dystrophy asthma colorectal_neoplasms" 0.05 comorbidity_1
# e.g. perl -w dismonet.pl -f gseadata.txt 1.00e-5 gseadata

# DATA
########################################################################################
# Main folder: "data" --> containing four subfolders
# Subfolder 1: "genes" --> default gene-disease information: Genopedia (May 2014) 
############# GENOPEDIA - http://www.hugenavigator.net/HuGENavigator/startPagePedia.do
############# genopedia.txt --> gene_name \t disease_name \t disease_name ...
############# complete disease name is joined by underscore, e.g. diabetes_mellitus
# Subfolder 2: "interactome" --> HI_022014_PQ_V_CC_simple.gr and HI_022014_PQ_V_CC_simple.clas
# Subfolder 3: "class_GO_enrich_10e-5" --> GO_BP_Class_10e5; GO_MF_Class_10e5; GO_CC_Class_10e5; GO_Id_Name.mapping
# Sulfolder 4: "R" --> Rscript_pvalueFisher.R

########################################################################################

#!/usr/bin/perl

use strict; 
use warnings;
use File::Copy;

print "\n\n####################################\n";
print "# DISMoNET: DISEASE MODULE NETWORK #\n";
print "####################################\n\n";

##########################################################################################
# OPTIONS. Processing a single disease (-s), multiple diseases (-m) or a given file (-f) #
##########################################################################################

#q-value and file name given as parameters
my $QValue = $ARGV[2];
my $filename = $ARGV[3];

#Creating folder "analysis" to store results
mkdir ("./analysis_$filename");

#Functions with the three options
if ($ARGV[0] eq "-s")
{
	print "Selection of human disease genes\n\n"; 	
	
	#Creating gene lists for the required disease 
	my $trait = $ARGV[1]; 
	
	open (GENES, ">./analysis_$filename/genes_$filename") or die "can't open file ./analysis_$filename/genes_$filename\n";
	open (DATA, "./data/genes/genopedia.txt") or die "can't open file ./data/genes/genopedia.txt\n";
	while (my $line = <DATA>) 
	{
		chomp ($line);	
		my @fields = split (/\t/, $line);
		if ($line =~ /$trait/i) 
		{ 
			print GENES "$fields[0]\t$trait\n";		
		}
	}
	close GENES;
	close DATA;

	print "list of genes of $trait done\n\n";
	
	print "Selection of genes implicated in the human interactome\n\n";
		
	ppi_interactions ("./analysis_$filename/genes_$filename");
	
	print "Selection of genes done\n\n";
		
	print "Performing Fisher´s exact test\n\n";	
	
	fisher_test ("./analysis_$filename/gene_ppi_$filename");	
}
elsif ($ARGV[0] eq "-m")
{	
	print "Selection of human disease genes\n\n";
	
	#Creating gene lists for the different diseases 
	my @pathologies = split (/\s/, $ARGV[1]);	
	open (GENES, ">./analysis_$filename/genes_$filename") or die "can't open file ./analysis_$filename/genes_$filename\n";
	open (DATA, "./data/genes/genopedia.txt") or die "can't open file ./data/genes/genopedia.txt\n";
	while (my $line = <DATA>) 
	{
		chomp ($line);	
		my @fields = split (/\t/, $line);
		foreach my $illness (@pathologies)
		{
			if ($line =~ /$illness/i) 
			{ 
					print GENES "$fields[0]\t$illness\n";	
			}
		}	
	}	
	close GENES;
	close DATA;

	print "list of genes done\n\n";

	duplicate ("./analysis_$filename/genes_$filename");		
	
	print "Selection of genes implicated in the human interactome\n\n";

	ppi_interactions ("./analysis_$filename/genesdupl_$filename");
	
	print "\n\nSelection of genes done\n\n";
	
	print "Performing Fisher´s exact test\n\n";
	
	fisher_test ("./analysis_$filename/gene_ppi_$filename");
}
elsif ($ARGV[0] eq "-f") 
{
	duplicate ("$ARGV[1]");
	
	print "Selection of genes implicated in the human interactome\n\n";		
	
	ppi_interactions ("./analysis_$filename/genesdupl_$filename");
	
	print "\n\nSelection of genes done\n\n";	

	print "Performing Fisher´s exact test\n\n";	
	
	fisher_test ("./analysis_$filename/gene_ppi_$filename");
}
else
{
	print "No valid parameter";
	exit;
}	
		
#Moving file.Rout to R folder
mv ("./Rscript_pvalueFisher.Rout", "./data/R/Rscript_pvalueFisher.Rout");

#Table Class Disease 
print "Table Class Disease\n\n";

#Creating and printing table with the relation class-disease
open (DISEASECLASS, "./analysis_$filename/class_$filename") or die "cannot open file ./analysis_$filename/class_$filename\n";
my %diseasome = ();
while (my $line = <DISEASECLASS>) 
{
	if ($line =~ /^(Class_[0-9]+\|[0-9]+)/)
	{	
		my @matches = ();
		while ($line =~ m/([a-zA-Z_a-z]+\t)+/g) 
		{
			push (@matches, $&);
		}
		foreach my $name (@matches)
		{
        		if (!exists ($diseasome{$name}))
        		{
        			$diseasome{$name} = $name;
    			}
		}
	}
}
close DISEASECLASS;

open (DC, ">./analysis_$filename/disease_class_$filename") or die "cannot open file ./analysis_$filename/disease_class_$filename";
open (NETWORK, ">./analysis_$filename/network_$filename") or die "cannot open file ./analysis_$filename/network_$filename";

my @list_disease = keys %diseasome;
foreach my $elem (@list_disease)
{
	#Eliminating NA from disease class table
	if ($elem ne "NA")
	{
		print $elem."\n";
		print DC "$elem\n";
	}	
	open (DISEASECLASS, "./analysis_$filename/class_$filename") or die "cannot open file ./analysis_$filename/class_$filename\n";
	while (my $line = <DISEASECLASS>) 
	{
		chomp $line;
		if ($line =~ /($elem)/)  
		{
		 	#Eliminating NA from disease class table
		 	if ($elem ne "NA")
		 	{
		 		$line =~ /^(Class_[0-9]+\|[0-9]+)/;
		 		print "\t".$1."\t";
		 		print DC "\t$1\t";		 	
		 		print NETWORK "$elem\t$1\n";
		 		$line =~ /($elem)(\d*\s\d*\s\d*\s\d*\s)(([0-9\-\.a-z]+)\s([0-9\-\.a-z]+))/;
		 		print $3."\n";
		 		print DC "$3\n";
		 	}		 	
		}
	}
	print "\n";
	print DC "\n";
}
close DISEASECLASS;
close DC;
close NETWORK;

#Annotation with Gene Ontology
go_annotation ("./analysis_$filename/disease_class_$filename");


################################################################################
################################## SUBFUNCTIONS ################################
################################################################################

################
# DUPLICATIONS #
################

sub duplicate
{
	#Creating a table gene associated disease(s)
	undef my %doublegene;
	undef my %synonyms;

	open (OUTPUT, ">./analysis_$filename/genesdupl_$filename") or die "cannot open file ./analysis_$filename/genesdupl_$filename\n";
	open (GENEDIS, "$_[0]") or die "cannot open file $_[0]\n";
	while (my $l = <GENEDIS>) 
	{
		chomp $l;
    	if ($l =~ /^(\w+)\t(\w+)$/) 
    	{
    		my $ref = $1;
			$doublegene{$ref}++;
			$synonyms{$ref}.="$2\t";
		}
	}
	while ((my $keys, my $values) = each (%synonyms)) 
	{
  		print "$keys\t$values\n";	
		print OUTPUT "$keys\t$values\n";
	}
	close GENEDIS;
	close OUTPUT;
}


##########################
# PPI selection of genes #
##########################

sub ppi_interactions
{
	#Creating a list of all genes from the human interactome
	open (INTERACTOME, ">./analysis_$filename/genes_interactome") or die "cannot open file ./analysis_$filename/genes_interactome\n";
	open (PPIDB, "./data/interactome/HI_022014_PQ_V_CC_simple.gr") or die "cannot open file ./data/interactome/HI_022014_PQ_V_CC_simple.gr\n";
	my %ppi = ();
	while (my $line = <PPIDB>)
	{
    	chomp $line;
    	if ($line =~ /^(.*?)\t(.*)$/)
    	{
	 		my $int1 = $1;
	 		if (not exists $ppi{$int1})
	 		{
	 			$ppi{$int1} = $int1;	 	
	 			print INTERACTOME $ppi{$1}."\n";
	 		}
	 		my $int2 = $2;
	 		if (not exists $ppi{$int2})
	 		{
	 			$ppi{$int2} = $int2;
	 			print INTERACTOME $ppi{$2}."\n";
	 		}
    	}
	}
	close PPIDB;
	close INTERACTOME;

	print "list of interactome genes done\n\n";

	#Creating a gene list from the overlapping interactome versus disease
	open (LISTFINAL, ">./analysis_$filename/gene_ppi_$filename") or die "cannot open file ./analysis_$filename/gene_ppi_$filename\n";
	open (INTER, "$_[0]") or die "cannot open file $_[0]";
	my %sel = ();
	while (my $line = <INTER>)
	{
    	chomp $line;
    	if ($line =~ /^(.*?)\t(.*?)$/)
    	{
			my $id = $1;
			my $disname = $2;
			$sel{$id} = $disname;
		}
	}
	open (DB, "./analysis_$filename/genes_interactome") or die "cannot open file ./analysis_$filename/genes_interactome";
	while (my $line = <DB>)
	{
		chomp $line;
		if (exists ($sel{$line}))
		{
			print $line."\t".$sel{$line}."\n";
			print LISTFINAL "$line\t$sel{$line}\n";	
		}
		else
		{
			print $line."\tNA\n";
			print LISTFINAL "$line\tNA\n";	
		}
	}
	close INTER;
	close LISTFINAL;
}


#######################
# Fisher´s exact test # 
#######################

sub fisher_test 
{
	print "Fisher test for each class/annotation\n\n";

	#Reading file with genes annotations and creating hash (genes = keys / value = annotations)
	my %asso=();
	open (GENELIST, "$_[0]") or die "cannot open file $_[0]";
	while (my $lin = <GENELIST>)
	{
		if ($lin =~ /^(.*?)\t(.*)$/)
		{
			chomp $lin;
			my $geneid = $1;
			my $endline = $2;
			$asso{$geneid} = $endline;
		}	
	}
	# Retrieving and counting annotations (%seen = all annotations from the annotation file)
	my %seen=();
	while ((my $keys, my $values) = each (%asso)) 
	{
		my @list= ();
		@list = split ("\t", $values);
		foreach my $item (@list)
		{
			chomp $item;
			unless ($seen{$item})
			{
				$seen{$item} = 1;
			}	 
			else 
			{
				$seen{$item}++;
			}	
		} 
	}
	# Reading class file (HI_022014_PQ_V_CC_simple.clas)
	my %ClassGenesAnnot = ();
	my %References = ();
	my %details = ();
	my @listvalues = ();

	my $j=0;
	my $f=1;

	open (FISHER, ">./analysis_$filename/class_$filename") or die "cannot open file ./analysis_$filename/class_$filename";
	open (CLASSFILE, "./data/interactome/HI_022014_PQ_V_CC_simple.clas") or die "cannot open file ./data/interactome/HI_022014_PQ_V_CC_simple.clas"; 
	while (my $l = <CLASSFILE>)
	{
		chomp $l;
		$j = 1;
		print "\n";
		print FISHER "\n";
		print "Class_$f|";
		print FISHER "Class_$f|";	
		
		$f++;
		my $i = 0;
		
		%ClassGenesAnnot = ();
		
		# Opening file to write the contingency table for R
		open (W, ">./data/R/contingencytable.tmp") or die "cannot open file ./data/R/contingencytable.tmp";
		# Creating table with all class genes
		my @genelist = split ("\t", $l);
		# Printing class size
		my $nbegene = $#genelist+1;
		print "$nbegene\t";
		print FISHER "$nbegene\t";
		
		# Searching annotations for class genes in the hash table created in the first step
		foreach my $gene (@genelist)
		{
			# For a given gene, retrieving its annotations
			if (exists ($asso{$gene})) 
			{
				# Creating a table listing all gene annotations
				my @ClassAnnot = split("\t", $asso{$gene});
				$i++;
				# For each annotation of each gene, creating a hash table with key = annotation and 
				# value = number of time the annotation is seen in the class
				foreach my $ClassAnnot (@ClassAnnot)
				{
					if (exists ($ClassGenesAnnot{$ClassAnnot})) 
					{
						$ClassGenesAnnot{$ClassAnnot}++;
					}	 
					else 
					{
						$ClassGenesAnnot{$ClassAnnot} = 1;
					}
				}
			}		 
		}
		#Reading the hash table containing annotation counts to create the contingency table for R
		#Contingency table contains DatasetPositive DatasetNegative GenomePositive GenomeNegative
		#(Interactome = Background)
		if (!keys %ClassGenesAnnot)			
		{ 
		 	next
		} 
		else 
		{
		 	while ((my $keys, my $values) = each (%ClassGenesAnnot)) 
		 	{
		 		# $values correspond to Dataset Positive			
				my $DatasetPositive = $values;
				# $i correspond to Dataset size
				my $DatasetSize = $i;
				# Dataset Negative
				my $DatasetNegative = $DatasetSize - $values;
				# Creating hash table with $j corresponding to annotation number
				# in order to read R output that is printed according to p-values
				$References{$j} = $keys;
				if (exists ($seen{$keys}))
				{
					$j++;
					my $GenomePositive = $seen{$keys};
					# Genome size correspond to %seen size
					my $GenomeSize = keys(%asso);
					my $GenomeNegative = $GenomeSize - $GenomePositive;
					@listvalues = ($DatasetPositive,$DatasetNegative,$GenomePositive, $GenomeNegative);
					@{$details{$j}} = @listvalues;
					print W "$DatasetPositive\t$DatasetNegative\t$GenomePositive\t$GenomeNegative\n";
				} 	
			}
		}
		# R calling to Rscript_pvalueFisher.R
		system("R CMD BATCH ./data/R/Rscript_pvalueFisher.R");
		# Parsing R results, and selecting them according to a q-value threshold
		open (H,"<./data/R/resultsR") or die "cannot open file ./data/R/resultsR\n";
		# Results for only one class
		while (my $lineFisher = <H>)
		{
			if ($lineFisher =~ /^(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)$/)
			{
				my @pvalues = ();
				my $num = $1;
				my $num2 = $2;
				my $rawp = $3;
				my $bonferro = $4;
				my $BH = $5;
				# Retrieving annotation name
				if (exists ($References{$num2}))
				{
				# Print only q-value under threshold
					if ($BH<=$QValue) 
					{
						print "$References{$num2}\t";
						print "@{$details{$num2+1}}\t";
						print "$rawp\t$BH\t";
						print FISHER "$References{$num2}\t";
						print FISHER "@{$details{$num2+1}}\t";
						print FISHER "$rawp\t$BH\t";					
					}		
				}
			} 
			elsif ($lineFisher =~ /^\"(.*?)\"\s(.*?)$/)
			{
				my @pvalues = ();
				my $numalone = $1;
				my $rawpalone = $2;
				# Retrieving annotation name
				if (exists ($References{$numalone}))
				{
					# Print only q-value under threshold
					if ($rawpalone<=$QValue) 
					{
						print "$References{$numalone}\t";
						print "@{$details{$numalone+1}}\t";
						print "$rawpalone\t";
						print FISHER "$References{$numalone}\t";
						print FISHER "@{$details{$numalone+1}}\t";
						print FISHER "$rawpalone\t";					
					}
				}
			}
		}
	}
	close CLASSFILE;
	close FISHER;
	close H;
	close W;

	print "Fisher test analysis_$filename done\n\n"; 
}

############################
# GENE ONTOLOGY ANNOTATION #
############################

sub go_annotation 
{
	#Creating hash with Gene Ontology, keys = GO identification qnd values = function
	my %gene_ontology = ();
	open (GONAMES, "./data/class_GO_enrich_10e-5/GO_Id_Name.mapping") or die "cannot open file ./data/class_GO_enrich_10e-5/GO_Id_Name.mapping";
	while (my $line =  <GONAMES>)
	{
		if ($line =~ /^(.+)*\t(.+)*$/)
		{
			my $go_id = $1;
			my $function = $2;
			$gene_ontology{$go_id} = $function;
		}
	}
	close GONAMES;

	#Creating three hashes with annotations from GO (BP, MF, CC)
	my %GOBP = ();
	open (BP, "./data/class_GO_enrich_10e-5/GO_BP_Class_10e5") or die "cannot open file ./data/class_GO_enrich_10e-5/GO_BP_Class_10e5";
	while (my $line = <BP>)
	{
		$line =~/(Class_[0-9]+\|[0-9]+)\t((.*)$)/;
		$GOBP{$1} = $2;
	}

	my %GOMF = ();
	open (MF, "./data/class_GO_enrich_10e-5/GO_MF_Class_10e5") or die "cannot open file ./data/class_GO_enrich_10e-5/GO_MF_Class_10e5";
	while (my $line = <MF>)
	{
		$line =~/(Class_[0-9]+\|[0-9]+)\t((.*)$)/;
		$GOMF{$1} = $2;
	}

	my %GOCC = ();
	open (CC, "./data/class_GO_enrich_10e-5/GO_CC_Class_10e5") or die "cannot open file ./data/class_GO_enrich_10e-5/GO_CC_Class_10e5";
	while (my $line = <CC>)
	{
		$line =~/(Class_[0-9]+\|[0-9]+)\t((.*)$)/;
		$GOCC{$1} = $2;
	}

	#Annotating file disease class with BP
	open (BPAN, ">./analysis_$filename/bpannotation.csv")  or die "cannot open file ./analysis_$filename/bpannotation";
	open (FILE1, "$_[0]") or die "cannot open file $_[0]";
	while (my $line = <FILE1>)
	{
		chomp $line;
		unless ($line =~ /(Class_[0-9]+\|[0-9]+)/)
		{
			print $line."\n";
			print BPAN "$line\n";
		}
		if ($line =~ /(Class_[0-9]+\|[0-9]+)/)	
		{
			if (exists ($GOBP{$1}))
			{
				my $string = $GOBP{$1};
				my @ma = ();
				while ($string =~ m/(GO:\d{7})+/g)
				{ 
					push (@ma, $&);
				}
				foreach my $mt (@ma)
				{
					print $line."\t".$mt."\t".$gene_ontology{$mt}."\n";
					print BPAN "$line\t$mt\t$gene_ontology{$mt}\n";						
				}
			}		
		}
	}			
	close BPAN;

	#Annotating file disease class with MF
	open (MFAN, ">./analysis_$filename/mfannotation.csv")  or die "cannot open file ./analysis_$filename/mfannotation";
	open (FILE1, "$_[0]") or die "cannot open file $_[0]";
	while (my $line = <FILE1>)
	{
		chomp $line;
		unless ($line =~ /(Class_[0-9]+\|[0-9]+)/)
		{
			print $line."\n";
			print MFAN "$line\n";
		}
		if ($line =~ /(Class_[0-9]+\|[0-9]+)/)	
		{
			if (exists ($GOMF{$1}))
			{
				my $string = $GOMF{$1};
				my @ma = ();
				while ($string =~ m/(GO:\d{7})+/g)
				{ 
					push (@ma, $&);
				}
				foreach my $mt (@ma)
				{
					print $line."\t".$mt."\t".$gene_ontology{$mt}."\n";
					print MFAN "$line\t$mt\t$gene_ontology{$mt}\n";						
				}
			}		
		}
	}	
	close MFAN;

	#Annotating file disease class with CC
	open (CCAN, ">./analysis_$filename/ccannotation.csv")  or die "cannot open file ./analysis_$filename/ccannotation";
	open (FILE1, "$_[0]") or die "cannot open file $_[0]";
	while (my $line = <FILE1>)
	{
		chomp $line;
		unless ($line =~ /(Class_[0-9]+\|[0-9]+)/)
		{
			print $line."\n";
			print CCAN "$line\n";
		}
		if ($line =~ /(Class_[0-9]+\|[0-9]+)/)	
		{
			if (exists ($GOCC{$1}))
			{
				my $string = $GOCC{$1};
				my @ma = ();
				while ($string =~ m/(GO:\d{7})+/g)
				{ 
					push (@ma, $&);
				}
				foreach my $mt (@ma)
				{
					print $line."\t".$mt."\t".$gene_ontology{$mt}."\n";
					print CCAN "$line\t$mt\t$gene_ontology{$mt}\n";						
				}
			}
		}		
	}	
	close CCAN;
}
