#########################################
#    DISMoNET: DISEASE MODULE NETWORK	#
#########################################

# GOAL: analysis of diseasome based on modules 
#########################################################################################

# PARAMETERS
#########################################################################################
# $ARGV[0]--> a file with your collection of gene-disease data tabulated in two columns
# gene name (HUGO Gene Nomenclature) \t disease name  
# $ARGV[1]--> q-value (e.g. 1.00e-5)
# $ARGV[2]--> homogeneity (number of disease susceptibility genes in a module divided by
# the total number of genes in the module); it is a cut-off selecting disease 
# enriched modules with homogeneity above the threshold established as parameter
# eg: 0.1 (10%), 0.25 (25%), 0.30 (30%) ... (0 in case of non-applicable)
# $ARGV[3]--> name to identify result files 

# e.g. perl -w dismonet.pl gseadata.txt 1.00e-5 0.50 gseadata

# DATA
########################################################################################
# Main folder: "data" --> containing FOUR subfolders
# Subfolder 1: "interactome" --> HI_022014_PQ_V_CC_simple.gr and HI_022014_PQ_V_CC_simple.clas
# Subfolder 2: "class_GO_enrich_10e-5" --> GO_BP_Class_10e5; GO_MF_Class_10e5; GO_CC_Class_10e5; GO_Id_Name.mapping
# Subfolder 3: "class_reactome_enrich_10e-5" --> Class_Reactome_Enrich_10e-5
# Sulfolder 4: "R" --> Rscript_pvalueFisher.R
########################################################################################

#!/usr/bin/perl

use strict; 
use warnings;
use File::Copy;

print "\n\n####################################\n";
print "# DISMoNET: DISEASE MODULE NETWORK #\n";
print "####################################\n\n";

########################################
# Processing a given gene-disease file #
########################################

#q-value, ratio of homogeneity and file name for results given as parameters
my $QValue = $ARGV[1];
my $homogeneity = $ARGV[2];
my $filename = $ARGV[3];

#Counting number of genes associated to each disease to summarize results
my %counting = ();

#Creating folder "analysis" to store results
mkdir("./analysis_$filename");

#Given file gene-disease data analysis
#Subfunction duplicate creates a tabulated file with gene name and associated diseases
duplicate ("$ARGV[0]");
print "Selection of genes implicated in the human interactome\n\n";		
	
ppi_interactions ("./analysis_$filename/genesdupl_$filename");
print "\n\nSelection of genes done\n\n";	

counting_genes ("./analysis_$filename/gene_ppi_$filename");
print "Performing Fisher´s exact test\n\n";	
	
fisher_test ("./analysis_$filename/gene_ppi_$filename");
		
#Moving file.Rout to R folder
move("./Rscript_pvalueFisher.Rout","./data/R/Rscript_pvalueFisher.Rout");

#Creating and printing table with the relation class-disease
print "Table Class Disease\n\n";
open (DISEASECLASS, "./analysis_$filename/class_$filename") or die "cannot open file ./analysis_$filename/class_$filename\n";
my %diseasome = ();
while (my $line = <DISEASECLASS>) 
{
	if ($line =~ /^(Class_[0-9]+\|[0-9]+)/)
	{	
		chomp $line;
		my @matches = ();
		while ($line =~ m/(\t\w+\t)/g) 
		{
			push (@matches, $&);
		}
		foreach my $tabname (@matches)
		{
        	chop $tabname;
        	my $name = substr $tabname,1;
        	if (!exists ($diseasome{$name}))
        	{
        		$diseasome{$name} = 1;
    		}
			else		
			{
				$diseasome{$name}++;			
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
		 		$line =~ /^(Class_\d+\|\d+)/;
		 		print "\t".$1."\t";
		 		print DC "\t$1\t";		 	
		 		print NETWORK "$elem\t$1\n";
		 		$line =~ /($elem)\t(\d*\s\d*\s\d*\s\d*)\t([\d\-\.e]+\t[\d\-\.e]+)/;
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

#Gene Ontology Annotation
go_annotation ("./analysis_$filename/disease_class_$filename");

#Reactome Annotation 
reactome_annotation ("./analysis_$filename/disease_class_$filename");

#Summarizing results (Disease - Number of genes - Number of modules)
open (SUMMARY, ">./analysis_$filename/results_experiment.txt") or die "cannot open file ./analysis_$filename/results_experiment";
print "Disease\tNumber_of_Genes\tNumber_of_Modules\n";
print SUMMARY "Disease\tNumber_of_Genes\tNumber_of_Modules\n";
while ((my $keys, my $values) = each (%diseasome))
{
	print $keys."\t".$counting{$keys}."\t".$values."\n";
	print SUMMARY "$keys\t$counting{$keys}\t$values\t\n";
}
close SUMMARY;

#Counting the number of modules associated to each disease
my %nosome = ();
my %counter = ();
open (FILENET, "./analysis_$filename/network_$filename") or die "cannot open file ./analysis_$filename/network_$filename";
while (my $line = <FILENET>)
{
	chomp $line;
	$line =~ /(.*)\t(.*)/;
	my $disease = $1;
	my $class = $2;
	$nosome{$class} = $disease;
	if (not exists ($nosome{$class}))
	{
		$counter{$disease} = 1;
	}
	else
	{
		$counter{$disease}++;
	}	
}
close FILENET;

#Creating hash of arrays with the diseases sharing the same module
my %HoA = ();
open (FILENET, "./analysis_$filename/network_$filename") or die "cannot open file ./analysis_$filename/network_$filename";
while (my $line = <FILENET>)
{
	chomp $line;
	$line =~ /(.*)\t(Class_\d+\|\d+)/; 
	my $disease = $1;
	my $class = $2;
	push (@{$HoA{$class}}, $disease);
}

foreach my $k (keys %HoA)
{
	foreach my $el (@{$HoA{$k}})
	{
		print $el."\t".$k."\n"
	}
}
close FILENET;

#Creating file disease-disease associated to shared module
my %diff = ();
open (NOTES, ">./analysis_$filename/Common_modules.txt") or die "cannot open ./analysis_$filename/Common_modules.txt";
open (FILE, "./analysis_$filename/network_$filename") or die "cannot open file ./analysis_$filename/network_$filename";
while (my $line = <FILE>)
{
	chomp $line;
	$line =~ /(.*)\t(.*)/;
	my $disease = $1;
	my $class = $2;
	foreach my $el (@{$HoA{$class}})
	{
		if ($el ne $disease)
		{
			if (exists ($diff{$el.$disease.$class}) or ($diff{$disease.$el.$class}))
			{
				next;
			}
			else
			{
				$diff{$el.$disease.$class} = $el.$disease.$class;
				$diff{$disease.$el.$class} = $disease.$el.$class;
				print "$disease\t$el\t$class\n";
				print NOTES "$disease\t$el\t$class\n";
			}
		}
	}
}
close FILE;
close NOTES;

#Ordering pair diseases
my %order = ();
open (NOTES2, ">./analysis_$filename/Pair_diseases_mod.txt") or die "cannot open ./analysis_$filename/Pair_disease_mod.txt";
open (FILE2, "./analysis_$filename/Common_modules.txt") or die "cannot open ./analysis_$filename/Common_modules.txt";
while (my $line = <FILE2>)
{
	chomp $line;
	$line =~ /(\w*\t\w*)\t(Class_\d+|\d+)/;
	my $pairdisease = $1;
	my $classes = $2;
	if (not exists ($order{$pairdisease}))
	{
		$order{$pairdisease} = $classes;
	}
	else
	{
		$order{$pairdisease} .= "\t$classes";
	}
}
close FILE2;

while (my ($key, $value) = each %order)
{
	print "$key\t$value\n";
	print NOTES2 "$key\t$value\n";
}
close NOTES2;

#Creating file with disease pair sharing modules and Jaccard and Dice-Czekanowski distances
open (NOTES3, ">./analysis_$filename/Module_distances.txt") or die "cannot open ./analysis_$filename/Module_distances.txt";
open (FILE3, "./analysis_$filename/Pair_diseases_mod.txt") or die "cannot open ./analysis_$filename/Pair_disease_mod.txt";
while (my $line = <FILE3>)
{
	chomp $line;
	$line =~ /(\w+)\t(\w+)\t(.*)/;
	my $disease1 = $1;
	my $disease2 = $2;
	my @intersec = split ("\t", $3);
	my $intersection = scalar@intersec;
	my $union = ($counter{$disease1}) + ($counter{$disease2});
	my $Jaccard = ($intersection) / ($union);
	my $Dice = (($union) - ($intersection)) / (($union) + ($intersection));
	print "$disease1\t$counter{$disease1}\t$disease2\t$counter{$disease2}\tUnion\t$union\tIntersection\t$intersection\tJac\t$Jaccard\tDice-Czekanowski\t$Dice\t$3\t\n";
	print NOTES3 "$disease1\t$counter{$disease1}\t$disease2\t$counter{$disease2}\tUnion\t$union\tIntersection\t$intersection\tJac\t$Jaccard\tDice-Czekanowski\t$Dice\t$3\t\n";
}
close NOTES3;
close FILE3;

system ("sort -k8nr,8 ./analysis_$filename/Module_distances.txt | cut -f1-8 > ./analysis_$filename/Ranking_diseases_modules_intersection.txt");		 
system ("sort -k10nr,10 ./analysis_$filename/Module_distances.txt | cut -f1-12 > ./analysis_$filename/Ranking_diseases_Jaccard_Dice.txt");		


################################################################################
################################## SUBFUNCTIONS ################################
################################################################################

#######################
# DUPLICATION MANAGER #
#######################

#To avoid data duplication merging all the diseases related to a gene in a single row
sub duplicate
{
	#Creating a table gene associated disease(s)
	undef my %doublegene;
	undef my %synonyms;

	open (OUTPUT, ">./analysis_$filename/genesdupl_$filename") or die "cannot open file ./analysis_$filename/genesdupl_$filename\n";
	open (GENEDIS, "$_[0]") or die "cannot open file $_[0]\n";
	while (my $line = <GENEDIS>) 
	{
		chomp $line;
    	if ($line =~ /^(\w+)\t(\w+)$/) 
    	{
    		my $ref = $1;
			$doublegene{$ref}++;
			$synonyms{$ref} .= "$2\t";
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

################
# GENE COUNTER #
################

#To count the number of genes associated to each disease
sub counting_genes
{
	open (COUNTS, "$_[0]") or die "cannot open file $_[0]\n";
	while (my $line = <COUNTS>)
	{
		chomp $line;
		$line =~ /^([A-Za-z0-9\-\@]+)\t(.*)$/;
		my @dises = split (/\t/, $2);
		foreach my $element (@dises)
		{
			if (!exists ($counting{$element}))
			{
				$counting{$element} = 1;
			}
			else
			{
				$counting{$element}++;
			}
		}
	}
}

######################
# PPI GENE SELECTION #
######################

#To identify protein-protein interactions in the products of our genes of interest
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
# FISHER´S EXACT TEST # 
#######################

#To identify the modules enriched in disease(s) through 2x2 matrices and multitesting
sub fisher_test 
{
	print "Fisher test for each class/annotation\n\n";

	#Reading file with genes annotations and creating hash (genes = keys / value = annotations)
	my %asso = ();
	open (GENELIST, "$_[0]") or die "cannot open file $_[0]";
	while (my $line = <GENELIST>)
	{
		if ($line =~ /^(.*?)\t(.*)$/)
		{
			chomp $line;
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
		@list = split (/\t/, $values);
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
	my $num2;

	my $j = 0;
	my $f = 1;

	open (FISHER, ">./analysis_$filename/class_$filename") or die "cannot open file ./analysis_$filename/class_$filename";
	open (CLASSFILE, "./data/interactome/HI_022014_PQ_V_CC_simple.clas") or die "cannot open file ./data/interactome/HI_022014_PQ_V_CC_simple.clas"; 
	while (my $line = <CLASSFILE>)
	{
		chomp $line;
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
		my @genelist = split ("\t", $line);
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
				my @ClassAnnot = split(/\t/, $asso{$gene});
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
					@listvalues = ($DatasetPositive, $DatasetNegative, $GenomePositive, $GenomeNegative);
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
				$num2 = $2;
				my $rawp = $3;
				my $bonferroni = $4;
				my $BH = $5;
				# Retrieving annotation name
				if (exists ($References{$num2}))
				{
					# Print only q-value under threshold
					if ($BH <= $QValue) 
					{
						if (@{$details{$num2+1}}[0]/(@{$details{$num2+1}}[0]+@{$details{$num2+1}}[1]) > $homogeneity)
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
					if ($rawpalone <= $QValue) 
					{
						if (@{$details{$num2+1}}[0]/(@{$details{$num2+1}}[0]+@{$details{$num2+1}}[1]) > $homogeneity)
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
	}
	close CLASSFILE;
	close FISHER;
	close H;
	close W;

	print "\n\nFisher test analysis_$filename done\n\n"; 
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
		chomp $line;
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
		chomp $line;
		$line =~/(Class_\d+\|\d+)\t((.*)$)/;
		$GOBP{$1} = $2;
	}

	my %GOMF = ();
	open (MF, "./data/class_GO_enrich_10e-5/GO_MF_Class_10e5") or die "cannot open file ./data/class_GO_enrich_10e-5/GO_MF_Class_10e5";
	while (my $line = <MF>)
	{
		chomp $line;
		$line =~/(Class_\d+\|\d+)\t((.*)$)/;
		$GOMF{$1} = $2;
	}

	my %GOCC = ();
	open (CC, "./data/class_GO_enrich_10e-5/GO_CC_Class_10e5") or die "cannot open file ./data/class_GO_enrich_10e-5/GO_CC_Class_10e5";
	while (my $line = <CC>)
	{
		chomp $line;
		$line =~/(Class_\d+\|\d+)\t((.*)$)/;
		$GOCC{$1} = $2;
	}

	#Annotating file disease class with BP
	open (BPAN, ">./analysis_$filename/bpannotation.csv")  or die "cannot open file ./analysis_$filename/bpannotation";
	open (FILE1, "$_[0]") or die "cannot open file $_[0]";
	while (my $line = <FILE1>)
	{
		chomp $line;
		unless ($line =~ /(Class_\d+\|\d+)/)
		{
			print $line."\n";
			print BPAN "$line\n";
		}
		if ($line =~ /(Class_\d+\|\d+)/)	
		{
			if (exists ($GOBP{$1}))
			{
				my $string = $GOBP{$1};				
				my @ma = ();
				while ($string =~ m/(GO:\d{7}\t\d+\s\d+\s\d+\s\d+\t[\d\-\.a-z]+\t[\d\-\.a-z]+)+/g)	
				{ 
					push (@ma, $&);
				}
				foreach my $mt (@ma)
				{
					if ($mt =~ /(GO:\d{7})(\t\d+\s\d+\s\d+\s\d+\t)([\d\-\.a-z]+\t[\d\-\.a-z]+)/)
					{										
						print $line."\t".$1."\t".$gene_ontology{$1}."\t".$3."\n";
						print BPAN "$line\t$1\t$gene_ontology{$1}\t$3\n";
					}							
				}
			}
			else
			{
				print $line."\tNA\n";
				print BPAN "$line\tNA\n";
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
		unless ($line =~ /(Class_\d+\|\d+)/)
		{
			print $line."\n";
			print MFAN "$line\n";
		}
		if ($line =~ /(Class_\d+\|\d+)/)	
		{
			if (exists ($GOMF{$1}))
			{
				my $string = $GOMF{$1};
				my @ma = ();
				while ($string =~ m/(GO:\d{7}\t\d+\s\d+\s\d+\s\d+\t[\d\-\.a-z]+\t[\d\-\.a-z]+)+/g)
				{ 
					push (@ma, $&);
				}
				foreach my $mt (@ma)
				{
					if ($mt =~ /(GO:\d{7})(\t\d+\s\d+\s\d+\s\d+\t)([\d\-\.a-z]+\t[\d\-\.a-z]+)/)
					{										
						print $line."\t".$1."\t".$gene_ontology{$1}."\t".$3."\n";
						print MFAN "$line\t$1\t$gene_ontology{$1}\t$3\n";
					}						
				}
			}		
			else
			{
				print $line."\tNA\n";
				print MFAN "$line\tNA\n";
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
		unless ($line =~ /(Class_\d+\|\d+)/)
		{
			print $line."\n";
			print CCAN "$line\n";
		}
		if ($line =~ /(Class_\d+\|\d+)/)	
		{
			if (exists ($GOCC{$1}))
			{
				my $string = $GOCC{$1};
				my @ma = ();
				while ($string =~ m/(GO:\d{7}\t\d+\s\d+\s\d+\s\d+\t[\d\-\.a-z]+\t[\d\-\.a-z]+)+/g)
				{ 
					push (@ma, $&);
				}
				foreach my $mt (@ma)
				{
					if ($mt =~ /(GO:\d{7})(\t\d+\s\d+\s\d+\s\d+\t)([\d\-\.a-z]+\t[\d\-\.a-z]+)/)
					{										
						print $line."\t".$1."\t".$gene_ontology{$1}."\t".$3."\n";
						print CCAN "$line\t$1\t$gene_ontology{$1}\t$3\n";
					}						
				}
			}
			else
			{
				print $line."\tNA\n";
				print CCAN "$line\tNA\n";
			}					
		}		
	}	
	close CCAN;
}

#######################
# REACTOME ANNOTATION #
#######################

sub reactome_annotation 
{
	my %reactome = ();
	open (REACTOME, "./data/class_reactome_enrich_10e-5/Class_Reactome_Enrich_10e-5") or die "cannot open file ./data/class_reactome_enrich_10e-5/Class_Reactome_Enrich_10e-5";
	while (my $line =  <REACTOME>)
	{
		chomp $line;
		$line =~ /^(Class_\d+\|\d+)\t((.)*)$/;
		my $class = $1;
		#print $class."\n";
		my $annotations = $2; 
		#print $annotations."\n";
		my @matches = ();		
		while ($annotations =~ m/((\w+)\t\d*\s\d*\s\d*\s\d*\s([0-9\-\.a-z]+\t[0-9\-\.a-z]+))+/g)
		{
			push (@matches, $&);			
		}				
		if ($annotations =~ m/((\w+)\t\d*\s\d*\s\d*\s\d*\s([0-9\-\.a-z]+\t[0-9\-\.a-z]+))+/g)		
		{		
			foreach my $mat (@matches)
			{
				$mat =~ /((\w+)\t\d*\s\d*\s\d*\s\d*\s([0-9\-\.a-z]+\t[0-9\-\.a-z]+))+/;
				$reactome{$class}.= "$2\t$3\;";
			}
		}
		else
		{
			$reactome{$class}.= "NA";
		}			
	}
	close REACTOME;
	
	open (ANNOT, ">./analysis_$filename/disclass_reactome") or die "cannot open file ./analysis_$filename/disclass_reactome\n";
	open (DISCLAS, "$_[0]") or die "cannot open file $_[0]\n";	
	while (my $line =  <DISCLAS>)
	{
		chomp $line;
		if ($line !~ /^\t(Class_\d+\|\d+)/)
		{
			print $line."\n";
			print ANNOT "$line\n";
		}
		else 
		{
			$line =~ /(Class_\d+\|\d+)/;
			if (exists ($reactome{$1}))
			{
				my @notes = split (/\;/, $reactome{$1});
				foreach my $pathway (@notes)
				{
					print $line."\t".$pathway."\n";
					print ANNOT "$line\t$pathway\n"
				}
			}
			else
			{
				print $line."\n";
				print ANNOT "$line\n";
			}
		}
	}	
	close DISCLAS;
	close ANNOT;
}	
		
