======================================================================
HOW TO USE DISMONET.PL
======================================================================

WHAT IS DISMONET.PL?

Dismonet.pl is a script developed in Perl (v 5.10.0) to analyze the genetic background of human diseasome. It is based on protein-protein interactions of genes related to diseases clustering into modules by an overlapping clustering generator (OCG; http://tagc.univ-mrs.fr/welcome/spip.php?rubrique197). This pipeline creates several result files, which allow to analyze outputs as network.

PROGRAM FOLDER HIERARCHY

It must be created a folder called “dismonet”, which must contained a subfolder called “data”.

# Programme folder: “dismonet”. The Perl script called “dismonet.pl” must be placed in this folder.
####Main subfolder: "data" --> containing four subfolders
######## Subfolder 1: "interactome" --> HI_022014_PQ_V_CC_simple.gr and HI_022014_PQ_V_CC_simple.clas
######## Subfolder 2: "class_GO_enrich_10e-5" --> GO_BP_Class_10e5; GO_MF_Class_10e5; GO_CC_Class_10e5; GO_Id_Name.mapping
######## Subfolder 3: "class_reactome_enrich_10e-5" --> Class_Reactome_Enrich_10e-5
######## Subfolder 4: "R" --> Rscript_pvalueFisher.R

COMMAND LINE PARAMETERS

# $ARGV[0]--> a file with your collection of gene-disease data tabulated in two columns
# gene name (HUGO Gene Nomenclature) \t disease name
# $ARGV[1]--> q-value (e.g. 1.00e-5)
# $ARGV[2]--> homogeneity (number of disease susceptibility genes in a module divided by
# the total number of genes in the module); it is a cut-off selecting disease 
# enriched modules with homogeneity above the threshold established as parameter
# eg: 0.1 (10%), 0.25 (25%), 0.30 (30%) ... (0 in case of non-applicable)
# $ARGV[3]--> name to identify result files 

Command line example: perl dismonet.pl gene-disease_file 1.00e-5 0.25 result_filename

PIPELINE STEPS

Gene-disease data
First step is to built a tabulated file with the gene name (HUGO Gene Nomenclature) and its related disease (each gene-disease in a single row). This information can be retrieved from several sources as OMIM, Genotator, Genopedia, etc., based on different studies (experimental biology, GWAS, etc.), which show various degrees of association. This gene-disease flat file is the first parameter in the command line when you execute the perl script.

Human interactome data
Regarding data about protein-protein interaction (PPI), an interactome can be built from experimental data considering human, direct, and binary interactions. An excellent alternative is the Proteomics Standard Initiative Common QUery InterfaCe (PSICQUIC) created by the Human Proteome Organization Proteomics Standards Initiative (HUPO-PSI) to enable computational access to molecular-interaction data resources by means of a standard Web Service and query language. This tool can be accessed in the following address: (http://www.ebi.ac.uk/Tools/webservices/psicquic/view/main.xhtml).

Interactome data are clustered into modules using OCG, which allow that a protein can be included in several modules. Only those modules with more than two proteins are considered.

It must be created a subfolder called "interactome" into the previously created folder called “data”. In the given example the human, direct and binary PPI data is stored in the file called “HI_022014_PQ_V_CC_simple.gr” and the clustered interactome in the file called “HI_022014_PQ_V_CC_simple.clas”

Genes related to diseases and human interactome
From previous files can be selected those human genes that produce proteins presenting binary interactions between them. This generates a tabulated file with gene-gene relations based on PPI called “gene-ppi_filename” and a file with all pair genes involved in the human interactome called “genes_interactome”. 

The results of the pipeline are identified through a filename added to the outputs. This filename is given as the fourth parameter in the command line after gene-disease file, q-value threshold, and homogeneity value.

Fisher´s Exact Test
A contingency table is created in order to determine if an interactome module is enriched with a specific disease. This table compares disease genes localized in a module versus their presence in the rest of the interactome performing a Fisher´s exact test. To perform the Fisher test is needed a R script called “pvalueFisher.R”. It must be created a subfolder called “R”, into the folder “data”, to store this script. 

A file is created with the interactome modules statistically significant enrich in diseases. This file is called “class-filename” and show the number of the module with the number of genes that contains this module following the number of genes of each cell of the contingency table and both, p and q-value (e.g.: Class_1|24 Lung_Neoplasms 3 21 25 11904 2.25459340560015e-05 0.00014866910616564). Further, it is generated a file called “disease_class_filename”, which shows all the modules related to a specific disease.

The program enables to chose the q-value threshold to determine the disease-module association statistically significant introducing a value (e.g. 1.00e-5) as second parameter. Interestedly, homogeneity (number of disease susceptibility genes in a module divided by the total number of genes in the module) can be also established as a threshold when its value is introduced as third parameter, e.g.: 0.1 (10%), 0.25 (25%), 0.30 (30%)... and 0 in case of non-applicable.

Biological annotation
Modules enriched in disease genes are annotated with their biological functions using information from Gene Ontology (GO) and Reactome. Four files are generated: disclass_reactome (Reactome), bpannotation.csv (GO Biological Process), ccannotation (GO Cellular Component) and mfannotation (GO Molecular Function). Previously, it is defined through a Fisher´s test the modules enrichment (q-value < 10e-5) in specific biological functions storing these data in the Subfolder 2 (“GO_BP_Class_10e5”; “GO_MF_Class_10e5”; “GO_CC_Class_10e5”; “GO_Id_Name.mapping”) and Subfolder 3 ("class_reactome_enrich_10e-5").

Network
A tabulated file (“network_filename”) including disease name and module identification is generated. This file can be uploaded as input file in a network generator program as Cytoscape for analysis and visualization.

Experiment results
A file called “results_experiment.txt” summarizes the number of genes and modules related to the different diseases.

Ranking of diseases
Diseases are paired taking into account their shared modules (“Common_modules.txt” and “Pair_diseases_mod.txt”). Two similarity coefficients (Jaccard and Czekanovski-Dice) are used to measure distances among paired diseases (“Module_distances.txt”). Finally, diseases are ranking by the number of shared modules (“Ranking_diseases_modules.txt”) and the Jaccard distance (“Ranking_diseases_Jaccard.txt”).
