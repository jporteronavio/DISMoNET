DISMoNET
========
Diseasome analysis by modules

# GOAL: analysis of diseasome based on modules
######################################################################################################################

# PARAMETERS
######################################################################################################################
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
#####################################################################################################################
# Main folder: "data" --> containing four subfolders
# Subfolder 1: "genes" --> default gene-disease information: Genopedia (May 2014)
############# GENOPEDIA - http://www.hugenavigator.net/HuGENavigator/startPagePedia.do
############# genopedia.txt --> gene_name \t disease_name \t disease_name ...
############# complete disease name is joined by underscore, e.g. diabetes_mellitus
# Subfolder 2: "interactome" --> HI_022014_PQ_V_CC_simple.gr and HI_022014_PQ_V_CC_simple.clas
# Subfolder 3: "class_GO_enrich_10e-5" --> GO_BP_Class_10e5; GO_MF_Class_10e5; GO_CC_Class_10e5; GO_Id_Name.mapping
# Sulfolder 4: "R" --> Rscript_pvalueFisher.R

#####################################################################################################################
