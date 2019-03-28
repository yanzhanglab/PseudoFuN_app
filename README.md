# PseudoFuN
PseudoFuN is a novel database and query tool for homologous PseudoGene and coding Gene (PGG) families.
It supports dynamic search, graphical visualization and functional analysis of pseudogenes and coding genes based on the PGG families. 
This work sets a start point for functional analysis of potentially regulatory pseudogenes.

This Shiny app supports queries of pseudogene/coding gene names. It is implemented by Travis S Johnson and Zhi Huang. 
The underelying pseudogene alignment database is generated by Travis S Johnson and Sihong Li using the Ohio Supercomputer Center (OSC). 
A second query tool which allows for direct sequence queries is also available through the OSC developed by Eric Franz and Travis S Johnson. 

The online Shiny app is still under improvement. The app has been tested on OSX El Captitan 10.11.6 and R version 3.4.2. 
The website version has been tested in Chrome, Firefox, and Safari.

Travis Johnson, Sihong Li, Eric Franz, Zhi Huang, Shuyu Dan Li, Moray J Campbell, Kun Huang, Yan Zhang. PseudoFuN: a resource to derive functional potentials in pseudogenes. Submitted.

PseudoFuN DB Search
https://integrativeomics.shinyapps.io/pseudofun_app/

## What does PseudoFuN do?
Pseudogene-gene families are sets of homologous genes and pseudogenes
based on sequence similarity. This app displays a set of the more than 26000
pseudogene-gene families that were generated using our pipeline. This database allows a
user to query gene symbols, Ensembl IDs and Entrez IDs in our database. The
pseudogene-gene family networks containing these queries are then displayed interactively
using networkD3. The genes in these networks can also be used to run Gene Ontology analysis
with topGO to assess the functional potentials of these pseudogenes. 
A more comprehensive query tool is also developed in association with the Ohio Supercomputer Center which
allows users to not only query gene names but also query sequences and then include the new query into
the assigned PGG network. The query will be performed on GPU-based clusters. This application also allows users to interactively view the GO
information for genes in the network individually, run multiple queries simultaneously and
download the results to a csv file.  For access to the more expansive tool please contact
Travis S Johnson (travis.johnson@osumc.edu) or Dr. Yan Zhang (yan.zhang@osumc.edu).

## Running PseudoFuN
### Network Analysis
After visiting our URL you must select the Database you wish to use and input a gene or
pseudogene to query. Once these are selected clicking on the tabs will display the networks
for pseudogene-gene families that contain the respective gene/pseudogene. Gene symbols,
Ensembl Gene IDs (for genes), Ensembl Transcript IDs (for pseudogenes), and Entrez IDs are
supported. The resulting graphs when nodes are selected will display gene/pseudogene, the
gene symbol if available, and the Ensembl ID.

### GO Analysis (Long runtime):
If the Run GO Analysis option is selected then all of the genes within the networks will be
aggregated into a single gene list. This gene list will then be used as input against all
other genes into the topGO R package. The resulting table displays GO terms with their
associated p-values and number of assigned genes. The GO term set must also be selected.
Options include Biological Process, Molecular Function or Cellular Component. The default
method is using a Fisher's exact test but other options include using Kolmogorov-Smirnov (KS
Classic), and Kolmogorov-Smirnov Elim (KS Elim, Alexa et. al. 2016). The table will always
be sorted by the Fisher's exact test p-value. The GO terms that have no assigned genes are removed
unless the "Include GO terms without any assigned genes" option is selected.

### Genomic Mapping analysis
All of the genes that are contained in the networks identified by the query can be mapped
back to the genome. These graphs are meant to show if the genes/pseudogenes localize to
specific regions of the genome and are displayed as a circular plot.

### TCGA mRNA and miRNA Expression analysis
For each network that is generated, the gene expression information can be displayed for each
cancer in the TCGA. The gene expression values were taken downloaded from Broad GDAC Firehose
and the pseudogene expression values were taken from the dreamBase pseudogene quantification
of TCGA reads (Zheng 2017). The miRNA correlation information was also used from Broad GDAC
Firebrowse. The  result is that for each network's nodes the tumor sample gene/PseudoGene
coexpression network, normal sample gene/pseudogene coexpression network, differential gene/pseudogene
expression between tumor and normal, and miRNA correlation with gene/pseudogenes are displayed
in the TCGA Expression tab. The gene/pseudogene expression matrix and miRNA correlations
are also available for download.

This application is still being tested. For any bugs please contact (travis.johnson@osumc.edu).
