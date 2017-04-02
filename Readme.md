# This directory contains scripts and results from analyses on tissue specific gene expression in GTEx.
 
## Introduction and Methods

GTEx dataset with FPKM normalized expression values from Tophat2 alignments and associated metadata listing for each GTEx sample its tissue name, were downloaded from synapse ([syn8105922](https://www.synapse.org/#!Synapse:syn8105922)) and [syn7596611](https://www.synapse.org/#!Synapse:syn7596611), respectively). Prior to the analysis testis samples were removed due to overrepresentation in the analysis. Tissues belonging to the digestive tract (Liver, Colon, Pancreas, Esophagus, Stomach, Biliary) and the female reproductive system (Breast, Cervix Uteri, Ovary, Uterus, Fallopian Tube, Vagina) were merged into meta-cohorts following guidelines in the PCAWG meta tumor group documentation. 11,148 of the total 57,820 genes having constantly FPKM values < 1 in all GTEx samples were filtered out from further analysis. The FPKM threshold was chosen on the basis of the expression of olfactory receptors, which are likely to have no function in almost all tissues [(Ezkurdia <i>et al.</i>)](http://pubs.acs.org/doi/abs/10.1021/pr500572z) and of which 99% had a smaller expression than 1 FPKM.To determine which of the remaining genes are expressed in a tissue specific manner, we followed the ideas from the [TiGER database](http://bioinfo.wilmer.jhu.edu/tiger/), where tissue specificity was defined by a minimum enrichment score between two tissues and an associated maximum p-value between their samples. The two tissues with the hightest (A) and 2nd highest (B) mean expression were selected for each gene and an enrichment score, equivalent to the quotient A/B was calculated. The enrichment score was supplemented with a p-value from a Wilcox-Rank Sum test comparing the set of expression values for the gene in the samples of the highest expressed tissue (A<sub>i</sub>) and 2nd highest expressed tissue (B<sub>i</sub>). Genes with an enrichment score > 5 and an p-value < 0.01 were regarded as tissue specific, while the rest were considered to be non-tissue specific. 

## This directory has following structure:
   * scripts/
     * *A set of PERL scripts for determining tissue specific gene expression*
   * input/
     * *Input data files*.
   * output/
     *  *Result data files* 
   * README.md
     * *This Readme file*.

## Requirements
`scripts/tissueSpecificGeneExpression.pl` requires the [Statistics::R](http://search.cpan.org/~gmpassos/Statistics-R-0.02/lib/Statistics/R.pm) PERL module. Try to install it via  

`perl -MCPAN -e 'install Statistics::R'`

## Running the analysis
First, remove genes with negligible expression:
`scripts/filterGeneExpression.pl -in GTEx.tophat2.gene.fpkm.tsv.gz -prefix ENSG -min 1 | gzip > input/GTEx.tophat2.gene.fpkm_filt.tsv.gz`

GTEx.tophat2.gene.fpkm.tsv.gz can be downloaded from [Synapse](https://www.synapse.org/#!Synapse:syn8105922).

Then, determine tissue specific gene expression:
`scripts/tissueSpecificGeneExpression.pl -meta GTEX_v4.metadata.tsv.gz -cohort input/metaCohorts.tsv -express GTEx.tophat2.gene.fpkm_filt.tsv.gz -mean -exclude testis | gzip > output/tissueSpecificGeneExpression_meta_mean.tsv.gz`

GTEX_v4.metadata.tsv.gz can be downloaded from [Synapse](https://www.synapse.org/#!Synapse:syn7596611).


Lastly, reformat output and put genes in order of a predefined gene list:
`scripts/assignTissueSpecificity2genes.pl -gene input/genes.txt -tissue output/tissueSpecificGeneExpression_meta_mean.tsv.gz -map input/ensg_ensp_enst_ense_geneName_v75.tsv.gz -min 5 -max 0.01 > output/genes_tissueEnrichmentPvalue.tsv`
