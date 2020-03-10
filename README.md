# waaler-et-al-2018


<br/>
<br/>
This repository contains R scripts and other scripts used to produce some of the figures and tables derived from the sequencing data in the article "Tankyrase inhibition sensitizes melanoma to PD-1 immune checkpoint blockade in syngeneic mouse models" by Waaler et al to appear in Communications Biology 2020.

This part of the paper is an exploratory analysis of the RNA expression of 18 human and one mouse melanoma cell line treated with G007LK. The samples were grouped in different constellation based on known biology and observed clustering. Based on the grouping, differences in gene expression were assessed.

The raw data is fastq files from Illumina RNASeq. The quantification was done using Kallisto. The r-scripts are r-markdown and provides some QC plots, explanations, exploratory plot and DEG test with the resulting gene lists. All scripts may not be reproducible since most of the raw and some of the processed data is not provided here. 

Short description of the important directories.


### humananalysis
Only related to the 18 human samples (treated / untreated)
Script used to quantify RNA abundance. 
r-script with some QC plots and a DEG analysis and some additional plots.


### mouseanalysis
Only for the mouse cell line, which had 4 different "treatments" sequenced in triplicates.
Script used to quantify RNA abundance. 
r-script with some QC plots and a DEG analysis and some additional plots.


### human_vs_mouse
Plots made to explore how the mouse correlates with the human cell lines. Mostly hierarchical clusters with heat-maps for subsets of genes and samples.


### mutationcall
The RNASeq data was also run with a restricted mutation call workflow. Known (from Cosmic) cancer mutations were searched for and reported. This was intersected with other known ( external DNA data) mutations.


## Figures

The figures used in the article are found here:

<pre>
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_6a.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_6d.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl21a.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl21b.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl21c.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl21d.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl21e.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl21f.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl24b.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl24c.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl24d.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl24e.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl25a.pdf
./human_vs_mouse/output_human_vs_mouse_2020-03-10/article_plots_2020-03-10/fig_suppl25b.pdf
./mutationcall/output_mutation_call_2020-03-10/fig_suppl22.pdf
</pre>
ouput from "find . -name fig*.pdf"

The figures as they appear in the article may be cosmetically altered.
Due to multiple version and submissions of the paper there are there are other figures and list that did not make it to the final paper, but still is in this repository.


## Software used

kallisto v0.44<br/>
Bray, N.L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 34, 525-527 (2016).
<br/>
<br/>
sleuth v0.29<br/>
Pimentel, H., Bray, N.L., Puente, S., Melsted, P. & Pachter, L. Differential analysis of RNA-seq incorporating quantification uncertainty. Nat Methods 14, 687-690 (2017).
<br/>
<br/>
NMF_0.23.6<br/>
Gaujoux, R. & Seoighe, C. A flexible R package for nonnegative matrix factorization. BMC Bioinformatics 11, 367 (2010).
<br/>
<br/>
HISAT2 v2.1.0<br/>
Kim, D., Langmead, B. & Salzberg, S.L. HISAT: a fast spliced aligner with low memory requirements. Nature Methods 12, 357 (2015).
<br/>
<br/>
VarDict v1.2<br/>
Lai, Z., et al. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Research 44, e108-e108 (2016).
<br/>
<br/>
ANNOVAR v2017-07-17<br/>
Wang, K., Li, M. & Hakonarson, H. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Research 38, e164-e164 (2010).
<br/>
<br/>
R version 3.4.2<br/>
R Core Team (2013). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
URL http://www.R-project.org/.
 <br/>
<br/>
## Databases

Ensembl <br/>
Kinsella, R.J., et al. Ensembl BioMarts: a hub for data retrieval across taxonomic space. Database (Oxford) 2011, bar030 (2011).
<br/>
<br/>
COSMIC v82<br/>
Forbes, S.A., et al. COSMIC: somatic cancer genetics at high-resolution. Nucleic Acids Research 45, D777-D783 (2017).
<br/>


<br/>
<br/>

Contact for this repository: vegards.email@gmail.com or 


