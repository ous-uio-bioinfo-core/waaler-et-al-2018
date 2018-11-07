# waaler-et-al-2018
R scripts and other scripts used to produce figures and tables related to RNASeq in the article "Tankyrase inhibition counteracts resistance to αPD-1 treatment in syngeneic mouse melanoma models" in preparation by Waaler et al.

This part of the paper is an exploratory analysis of the RNA expression of 18 human and one mouse melanoma cell line treated with G007LK. The samples were grouped in different constellation based on known biology and observed clustering. Based on the grouping, differences in gene expression was assessed.

The raw data is fastq files from illuminas RNASeq. The quantification was done using Kallisto. The r-scripts are r-markdown and provides some QC, explanations, exploratory plot and DEG test with the resulting gene lists. All scripts may not be reproducible since most of the raw and some of the processed data is not provided here.

Short description of the important directories


### humananalysis
Only related to the 18 human samples (treated / untreated)
Script used to quantify RNA abundance. 
r-script with some QC plots and a DEG analysis and some additional plots.


### mouseanalysis
Only for the mouse cell line, which had 4 different "treatments" sequenced in triplicates.
Script used to quantify RNA abundance. 
r-script with some QC plots and a DEG analysis and some additional plots.


### human_vs_mouse
Plots made to explore how the mouse correlates with the human cell lines. Mostly hierarchical clusters with heatmaps for subsets of genes and samples.


### mutationcall
The RNASeq data was also run with a restricted mutation call workflow. Known (from Cosmic) cancer mutations were searched for and reported. This was intersected with other known ( external DNA data) mutations.


## Figures

The figures used in the article are found here:

<pre>
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_4a.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_4d.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl13a.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl13b.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl13c.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl13d.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl13e.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl13f.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl17b.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl17c.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl17d.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl17e.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl18a.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl18b.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl19b.pdf
./human_vs_mouse/output_human_vs_mouse_2018-11-07/article_plots_2018-11-07/fig_suppl19c.pdf
./mutationcall/output_mutation_call_2018-11-07/fig_suppl14.pdf
</pre>
ouput from "find . -name fig*.pdf"

There are figures and list that did not make it to the paper, but still is in this repository.




