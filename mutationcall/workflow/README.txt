
"Finding Common Mutations In Common Cancers From RNASeq Without Germ line"

A short description of the files used by the workflow.

variantcall_RNASeq_workflow.sh: Main script, uses the other files. Calls vardict and annovar.

cosmicfilter.r: Script to read in mutations reporter in Cosmic and make a bed file of the positions with at least two mutations reported.

CosmicMutantExportCensus_v82_unique_2ormorehits.bed: The output from cosmicfilter used as input to vardict to restrict search to those ca 25k positions.

mergeandfilter.r: Merging output from vardict and annovar, and filter out known exac positions.

vardictheaders.txt: Hack to mitigate missing or erroneous headers in the output from Vardict.



The intended use setting is when you have RNASeq for a tumor, but no mutation information from DNA, exome etc.

Using Cosmic-records as filter I managed to get a set of mutations from RNASeq samples that I think is comparable with what one could get from targeted DNA sequencing for instance the cancer hotspot panel2 from IonTorrent.
