


library(dplyr)
library(tidyr)

version = 82
cosmicfn = paste( "/Volumes/nobackup/annotation/Homo_sapiens/cosmic/v",version,"/CosmicMutantExportCensus.tsv.gz", sep="")
rawcosmictab = read.table(cosmicfn, header=TRUE, stringsAsFactors = FALSE, sep="\t")

modcosmictab = rawcosmictab %>% 
	filter(grepl( "somatic",  Mutation.somatic.status) ) %>%
	filter(Mutation.genome.position!="") %>%
	group_by(Mutation.genome.position, Gene.name) %>%
	summarise( count=n()) %>%
	separate(Mutation.genome.position, c("chr", "start", "end"), convert=TRUE) %>%
	filter(!is.na(start)) %>%
	filter( (end-start) < 4) %>%
	filter( count > 1)

modcosmictab[modcosmictab$chr==23, "chr"]="X" # cosmic uses 23 for X . Will not work in VarDict.
modcosmictab$start = modcosmictab$start-1 # Just to be sure, if a deletion is described to start earlier in the bam file 
write.table(modcosmictab, file=paste("CosmicMutantExportCensus_v", version, "_unique_2ormorehits.bed", sep=""), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

#cosmictab = read.table("dummy.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")



# just a count of the genes and coverage if one uses the filtered cosmic regions.
tmp = modcosmictab  %>%
	group_by(Gene.name) %>%
	summarise( bpcover=n(), cosmicobscount=sum(count)) %>%
	arrange(desc(cosmicobscount, bpcover))
# write.table(tmp, file=paste("CosmicMutantExportCensus_v", version, "_unique_2ormorehits_summary.txt", sep=""), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)




