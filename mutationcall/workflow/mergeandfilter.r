
# Results from vardict need to be annotated and further filtered for non-somatic mutations and unreliable calls.
# 1. Merge the vardict and annovar results.
# 2. Filter results based on exac, read depth and alt-allele freq.

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
	stop("Need at least 2 arguments, vardict result file, table_annovar result file. Sample name is optional.,\n", call.=FALSE)
}

vardictfn = args[1]
annovarfn = args[2]
if(length(args)==3){
	samplename=args[3]
}else{
	samplename = tail(strsplit(vardictfn, "/")[[1]], 1)
	samplename = strsplit(samplename, "\\.")[[1]][1]
	cat("Infering sample name from file name: ", samplename, "\n")
}
#scriptdir <- dirname(sys.frame(1)$ofile)#not working
#scriptdir=getSrcDirectory(function(x) {x})#also not working
cmd.args <- commandArgs()
m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
scriptdir <- dirname(regmatches(cmd.args, m))
if(length(scriptdir) == 0) stop("can't determine script dir: please call the script with Rscript")
if(length(scriptdir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
cat("scriptdir=", scriptdir, "=\n")


# vardictfn="160622_UNC17-D00216_0399_BC8L3EANXX_TTAGGC_L001.fixed-rg.sorted_vardict.txt"
# annovarfn="160622_UNC17-D00216_0399_BC8L3EANXX_TTAGGC_L001.fixed-rg.sorted.hg19_multianno.txt"

vardictrestab=read.table(vardictfn, header=FALSE, sep="\t", stringsAsFactors = FALSE)
# vardictjava does not output correct number of columnames by default so this is a retrofit-hack
vdhead=read.table(paste(scriptdir, "/vardictheaders.txt",sep=""), sep="\t", stringsAsFactors = FALSE)
colnames(vardictrestab)=vdhead[,1]

annovarrestab=read.table(annovarfn, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# make match-strings
annovarrestab$matchstring = paste(annovarrestab$Chr, annovarrestab$Start, annovarrestab$Ref, annovarrestab$Alt, sep="_")
vardictrestab$matchstring = paste(vardictrestab$Chr, vardictrestab$Start, vardictrestab$Ref, vardictrestab$Alt, sep="_")
colnames(annovarrestab) = paste("AN_", colnames(annovarrestab), sep="")
if(any(!annovarrestab$AN_matchstring %in% vardictrestab$matchstring) )
	warning("Some results from annovar not found in input. Something wrong? ", annovarrestab$AN_matchstring[!annovarrestab$AN_matchstring %in% vardictrestab$matchstring ])

# start making combined result table
# bigtab is gonna be the combined vardict + annovar table
bigtab = vardictrestab[!duplicated(vardictrestab$matchstring),] # for some reason vardict returns almost a few almost duplicated lines.
a=match(bigtab$matchstring, annovarrestab$AN_matchstring)
bigtab = cbind(bigtab, annovarrestab[a,])


# Filter

# Filter out probable germline variants, those in exac
# I am not sure if dbSNP contains only germline mutations. If it also contains somatic mutations it is wrong to filter out.
# bigtab[, "AN_avsnp144"]

# low in exac. Using absent could filter out real driver mutations found in a low freq in exac, for instanse KRAS12 (25398285 c->a)
bigtab[bigtab$AN_ExAC_ALL==".","AN_ExAC_ALL"]="0"
bigtab[,"AN_ExAC_ALL"] = as.numeric(bigtab[,"AN_ExAC_ALL"])
lowinexacg = bigtab[,"AN_ExAC_ALL"] < 0.0001
# notinexacg = bigtab[,"AN_X1000g2015aug_all"]=="."
cat( "Number of vardict calls with less than one in 10000 in ExAC_ALL: ", sum(lowinexacg), "\n")

# 
# trunc <- function(x) data.frame(lapply(x, substr, 1, 30))
# a = bigtab$AN_Gene.refGene=="KRAS"
# trunc(bigtab[a,c("matchstring", "AN_Gene.refGene", "AN_X1000g2015aug_all", "AN_ExAC_ALL", "AN_cosmic82")])
# are they in dbSNP?
# bigtab[notin1000g, "AN_avsnp144"]  # some sure are

# filter by read-depth and alternat allele freq. Somewhat arbitrary values
mindepth = 0
minaltdepth = 0
reliable = bigtab$Depth > mindepth & bigtab$AltDepth > minaltdepth
cat("Filter using mindepth:",mindepth, ", minaltdepth:", minaltdepth, ". Number of vardict calls above minimums:",sum(reliable) , "\n")

##
ok = lowinexacg & reliable
cat("Total:", sum(ok), ".  Probable somatic mutations for ", samplename, "\n")

bigtab = bigtab[ok,]

bigtabfn = paste(samplename, "_mutations_bigtab.csv", sep="")
write.csv(bigtab, file=bigtabfn, row.names=FALSE)

nicecols = c("Sample", "Chr", "Start", "End", "AN_Gene.refGene",
						"Ref", "Alt", "AF", "Depth", "AltDepth",
						"AN_cosmic82", "AN_ExAC_ALL", "AN_avsnp144", "AN_AAChange.refGene", "AN_Func.refGene")

smallertabfn = paste(samplename, "_mutations_smallertab.csv", sep="")
write.csv(bigtab[,nicecols], file=smallertabfn, row.names=FALSE)

cat("Wrote output tables; ", bigtabfn, " and ", smallertabfn, " for sample:",samplename, "\n")


