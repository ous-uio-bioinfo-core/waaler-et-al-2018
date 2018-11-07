#!/bin/bash

# usage:
# ./variantcall_RNASeq_workflow.sh ../testdata/160622_UNC17-D00216_0399_BC8L3EANXX_TTAGGC_L001.fixed-rg.sorted.bam samplename
# ./variantcall_RNASeq_workflow.sh /Users/vegardnygaard/prosjekter/rprojects/MetActionExplore/div/not_in_github/RNASeq_variantcalling/testdata/160622_UNC17-D00216_0399_BC8L3EANXX_TTAGGC_L001.fixed-rg.sorted.bam samplename
# ./variantcall_RNASeq_workflow.sh /Users/vegardnygaard/prosjekter/rprojects/krauss/mutationcall/LOX-DMSO_sorted.bam LOX-DMSO


# This is an ad-hoc script to get the more reliable driver mutations from a RNASeq-bam file.
# Only meaningful for tumor samples and not exhaustive, that is; many of the real
# mutations will probably be filtered out (or not search for) along with the many false
# ones. Net result wil be a handful (0-10) of likely somatic mutations.
# Made by Vegard Nygaard at the core facility for bioinformatics at IKF-OUS

 
# input: a bam file from tumor RNASeq. and an optional samplename to be used in the output
# output: a short list of probable driver mutations, with annotation from annovar

# programs needed in this script: java, vardictjava, annovar, R
# Several paths are hardcoded and need to be set!
#module load java

SCRIPTFOLDER=$(dirname "$0")

# Vardict JAVA version
VARDICT=/Users/vegardnygaard/programmer/VarDictJava/build/install/VarDict/bin/VarDict
# VARDICT=/net/tsd-evs.tsd.usit.no/p65/data/durable/programs/VarDictJava/build/install/VarDict/bin/VarDict
#using a hard cutoff in vardict for minimum alt-allele frequency
AF_THR=0.1
THREADS=3

# Annovar with at least refGene,1000g2015aug_all, cosmic82
# ANNOVAR=/net/tsd-evs.tsd.usit.no/p65/data/durable/programs/annovar/table_annovar.pl
# ANNOVARDBDIR=/net/tsd-evs.tsd.usit.no/p65/data/no-backup/annovardb/humandb/
ANNOVAR=table_annovar.pl
ANNOVARDBDIR=/Users/vegardnygaard/programmer/annovar/humandb/

# to downalod a annovar database:
# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
# but cosmic need a special script.

# Small R-script to do some filtering and output the final tables.
FILTERSCRIPT="${SCRIPTFOLDER}/mergeandfilter.r"

# Prepared bed file with positions to interrogate by vardict. Reported to have been 
# observed in tumors two or more times.
COSMBED="${SCRIPTFOLDER}/CosmicMutantExportCensus_v82_unique_2ormorehits.bed"

# ref to genome fastafile neede as input for vardict
#REF19=/tsd/p65/data/durable/annotation/ucsc/hg19/hg19.fa
REF19=/Users/vegardnygaard/annotation/hg19/hg19_ucsc/hg19.fa


TUMORBAM=$1
SAMPLENAME=$2
if [ -z "$2" ]
  then
    echo "No SAMPLENAME supplied, creating from bamfile"
    SAMPLENAME=$(basename $TUMORBAM)
	SAMPLENAME="${SAMPLENAME%.*}"
fi

printf "Finding mutations for ${SAMPLENAME}\n"
VARDICTOUTPUT="${SAMPLENAME}_vardict.txt"

COMMAND="$VARDICT -G $REF19 -f $AF_THR -th $THREADS -VS SILENT -N $SAMPLENAME -b $TUMORBAM -c 1 -S 2 -E 3 -g 4 $COSMBED > $VARDICTOUTPUT"
printf "Running command: ${COMMAND}\n" && eval "$COMMAND"
printf "Outputfile from Vardict: ${VARDICTOUTPUT}\n"


# format Vardict-output to annovar-input, i.e use some of the columns. In addition vardicts end-formatting did not work in annovar for insertions
# vardict writes:   chr6	157100413	157100416	G	GCCG
# but annovar needs:chr6	157100413	157100413	G	GCCG
ANNOVARINPUT="${SAMPLENAME}_annovarinput.txt"
awk -F $'\t' '{end=$4+length($6)-1}; {print $3 "\t" $4 "\t" end "\t" $6 "\t" $7}'  $VARDICTOUTPUT > $ANNOVARINPUT
printf "File used as input to annovar: ${ANNOVARINPUT}\n"

# the annotationtables need to be installed in annovar.
COMMAND="$ANNOVAR ${ANNOVARINPUT} ${ANNOVARDBDIR} -buildver hg19 -nastring . -remove -protocol refGene,avsnp144,1000g2015aug_all,exac03,cosmic82,ljb26_all -operation g,f,f,f,f,f --outfile ${SAMPLENAME}"
printf "Running command: ${COMMAND}\n" && eval "$COMMAND"

ANNOVAROUTPUT="${SAMPLENAME}.hg19_multianno.txt"
printf "Outputfile from annovar: ${ANNOVAROUTPUT}\n"

# Merge Vardict and Annovar and filter. Print table with all annotation and table with less but usefull annotation.

COMMAND="Rscript --vanilla $FILTERSCRIPT $VARDICTOUTPUT $ANNOVAROUTPUT $SAMPLENAME"
printf "Running command: ${COMMAND}\n" && eval "$COMMAND"


printf "Finsihed for sample ${SAMPLENAME}\n"

