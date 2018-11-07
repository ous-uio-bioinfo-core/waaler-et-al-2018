#!/bin/bash




DUMMY=/Users/vegardnygaard/prosjekter/rprojects/krauss/mutationcall/test10k.fastq
BAMDIR=/Volumes/nobackup/largedata/krauss/celllinedata/bam
SCRIPTFILE=/Users/vegardnygaard/prosjekter/rprojects/krauss/mutationcall/workflow/variantcall_RNASeq_workflow.sh


####NB,  SETTING IFS to , makes , not work later !!!!!

SAFILE=/Users/vegardnygaard/prosjekter/rprojects/krauss/annotation/cellline_annotation.csv
IFS=,
sed -n '2,$p' $SAFILE | while read cell_line hippo_group_untreated hippo_group_response group_response
do
	unset IFS
	echo $cell_line
	F1="${BAMDIR}/${cell_line}-DMSO.bam"
	F2="${BAMDIR}/${cell_line}-G007-LK.bam"
	MERGED="${BAMDIR}/${cell_line}-merged.bam"
	#COMMAND="samtools merge $MERGED $F1 $F2"
	#COMMAND="samtools index $MERGED"
	COMMAND="$SCRIPTFILE $MERGED $cell_line"
	#echo $COMMAND
	printf "\n\ndoing: $COMMAND" && eval "$COMMAND"
	#COMMAND2="samtools index $BAM"
	#printf "\n\ndoing: $COMMAND2" && eval "$COMMAND2"
	#COMMAND3="ls -l $F2"
	#printf "\n\ndoing: $COMMAND2" && eval "$COMMAND2"
	IFS=,
done

