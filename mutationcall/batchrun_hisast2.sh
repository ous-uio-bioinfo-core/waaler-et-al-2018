#!/bin/bash




HISAT2_HOME=/Users/vegardnygaard/programmer/hisat2-2.1.0
FASTQDIR=/Volumes/nobackup/largedata/krauss
INDEXDIR=/Volumes/nobackup/annotation/Homo_sapiens/hisat2index/hg19/genome
THREADS=3


#${HISAT2_HOME}/hisat2 -f -x ${HISAT2_HOME}/example/index/22_20-21M_snp -1 ${HISAT2_HOME}/example/reads/reads_1.fa -2 ${HISAT2_HOME}/example/reads/reads_2.fa -S eg2.sam
#${HISAT2_HOME}/hisat2 -f -x ${HISAT2_HOME}/example/index/22_20-21M_snp -U ${HISAT2_HOME}/example/reads/reads_1.fa,${HISAT2_HOME}/example/reads/reads_2.fa -S eg3.sam
#${HISAT2_HOME}/hisat2 -f -x $INDEXDIR -U ${HISAT2_HOME}/example/reads/reads_1.fa,${HISAT2_HOME}/example/reads/reads_2.fa -S eg3.sam

DUMMY=/Users/vegardnygaard/prosjekter/rprojects/krauss/mutationcall/test10k.fastq
OUTPUTDIR=/Volumes/nobackup/largedata/krauss/celllinedata/bam

####NB,  SETTING IFS to , makes , not work later !!!!!

SAFILE=/Users/vegardnygaard/prosjekter/rprojects/krauss/humananalysis/abundance_scripts/sa_first_run.csv
IFS=,
sed -n '2,$p' $SAFILE | while read run_sample_name sample_name run_id cell_line treatment filename
do
	unset IFS
	#echo $run_sample_name
	F1="${FASTQDIR}/celllinedata/demultiplex-H5G25BGX5/GCF-0269S-Krauss-mRNA/${filename}"
	F2="${FASTQDIR}/celllinedata/demultiplex-HCCF2BGX5/GCF-0269S-Krauss-mRNA/${filename}"
	F3="${FASTQDIR}/celllinedata/demultiplex-HG7K5BGX5/GCF-0269S-Krauss-mRNA/${filename}"
	BAM="${OUTPUTDIR}/${sample_name}.bam"
	COMMAND="${HISAT2_HOME}/hisat2 -x $INDEXDIR -U ${F1},${F2},${F3} | samtools sort -o $BAM"
	#COMMAND="${HISAT2_HOME}/hisat2 -x $INDEXDIR -U $DUMMY | samtools sort -o $BAM"
	#echo $COMMAND
	#printf "\n\ndoing: $COMMAND" && eval "$COMMAND"
	COMMAND2="samtools index $BAM"
	#printf "\n\ndoing: $COMMAND2" && eval "$COMMAND2"
	#COMMAND3="ls -l $F2"
	#printf "\n\ndoing: $COMMAND2" && eval "$COMMAND2"
	IFS=,
done


SAFILE=/Users/vegardnygaard/prosjekter/rprojects/krauss/humananalysis/abundance_scripts/sa_second_run.csv
IFS=,
sed -n '2,$p' $SAFILE | while read run_sample_name sample_name run_id cell_line treatment filename
do
	unset IFS
	#echo $run_sample_name
	F1="${FASTQDIR}/celllinedata/demultiplex-H5G22BGX5/GCF-0269S-Krauss-mRNA/${filename}"
	F2="${FASTQDIR}/celllinedata/demultiplex-H5HYJBGX5/GCF-0269S-Krauss-mRNA/${filename}"
	BAM="${OUTPUTDIR}/${sample_name}.bam"
	COMMAND="${HISAT2_HOME}/hisat2 -x $INDEXDIR -U ${F1},${F2} | samtools sort -o $BAM"
	#COMMAND="${HISAT2_HOME}/hisat2 -x $INDEXDIR -U $DUMMY | samtools sort -o $BAM"
	#echo $COMMAND
	printf "\n\ndoing: $COMMAND" && eval "$COMMAND"
	COMMAND2="samtools index $BAM"
	printf "\n\ndoing: $COMMAND2" && eval "$COMMAND2"
	IFS=,
done


