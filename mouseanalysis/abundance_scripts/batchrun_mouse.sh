#!/bin/bash


# use http://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz from release 92.
# $KALLISTO index -i Mus_musculus.GRCm38.rel92.cdna.all.idx Mus_musculus.GRCm38.rel92.cdna.all.fa.gz 

BASEDIR=/Users/vegardnygaard/prosjekter/rprojects/krauss
KALLISTO=/Users/vegardnygaard/programmer/kallisto/kallisto
FASTQDIR=/Volumes/nobackup/largedata/krauss/mousedata
INDEXDIR=/Users/vegardnygaard/annotation/GRCm38
INDEXNAME=Mus_musculus.GRCm38.rel92.cdna.all
INDEX=${INDEXDIR}/${INDEXNAME}.idx
BOOTSTRAP=30
THREADS=3
OUTPUTDIR=results2
SAFILE=sa_mouse.csv
# building index
# cd $INDEXDIR
# time $KALLISTO index -i ${INDEXNAME}.idx ${INDEXNAME}.fa.gz


# cd $BASEDIR
mkdir -p $OUTPUTDIR

IFS=,
sed -n '2,$p' $SAFILE | while read run_sample_name sample_name run_id cell_line treatment rep filename
do
	#echo $run_sample_name
	FN="${FASTQDIR}/${filename}"
	R1="${FN/_RX_/_R1_}"
	R2="${FN/_RX_/_R2_}"	
	COMMAND="$KALLISTO quant -i $INDEX -b $BOOTSTRAP -t $THREADS -o ${OUTPUTDIR}/${sample_name} $R1 $R2"
	#echo $COMMAND
	printf "\n\ndoing: $COMMAND" && eval "$COMMAND"
done

# time /Users/vegardnygaard/programmer/kallisto/kallisto quant -i /Users/vegardnygaard/annotation/GRCm38/Mus_musculus.GRCm38.rel92.cdna.all.idx -b 30 -t 3 -o results/WNT3 /Volumes/nobackup/largedata/krauss/mousedata/GCF251S-WNT3_S39_R1_001.fastq.gz /Volumes/nobackup/largedata/krauss/mousedata/GCF251S-WNT3_S39_R2_001.fastq.gz





