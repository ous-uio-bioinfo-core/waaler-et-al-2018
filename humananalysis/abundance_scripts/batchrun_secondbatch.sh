#!/bin/bash




BASEDIR=/Users/vegardnygaard/prosjekter/rprojects/krauss
KALLISTO=/Users/vegardnygaard/programmer/kallisto/kallisto
FASTQDIR=${BASEDIR}/data
INDEXDIR=/Users/vegardnygaard/annotation/GRCh38
INDEXNAME=Homo_sapiens.GRCh38.rel91.cdna.all
INDEX=${INDEXDIR}/${INDEXNAME}.idx
BOOTSTRAP=30
THREADS=3
OUTPUTDIR=results_second_run
SAFILE=sa_second_run.csv
# building index
# cd $INDEXDIR
# time $KALLISTO index -i ${INDEXNAME}.idx ${INDEXNAME}.fa.gz


cd $BASEDIR
mkdir -p $OUTPUTDIR

IFS=,
sed -n '2,$p' $SAFILE | while read run_sample_name sample_name run_id cell_line treatment filename
do
	#echo $run_sample_name
	F1="${BASEDIR}/data/demultiplex-H5G22BGX5/GCF-0269S-Krauss-mRNA/${filename}"
	F2="${BASEDIR}/data/demultiplex-H5HYJBGX5/GCF-0269S-Krauss-mRNA/${filename}"
	#F3="${BASEDIR}/data/demultiplex-HG7K5BGX5/GCF-0269S-Krauss-mRNA/${filename}"
	# demultiplex-HCCF2BGX5
	# demultiplex-HG7K5BGX5
	#R2="${FASTQDIR}/${p}*R2_001.fastq.gz"
	COMMAND="$KALLISTO quant -i $INDEX -b $BOOTSTRAP -t $THREADS --single -l 200 -s 20 -o ${OUTPUTDIR}/${sample_name} $F1 $F2"
	#echo $COMMAND
	printf "\n\ndoing: $COMMAND" && eval "$COMMAND"
done

