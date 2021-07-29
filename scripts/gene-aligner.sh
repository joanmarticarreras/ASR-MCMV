#!/usr/bin/env bash


# Author:       Joan Martí-Carreras
# E-mail:       joan.marti.carreras@gmail.com
# Info:         www.joanmarticarreras.com


# Released under GPL3

# OPTIONS

## Missing options execption case
if [ $# -eq 0 ]
then
        echo "Missing options!"
        echo "Run $0 -h for help"
        echo ""
        exit 0
fi

OPTIND=1

## Option assignation
while getopts f:g:a:i:o:h OPTIONS
do
        case ${OPTIONS} in
		f)
				f=${OPTARG};;
		g)
				g=${OPTARG};;
		a)
				a=${OPTARG};;
		i)
				i=${OPTARG};;
		o)
				o=${OPTARG};;
		h)
				echo "Usage ./dNdS.sh [OPTIONS]"
				echo ""
				echo "-f	Reference fasta file"
				echo "-g	Reference genbank file"
				echo "-a	Reference gff file"
				echo "-i	Input (multi)fasta"
				echo "-o	Output string"
				echo "-h	Visualize this help"
				exit
				;;

	esac
done


mkdir -p ${o}

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $i > ${o}/$i.tsv

while read line
do
	name=$(echo $line | grep ">" | awk '{print $1}' | sed 's/>//g')
	echo $name
	echo $line | tr " " "\n" > ${o}/$name.fasta

	randomreads.sh ref=${o}/$name.fasta coverage=50 minlength=100 maxlength=150 overwrite=t paired=t out=${o}/${name}_1.fastq out2=${o}/${name}_2.fastq

	bwa index $f
	bwa mem -t 96 $f ${o}/${name}_1.fastq ${o}/${name}_2.fastq | samtools view -@ 96 -b - -T $f | samtools sort -@ 96 - > ${o}/${name}.sort.bam

	samtools index ${o}/${name}.sort.bam

	lofreq faidx $f
	lofreq viterbi -f $f ${o}/${name}.sort.bam | samtools sort -@ 96  - > ${o}/${name}.sort.vit.bam
	lofreq alnqual -b -r ${o}/${name}.sort.vit.bam $f | samtools sort -@ 96 - > ${o}/${name}.sort.vit.alnqual.bam
	lofreq index ${o}/${name}.sort.vit.alnqual.bam
	lofreq call -f $f -o ${o}/${name}.vcf --call-indels -b dynamic -s -a 0.001 --use-orphan ${o}/${name}.sort.vit.alnqual.bam

	vcf-annotator.py --output ${o}/${name}_annotated.vcf ${o}/${name}.vcf ${g}

	grep "IsSynonymous=0" ${o}/${name}_annotated.vcf > ${o}/${name}_annotated_nsyn.vcf
	grep "IsSynonymous=1" ${o}/${name}_annotated.vcf > ${o}/${name}_annotated_syn.vcf

	grep -v "#" $a  | awk '{print $3}' | sort -u > ${o}/features.list
	liftoff -p 96 -f ${o}/features.list -copies -g $a -u ${name}_untransferred.gff -o ${o}/${name}.gff ${o}/${name}.fasta $f
	bedtools getfasta -fi ${o}/${name}.fasta -bed ${o}/${name}.gff -s -name -fo ${o}/${name}_genes.fasta
	sed -i "s/\(>.*\)/\1_$name/" ${o}/${name}_genes.fasta

	awk '{print $3}' ${o}/${name}.gff > ${o}/gene.list

	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ${o}/${name}_genes.fasta > ${o}/${name}_genes.tsv

done < ${o}/$i.tsv

bedtools getfasta -fi $f -bed $a -s -name -fo ${o}/reference_genes.fasta

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ${o}/reference_genes.fasta > ${o}/reference_genes.tsv

while read i
do
	echo $i
	grep $i ${o}/*_genes.tsv | sed 's/.*://g' | tr "\t" "\n" > ${o}/${i}_gene.fa
	transeq -sequence ${o}/${i}_gene.fa -frame 1 -table 1 -outseq ${o}/${i}_gene.faa
	mafft --thread 96 ${o}/${i}_gene.faa > ${o}/${i}_gene_aln.faa
done < ${o}/gene.list

