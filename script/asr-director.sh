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

while read i
do
	echo $i
	mkdir -p ${o}/${i}/
	cd ${o}/${i}

	grep "GU305914\|KY348373\|AM886412" ${o}/tmp_genes.tsv | sed 's/\t/\n/g' |  sed 's/(-)//g' |  sed 's/(+)//g' |awk '/^>/ {$0=$0"Lab"}1' > ${o}/${i}/${i}_lab_strains.fa
	grep "HE610451\|HE610452\|HE610453\|HE610454\|HE610455\|HE610456\|MG957497\|MH118555|\MH118556|\MH118557|\MH118558\|MCMV_HaNa1_231348\|EU579859\|EU579861\|EU579860" ${o}/tmp_genes.tsv | sed 's/\t/\n/g' | sed 's/(-)//g' |  sed 's/(+)//g' | awk '/^>/ {$0=$0"WT"}1' > ${o}/${i}/${i}_WT_strains.fa
	cd ${o}/${i}
	cat ${i}_lab_strains.fa ${i}_WT_strains.fa > ${i}_analysis_set.fa
	iqtree -keep-ident -s ${i}_analysis_set.fa -bb 10000 -redo -m TEST -pre ${i}_iqtree
	bioseq --remove-stop ${i}_analysis_set.fa | sed 's/Lab//g' | sed 's/WT//g' | sed 's/wt-//g' | sed 's/lab-//g' > ${i}_analysis_set_nonSTOP.fa
	sed 's/Lab/{Lab}/g' ${i}_iqtree.contree | sed 's/WT/{WT}/g' | sed 's/wt_//g' | sed 's/lab_//g' > ${i}_iqtree_annotated.contree

	hyphy contrast-fel --code Universal --alignment  ${i}_analysis_set_nonSTOP.fa --tree ${i}_iqtree_annotated.contree --branch-set "Lab" --branch-set "WT" --srv Yes --permutations Yes --p-value 0.05 --q-value 0.1 --output ${i}.json > ${i}_report.txt

done < ${o}/gene.list
