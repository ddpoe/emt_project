#!/bin/bash

case "$1" in 
	-p)
		if [[ $2 == '*.txt' ]]
		then
			for accession in $(cat $2)
			do
				fasta-dump --split-files $accession
			done
		elif [[ ! -z $2 ]]	
		then
			accession=$2
			fastq-dump --split-files $accession
		else
			echo No file or accession given
		fi
		;;
	-s) 
		if [[ $2 == '*.txt' ]]
		then
			for accession in $(cat $2)
			do
				fasta-dump  $accession
			done
		elif [[ ! -z $2 ]]	
		then
			accession=$2
			fastq-dump  $accession
		else
			echo No file or accession given
		fi
		;;
	"")
		echo
		echo Use $0 -p \<accession file\> for paired end scRNAseq data. This accession file is a text file of all the accession numbers for the samples that will be processed \(one per line\). A single accession number can also be replaced with file. If the samples you want to process are single end reads use the -s flag with the same arguments as the paired end. Files will be outputted in the file where this is command is ran.
		echo
		;;
	*)
		echo
		echo single or paired read flag not defined or defined incorrectly. Use $0 for more information.
		echo
		;;
esac
