# rtPseudoCall
Readthrough Transcribed Pseudogenes

Material for bookchapter XXX

Project Name

	Acronym??? - Readthrough Transcribed Pseudogenes 


# Description

	Framework to detect pseudogenes transcribed due to transcription readthrough of upstream genes, employing already existent bioinformatics tools with parameters adjusted for pseudogenes inspection.  


# Table of Contents

Installation
Usage
Prerequisites

#Prerequisites
		Snakemake pipeline
		- Snakemake
		- python3
		- pandas
		- docker images: trinityctat/starfusion:1.8.1, comicspt/dogfinder:3.0

		filter_pseudogenes.py
		- python3
		- pandas
		- bedtools


# Installation
	


# Usage
##Snakemake pipeline
	
The steps of the pipeline that consist of running STAR-Fusion, preprocessing alignment files and detecting readthrough regions/genes with DoGFinder (except processing the annotation file, this must be run independently) can be automatically run using the snakemake file TRTdogfinder2.py and providing a configuration file in yaml format, as shown in example TRT_configfile.yaml. This pipeline takes the [sample]_1.fastq and [sample]_2.fastq files as input and outputs for every sample a star_fusion output folder [sample_starout_fusion], a bedfile with DoGs coordinates called Final_Dog_annotation_[sample].bed and a csv table called DoGs_rpkm_table_[sample]. 




##Detect pseudogenes in readthrough regions

To cross readthrough regions’ genomic coordinates with pseudogenes, you can use the following script. This script has the following dependencies: pandas, bedtools and python3.
	usage: filter_pseudogenes.py [-h] [-f FILE [FILE ...]]
                             [-o OUTPUT_PREFIX [OUTPUT_PREFIX ...]]
                             [-p PSEUDO_ANNOT]
                             [-r PSEUDOS_REMOVED_FROM_ANNOTATION]
                             [-w WINDOW_SIZE] [--rpkm RPKM [RPKM ...]]

	Cross readthrough genomic coordinates with pseudogene coordinates and rpkm information

	optional arguments:
	  -h, --help            show this help message and exit
	  -f , --file FILE 
							bed file or files with readthrough coordinates, output
							of dogfinder
	  -o , --output_prefix OUTPUT_PREFIX 
							output prefix for the output each file, must be as
							many as number of files and in same order
	  -p, --pseudo_annot PSEUDO_ANNOT
							gtf with pseudogenes which you want to know if are
							inside readthrough region or not
	  -r, --pseudos_removed_from_annotation PSEUDOS_REMOVED_FROM_ANNOTATION
							file with the pseudogene lines removed from gtf, you
							want to know if pseudogenes are inside pseudogenes on
							the opposite strand
	  -w, --window_size WINDOW_SIZE
							size of window for bedtools window in first step, default is 30bp
	  --rpkm RPKM 
							output of dogfinder with rpkm values


This python script works as a user-friendly wrapper to bash script filter_pseudogenes.sh, so it is necessary to have both in the same folder. If you have readthrough region coordinates from Get_DoGs and the corresponding output from Get_DoGs_rpkm, you can do as follows:

	python filter_pseudogenes.py -f Final_DoG_annotation_SRR7537190  -p gencode.v31.Pseudogenes.gtf --rpkm DoGs_rpkm_table_SRR7537190.csv --output_prefix pseudogenes_SRR7537190

Depending on the annotation with pseudogenes, you might need to remove some columns so that If you only want to cross readthrough regions and pseudogene coordinates, you can use the script as follows:

	python filter_pseudogenes.py -f Final_DoG_annotation_SRR7537190  --rpkm DoGs_rpkm_table_SRR7537190.csv --output_prefix pseudogenes_SRR7537190 -p gencode.v31.Pseudogenes.gtf 

And later, you can cross the information with rpkm information. This will output a file called {output_prefix}_final_results.csv:

	python filter_pseudogenes.py -f Final_DoG_annotation_SRR7537190  --rpkm DoGs_rpkm_table_SRR7537190.csv --output_prefix pseudogenes_SRR7537190

You can also use a different annotation with pseudogenes than the one with the pseudogenes you removed from the original annotation file as long as you supply these as well:

	python filter_pseudogenes.py -f Final_DoG_annotation_SRR7537190  -p gencode.v31.AlternativePseudogenes.gtf --rpkm DoGs_rpkm_table_SRR7537190.csv --output_prefix pseudogenes_SRR7537190 -r gencode.v31.Pseudogenes.gtf


This script first crosses the coordinates of readthrough regions with the coordinates of pseudogenes and reports readthrough regions that intersect with a pseudogene, along with the pseudogene with which it intersects. It uses bedtools window, outputting a file called [output_prefix]_RTs_ovlap_pseudos.bed.This is then filtered for readthrough regions that overlap pseudogenes on the same strand, outputted to a file called [output_prefix]_RTs_same_strand_ovlap_pseudos.bed.
After that, it detects pseudogenes which are between two opposite strand readthrough genes in unstranded libraries and thus it is difficult to understand if the reads align to the strand of the pseudogene or to the opposite strand.  In the end, these are present in a file called [output_prefix]_between2RTs.bed and are removed from the final filtered output called [output_prefix]_RTs_samestrand_ovlap_pseudos_filtered.bed.
It then gets pseudogenes which overlap with pseudogenes on the opposite strand. This is because, on data that is not strand-specific, DogFinder stops elongating readthrough regions when it reaches another gene downstream, even if that gene is on the opposite strand. By removing pseudogenes from the annotation, we are allowing the elongation of readthrough regions to continue past a pseudogene. A pseudogene that overlaps with a readthrough region and with another pseudogene on the opposite strand typically would not be reported. However, since we remove the pseudogenes from the annotation and we have no indication to which strand the reads correspond, this script outputs these cases to a different file ([output_prefix]_inside_opposite_strand_pseudogenes.gtf) and removes them from the filtered output. 






# Contributions

	sadsadsad


# Credits

	sadsadsad


# License

	sadsadsad

