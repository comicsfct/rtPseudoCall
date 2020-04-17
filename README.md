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

# Prerequisites

### Snakemake pipeline
		- Snakemake
		- python3
		- pandas
		- docker images: trinityctat/starfusion:1.8.1, comicspt/dogfinder:3.0

### Pseudogenes detection
		- python3
		- pandas
		- bedtools

### Additional tools

#### ARTDeco

	- docker images: comicspt/artdeco:1.0
	- samtools
	
#### Dogcatcher
	- docker images: Dogcatcher
	- bedtools
	
# Installation
	
We suggest pulling the docker images mentioned in the prerequisites. To run the pipeline we suggest installing Anaconda/Miniconda and creating a conda environment on which to install
python3 and Snakemake:
		conda create -n TRTpipeline
		conda activate TRTpipeline
		conda install -c bioconda snakemake


# Usage


## Automated pipeline
	
In order to detect pseudogenes transcribed due to readthrough transcription of upstream genes, we propose an automatic pipeline that employs already existent tools to detect readthrough transcripts spliced/unspliced and depicts subsequently the overlap with pseudogenes.  Before the automatic pipeline can be run, it is necessary to create a global annotation file to ensure that only non-genic regions are considered, where all genes/transcripts that physically overlap will be fused and the most upstream and downstream coordinates will be set as the start and end, respectively. Since DoGFinder elongates readthrough regions until it finds another gene on the annotation, we will first remove the pseudogenes from the annotation file to improve the chances of finding readthrough regions that overlap pseudogenes and then run “Get_loci_annotation” to preprocess the annotation for readthrough detection

	grep -v pseudogene gencode.v31.annotation.gtf  > gencode.v31.annotation.noPseudogenes.gtf
	Get_loci_annotation -out /annotation -gtf gencode.v31.annotation.noPseudogenes.gtf

The steps of the pipeline that consist of running STAR-Fusion, preprocessing alignment files and detecting readthrough regions/genes with DoGFinder (except processing the annotation file, this must be run independently) can be automatically run through a pipeline built in Snakemake, present in TRTdogfinder2.py file, and providing a configuration file in yaml format, as shown in example TRT_configfile.yaml. This pipeline takes the "[sample]_1.fastq" and "[sample]_2.fastq" files as input and outputs for every sample a STAR-Fusion output folder "[sample_starout_fusion]", a bedfile with DoGs coordinates called "Final_Dog_annotation_[sample].bed" and a csv table called "DoGs_rpkm_table_[sample]".  
The configuration file must have the annotation previously built using “Get_loci_annotation”, parameters for STAR-Fusion, DoGFinder’s “Get_DoGs” and “Get_DoGs_rpkm”, docker containers/images for STAR-Fusion and DoGFinder, names of samples, file with chromosome sizes and output directory. Certain fields may be empty (e.g. "Get_DoGs_rpkm”) but must be present nonetheless (e.g “Get_DoGs_rpkm:  “). 
It is recommended to read Snakemake’s documentation but, in summary, the pipeline can be run as:

	snakemake -s TRTdogfinder2.py --configfile TRT_configfile -p  --cores 8 --resources load=50

The flag ‘-s’ requires the snakemake file “TRTdogfinder2.py”, --configfile determines the configuration file, ‘-p’ prints the pipeline’s bash commands to the screen and ‘--cores’ allows the pipeline to do more than one job at the same time (this does not mean, however, that it will only use the established number of processors). “--resources load” determines the number of STAR-Fusion jobs that are run in parallel, this option should be modified according to the available memory (e.g. if you want to run three jobs in parallel, do --resources load=150).

## Additional tools

There is other software that are adequate for detection of unspliced transcripts called Dogcatcher and ARTDeco but are better suited for strand-specific data. As the difference between the tools’ algorithms and criteria leads to a considerable disparity in the results, it is advised to test different programs and, possibly, use a “wisdom of crowds” approach.
 
### Dogcatcher

To run Dogcatcher, we suggest pulling the following docker image: comicspt/Dogcatcher. Then, use the docker container for the following steps.
First, create a file called bedgraphs with bedtools genomecov:

	bedtools genomecov -bg -split -strand + -ibam SRR7537190_plu.BedGraph -g chromsizes.genome >  SRR7537190_plu.BedGraph
	bedtools genomecov -bg -split -strand - -ibam SRR7537190_min,BedGraph -g chromsizes.genome > SRR7537190_min.BedGraph
 

Now, it is necessary to use an ensembl annotation file, which must be preprocessed before further analyses; this works in a similar way to DoGFinder “Get_loci_annotation” step, removing genes inside other genes and keeping track of inside genes and any overlap on either strand. It will create several files in the ensembl annotation folder.

	python Dogcatcher/1.0_Dogcatcher_flatten_gtf.py --annotation_file_with_path annotation/ensembl_annotation

Now, this next step discovers readthrough genes/regions:

	python Dogcatcher/2.0_Dogcatcher.py [--window_size 100 --coverage_percentage 0.95 --BedGraph_input_path --cpus --get_biotypes True] --BedGraph_input_min_strand bedgraphs/SRR7537190_bedGraph_min --BedGraph_input_plu_strand bedgraphs/SRR7537190_bedGraph_plu --output_prefix  CaughtDogs/ --annotation_file_with_path ensembl_annotation/ensembl_gtf

		--window_size
				This is the sliding window size. Default set to 100

		--coverage_percentage
				This is the coverage percent calculated at each window. Default 0.95

		--BedGraph_input_path
				Example: BedGraph_input_files/
				This script will take the plu and min BedGraph files as input and create BedGraph files for each chromosome in a temporary folder.
				You can delete this file when you are finished calculating data. If you run multiple sessions it will find the temporary folder so you do not have to
				remake the BedGraphs for each chromosome every time.

		--get_biotypes
				Default True. Set to False if you dont care about biotype overlap. Although I would leave it on if doing DE
		--output_prefix
				This is the name of the output file
				Example: Dogcatcher_output/
		--cpus
				This is the amount of cores the program will use.

After, you need to filter the resulting readthrough regions to get the longest ones. It can also filter by the shortest but, in this case, we are interested in the longest. Dont forget to add that thing about independent samples

	python Dogcatcher/2.5_Dogcatcher_filter.py --input_prefix CaughtDogs/ --Dogcatcher_plu_strand_list SRR7537190_plu.BedGraph --Dogcatcher_min_strand_list SRR7537190_min.BedGraph --output_prefix CaughtDogs/filtered/

In the final folder, in this case “CaughtDogs/filtered/”, two sub-directories will be present, corresponding to antisense and sense transcripts, with three more folders in each: bed/, csv/ and gtf/. Dogcatcher additionally presents downstream steps for differential expression.



### ARTDeco 

ARTDeco has several steps but can be easily run like this if you don’t need differential expression information:

	ARTDeco -home-dir ARTDECO_DIR -bam-files-dir BAM_FILES_DIR -gtf-file GTF_FILE -cpu NUM_CPU -chrom-sizes-file CHROM_SIZES_FILE -layout [PE/SE] -orientation [Forward/Reverse] -stranded [True/False]

Or like this if you want differential expression information:

	ARTDeco -home-dir ARTDECO_DIR -bam-files-dir BAM_FILES_DIR -gtf-file GTF_FILE -cpu NUM_CPU -chrom-sizes-file CHROM_SIZES_FILE
	
This will output a folder like Examples/ARTDeco_output.


## Detect pseudogenes in readthrough regions

To cross genomic coordinates of readthrough coordinates with pseudogene coordinates in strand-specific data, you can use bedtools window: 

	bedtools window -l 0 -r 30 -sw -sm -a Final_Dog_annotation_SRR7537190.bed -b gencode.v31.annotation.Pseudogenes.gtf > readthroughPseudogenes.bed



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

