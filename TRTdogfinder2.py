#for now this only works for paired-end
import os
import pandas as pd
#configfile: 'RTsAndChimeras_configfile.yaml'   change this to automatically recognize config file  
#config file must have:
#samples
#dir
#star fusion container
#star fusion parameters, even if empty
#DogFinder options, even if empty
#dogfinder container
#annotation_ref for dogfinder
#Get_DoGs_rpkm_options
dir=config['dir']
samples=config['samples']
out_dir=config['out_dir']
#samples_dir=config['samples_dir']

rule all:
    input:
        #expand('{out_dir}starout_fusion_{samples}/Aligned.out.bam', samples=config['samples'], out_dir=out_dir),
        #expand('{out_dir}bams/{samples}_sorted.bam', samples=config['samples'], out_dir=out_dir),
        expand('{out_dir}Final_Dog_annotation_{sample}.bed', sample=config['samples'], out_dir=out_dir),   
        #expand('{out_dir}bams/{sample}.bam', sample = config["samples"], out_dir=out_dir),
        expand('{out_dir}DoGs_rpkm_table_{sample}.csv', sample = config["samples"], out_dir=out_dir)
 
rule star_fusion:
    input:
        left='{out_dir}{sample}_1.fastq',
        right='{out_dir}{sample}_2.fastq'
    params:
        star_fusion_args = config['star_fusion_parameters'], #necessary to define --genome_lib_dir [path/to/genome_ctat_lib] and --CPU  or chimeric_junction if you want to use output from previous STAR run
        star_fusion_container = config['star_fusion_container'],
        dir='{out_dir}starout_fusion_{sample}'
    resources:
        load=50
    output:
        bam='{out_dir}starout_fusion_{sample}/Aligned.out.bam',
        #dir='{out_dir}starout_fusion_{sample}'
    priority:10
    run:
        shell('udocker run -v {dir}:{dir} {params.star_fusion_container} /usr/local/src/STAR-Fusion/STAR-Fusion {params.star_fusion_args} --left_fq {dir}{input.left} --right_fq {dir}{input.right} --output_dir {dir}{params.dir}')
        #shell('mv starout_fusion_{sample}/Aligned.out.bam starout_fusion_{sample}/{sample}.bam')

rule prepare_bams:
    input:
        bams="{out_dir}starout_fusion_{sample}/Aligned.out.bam"
    output:
        bams_out = '{out_dir}bams/{sample}.bam'
    priority:9
    run:
        shell('cp {input} {output}')

rule unique_reads:
    #maintain only uniquely mapped reads for downstream analyses
    input:
        '{out_dir}bams/{sample}.bam'
    output:
        '{out_dir}bams/{sample}_unique.bam'
    params:
        config['DoGFindercontainer']
    shell:     
        'udocker run -v {dir}:{dir} {params} samtools view -h -b -q 255 {dir}{input} > {dir}{output}'
        
        


     #sorting and indexing bam files, needed for dogfinder
rule samtools_sort_index:
    input:
        '{out_dir}bams/{sample}_unique.bam'
    params:
        config['DoGFindercontainer']
    output:
        sorted='{out_dir}bams/{sample}_sorted.bam',
       # idx='{out_dir}bams/{sample}_sorted.bam.bai'
    priority:8
    run:
        shell('udocker run -v {dir}:{dir} {params} samtools sort -o {dir}{output.sorted} {dir}{input}')
        shell('udocker run -v {dir}:{dir} {params} samtools index {dir}{output.sorted}')

#all files must be given together to Pre_Process so that it downsamples to the number of reads of the smallest file, enables comparison between samples
pre_process_input_bams= '' 
pre_process_output_bams= '' 


for sample in samples:
    if len(pre_process_input_bams) == 0:
        pre_process_input_bams = pre_process_input_bams + '{}{}bams/{}_sorted.bam'.format(dir,out_dir, sample)
        pre_process_output_bams = pre_process_output_bams + dir + out_dir
    else:
         pre_process_input_bams = pre_process_input_bams + ',' + '{}{}bams/{}_sorted.bam'.format(dir,out_dir, sample)
         pre_process_output_bams = pre_process_output_bams +  ',' + dir + out_dir

assert (len(pre_process_input_bams)) > 1

Q = len(config['samples'])

rule DogFinder_pre_process:            #see if you can define the output to send bams to<
    input:
        expand('{out_dir}bams/{sample}_sorted.bam', sample=samples, out_dir=out_dir)
    params:
        ref = config['annotation_ref'],
        container = config['DoGFindercontainer'],
        input = pre_process_input_bams,
        output = pre_process_output_bams
    priority:7
    output:
        expand('{out_dir}{sample}_sorted.sorted_PE.sorted_DS.bam', sample=samples, out_dir=out_dir)
    shell:
        'udocker run -v {dir}:{dir} {params.container} Pre_Process -Q {Q} -bam {params.input} -ref {dir}{params.ref} -out {params.output}'

print(config['DogFinderOptions'])

if config['DogFinderOptions'] == None:
    DogFinderOptions = ' '
else:
    DogFinderOptions = config['DogFinderOptions']

if config['Get_DoGs_rpkm_options'] == None:
    Get_DoGs_rpkm_options = ' '
else:
    Get_DoGs_rpkm_options = config['Get_DoGs_rpkm_options']


rule Get_Dogs:
    input:
        '{out_dir}{sample}_sorted.sorted_PE.sorted_DS.bam'
    params:
        ref = config['annotation_ref'],
        options = DogFinderOptions,
        container = config['DoGFindercontainer'],
       # suff={sample}
    output:
        '{out_dir}Final_Dog_annotation_{sample}.bed'
    priority:6
    shell:
        #shell('mkdir {out_dir}')
        'udocker run -v {dir}:{dir} {params.container} Get_DoGs {params.options} -out {dir}{out_dir} -bam {dir}{input} -suff {wildcards.sample} -a {dir}{params.ref}'

rule Get_DoGs_rpkm:
    input:
        bam = '{out_dir}{sample}_sorted.sorted_PE.sorted_DS.bam',
        dogs= '{out_dir}Final_Dog_annotation_{sample}.bed'
    params:
        options = Get_DoGs_rpkm_options,
        container = config['DoGFindercontainer'],
        genome_chromsizes = config['chromsizes']
    output:
        '{out_dir}DoGs_rpkm_table_{sample}.csv'
    priority:5
    shell:
        'udocker run -v {dir}:{dir} {params.container} Get_DoGs_rpkm {params.options} -out {dir}{out_dir} -bam {dir}{input.bam} -dog {dir}{input.dogs} -g {dir}{params.genome_chromsizes} -suff {wildcards.sample}' 

