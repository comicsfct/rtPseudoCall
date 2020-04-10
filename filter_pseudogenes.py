import os
import argparse
import pandas as pd
import subprocess 

parser= argparse.ArgumentParser(description="Cross readthrough genomic coordinates with pseudogene coordinates and rpkm information")

parser.add_argument('-f','--file', nargs="+", default=None, type=str, help="bed file or files with readthrough coordinates, output of dogfinder")
parser.add_argument('-o', '--output_prefix', nargs="+", type=str, help="output prefix for the output each file, must be as many as number of files and in same order")
parser.add_argument('-p', '--pseudo_annot', default=None, help="gtf with pseudogenes which you want to know if are inside readthrough region or not")
parser.add_argument('-r', '--pseudos_removed_from_annotation', default=None, help="file with the pseudogene lines removed from gtf, you want to know if pseudogenes are inside pseudogenes on the opposite strand")
parser.add_argument('-w', '--window_size', default= '30', help="size of window for bedtools window in first step, default is 30bp")
parser.add_argument('--rpkm', nargs="+", default=None, help="output of dogfinder with rpkm values")

args = parser.parse_args()

def remove_parenthesis(string):
    
    new_string = ''
    for idx, i in enumerate(string):
        if i != "\"":
            new_string += i

    return(new_string)


def get_gene_id_name(info):
#transform all the annotation info for the genes into only id and name 
    sep_info = info.split(";")

    fields_idx = []

    for idx, field  in enumerate(sep_info):
        if 'gene_id' in field or 'gene_name' in field:
            fields_idx.append(idx)

    final_info = [remove_parenthesis(sep_info[x].split()[1]) for x in fields_idx]

    return(final_info)


def merge_outputs(filtered, rpkm, same_strand):
    #load files
    df_filtered = pd.read_table(filtered, header=None, sep='\t')
    df_rpkm = pd.read_csv(rpkm)
    df_same_strand = pd.read_csv(same_strand, header=None, sep='\t')

    df_same_strand.iloc[:,10] = df_same_strand.iloc[:,10].apply(get_gene_id_name)
    
    # new data frame with split value columns 
    df_same_strand.rename(columns={df_same_strand.columns[10]:'info'}, inplace=True)
    
    new = pd.DataFrame(df_same_strand['info'].tolist(), columns=['gene_id', 'gene_name'])

    # making separate columns from new data frame 
    df_same_strand['gene_id'] = new['gene_id']
    df_same_strand['gene_name'] = new['gene_name']

    # Dropping old info column 
    df_same_strand.drop(columns = ["info"], inplace = True) 


    df_filtered.iloc[:,4] = df_filtered.iloc[:,4].apply(get_gene_id_name)
    

    df_filtered.rename(columns={df_filtered.columns[4]:'info'}, inplace=True)

    new2 = pd.DataFrame(df_filtered['info'].tolist(), columns=['gene_id', 'gene_name'])

    df_filtered['gene_id'] = new2['gene_id'] 
    df_filtered['gene_name'] = new2['gene_name']

    df_filtered.drop(columns = ["info"], inplace = True) 
    
    #will hold existent pseudos for each RT
    dict_RT_pseudos = {}

    df_same_strand = df_same_strand.drop_duplicates(subset = "gene_name")
    df_filtered = df_filtered.drop_duplicates(subset="gene_name")
    

    df_final = pd.DataFrame(columns = ['RT_gene_id', 'chr_RT', 'start_RT', 'end_RT', 'RT_length', 'RT_strand', 'RT_rpkm',  'chr_pseudo', 'start_pseudo', 'end_pseudo', 'strand_pseudo', 'gene_id_pseudo', 'gene_name_pseudo'])

   #append to dict_RT_pseudos
    for idx, row in df_same_strand.iterrows():
        
        if row['gene_name'] in df_filtered['gene_name'].tolist():

            if row[3] not in dict_RT_pseudos.keys():
                dict_RT_pseudos[row[3]] = list()
                dict_RT_pseudos[row[3]].append(row['gene_name'])
            else:
                dict_RT_pseudos[row[3]].append(row['gene_name'])

    #create final df with the combined information from the 3 datasets
    for gene_id in dict_RT_pseudos.keys():

        for pseudo in dict_RT_pseudos[gene_id]:

            
            row_rpkm = df_rpkm.loc[df_rpkm['DoG_name'] == gene_id]
            
            row_pseudo = df_filtered.loc[df_filtered['gene_name'] == pseudo]
           
            row_rpkm = row_rpkm.squeeze().tolist()
            row_pseudo = row_pseudo.squeeze().tolist()
            
            row_to_append = row_rpkm + row_pseudo
            row_to_append = row_to_append[1:]
            
            
            df_final.loc[len(df_final)] = row_to_append
            
    
    return(df_final)

def commands(output_prefix, file, annot, pseudos_removed_from_annot , window):
    #runs bash script
    com = f"bash ./filter_pseudogenes.sh {annot} {file} {output_prefix} {pseudos_removed_from_annot} {window}"
    os.system(com)

annot = args.pseudo_annot
window_size = args.window_size

#if no annotation is provided for pseudogenes removed from the annotation
#  it will assume it is the same annotation as provided to intersect for pseudogenes

if args.pseudos_removed_from_annotation is None:
    pseudos_removed_from_annotation = annot
else:
    pseudos_removed_from_annotation = args.pseudos_removed_from_annotation

#if last file from dogfinder with the rpkm is not provided it just runs with the Final_Dog_annotation(RT coordinates file)
if args.rpkm == None:
    for idx, file in enumerate(args.file):

        commands(args.output_prefix[idx], file, annot, pseudos_removed_from_annotation, window = window_size)

#when bash script was already run cross the information and output final_results.csv with the same prefix
elif args.file == None:
    for idx, file in enumerate(args.rpkm):
        df = merge_outputs(f"{args.output_prefix[idx]}_samestrand_ovlap_pseudos_filtered.bed", \
        args.rpkm[idx], f"{args.output_prefix[idx]}_RTs_same_strand_ovlap_pseudos.bed")

        df.to_csv(f"{args.output_prefix[idx]}_final_results.csv")

#when a rpkm file is provided for each RT coordinates file
elif len(args.rpkm) == len(args.file):
    for idx, file in enumerate(args.file):
    
        commands(args.output_prefix[idx], file, annot, pseudos_removed_from_annotation, window_size)

        #check that commands worked
        try:
            assert os.path.isfile(f'{args.output_prefix[idx]}_RTs_same_strand_ovlap_pseudos.bed')
            assert os.path.isfile(f'{args.output_prefix[idx]}_samestrand_ovlap_pseudos_filtered.bed')
            
        except AssertionError as error:
            print("Output files weren't created")
            break

        df = merge_outputs(f"{args.output_prefix[idx]}_samestrand_ovlap_pseudos_filtered.bed", \
        args.rpkm[idx], f"{args.output_prefix[idx]}_RTs_same_strand_ovlap_pseudos.bed")

        df.to_csv(f"{args.output_prefix[idx]}_final_results.csv")

else:
    print("Either supply as many rpkm files as readthrough coordinates files or don't supply any or supply only rpkm files and output prefixes.")

