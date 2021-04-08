from Bio import SeqIO
import pandas as pd
import os
import argparse
import findcirc_gmap,findcirc_blast,findcirc_minimap,findcirc_blat 
import multiprocessing
parser = argparse.ArgumentParser(description='Predict circRNA')

parser.add_argument("-g",
                    "--genome",
                    help='reference genome files are Fasta')
parser.add_argument("-q",
                    "--query",
                    help='reference query files are Fasta')
parser.add_argument("-o",
                    "--output_file",
                    help='Output_file saved reasult')
parser.add_argument("-t",
                    "--num_threads",
                    default=4)


args = parser.parse_args()

num_threads = args.num_threads
output_file = args.output_file
query_filename = args.query
genome_filename = args.genome



def map(map_method):
    if map_method == "blast":
        os.system("makeblastdb -in %s -dbtype nucl -out %s_blastgenome" %
          (genome_filename, output_file))
        os.system("blastn -db %s_blastgenome -query %s -out %s_blast -outfmt 6 -evalue 1e-5 -num_threads %d" %
          (output_file, query_filename, output_file, num_threads))

    elif map_method == "blat":
        os.system("blat  %s  %s -trimT -trimHardA -maxIntron=1  %s_blat  "%
          (genome_filename, query_filename, output_file))

    elif map_method == "gmap":
        os.system("sudo gmap_build  -d genome_gmap  %s  " %
          (genome_filename))
        os.system("gmap -t %d -d genome_gmap  %s --split-large-introns --max-intronlength-middle=1 --max-intronlength-ends=1 -w 2 -L 2 -n 10 -f psl > %s_gmap.psl   "%
          (num_threads, query_filename,output_file))

    elif map_method == "minimap2":
        os.system("minimap2 -x splice -G 50 -N 10 %s  %s -Y -o  %s_minimap " %
          (genome_filename, query_filename, output_file))

if __name__ == '__main__':
    map_methods = ["blast", "blat","gmap","minimap2"]
    pool = multiprocessing.Pool(len(map_methods))
    for map_method in map_methods:
        pool.apply_async(map, (map_method,))
    pool.close()
    pool.join()

    df_blast = findcirc_blast.exactReads("%s_blast" % output_file)
    print(df_blast.shape)
    df_list_blast = findcirc_blast.cutDatafram(df_blast)
    print(len(df_list_blast))
    pool_blast = multiprocessing.Pool(len(df_list_blast))
    for datafram in df_list_blast:
        print(datafram.shape)
        pool_blast.apply_async(findcirc_blast.getCircInformation, (datafram,output_file,))   
    pool_blast.close()
    pool_blast.join()

    df_blat = findcirc_blat.exactReads("%s_blat" % output_file)
    print(df_blat.shape)
    df_list_blat = findcirc_blat.cutDatafram(df_blat)
    pool_blat = multiprocessing.Pool(len(df_list_blat))
    for datafram in df_list_blat:
        print(datafram.shape)
        pool_blat.apply_async(findcirc_blat.getCircInformation, (datafram,output_file,))   
    pool_blat.close()
    pool_blat.join()

    df_gmap = findcirc_gmap.exactReads("%s_gmap.psl" % output_file)
    print(df_gmap.shape)
    df_list_gmap = findcirc_gmap.cutDatafram(df_gmap)
    pool_gmap = multiprocessing.Pool(len(df_list_gmap))
    for datafram in df_list_gmap:
        print(datafram.shape)
        pool_gmap.apply_async(findcirc_gmap.getCircInformation, (datafram,output_file,))   
    pool_gmap.close()
    pool_gmap.join()

    df_minimap = findcirc_minimap.exactReads("%s_minimap" % output_file)
    print(df_minimap.shape)
    df_list_minimap = findcirc_minimap.cutDatafram(df_minimap)
    pool_minimap = multiprocessing.Pool(len(df_list_minimap))
    for datafram in df_list_minimap:
        print(datafram.shape)
        pool_minimap.apply_async(findcirc_minimap.getCircInformation, (datafram,output_file,))   
    pool_minimap.close()
    pool_minimap.join()

os.system("cat %s_*_* > %s_all" %(output_file,output_file))
infile = '%s_all' % (output_file)
outfile = '%s_circRNA' % (output_file)
df = pd.read_csv(infile, sep=',', header=None)
print(df.shape)
df = df.drop_duplicates(subset=[1 ,2, 3, 4], keep='first')
df.to_csv(outfile, header=None, index=None)

