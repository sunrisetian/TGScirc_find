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
        os.system("blastn -db %s_blastgenome -query %s -out %s_blast -outfmt 6 -evalue 1e-5 -num_thread %d" %
          (output_file, query_filename, output_file, num_threads))
        df = findcirc_balst.exactReads("%s_blast" % output_file)
        print(df.shape)
        df_list = findcirc_balst.cutDatafram(df)
        pool = multiprocessing.Pool(len(df_list))
        for datafram in df_list:
            print(datafram.shape)
            pool.apply_async(findcirc_balst.getCircInformation, (datafram,ouput_file,))    
        pool.close()
        pool.join()
    elif map_method == "blat":
        os.system("blat  %s  %s -trimT -trimHardA -maxIntron=1  %s_blat  "%
          (genome_database, query_filename, output_file))
        df = findcirc_blat.exactReads("%s_blat" % output_file)
        print(df.shape)
        df_list = findcirc_blat.cutDatafram(df)
        pool = multiprocessing.Pool(len(df_list))
        for datafram in df_list:
            print(datafram.shape)
            pool.apply_async(findcirc_blat.getCircInformation, (datafram,ouput_file,))    
        pool.close()
        pool.join()
    elif map_method == "gmap":
        os.system("sudo gmap_build  -d %s_gmap  %s  " %
          (genome_database, genome_filename))
        os.system("gmap -t %s -d %s_gmap  %s --split-large-introns --max-intronlength-middle=1 --max-intronlength-ends=1 -w 2 -L 2 -n 10 -f psl > %s_gmap   "%
          (genome_database, query_filename,num_threads,output_file))
        df = findcirc_gmap.exactReads("%s_gmap" % output_file)
        print(df.shape)
        df_list = findcirc_gmap.cutDatafram(df)
        pool = multiprocessing.Pool(len(df_list))
        for datafram in df_list:
            print(datafram.shape)
            pool.apply_async(findcirc_gmap.getCircInformation, (datafram,ouput_file,))    
        pool.close()
        pool.join()
    elif map_method == "minimap2":
        os.system("minimap2 -x splice -G 50 -N 10 %s  %s -Y -o  %s_minimap " %
          (genome_filename, query_filename, output_file))
        df = findcirc_minimap.exactReads("%s_minimap" % output_file)
        print(df.shape)
        df_list = findcirc_minimap.cutDatafram(df)
        pool = multiprocessing.Pool(len(df_list))
        for datafram in df_list:
            print(datafram.shape)
            pool.apply_async(findcirc_minimap.getCircInformation, (datafram,ouput_file,))    
        pool.close()
        pool.join()
if __name__ == '__main__':
    map_methods = ["blast", "blat","gmap","minimap2"]
    pool = multiprocessing.Pool(len(map_methods))
    for map_method in map_methods:
        pool.apply_async(map, (map_method,))
    pool.close()
    pool.join()
'''
os.system("cat %s_*_* > %s_all" %(output_file,output_file))
infile = '%s_all' % (output_file)
outfile = '%s_circRNA' % (output_file)
df = pd.read_csv(infile, sep=',', header=None)
print(df.shape)
df = df.drop_duplicates(subset=[1 ,2, 3, 4], keep='first')
df.to_csv(outfile, header=None, index=None)
'''




