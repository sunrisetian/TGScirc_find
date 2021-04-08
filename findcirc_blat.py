import pandas as pd
import numpy as np
import multiprocessing

#use gmap
def exactReads(filename):
    df = pd.read_csv(filename, sep='\t', skiprows =5 ,header=None)
    df = df[df.duplicated(subset=9, keep=False)]  
    df = df.drop_duplicates(subset=[9, 13, 0, 11, 12], keep='first') 
    return df


def cutDatafram(df):  # cut df based on chromnome
    set_chr, df_list = set(df[13].values), []  # get all chromnome
    for chr_num in set_chr:
        index_list = df[df[13] == chr_num].index.tolist()
        df_cut = df.loc[index_list, :17]
        df_list.append(df_cut)
    return df_list


def getCircInformation(df,output_path):
    set_name = set(df[9].values)
    circ_start_end = []
    for name in set_name:
        index_list = df[df[9] == name].index.tolist()
        little_df = df.loc[index_list, :17]
        little_df.sort_values(11, inplace=True)
        for i in range(len(index_list) - 1):
            res_list = []
            if little_df.iloc[i, 13] == little_df.iloc[i + 1, 13]:
                chr_name, querry_name = little_df.iloc[i, 13], name
                spacing = little_df.iloc[i + 1, 11] - little_df.iloc[i, 12]
                reads_po = sorted([little_df.iloc[i, 11], little_df.iloc[i, 12], little_df.iloc[i + 1, 11],
                                   little_df.iloc[i + 1, 12]])
                if (little_df.iloc[i, 8] == '+') and (little_df.iloc[i + 1, 8] == '+') and \
                        (-20 < little_df.iloc[i, 15] - little_df.iloc[i + 1, 16] < 1000000)   and (spacing in range(-5, 10)):
                    res_list.append(querry_name)
                    res_list.append(chr_name)
                    res_list.append(little_df.iloc[i, 8])
                    res_list.append(little_df.iloc[i+1, 15])
                    res_list.append(little_df.iloc[i, 16])
                    circ_start_end.append(res_list)
                    res_list = []
                elif (little_df.iloc[i, 8] == '-') and (little_df.iloc[i + 1, 8] == '-') and \
                        (-20 < little_df.iloc[i + 1, 15] - little_df.iloc[i, 16] < 1000000)  and (spacing in range(-5, 10)):
                    res_list.append(querry_name)
                    res_list.append(chr_name)
                    res_list.append(little_df.iloc[i, 8])
                    res_list.append(little_df.iloc[i, 15])
                    res_list.append(little_df.iloc[i+1, 16])
                    circ_start_end.append(res_list)
                    res_list = []
            else:
                pass
    res_df = pd.DataFrame(circ_start_end)
    filepath = '%s_%s_blat' % (output_path,chr_name)
    res_df.to_csv(filepath, header=None, index=None)
    return res_df


