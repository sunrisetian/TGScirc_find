import pandas as pd
import numpy as np
import multiprocessing

#use gmap
def exactReads(filename):
    df = pd.read_csv(filename, sep='\t', header=None,usecols=[0,1,2,3,4,5,6,7,8,9,10,11])
    df = df[df.duplicated(subset=0, keep=False)]  # 去除只有单一比对结果的reads
    df = df.drop_duplicates(subset=[0, 1, 2, 3, 5], keep='first')  # 去除完全重复的全长比对reads
    return df


def cutDatafram(df):  # cut df based on chromnome
    set_chr, df_list = set(df[5].values), []  # get all chromnome
    for chr_num in set_chr:
        index_list = df[df[5] == chr_num].index.tolist()
        df_cut = df.loc[index_list, :11]
        df_list.append(df_cut)
    return df_list



def getCircInformation(df,output_path):
    set_name = set(df[0].values)
    circ_start_end = []
    for name in set_name:
        index_list = df[df[0] == name].index.tolist()
        little_df = df.loc[index_list, :11]
        little_df.sort_values(2, inplace=True)
        for i in range(len(index_list) - 1):
            res_list = []
            if little_df.iloc[i, 5] == little_df.iloc[i + 1, 5]:
                chr_name, querry_name = little_df.iloc[i, 5], name
                spacing = little_df.iloc[i + 1, 2] - little_df.iloc[i, 3]
                reads_po = sorted([little_df.iloc[i, 2], little_df.iloc[i, 3], little_df.iloc[i + 1, 2],
                                   little_df.iloc[i + 1, 3]])
                if (little_df.iloc[i, 4] == '+') and (little_df.iloc[i + 1, 4] == '+') and \
                        (-20 < little_df.iloc[i, 7] - little_df.iloc[i+ 1 , 8] < 1000000) and (spacing in range(-5, 10)):
                    res_list.append(querry_name)
                    res_list.append(chr_name)
                    res_list.append(little_df.iloc[i, 4])
                    res_list.append(little_df.iloc[i+1, 7])
                    res_list.append(little_df.iloc[i, 8])
                    circ_start_end.append(res_list)
                    res_list = []
                elif (little_df.iloc[i, 4] == '-') and (little_df.iloc[i + 1, 4] == '-') and \
                        (-20 < little_df.iloc[i + 1, 7] - little_df.iloc[i, 8] < 1000000) and (spacing in range(-5, 10)):
                    res_list.append(querry_name)
                    res_list.append(chr_name)
                    res_list.append(little_df.iloc[i, 4])
                    res_list.append(little_df.iloc[i, 7])
                    res_list.append(little_df.iloc[i+1, 8])
                    circ_start_end.append(res_list)
                    res_list = []
            else:
                pass
    res_df = pd.DataFrame(circ_start_end)
    filepath = '%s_%s_minimap '% (output_path,chr_name)
    res_df.to_csv(filepath, header=None, index=None)
    return res_df

'''
if __name__ == '__main__':
    filename = ''
    df = exactReads(filename)
    print(df.shape)
    df_list = cutDatafram(df)
    pool = multiprocessing.Pool(len(df_list))
    for datafram in df_list:
        print(datafram.shape)
        pool.apply_async(getCircInformation, (datafram,))

    pool.close()
    pool.join()
'''
