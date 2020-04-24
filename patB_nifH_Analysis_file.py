#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import scipy as sp
import numpy as np
import multiprocessing as mp
from sklearn.feature_selection import mutual_info_regression
import seaborn as sns

# Stockel Dataset

df_stockel = pd.read_csv('MicroarrayData/StockelProcessed.csv')
duplicates = df_stockel.loc[df_stockel.duplicated(subset='ORF',keep=False)].sort_values(by=['ORF'])
mean_columns = list(duplicates.columns)
mean_columns.remove('Contig')
mean_columns.remove('ORF')
df_expression = df_stockel.groupby('ORF')[mean_columns].mean().reset_index()

# Genome Dataset

GenCyanoDB = pd.read_excel('GenCyanoDB.xlsx',index_col=0,usecols=[0,1,2,3])

# Merge the two DataFrames by their ORF column

df = df_expression.merge(GenCyanoDB,on='ORF',how='inner')

# Get the 2-fold change data

filename = 'Pakrasi PatN Nif H analysis file.xlsx'
filepath = '/Users/dzb5732/Box Sync/PhD-dzb5732@psu.edu/Cyanothece/Gene Expression Dataset/UWash-Raw/'
file = filepath + filename
df_2_fold = pd.read_excel(file,sheet_name='d4',usecols=[0,1,2,3,4,5])


df_2_fold_nifH = df_2_fold.iloc[:,0:2]
df_2_fold_nifH.columns = ['orf','fold_change']

df_2_fold_patB = df_2_fold.iloc[:,[0,2]]
df_2_fold_patB.columns = ['orf','fold_change']

# Get Mutual Information

class Interaction:
    def __init__(self,Exp_data,gene='all',mi_thresh=0):
        self.Exp_data = Exp_data
        if self.Exp_data.isnull().values.any():
            self.Exp_df = self.Exp_data.iloc[:,:-2].set_index('ORF').interpolate(method='linear',axis=1,limit_direction='both').T
        else:
            self.Exp_df = self.Exp_data.iloc[:,:-2].set_index('ORF').T
        if gene=='all':
            self.mi_dict = self._get_dict()
        else:
            self.gene_orf = gene
            self.mi_list = self._miscorelist(self.gene_orf)
            self.mi_thresh = mi_thresh
            self.df = self._get_df(self.mi_list,self.mi_thresh)
           
    
    def _get_dict(self):
        all_genes = list(self.Exp_df.columns)
        pool = mp.Pool(mp.cpu_count())
        results = pool.map(self._miscorelist,all_genes)
        fast_dict= dict(zip(all_genes,results))
        return fast_dict

    
    def _miscorelist(self,gene):
        all_other_genes_df = self.Exp_df.loc[:,self.Exp_df.columns!=gene]
        all_other_genes = np.array(all_other_genes_df.columns)
        this_gene_df = self.Exp_df[gene]
        mi_score = mutual_info_regression(all_other_genes_df,this_gene_df,discrete_features=False,random_state=7)
        miscore_genes = list(zip(all_other_genes,mi_score))
        sorted_miscore = sorted(miscore_genes,key = lambda x:x[1],reverse=True)
        return sorted_miscore
    
    def _get_df(self,mi_list,mi_thresh):
        my_dict = {'orf':[],'function':[],'CommonName':[],'mi':[]}
        for orf,mi in mi_list:
            if mi<=mi_thresh:
                break

            my_dict['orf'].append(orf)
            my_dict['function'].append(self.Exp_data.loc[self.Exp_data.ORF==orf].Function.values[0])
            my_dict['CommonName'].append(self.Exp_data.loc[self.Exp_data.ORF==orf].CommonName.values[0])
            my_dict['mi'].append(mi)

        return pd.DataFrame(my_dict)

## patB

patB = Interaction(df,'cce_1898',-1)

patB_info = df_2_fold_patB.merge(patB.df,how='outer',on='orf')

## nifH

nifH = Interaction(df,'cce_0559',-1)

nifH_info = df_2_fold_nifH.merge(nifH.df,how='outer',on='orf')

# Helper Functions

def top_fold_changed_genes(df,direction='negative',top=25):
    df.dropna(inplace=True)
    return df.sort_values(by='fold_change',ascending=True).iloc[0:top] if direction=='negative' else df.sort_values(by='fold_change',ascending=False).iloc[0:top]

def get_regulators(df,direction='negative'):
    df.dropna(inplace=True)
    return df.loc[abs(df.function.str.contains('regulator')) & (df.fold_change<-1)].sort_values(by='fold_change',ascending=True) if direction=='negative' else df.loc[abs(df.function.str.contains('regulator')) & (df.fold_change>1)].sort_values(by='fold_change',ascending=False) 

def get_sensors(df,direction='negative'):
    df.dropna(inplace=True)
    return df.loc[abs(df.function.str.contains('sensor')) & (df.fold_change<-1)].sort_values(by='fold_change',ascending=True) if direction=='negative' else df.loc[abs(df.function.str.contains('sensor')) & (df.fold_change>1)].sort_values(by='fold_change',ascending=False) 

def get_circadian(df,direction='negative'):
    df.dropna(inplace=True)
    return df.loc[abs(df.function.str.contains('circadian')) & (df.fold_change<-1)].sort_values(by='fold_change',ascending=True) if direction=='negative' else df.loc[abs(df.function.str.contains('circadian')) & (df.fold_change>1)].sort_values(by='fold_change',ascending=False) 

def get_sigma_factors(df,direction='negative'):
    df.dropna(inplace=True)
    return df.loc[abs(df.function.str.contains('sigma')) & (df.fold_change<-1)].sort_values(by='fold_change',ascending=True) if direction=='negative' else df.loc[abs(df.function.str.contains('sigma')) & (df.fold_change>1)].sort_values(by='fold_change',ascending=False) 

