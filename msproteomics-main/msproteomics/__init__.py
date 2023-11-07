# coding=utf-8
# Copyright 2023 Thang V Pham
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# With contributions from Truong Xuan Nam and others at Thuy Loi University.

__version__ = '0.0.2'
__author__ = 'Thang V Pham'

import sys
import numpy as np
import pandas as pd
import pyarrow
from pandas.core.dtypes.common import is_numeric_dtype
from scipy.linalg import lstsq
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm

import re

def read(file_path,
         primary_id="PG.ProteinGroups",
         secondary_id=np.array(["EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"]),
         sample_id="R.Condition",
         intensity_col="F.PeakArea",
         other_col=np.array(["F.ExcludedFromQuantification", "PG.Qvalue", "EG.Qvalue"])):
    cols = np.concatenate(([primary_id],
                              secondary_id,
                              [sample_id],
                              [intensity_col],
                              other_col),
                             axis=0)
    return pd.read_csv(file_path,
                       delimiter="\t",
                       #usecols=cols,
                       engine="pyarrow")


def preprocess(quant_table,
               primary_id="PG.ProteinGroups",
               secondary_id=np.array(["EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"]),
               sample_id="R.Condition",
               intensity_col="F.PeakArea",
               median_normalization=True,
               log2_intensity_cutoff=0,
               pdf_out="qc-plots.pdf",
               pdf_width=12,
               pdf_height=8):
    if isinstance(quant_table, pd.DataFrame):
        if not is_numeric_dtype(quant_table[intensity_col]):
            raise TypeError("Intensity column must be numeric")

        print("Concatenating secondary ids...")
        second_id = quant_table[secondary_id[0]]
        for col in range(1, len(secondary_id)):
            second_id += quant_table[secondary_id[col]].astype(str)
        df = pd.DataFrame({'protein_list': quant_table[primary_id],
                           'sample_list': quant_table[sample_id],
                           'quant': np.log2(quant_table[intensity_col]),
                           'id': second_id})
        df.dropna(axis=0)
        # intensity cut off
        figs = []
        if log2_intensity_cutoff is not None:
            print('Removing low intensities...')
            if pdf_out is not None:
                # histogram
                fig1 = plt.figure(figsize=(pdf_width, pdf_height))
                n, bins, patch = plt.hist(df['quant'], bins=50, density=True)
                plt.ylabel('Density')
                plt.xlabel('log2 intensity')
                plt.title('Histogram of log2 intensities')
                #plt.arrow(log2_intensity_cutoff, 0, log2_intensity_cutoff, max(n)/2, color='r')
                plt.annotate('Cutoff',
                              xy =(log2_intensity_cutoff, 0),
                              xytext =(log2_intensity_cutoff, max(n)/2), 
                              arrowprops = dict(width = 0.001,
                              color = 'r',
                              facecolor ='red',
                              shrink = 0.0),)
                figs.append(fig1)
            df = df[df["quant"] > log2_intensity_cutoff]
        samples = df['sample_list'].unique()

        if pdf_out is not None:
            dl = []
            m = []
            for index, sample in enumerate(samples):
                dl.append(df.loc[df['sample_list'] == sample, 'quant'])
                m.append(np.median(dl[index]))
            print("Barplotting raw data ...")
            fig2 = plt.figure(figsize=(pdf_width, pdf_height))
            y_pos = np.arange(1, len(samples) + 1)
            plt.boxplot(dl,
                        flierprops=dict(marker='o', markerfacecolor='blue', markersize=1, markeredgecolor='none'))
            plt.xticks(y_pos, samples, rotation=90)
            figs.append(fig2)

        if median_normalization is True:
            print("Median normalization ...")
            for sample, f in zip(samples, np.mean(m) - m):
                df.loc[df['sample_list'] == sample, 'quant'] += f
            if pdf_out is not None:
                dl = []
                m = []
                for index, sample in enumerate(samples):
                    dl.append(df.loc[df['sample_list'] == sample, 'quant'])
                    m.append(np.median(dl[index]))
                fig3 = plt.figure(figsize=(pdf_width, pdf_height))
                y_pos = np.arange(1, len(samples) + 1)
                plt.boxplot(dl,
                            flierprops=dict(marker='o', markerfacecolor='green', markersize=1, markeredgecolor='none'))
                plt.xticks(y_pos, samples, rotation=90)
                figs.append(fig3)

        with PdfPages(pdf_out) as pdf:
            for fig in figs:
                plt.figure(fig)
                pdf.savefig()
        return df
    else:
        raise TypeError("quant_table isn't pd.Dataframe")


def create_protein_list(preprocessed_data):
    if isinstance(preprocessed_data, pd.DataFrame):
        if any(pd.isna(preprocessed_data['protein_list'])):
            raise Exception("NA value in protein_list")
            
        if any(pd.isna(preprocessed_data['sample_list'])):
            raise Exception("NA value in sample_list")
            
        if any(pd.isna(preprocessed_data['id'])):
            raise Exception("NA value in id")
            
        if any(pd.isna(preprocessed_data['quant'])):
            raise Exception("NA value in quant") 
            
        proteins = preprocessed_data['protein_list'].unique()
        samples = preprocessed_data['sample_list'].unique()
        print("Create protein list..")
        print("#proteins = {0}, #samples = {1}".format(proteins.shape[0], samples.shape[0]))
        
        p_list = {}
        # progress display
        threes_display = 0
        filled = 0
        step = proteins.shape[0] / 20
        for i in range(proteins.shape[0]):
            if i >= threes_display-1:
                print('\r[{:}] {:.0%}'.format('#'*filled + ' '*(20-filled), i/proteins.shape[0]), end = '')
                threes_display+=step
                filled+=1
        # progress display
            tmp = preprocessed_data[preprocessed_data["protein_list"] == proteins[i]]
            if tmp.shape[0] > 0:
                dupl = tmp[['id', 'sample_list']].duplicated()
                if dupl.any():
                    for _sample, _id in tmp.loc[dupl, ['id', 'sample_list']]:
                        print("sample {0}; id {1} not unique.".format(_sample, _id))
                    raise Exception("Duplicate entry")
                else:
                    m = pd.DataFrame(columns=samples, index=tmp['id'].unique())
                    for j in tmp.index:
                        m.loc[tmp.at[j, 'id'], tmp.at[j, 'sample_list']] = tmp.at[j, 'quant']
                    p_list[proteins[i]] = m
        print("Complete!!")
        return p_list
    else:
        raise TypeError("preprocessed_data isn't pd.Dataframe")



def maxLFQ(X):
    X = np.array(X, dtype=np.float64)
    
    # check value None in data
    if np.all(np.isnan(X)):
        return dict({"estimate": None, "annotation": "NA"})
    # check number row
    if X.shape[0] == 1:
        return dict({"estimate": X[0, ], "annotation": ""})
    
    N = X.shape[1]  # [row, col]
    cc = 0
    g = np.full(N, np.nan)
    
    def spread(i):
        g[i] = cc
        for r in range(X.shape[0]):
            if not np.isnan(X[r, i]):
                for k in range(X.shape[1]):
                    if (not np.isnan(X[r, k])) and np.isnan(g[k]):
                        spread(k)
    # MaxLFQ
    def maxLFQdo(X):
        X = np.array(X)
        Ncol = X.shape[1]
        AtA = np.zeros((Ncol, Ncol))
        Atb = np.zeros(Ncol)
        for i in range(Ncol - 1):
            for j in range(i + 1, Ncol):
                r_i_j = np.nanmedian(np.array(-X[:, i] + X[:, j]))
                if not np.isnan(r_i_j):
                    AtA[i, j] = AtA[j, i] = -1
                    AtA[i, i] = AtA[i, i] + 1
                    AtA[j, j] = AtA[j, j] + 1

                    Atb[i] = Atb[i] - r_i_j
                    Atb[j] = Atb[j] + r_i_j

        l = np.append(np.append(2*AtA, np.ones((Ncol, 1)), axis=1),
                      np.append(np.ones(Ncol), 0).reshape(1, Ncol+1),
                      axis=0)
        r = np.append(2*Atb,
                      [np.nanmean(X)*Ncol],
                      axis=0).reshape((Ncol+1, 1))
        x = np.linalg.solve(l, r)
        return x.flatten()[:Ncol]
       
    for i in range(N):
        if np.isnan(g[i]):
            cc += 1
            spread(i)
    
    w = np.full(N, np.nan)
    for i in range(cc):
        ind = np.array(g == i + 1)
        if sum(ind) == 1:
            w[ind] = np.nanmedian(np.array(X[:,ind]))
        else:
            w[ind] = maxLFQdo(X[:,ind])
    
    if np.all(np.isnan(w)):
        return dict({"estimate": w, "annotation": "NA"})
    else:
        if np.all(g[~np.isnan(w)]):
            return dict({"estimate": w, "annotation": ""})
        else:
            return dict({"estimate": w, "annotation": ";".join(np.put(g, np.nan, "NA"))})

          
def create_protein_table(protein_list, method = "maxLFQ"):
    if not isinstance(protein_list, dict):
        raise TypeError("Only dict are allowed")

    if len(protein_list) == 0:
        return None

    tab = pd.DataFrame(None, columns=list(protein_list.values())[0].columns, index=list(protein_list))
    annotation = pd.Series(np.full(len(protein_list), np.nan))
    

    # progress display
    nrow = tab.shape[0]
    threes_display = 0
    filled = 0
    step = nrow / 20
    for i in range(nrow):
        if i >= threes_display-1:
            print('\r[{:}] {:.0%}'.format('#'*filled + ' '*(20-filled), i/nrow), end = '')
            threes_display += step
            filled += 1
    
        if method == "maxLFQ":
            out = maxLFQ(list(protein_list.values())[i])
    
        else:
            raise Exception("Unknown method: ", method)

        tab.iloc[i, :] = out['estimate']
        annotation[i] = out['annotation']

    print(" Complete!!")
    return dict({"estimate": tab, "annotation": annotation})

def plot_protein(X, 
                 main = "", 
                 col = [], 
                 split = 0.6):
    #Input X is a matrix
    if col == []:
        for i in range(1, X.shape[0]):
            col.append(i)

    # if split is not None:
    #     old_par =    

    plt.close('all')
    X = X.to_numpy()
    a = []
    j = 0
    x = np.arange(len(X))
    ys = [i+x+(i*x)**2 for i in range(len(X))]
    colors = cm.rainbow(np.linspace(0, 1, len(ys)))
    for y, c in zip(ys, colors):
        a.append(c)
    fig4 = plt.figure()
    plt.ylabel('Intensity')
    plt.xlabel('Sample')
    plt.title(main)
    plt.xticks(np.linspace(1,24,24))
    plt.yticks(np.linspace(0,int(np.nanmax(X)),9))
    for i in range(len(X)):
        plt.plot(np.linspace(1,24,24),X[i],color = a[j],marker = 'o')
        j += 1
    plt.show() 

def extract_annotation(protein_ids, 
                       quant_table,
                       primary_id = "PG.ProteinGroups",
                       annotation_columns = None):
    
    #concatenate input columns
    all_columns = np.concatenate(([primary_id], annotation_columns), axis = 0) 

    #get column name of input data 
    colnames, index = [], []
    for a, b in quant_table.items():
        colnames.append(a)

    #check all_columns in data ?
    for i in range(len(all_columns)):
        if(all_columns[i] in colnames):
            index.append(colnames.index(all_columns[i]))
        else:
            index.append(None)
    
    #if exist value NaN then break
    if -1 in index:
        raise Exception("The input table has no column")

    for i in range(len(protein_ids)):
        if(protein_ids[i] in quant_table[primary_id]):
            index.append(quant_table[primary_id].index(protein_ids[i]))
        else:
            index.append(None)

    if -1 in index:
        raise Exception("Cannot find")

    #get rows index of column in all_columns    
    ind = [quant_table[primary_id].index(x) for x in protein_ids]
    tab = quant_table[np.concatenate(([primary_id], annotation_columns), axis = 0)].iloc[ind]
    
    return(tab)


def create_site_key(protein_ids, ptm_locations):
    n = len(protein_ids)
    
    all_sites = [''] * n
    
    secondary_keys = [''] * n
    
    first_sites = [''] * n
    
    for i in range(n):
        
        pp = protein_ids[i].split(';')

        ss = ptm_locations[i].split(';')
        
        if len(pp) > 0 :
            a = [''] * len(pp)
            for j in range(len(pp)):
                p = pp[j]
                s = ss[j]
                
                s2 = re.findall( '\(.*?\)', s) 
                
                b = [''] * len(s2) 
                
                for k in range(len(s2)):
                    r = re.findall( '[STY][0-9]+', s2[k]) 
                    if len(r) > 0:
                        tmp = [p + '_' + rr + '_M' + str(min(len(r), 3)) for rr in r]
                        b[k] = ';'.join(tmp)                                    
                a[j] = ';'.join(b)
                        
            
            first = True
            for j in range(len(pp)):
                if pp[j] != '' and a[j] != '':
                    if first:
                        all_sites[i] = a[j]
                        first_sites[i] = a[j]
                        first = False
                    else:
                        all_sites[i] = all_sites[i] + ';' + a[j]

        
        if i % 100000 == 0:
            print(int(i * 100.0 / n), '%\n')        
    
    unique_sites = {y for x in all_sites for y in x.split(';')}
    
    unique_first_sites = {y for x in first_sites for y in x.split(';')}

    return(all_sites, unique_sites, first_sites, unique_first_sites)

def create_site_report_longformat(tab, keys):
    
    ids = [y for x in keys for y in x.split(';')]
    
    tmp2 = [len(x.split(';')) for x in keys]
    r = [i for i in range(len(tmp2)) for j in range(tmp2[i])]
        
    tab2 = tab.reset_index(drop = True)
    return(tab2.iloc[r].assign(site = ids))

def create_site_report_wideformat(report_lf,
                                  sample_id = "R.FileName",
                                  intensity_col = "log2_intensity",
                                  primary_id = 'site',
                                  secondary_id = ["EG.PrecursorId", "EG.Library", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType"],
                                  annotation_cols = ["PG.Organisms"],
                                  method = "maxLFQ",
                                  check_uniqueness = False):
    
    print("Concatenating secondary ids...")
    second_id = report_lf[secondary_id[0]]
    for col in range(1, len(secondary_id)):
        second_id += ':' + report_lf[secondary_id[col]].astype(str)

    aa = pd.DataFrame({'protein_list': report_lf[primary_id].copy(),
                       'sample_list': report_lf[sample_id].copy(),
                       'quant': report_lf[intensity_col].copy(),
                       'id': second_id.copy()})
    aa.dropna(axis=0)
    

    if method == 'maxLFQ':
        #res = fast_MaxLFQ(aa)
    
        p_list = create_protein_list(aa)        
        result = create_protein_table(p_list, method = 'maxLFQ')
        
        # check difference    

    elif method == 'sum':        
        p_list = create_protein_list(aa)
        result = {'estimate': pd.DataFrame(np.nan, columns=list(p_list.values())[0].columns, index = list(p_list))}
        for i in p_list:
            tmp = np.exp2(p_list[i].to_numpy().astype(float))
            tmp2 = np.nansum(tmp, axis = 0)
            result['estimate'].loc[i] = np.log2(tmp2, out = np.full_like(tmp2, np.nan), where = tmp2 > 0)

    else:
        raise Exception('Unknown method.')
         
    a = list(report_lf[primary_id])
    ind = [a.index(x) for x in list(result['estimate'].index)]
    tmp = report_lf[[primary_id] + ['all_sites'] + annotation_cols].iloc[ind]    
    ret = tmp.set_axis(result['estimate'].index).join(result['estimate'])

    return ret


def normalize(report_lf, sample_id, intensity_col):
    
    s = set(report_lf[sample_id])
    print(len(s) , " samples:")
    for si in s:
        print(si)

    int_log2 = np.log2(report_lf[intensity_col],  out = np.zeros_like(report_lf[intensity_col]), where = (report_lf[intensity_col] != 0))
    print(sum(np.isnan(int_log2)), " NA(s).")
    
    m = {si: np.median(int_log2[report_lf[sample_id] == si]) for si in s}
        
    sm = 0
    for i in m:
        sm += m[i]
    sm /= len(m)  

    for i in m:
        m[i] = sm - m[i]
        
    for i in m:
        idx = report_lf[sample_id] == i
        int_log2[idx] = int_log2[idx] + m[i]

    return int_log2

