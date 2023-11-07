import numpy as np
import pandas as pd
import pyarrow
from pandas.core.dtypes.common import is_numeric_dtype
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


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
    engine = 'pyarrow' if pd.__version__>='1.4.0' else 'c'
    return pd.read_csv(file_path,
                       delimiter="\t",
                       usecols=cols,
                       engine=engine)


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
                plt.xlabel('Density')
                plt.xlabel('log2 intensity')
                plt.title('Histogram of log2 intensities')
                plt.annotate('Cutoff',
                             xy =(log2_intensity_cutoff, 0),
                             xytext =(log2_intensity_cutoff, max(n)/2), 
                             arrowprops = dict(width = 0.001,
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
                m.append(np.nanmedian(dl[index]))
            print("Barplotting raw data ...")
            fig2 = plt.figure(figsize=(pdf_width, pdf_height))
            y_pos = np.arange(1, len(samples) + 1)
            plt.boxplot(dl,
                        flierprops=dict(marker='o', 
                                        markerfacecolor='blue', 
                                        markersize=1, 
                                        markeredgecolor='none'))
            plt.xticks(y_pos, samples, rotation=90)
            figs.append(fig2)

        if median_normalization is True:
            print("Median normalization ...")
            for sample, f in zip(samples, np.nanmean(m) - m):
                df.loc[df['sample_list'] == sample, 'quant'] += f
            if pdf_out is not None:
                dl = []
                m = []
                for index, sample in enumerate(samples):
                    dl.append(df.loc[df['sample_list'] == sample, 'quant'])
                    m.append(np.nanmedian(dl[index]))
                fig3 = plt.figure(figsize=(pdf_width, pdf_height))
                y_pos = np.arange(1, len(samples) + 1)
                plt.boxplot(dl,
                            flierprops=dict(marker='o',
                                            markerfacecolor='green',
                                            markersize=1,
                                            markeredgecolor='none'))
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
        print(" Complete!!")
        return p_list
    else:
        raise TypeError("preprocessed_data isn't pd.Dataframe")


def maxLFQ(X):
    X = np.array(X, dtype=np.float64)
    
    # kiểm tra nan toàn bộ dữ liệu
    if np.all(np.isnan(X)):
        return dict({"estimate": None, "annotation": "NA"})
    # kiểm tra số hàng
    if X.shape[0] == 1:
        return dict({"estimate": X[1, ], "annotation": ""})
    
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

          
def create_protein_table(protein_list, method="maxLFQ"):
    
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
            threes_display+=step
            filled+=1
    # progress display 
    
        if method == "maxLFQ":
            out = maxLFQ(list(protein_list.values())[i])
        # elif method == "medpolish":
        #     out = median_polish(protein_list[[i]], ...)
        # elif method == "topN":
        #     out = topN(protein_list[[i]], ...)
        # elif method == "meanInt":
        #     out = meanInt(protein_list[[i]])
        else:
            raise Exception("Sorry, Unknown method: ", method)

        tab.iloc[i, :] = out['estimate']
        annotation[i] = out['annotation']

    print(" Complete!!")
    return dict({"estimate": tab, "annotation": annotation})