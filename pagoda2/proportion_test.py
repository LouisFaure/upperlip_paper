import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import matplotlib.patches as mpatches

def proportion_test(adata,
                sample_identity,
                cluster_identity,
                sample_1,
                sample_2,
                n_permutations=1000,
                log2FD_cutoff=round(np.log2(1.5),2),
                n_jobs=1):
    
    meta = adata.obs[[sample_identity,cluster_identity]]

    meta = meta[meta[sample_identity].isin([sample_1,sample_2])]
    meta[sample_identity]=meta[sample_identity].astype(str)
    cluster_cases=meta[cluster_identity].unique()

    obs_diff=meta.groupby([sample_identity,cluster_identity]).size()
    obs_diff=obs_diff/obs_diff.groupby(sample_identity).sum()
    obs_diff=obs_diff.unstack().T
    obs_diff.columns=obs_diff.columns.astype(str)
    obs_diff["obs_log2FD"]=(np.log2(obs_diff[sample_2])-np.log2(obs_diff[sample_1])).values

    def permut(meta):
        perm = meta.copy()
        perm.loc[:,sample_identity]=np.random.permutation(perm[sample_identity])
        obs_diff=perm.groupby([sample_identity,cluster_identity]).size()
        obs_diff=obs_diff/obs_diff.groupby(sample_identity).sum()
        obs_diff=obs_diff.unstack().T
        obs_diff.columns=obs_diff.columns.astype(str)
        return (np.log2(obs_diff[sample_2])-np.log2(obs_diff[sample_1])).values

    res=Parallel(n_jobs=n_jobs)(delayed(permut)(meta) for i in tqdm(range(n_permutations)))

    perm_results=pd.DataFrame(np.vstack(res).T,index=obs_diff.index)

    increased = (perm_results-obs_diff.obs_log2FD.values.reshape(-1,1)>=0).sum(axis=1)
    increased = (increased + 1) / (n_permutations + 1)

    decreased=(perm_results-obs_diff.obs_log2FD.values.reshape(-1,1)<=0).sum(axis=1)
    decreased = (decreased + 1) / (n_permutations + 1)


    obs_diff["pval"] = list(map(lambda i: increased[i] if obs_diff.obs_log2FD.values[i]>0 else decreased[i],range(obs_diff.shape[0])))
    obs_diff.loc[obs_diff.obs_log2FD==np.inf,"pval"]=np.nan
    obs_diff.loc[np.isnan(obs_diff.obs_log2FD),"pval"]=np.nan

    obs_diff.loc[~np.isnan(obs_diff.pval),"fdr"]=\
        fdrcorrection(obs_diff.loc[~np.isnan(obs_diff.pval),"pval"])[1]

    def bootstrap(meta):
        boot = meta.copy()
        boot[cluster_identity]= boot.groupby(sample_identity)[cluster_identity].sample(frac=1,replace=True).values
        obs_diff=boot.groupby([sample_identity,cluster_identity]).size()
        obs_diff=obs_diff/obs_diff.groupby(sample_identity).sum()
        obs_diff=obs_diff.unstack().T
        obs_diff.columns=obs_diff.columns.astype(str)
        return (np.log2(obs_diff[sample_2])-np.log2(obs_diff[sample_1])).values

    res_boot=Parallel(n_jobs=n_jobs)(delayed(bootstrap)(meta) for i in tqdm(range(n_permutations)))

    boot_results=pd.DataFrame(np.vstack(res_boot).T,index=obs_diff.index)

    boot_results=pd.concat([boot_results.mean(axis=1),
               boot_results.quantile(q=[0.025,0.975],axis=1).T],axis=1)
    boot_results.columns=["boot_mean_log2FD","boot_CI_2.5","boot_CI_97.5"]
    
    prop_res = pd.concat([obs_diff,boot_results],axis=1)
    
    prop_res["significance"]="n.s."
    prop_res.loc[(prop_res.obs_log2FD.abs()>log2FD_cutoff) & (prop_res.fdr<0.05),
              "significance"]="FDR < 0.05 & abs(Log2FD) > "+\
              str(log2FD_cutoff)
    
    prop_res.index.name="clusters"
    return prop_res


import matplotlib.pyplot as plt

def plot_proportion_test(prop_res,sort=True):
    prop_res = prop_res.copy()
    
    exclu=(prop_res.iloc[:,:2]==0)
    exclu_cl=exclu[exclu.apply(lambda x: x.sum(),axis=1)==1].idxmin(axis=1)
    exclu_all=exclu[exclu.apply(lambda x: x.sum(),axis=1)==2].index
    txt=list(map(lambda i: "cells from cluster "+str(exclu_cl.index[i])+" only present in "+str(exclu_cl.values[i]),range(len(exclu_cl))))
    txt2=list(map(lambda i: "no cells from cluster "+exclu_all[i],range(len(exclu_all))))
    toadd = "\n".join(txt+txt2)
    
    prop_res.reset_index(inplace=True)
    
    prop_res = prop_res.loc[np.isfinite(prop_res[["pval",
                                                   "fdr",
                                                   "boot_mean_log2FD",
                                                   "boot_CI_2.5",
                                                   "boot_CI_97.5"]]).apply(all,axis=1)]
    if sort:
        prop_res=prop_res.sort_values("obs_log2FD",ascending=False)

    fig, ax = plt.subplots()

    ax.scatter(prop_res.obs_log2FD,
               prop_res.clusters,c="grey")

    ax.scatter(prop_res.obs_log2FD[prop_res.significance=="n.s."],
               prop_res.clusters[prop_res.significance=="n.s."],
               label="n.s.",c="darkgrey")

    for c in prop_res.clusters[prop_res.significance=="n.s."].index:
        ax.plot([prop_res.loc[c,"boot_CI_2.5"],prop_res.loc[c,"boot_CI_97.5"]],
                [prop_res.loc[c,"clusters"],prop_res.loc[c,"clusters"]],c="darkgrey")
        
    if (prop_res.significance!="n.s.").sum()>0:
        signi=prop_res.significance[prop_res.significance!="n.s."].unique()[0]

        ax.scatter(prop_res.obs_log2FD[prop_res.significance!="n.s."],
                   prop_res.clusters[prop_res.significance!="n.s."],
                   label=signi,c="tomato")

        for c in prop_res.clusters[prop_res.significance!="n.s."].index:
            ax.plot([prop_res.loc[c,"boot_CI_2.5"],prop_res.loc[c,"boot_CI_97.5"]],
                    [prop_res.loc[c,"clusters"],prop_res.loc[c,"clusters"]],c="tomato")
        ax.axvline(-np.float32(signi.split("> ")[1]),c="k",linestyle="--")
        ax.axvline(np.float32(signi.split("> ")[1]),c="k",linestyle="--")

    ax.axvline(0,c="k")

    
    
    ax.set_xlabel("obs_log2FD")
    ax.set_ylabel("clusters")
    
    

    # where some data has already been plotted to ax
    handles, labels = ax.get_legend_handles_labels()

    # manually define a new patch 
    patch = mpatches.Patch(color='grey', label=toadd,alpha=0)

    # handles is a list, so append manual patch
    handles.append(patch) 
    ax.legend(handles=handles,bbox_to_anchor=(1.05, 1), loc='upper left')
    
    samples=exclu.columns
    
    ax.set_title(samples[0]+r"$\leftrightarrow$"+samples[1])
    