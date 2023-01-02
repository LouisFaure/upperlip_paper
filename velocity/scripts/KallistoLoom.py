print("Loading libraries")
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as scp
import anndata
import scvelo as scv
import sys
import mygene


pl=sys.argv[1]
print(pl)
print("reading spliced matrix")
s=sc.read_mtx(pl+"/spliced.mtx")
print("reading unspliced matrix")
u=sc.read_mtx(pl+"/unspliced.mtx")
print("reading barcodes and genes")
s_bcs=pd.read_csv(pl+"/spliced.barcodes.txt",header=None)+"-1"
u_bcs=pd.read_csv(pl+"/unspliced.barcodes.txt",header=None)+"-1"
s_genes=pd.read_csv(pl+"/spliced.genes.txt",header=None)
u_genes=pd.read_csv(pl+"/unspliced.genes.txt",header=None)
s.obs.index = s_bcs[0].values
u.obs.index = u_bcs[0].values
s.var.index = s_genes[0].values
u.var.index = u_genes[0].values
genes = s_genes
genes.columns=["gid"]
s_bcs.columns = ["bcs"]
u_bcs.columns = ["bcs"]


t2g=pd.read_table("/home/lfaure/tools/Kallisto/transcripts_to_genes.txt",header=None)
del t2g[0]
t2g=t2g.loc[~t2g[1].duplicated(),:]
t2g.index=t2g[1]

t2g=t2g.loc[genes["gid"],:]
t2g.columns=["EnsemblID","symbol"]

print("initialising spliced and unspliced matrices")
sadata = anndata.AnnData(X=s.X, obs=s_bcs, var=t2g)
uadata = anndata.AnnData(X=u.X, obs=u_bcs, var=t2g)
print("get filtered barcodes")

pathtobcs=sys.argv[2]
bcs=pd.read_csv(pathtobcs,header=None)
sadata.obs.index=s_bcs["bcs"]
uadata.obs.index=u_bcs["bcs"]
sadata = sadata[list(set(bcs[0]) & set(s_bcs["bcs"])),:]
uadata = uadata[list(set(bcs[0]) & set(s_bcs["bcs"])),:]
adata = anndata.AnnData(X=sadata.X+uadata.X,obs=bcs[0], var=t2g)
adata.layers["spliced"] = sadata.X
adata.layers["unspliced"] = uadata.X
#adata.layers["ambiguous"] = scp.sparse.csr_matrix(np.zeros(adata.X.shape))
adata.obs = sadata.obs
adata.obs["CellID"] = adata.obs.index
adata.var.index=adata.var["symbol"]
print("writing loom")
adata.write_loom(pl+"/velocyted.loom")