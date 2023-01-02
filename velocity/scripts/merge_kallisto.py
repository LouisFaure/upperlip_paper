import glob
import scvelo as scv
import panda as pd

palette=["#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F","#8CD17D","#B6992D",
         "#F1CE63","#499894","#86BCB6","#E15759","#FF9D9A","#79706E","#BAB0AC"]

files = ["ML8/velo/velocyted.loom","ML9/velo/velocyted.loom",
         "ML10/velo/velocyted.loom","ML11/velo/velocyted.loom"]
adata=scv.read(files[0])
adata.obs.index=adata.obs.index.map(lambda x: files[0].split("/")[0]+"_"+x)
for f in files[1:]:
    bdata=scv.read(f)
    bdata.obs.index=bdata.obs.index.map(lambda x: f.split("/")[0]+"_"+x)
    adata=adata.concatenate(bdata,index_unique=None)
    
adata.write_h5ad("ML8-11_ALL.h5ad")

umap=pd.read_csv("../_Output/UMAP.csv",index_col=0) # general UMAP
pca=pd.read_csv("../_Output/PCA_ML8-11.csv",index_col=0)
clusters=pd.read_csv("../_Output/clusters_ML8-11.csv",index_col=0)

adata=adata[umap.index,:]
adata.obs["batch"] = adata.obs.index.map(lambda x: x.split("_")[0])
adata.obs["devtime"]=adata.obs["batch"]
adata.obs["devtime"][adata.obs["batch"].isin(["ML11"])]="E9.5"
adata.obs["devtime"][adata.obs["batch"].isin(["ML10"])]="E10.5"
adata.obs["devtime"][adata.obs["batch"].isin(["ML8","ML9"])]="E11.5"

adata.obsm["X_umap"]=umap.values
adata.obsm["X_pca"]=pca.values

adata.write_h5ad("ML8-11_F.h5ad")


files = ["ML6/velo/velocyted.loom","ML7/velo/velocyted.loom"]
adata=scv.read(files[0])
adata.obs.index=adata.obs.index.map(lambda x: files[0].split("/")[0]+"_"+x)
for f in files[1:]:
    bdata=scv.read(f)
    bdata.obs.index=bdata.obs.index.map(lambda x: f.split("/")[0]+"_"+x)
    adata=adata.concatenate(bdata,index_unique=None)
    
adata.write_h5ad("ML6-7_ALL.h5ad")

umap=pd.read_csv("../_Output/UMAP_ML6-7_epcam.csv",index_col=0) # general UMAP
clusters=pd.read_csv("../_Output/clusters_ML6-7_epcam.csv",index_col=0)

adata=adata[umap.index,:]
adata.obs["batch"] = adata.obs.index.map(lambda x: x.split("_")[0])
adata.obsm["X_umap"]=umap.values

adata.write_h5ad("ML6-7_F.h5ad")

adata.obs.leiden[adata.obs.leiden=="4"]="2"
adata.obs['leiden']=adata.obs['leiden'].astype('category')
palette=[palette[x] for x in (np.array(list(map(float,adata.obs["leiden"].cat.categories))).astype(int)-1).tolist()]

adata.write("adata_ML6-7_velo.h5ad")

files = ["ML6/velo/velocyted.loom","ML7/velo/velocyted.loom"]
adata=scv.read(files[0])
adata.obs.index=adata.obs.index.map(lambda x: files[0].split("/")[0]+"_"+x)
for f in files[1:]:
    bdata=scv.read(f)
    bdata.obs.index=bdata.obs.index.map(lambda x: f.split("/")[0]+"_"+x)
    adata=adata.concatenate(bdata,index_unique=None)
    
adata.write_h5ad("ML6-7_ALL.h5ad")

umap=pd.read_csv("ML6-7/UMAP.csv",index_col=0)
adata=adata[umap.index]
adata.obs["leiden"]=pd.read_csv("ML6-7/leiden.tsv",header=None)[0].values.astype(str)
adata.obsm["X_umap"]=umap.values

adata.obs.leiden[adata.obs.leiden=="4"]="2"
adata.obs['leiden']=adata.obs['leiden'].astype('category')
palette=[palette[x] for x in (np.array(list(map(float,adata.obs["leiden"].cat.categories))).astype(int)-1).tolist()]

adata.write("adata_ML6-7_velo.h5ad")

adata_ML4=scv.read("ML4/velo/counts_unfiltered/adata.h5ad")
adata_ML5=scv.read("ML5/velo/counts_unfiltered/adata.h5ad")
adata_ML4.obs.index=adata_ML4.obs.index.map(lambda x: "ML4"+"_"+x+"-1")
adata_ML5.obs.index=adata_ML5.obs.index.map(lambda x: "ML5"+"_"+x+"-1")
adata=adata_ML4.concatenate(adata_ML5,index_unique=None)
adata.obs["batch"] = adata.obs.index.map(lambda x: x.split("_")[0])
adata=adata[umap.index,:]
adata.obsm["X_umap"]=umap.values
adata.write_h5ad("ML4-5_F.h5ad")