import numpy as np
import pandas as pd
import scipy.sparse as sp
import cupyx.scipy.sparse as csp
import scanpy as sc
import cudf 
import gc
import nvtx


from tqdm import tqdm
from joblib import Parallel, delayed
class ProgressParallel(Parallel):
    def __init__(
        self, use_tqdm=True, total=None, file=None, desc=None, *args, **kwargs
    ):
        self._use_tqdm = use_tqdm
        self._total = total
        self._desc = desc
        self._file = file
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        with tqdm(
            disable=not self._use_tqdm,
            total=self._total,
            desc=self._desc,
            file=self._file,
        ) as self._pbar:
            return Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        if self._total is None:
            self._pbar.total = self.n_dispatched_tasks
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()

@nvtx.annotate(color="blue")
def read_mtx_gpu(filename,nchunks=None,transpose=False):
    print("Reading mtx")
    mtxinfo=pd.read_csv(filename,nrows=1,sep=" ",comment="%",header=None).values[0]
    shape = mtxinfo[[1,0]] if transpose else mtxinfo[[0,1]]
    shape=tuple((shape).astype(int))
    
    if nchunks is not None:
        print("    loading to device chunkwise")
        chunks=np.round(np.arange(0, mtxinfo[2]+mtxinfo[2]/nchunks, mtxinfo[2]/nchunks)).astype(int)
        toadata = np.zeros((mtxinfo[2],3),dtype=np.float32)
        for i in range(nchunks):
            with nvtx.annotate("cudf loop", color="green"):
                mtx_data=cudf.read_csv(filename,sep=" ",dtype=['float32' for i in range(3)],
                                       comment="%",header=None,nrows=chunks[i+1]-chunks[i],
                                       skiprows=chunks[i]+2)
                # offseting row and column indices to fit python indexing
                mtx_data["0"]=mtx_data["0"]-1
                mtx_data["1"]=mtx_data["1"]-1

                toadata[chunks[i]:chunks[i+1],:]=mtx_data.to_numpy()
            

        print("    generating csr sparse")
        if transpose:
            toadata=sp.csr_matrix((toadata[:,2], (toadata[:,1], toadata[:,0])),
                                shape=shape,dtype=np.float32)
        else:
            toadata=sp.csr_matrix((toadata[:,2], (toadata[:,0], toadata[:,1])),
                                shape=shape,dtype=np.float32)
    else:
        print("    loading to device")
        mtx_data=cudf.read_csv(filename,sep=" ",dtype=['float32' for i in range(3)],
                               comment="%",header=None,skiprows=2)
        # offseting row and column indices to fit python indexing
        mtx_data["0"]=mtx_data["0"]-1
        mtx_data["1"]=mtx_data["1"]-1
        
        mtx_data=mtx_data.to_cupy()
        
        print("    generating csr sparse")
        if transpose:
            mtx_data=csp.coo_matrix((mtx_data[:,2], (mtx_data[:,1], mtx_data[:,0])),
                            shape=shape,dtype=np.float32)
        else:
            mtx_data=csp.coo_matrix((mtx_data[:,2], (mtx_data[:,0], mtx_data[:,1])),
                            shape=shape,dtype=np.float32)
        toadata=mtx_data.get().tocsr()

    return sc.AnnData(toadata)

def read_mtx_parallel(filename,n_jobs=1,nchunks=None):
    print("Reading mtx")
    nchunks = n_jobs if nchunks is None else nchunks
    strchunk = " chunkwise" if nchunks>1 else ""
    print(f"    loading to host{strchunk}")
    mtxinfo=pd.read_csv(filename,nrows=1,sep=" ",comment="%",header=None).values[0]
    shape=tuple((mtxinfo[:2]).astype(int))
    
    chunks = np.round(np.arange(0, mtxinfo[2]+mtxinfo[2]/nchunks, mtxinfo[2]/nchunks)).astype(int)
    toadata = np.zeros((mtxinfo[2],3),dtype=np.float32)
    def load_chunk(i):
        nr=chunks[i+1]-chunks[i]
        sk=chunks[i]+2
        toadata[chunks[i]:chunks[i+1],:] = pd.read_csv(filename,sep=" ",comment="%",
                                                       header=None,nrows=nr,skiprows=sk,
                                                       dtype=np.float32).values
        # offseting row and column indices to fit python indexing
        toadata[chunks[i]:chunks[i+1],0] = toadata[chunks[i]:chunks[i+1],0]-1
        toadata[chunks[i]:chunks[i+1],1] = toadata[chunks[i]:chunks[i+1],1]-1
        gc.collect()
   
    ProgressParallel(
        total=nchunks,
        n_jobs=n_jobs,
        require='sharedmem'
    )(
        delayed(load_chunk)(i) for i in range(nchunks)
    )
    
    print("    generating csr sparse")
    gc.collect()

    toadata=sp.csr_matrix((toadata[:,2], (toadata[:,0], toadata[:,1])),
                            shape=shape,dtype=np.float32)
    return sc.AnnData(toadata)


if __name__=="__main__":
    adata=read_mtx_gpu("./matrix.mtx",nchunks=10)
