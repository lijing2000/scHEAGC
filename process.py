import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
from scipy.sparse import csr_matrix
import numpy as np


rna = ad.read_h5ad("./input/Ma/Ma-2020-RNA.h5ad")
#atac = ad.read_h5ad("./input/Ma/Ma-2020-ATAC.h5ad")

atac.X = csr_matrix(atac.X)

sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")


#scglue.data.lsi(atac, n_components=100, n_iter=15)

np.save('./input/Ma/RNA_fea.npy',rna.obsm['X_pca'])
np.save('./input/Ma/ATAC_fea.npy',atac.obsm['X_lsi'])

print(rna.obs["batch"])
