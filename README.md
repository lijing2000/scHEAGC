scHEAGC

This toolkit is built upon the following dependencies:
- scanpy
- anndata
- pandas
- numpy
- torch
- scipy

The encoder module supports two models: VQVAE and GAE.

Main functionalities consist of two stages: pretraining and training.

Usage Steps:
1. Pretrain: Run `main.py --pretrain True`
2. Train: Run `main.py`


dataset:
PBMC-10k: https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-10-k-1-standard-2-0-0
PBMC-3k：https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0}
GSE126074：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126074
GSE100866：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866
CITE_PBMC：Clustering of single-cell multi-omics data with a multimodal deep learning method.
GSE140203：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140203
mouse_retina：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201402
