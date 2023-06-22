#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## Load all required packages

import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import pandas as pd
import scvelo as scv
import anndata as ad
import os
import glob


# In[ ]:


os.chdir("<path/to/files>")


# ## R script to convert Seurat object to anndata
# ### Refer SeuratObj_to_anndata.R and then run the code below

# In[ ]:



# load anndata object
adata = sc.read_mtx("counts.mtx").transpose()

# load cell metadata
cell_meta = pd.read_csv("metadata.csv", index_col=0)

# load gene names
gene_names = pd.read_csv("gene_names.csv", header=None, index_col=0).index.tolist()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.var_names = gene_names
adata.obs.index = adata.obs['barcode']

# load dimensional reduction
pca = pd.read_csv("pca.csv", index_col=0)

# set pca and umap
adata.obsm['X_pca'] = pca.values
adata.obsm['X_umap'] = adata.obs[['UMAP_1', 'UMAP_2']].values


# In[ ]:


adata


# In[ ]:


adata.obs.head()


# # plot a UMAP colored by clusters to test:
# sc.pl.umap(adata, color=['cluster'], frameon=False,title='UMAP based on Cell Types',legend_fontsize='large')
# sc.pl.umap(adata, color=['cell_type'], frameon=False,title='UMAP based on Cell Types',legend_fontsize='large')
# sc.pl.umap(adata, color=['object'], frameon=False,title='UMAP based on Cell Types',legend_fontsize='large')
# 

# In[ ]:


# load loom files for spliced/unspliced matrices for each sample:

ldata1 = scv.read('loom_files/sample1.loom', cache=False)
ldata2 = scv.read('loom_files/sample2.loom', cache=False)
ldata3 = scv.read('loom_files/sample3.loom', cache=False)
ldata4 = scv.read('loom_files/sample4.loom', cache=False)


# In[ ]:


## create a list depending on the number on samples
ldata_list = [ldata1, ldata2, ldata3, ldata4]

## add the sample names either based on loom file names or from metadata
sample_names = ["sample1", "sample2", "sample3", "sample4"]
df = cell_meta

## generate new column in metadata for the barcode suffixes that will be used to 
df['suffix'] = df['barcode'].str.split('-', n=1).str[1]
df['suffix'] = '-' + df['suffix']

new_df = df[['object', 'suffix']].drop_duplicates()


# In[ ]:


## Rename barcodes in ldata objects
## This allows renaming the barcodes for loom files to match with the metadata
## Here we modify the observation names of each object in ldata_list based on the associated suffix value 
## from the metadata

for i in range(len(ldata_list)):
    ldata = ldata_list[i]
    sample_name = sample_names[i]
    suffix = new_df[new_df['object'] == sample_name]['suffix'].values[0]

    barcodes = [bc.split(':')[1] if len(bc.split(':')) > 1 else bc for bc in ldata.obs_names]
    barcodes = [bc[:-1] + suffix for bc in barcodes]
    ldata.obs_names = barcodes

    ldata_list[i] = ldata


# In[ ]:


## Making var names unique

for data in ldata_list:
    data.var_names_make_unique()


# In[ ]:


# make obs names unique

for data in ldata_list:
    data.obs_names_make_unique()
    


# In[ ]:


## Remove the last two characters(normally we we have more than 10 samples) 
## from index to match with anndata barcodes

ldata.obs.index = ldata.obs.index.map(lambda x: str(x)[:-2])
ldata.obs.index


# In[ ]:


ldata.obs.index.name = 'barcode'
#ldata.obs.head()


# In[ ]:


## concatenate all the loom files

ldata = ldata1.concatenate([ldata2,ldata3,ldata4])


# In[ ]:


## Check dimensions for adata

adata


# In[ ]:


## Check dimensions for ldata

ldata


# ## Merge the files

# In[ ]:


## merge matrices into the original adata object
## may need to change the id_length depending on the lenght of cell_ID
## id_length parameter specifies the length of the unique identifier (ID) that is assigned to 
## each merged observation in the resulting AnnData object

adata2 = scv.utils.merge(adata, ldata, id_length=16)


# In[ ]:


adata2


# In[ ]:


## save intermediate merged object

adata2.write_h5ad("<path/to/dir>/merged_adata.h5ad")


# ### Inspect the proportion of spliced and unspliced transcripts in each of our cell clusters
# ##### pie chart of spliced/unspliced proprtions.

# In[ ]:


## Need to run this step if scv.pl.proportions fails.
## here we make the annotation column type as category

#adata2.obs['<annotation/column/name>'] = adata2.obs['<annotation/column/name>'].astype('category')


# In[ ]:


# merge matrices into the original adata object

scv.pl.proportions(adata2, groupby='<annotation/column/name>',figsize=(27, 7),fontsize=11)


# In[ ]:


## Preprocess the Data

scv.pp.filter_and_normalize(adata2)
scv.pp.moments(adata2, n_pcs=30, n_neighbors=30)


# In[ ]:


## compute velocity
## default mode is stochastic. You can opt for deterministic or dynamical modeling as well
## refer to scvelo documentation for more details.
## Downstreaming may change depending on the model user selects

scv.tl.velocity(adata2, vkey='velocity', mode='stochastic')


# In[ ]:


scv.tl.velocity_graph(adata2)


# In[ ]:


adata2


# In[ ]:


## save the object after comouting velocity

adata2.write_h5ad("<path/to/dir>/velocity_adata.h5ad")


# In[ ]:


## Read the .h5ad object

#adata2 = anndata.read_h5ad("<path/to/dir>/velocity_adata.h5ad")


# ## User can either decide to run the visualization as part of the same script and save the results
# ## Or run this locally by taking the above velocity_adata.h5ad object to explore more 
# ## Refer scvelo documentation for more - https://scvelo.readthedocs.io/en/stable/VelocityBasics/

# In[ ]:


## Visualize velocity fields

scv.pl.velocity_embedding(adata2, basis='umap', frameon=False)


# In[ ]:


## make sure to add annotation/column/name that has the seurat cluster annotation information in the metadata

scv.pl.velocity_embedding_grid(adata2, basis='umap', color='<annotation/column/name>',title='', scale=0.25, 
                               size=10,figsize=(9,9), dpi=300)


# In[ ]:


scv.pl.velocity_embedding_stream(adata2, basis='umap', color=['<annotation/column/name>'], fontsize=8,
                                 title='Stream plot of velocities on the embedding',
                                 figsize=(10,10),save="Stream_plot.png")



# In[ ]:


scv.tl.paga(adata2, groups='<annotation/column/name>')

## If you want to view this as a dataframe
#df = scv.get_df(adata2, 'paga/transitions_confidence', precision=2).T


# In[ ]:


scv.pl.paga(adata2, basis='umap', size=75, alpha=1,figsize=(15,12),legend_fontsize=18,fontsize=22,
            min_edge_width=2, node_size_scale=1.5,title ='PAGA graph', save="PAGA_graph.pdf")





