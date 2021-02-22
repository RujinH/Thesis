"""

Pack the scRNA-seq data using scanpy, prep for scran normalisation

"""

import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from rpy2.robjects.packages import importr
plt.rcParams['figure.figsize'] = (8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 10
sc.settings.autoshow = False

sam1 = sc.read("../scte_data/ss.gastrulation_E6.5_Sam1.h5ad")    ; sam1.obs['stage'] = "E6.5"   ; sam1.obs['replicate'] = "E6.5-1"
sam2 = sc.read("../scte_data/ss.gastrulation_E6.5_Sam5.h5ad")    ; sam2.obs['stage'] = "E6.5"   ; sam2.obs['replicate'] = "E6.5-2"
#sam3 = sc.read("../scte_data/ss.gastrulation_E6.5_Sam18.h5ad")   ; sam3.obs['stage'] = "E6.5"   ; sam3.obs['replicate'] = "E6.5-3"
#sam4 = sc.read("../scte_data/ss.gastrulation_E6.75_Sam7.h5ad")   ; sam4.obs['stage'] = "E6.75"  ; sam4.obs['replicate'] = "E6.75-1"
sam5 = sc.read("../scte_data/ss.gastrulation_E7.0_Sam10.h5ad")   ; sam5.obs['stage'] = "E7.0"   ; sam5.obs['replicate'] = "E7.0-1"
#sam6 = sc.read("../scte_data/ss.gastrulation_E7.0_Sam15.h5ad")   ; sam6.obs['stage'] = "E7.0"   ; sam6.obs['replicate'] = "E7.0-3"
sam7 = sc.read("../scte_data/ss.gastrulation_E7.0_Sam30.h5ad")   ; sam7.obs['stage'] = "E7.0"   ; sam7.obs['replicate'] = "E7.0-4"
sam8 = sc.read("../scte_data/ss.gastrulation_E7.0_Sam31.h5ad")   ; sam8.obs['stage'] = "E7.0"   ; sam8.obs['replicate'] = "E7.0-5"
sam9 = sc.read("../scte_data/ss.gastrulation_E7.0_Sam32.h5ad")   ; sam9.obs['stage'] = "E7.0"   ; sam9.obs['replicate'] = "E7.0-6"
sam10 = sc.read("../scte_data/ss.gastrulation_E7.25_Sam23.h5ad") ; sam10.obs['stage'] = "E7.25" ; sam10.obs['replicate'] = "E7.25-2"
sam11 = sc.read("../scte_data/ss.gastrulation_E7.25_Sam26.h5ad") ; sam11.obs['stage'] = "E7.25" ; sam11.obs['replicate'] = "E7.25-3"
sam12 = sc.read("../scte_data/ss.gastrulation_E7.25_Sam27.h5ad") ; sam12.obs['stage'] = "E7.25" ; sam12.obs['replicate'] = "E7.25-4"
sam13 = sc.read("../scte_data/ss.gastrulation_E7.5_Sam2.h5ad")   ; sam13.obs['stage'] = "E7.5"  ; sam13.obs['replicate'] = "E7.5-1"
sam14 = sc.read("../scte_data/ss.gastrulation_E7.5_Sam3.h5ad")   ; sam14.obs['stage'] = "E7.5"  ; sam14.obs['replicate'] = "E7.5-2"
sam15 = sc.read("../scte_data/ss.gastrulation_E7.5_Sam4.h5ad")   ; sam15.obs['stage'] = "E7.5"  ; sam15.obs['replicate'] = "E7.5-3"
sam16 = sc.read("../scte_data/ss.gastrulation_E7.5_Sam6.h5ad")   ; sam16.obs['stage'] = "E7.5"  ; sam16.obs['replicate'] = "E7.5-4"
sam17 = sc.read("../scte_data/ss.gastrulation_E7.5_Sam19.h5ad")  ; sam17.obs['stage'] = "E7.5"  ; sam17.obs['replicate'] = "E7.5-5"
sam18 = sc.read("../scte_data/ss.gastrulation_E7.5_Sam20.h5ad")  ; sam18.obs['stage'] = "E7.5"  ; sam18.obs['replicate'] = "E7.5-6"
sam19 = sc.read("../scte_data/ss.gastrulation_E7.75_Sam8.h5ad")  ; sam19.obs['stage'] = "E7.75" ; sam19.obs['replicate'] = "E7.75-1"
sam20 = sc.read("../scte_data/ss.gastrulation_E7.75_Sam9.h5ad")  ; sam20.obs['stage'] = "E7.75" ; sam20.obs['replicate'] = "E7.75-2"
sam21 = sc.read("../scte_data/ss.gastrulation_E7.75_Sam12.h5ad") ; sam21.obs['stage'] = "E7.75" ; sam21.obs['replicate'] = "E7.75-3"
sam22 = sc.read("../scte_data/ss.gastrulation_E7.75_Sam13.h5ad") ; sam22.obs['stage'] = "E7.75" ; sam22.obs['replicate'] = "E7.75-4"
sam23 = sc.read("../scte_data/ss.gastrulation_E8.0_Sam16.h5ad")  ; sam23.obs['stage'] = "E8.0"  ; sam23.obs['replicate'] = "E8.0-1"
sam24 = sc.read("../scte_data/ss.gastrulation_E8.0_Sam33.h5ad")  ; sam24.obs['stage'] = "E8.0"  ; sam24.obs['replicate'] = "E8.0-2"
sam25 = sc.read("../scte_data/ss.gastrulation_E8.0_Sam34.h5ad")  ; sam25.obs['stage'] = "E8.0"  ; sam25.obs['replicate'] = "E8.0-3"
sam26 = sc.read("../scte_data/ss.gastrulation_E8.0_Sam35.h5ad")  ; sam26.obs['stage'] = "E8.0"  ; sam26.obs['replicate'] = "E8.0-4"
sam27 = sc.read("../scte_data/ss.gastrulation_E8.25_Sam24.h5ad") ; sam27.obs['stage'] = "E8.25" ; sam27.obs['replicate'] = "E8.25-1"
sam28 = sc.read("../scte_data/ss.gastrulation_E8.25_Sam25.h5ad") ; sam28.obs['stage'] = "E8.25" ; sam28.obs['replicate'] = "E8.25-2"
sam29 = sc.read("../scte_data/ss.gastrulation_E8.25_Sam28.h5ad") ; sam29.obs['stage'] = "E8.25" ; sam29.obs['replicate'] = "E8.25-3"
sam30 = sc.read("../scte_data/ss.gastrulation_E8.5_Sam17.h5ad")  ; sam30.obs['stage'] = "E8.5"  ; sam30.obs['replicate'] = "E8.5-1"
sam31 = sc.read("../scte_data/ss.gastrulation_E8.5_Sam29.h5ad")  ; sam31.obs['stage'] = "E8.5"  ; sam31.obs['replicate'] = "E8.5-2"
sam32 = sc.read("../scte_data/ss.gastrulation_E8.5_Sam36.h5ad")  ; sam32.obs['stage'] = "E8.5"  ; sam32.obs['replicate'] = "E8.5-3"
sam33 = sc.read("../scte_data/ss.gastrulation_E8.5_Sam37.h5ad")  ; sam33.obs['stage'] = "E8.5"  ; sam33.obs['replicate'] = "E8.5-4"
sam34 = sc.read("../scte_data/ss.gastrulation_mixed_Sam21.h5ad") ; sam34.obs['stage'] = "mixed" ; sam34.obs['replicate'] = "mixed-1"
sam35 = sc.read("../scte_data/ss.gastrulation_mixed_Sam22.h5ad") ; sam35.obs['stage'] = "mixed" ; sam35.obs['replicate'] = "mixed-2"

print('Loaded Samples...')

# Do very simple prefiltering:
samples = [sam1, sam2, #sam3, sam4,
            sam5, #sam6,
            sam7, sam8, sam9, sam10,
            sam11, sam12, sam13, sam14, sam15,
            sam16, sam17, sam18, sam19, sam20,
            sam21, sam22, sam23, sam24, sam25,
            sam26, sam27, sam28, sam29, sam30,
            sam31, sam32, sam33, sam34, sam35]

# Quick pre-filtering, these should be low, otherwise it can mess up downstream analysis, but also can get rid of trivial uninteresting things
[sc.pp.filter_cells(sam, min_genes=2000) for sam in samples]
[sc.pp.filter_cells(sam, max_counts=100000) for sam in samples]
[sc.pp.filter_cells(sam, min_counts=5000) for sam in samples]
# Do not filter gene here; concatenate joins on the union, so if a gene fails in a single sample, it will also be deleted from all other samples;

print('Concatenating')
adata = sam1.concatenate(samples[1:])

del samples

adata.X = adata.X.astype('float32')

print(adata)

sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-pre-norm-replicates.pdf')

# Base filtering for trivial QC failures:
sc.pp.filter_cells(adata, min_genes=3000)
sc.pp.filter_cells(adata, min_counts=8000)
sc.pp.filter_cells(adata, max_counts=100000)
sc.pp.filter_genes(adata, min_cells=50) # Only filter genes here;

print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

#sc.pl.violin(adata, ['n_genes','n_counts'], groupby='stage', size=0, log=False, cut=0, show=False, save='qc1.pdf')
sc.pl.violin(adata, ['n_genes','n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-replicates.pdf')

p = sb.distplot(adata.obs['n_counts'], kde=False)
p.get_figure().savefig('figures/distplot_ncounts1.pdf')
p = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']<4000], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ncounts2.pdf')
p = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']>10000], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ncounts3.pdf')
#Thresholding decision: genes
p = sb.distplot(adata.obs['n_genes'], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ngenes1.pdf')
p = sb.distplot(adata.obs['n_genes'][adata.obs['n_genes']<2000], kde=False, bins=60)
p.get_figure().savefig('figures/distplot_ngenes2.pdf')

print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

adata.write('./raw_data.h5ad')
