import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
from matplotlib import rcParams
import os
import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run C2L on bash")
    parser.add_argument("--outs_path", type=str, required=True, help="Visium outs folder")
    parser.add_argument("--infer_csv", type=str, required=True, help="The infer C2L csv output")
    parser.add_argument("--N_cells_per_location", type=float, required=True, help="N_cells_per_location")
    parser.add_argument("--max_epochs", type=float, required=True, help="max_epochs")
    parser.add_argument("--output", type=str, required=True, help="Return csv file")
    args = parser.parse_args()
    
    # Import file
    adata_vis = sc.read_visium(args.outs_path)
    print("No Processed file found: Start Pipeline")
    adata_vis.var_names_make_unique()
    adata_vis.var['SYMBOL'] = adata_vis.var_names
    adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()

    # Read inf_aver
    inf_aver = pd.read_csv(args.infer_scv, index_col=0)
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)
    mod = cell2location.models.Cell2location(adata_vis, cell_state_df=inf_aver,N_cells_per_location=N_cells_per_location,detection_alpha=detection_alpha)
    mod.train(max_epochs=max_epochs,batch_size=None,train_size=1,use_gpu=False)
    adata_vis = mod.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False})
    
    cell_types = adata_vis.obsm["q05_cell_abundance_w_sf"]
    cell_types.to_csv(args.output)
    
    
    
    
    
    
    
