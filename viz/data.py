import os
from typing import Set

import anndata as ad
import numpy as np
import scanpy as sc
import pandas as pd


def read_ct_data(path: str) -> ad.AnnData:
    ct = sc.read(path, backed='r')  # Low memory mode
    return ct


def prune_ct_data(path: str, inplace=False):
    """
    Strip out all unused data from the data file to improve performance.
    """
    if inplace:
        pruned_path = path
    else:
        pruned_path = path.replace('.h5ad', '_pruned.h5ad')
        if os.path.exists(pruned_path):
            return pruned_path

    ct = sc.read(path, backed='r')  # Low memory mode
    new_ct = ad.AnnData(np.zeros_like(ct.X, dtype=np.float32))
    # Select only the data we need
    new_ct.obs['cell type'] = ct.obs['cell type']
    new_ct.obs['target'] = ct.obs['target']
    new_ct.obs['ligand'] = ct.obs['ligand'].astype(int)
    new_ct.obs['receptor'] = ct.obs['receptor'].astype(int)
    #new_ct.obs['DA_score'] = ct.obs['DA_score']
    new_ct.var.index = ct.var.index
    new_ct.layers['fdr'] = ct.layers['fdr']
    new_ct.layers['log2FC'] = ct.layers['log2FC']
    new_ct.layers['p.bonf'] = ct.layers['p.bonf']
    new_ct.layers['fdr.i1'] = ct.layers['fdr.i1']
    del ct
    gene_is_ligand = dict()
    for gene in new_ct.var.index:
        gene_is_ligand[gene] = new_ct[new_ct.obs['target'] == gene].obs.head(1)['ligand'].get(0, False)
    new_ct.uns['gene_is_ligand'] = gene_is_ligand
    if inplace:
        os.remove(path)
    new_ct.write_h5ad(pruned_path)
    return pruned_path


def read_interactions(path: str, condition_name='highCIN_vs_lowCIN') -> pd.DataFrame:
    df = pd.read_table(path)
    df["MAST_fdr_ligand"] = df[f'MAST_fdr_{condition_name}_ligand']
    df["MAST_log2FC_ligand"] = df[f'MAST_log2FC_{condition_name}_ligand']
    df["MAST_fdr_receptor"] = df[f'MAST_fdr_{condition_name}_receptor']
    df["MAST_log2FC_receptor"] = df[f'MAST_log2FC_{condition_name}_receptor']
    df["numSigI1_fdr_receptor"] = df[f'numSigI1_fdr25_receptor']
    df["numDEG_fdr_receptor"] = df[f'numDEG_fdr25_receptor']
    return df


def get_downstream_ligands(ct_res, receptor, receptor_celltype, min_log2fc=0.01, max_fdr=0.05) -> Set[str]:
    """
    Get activated ligands downstream of a given receptor.
    :param ct_res: The results object.
    :param receptor: The receptor to get downstream ligands for.
    :param receptor_celltype: The cell type of the receptor.
    :param min_log2fc: The minimum log2FC to consider.
    :param max_fdr: The maximum FDR to consider.
    :return: Set of downstream ligands.
    """
    gene2ligand = ct_res.uns['gene_is_ligand']
    receptor_res = ct_res[(ct_res.obs['target'] == receptor) & (ct_res.obs['cell type'] == receptor_celltype)]
    if receptor_res.shape[0] < 1:
        return set()
    downstream_genes = pd.DataFrame({'gene': [g for g in receptor_res.var.index],
                                     'log2FC': receptor_res.layers['log2FC'].flatten(),
                                     'fdr': receptor_res.layers['fdr'].flatten(),
                                     'p.bonf': receptor_res.layers['p.bonf'].flatten(),
                                     'is_ligand': [bool(gene2ligand[g]) for g in receptor_res.var.index]})
    downstream_ligands = downstream_genes[downstream_genes.is_ligand &
                                          (abs(downstream_genes.log2FC) > min_log2fc) &
                                          (downstream_genes.fdr <= max_fdr)].gene

    return set(downstream_ligands)

#
# def get_effective_logfc(from_logfc: float, to_logfc: float) -> float:
#     """
#     Get effective logFC, correcting for the ordering of the ligand and receptor.
#     :param from_logfc: LogFC of the ligand.
#     :param to_logfc: LogFC of the receptor.
#     :return: Direction corrected logFC.
#     """
#     # logfc from mast
#     if from_logfc < 0 and to_logfc < 0:
#         # Same direction, meaning positive correlation
#         return abs(to_logfc)
#     elif from_logfc > 0 and to_logfc > 0:
#         # Same direction, meaning positive correlation
#         return to_logfc
#     elif from_logfc < 0 and to_logfc > 0:
#         # Differing directions, negative correlation
#         return -to_logfc
#     elif from_logfc > 0 and to_logfc < 0:
#         # Differing directions, negative correlation
#         return to_logfc
#     else:
#         # Both 0
#         return 0


def get_diff_abundance(ct_res, celltype, target):
    return ct_res[(ct_res.obs['cell type'] == celltype) & (ct_res.obs['target'] == target)].obs['DA_score'].iloc[0]


def get_interaction_fdr(ct_res, celltype, target, gene):
    subset_ct = ct_res[(ct_res.obs['cell type'] == celltype) & (ct_res.obs['target'] == target)].layers['fdr.i1']
    return subset_ct[0, ct_res.var.index == gene][0]

