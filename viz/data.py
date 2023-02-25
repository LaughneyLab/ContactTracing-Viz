import itertools
import os
from glob import glob
from typing import Optional

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd


# TODO: Remove
def _get_original_data():
    df = pd.read_csv('old_data/circle_plot_tabular.tsv', sep='\t')
    df['cell_type'] = df['cell type']
    df['MAST_fdr_max'] = df['MAST_fdr']
    df['numSigI1_max'] = df['numSigI1_fdr25_max']
    df['cell_type_dc1_max'] = df['cell_type_dc1']
    return df


def _filename(prefix, condition, fdr, extension):
    return os.path.join(prefix, f'{condition}_{fdr}.{extension}')


def _normalize_condition_name(condition: str) -> str:
    if condition == 'cin':
        return 'highCIN_vs_lowCIN'
    elif condition == 'sting':
        return 'highCIN_vs_noSTING'
    elif condition == 'max':
        return condition
    elif condition == 'highCIN_vs_lowCIN' or condition == 'highCIN_vs_noSTING':
        return condition
    else:
        raise ValueError(f'Unknown condition {condition}')


def _normalize_fdr_name(fdr: str) -> str:
    return fdr.split('.')[-1].replace("fdr", "")


def read_circos_file(condition: str, fdr: str) -> pd.DataFrame:
    prefix = 'data/compiled/circos'
    df = pd.read_csv(_filename(prefix, _normalize_condition_name(condition), _normalize_fdr_name(fdr), 'csv'))
    return df


def read_ligand_effect_for(condition: str, celltype: str, target: str) -> pd.DataFrame:
    prefix = 'data/compiled/ligand_effects/'
    condition = _normalize_condition_name(condition)
    prefix += condition + '_' + celltype.replace('/', '') + '_' + target + '.csv'
    if not os.path.exists(prefix):
        return None
    df = pd.read_csv(prefix)
    return df


def read_all_effects_for(condition: str, target: str) -> pd.DataFrame:
    prefix = 'data/compiled/ligand_effects/'
    condition = _normalize_condition_name(condition)
    prefix += condition + '_'
    df = pd.DataFrame()
    for filename in glob(prefix + '*' + target + '.csv'):
        df = pd.concat([df, pd.read_csv(filename)])
    if df.shape[0] == 0:
        return None
    return df


def read_interactions_file(condition: str, fdr: str) -> pd.DataFrame:
    prefix = 'data/compiled/interactions'
    df = pd.read_csv(_filename(prefix, _normalize_condition_name(condition), _normalize_fdr_name(fdr), 'csv'))
    return df


def read_ligand_receptor_file(filename: str = 'data/allgenes/all_interactions.tsv') -> pd.DataFrame:
    df = pd.read_csv(filename, sep='\t')[['ligand', 'receptor']].drop_duplicates()
    return df


def _combine_obs_for_interactions(interactions: pd.DataFrame,
                                  cin_adata: ad.AnnData,
                                  sting_adata: ad.AnnData,
                                  condition_name: str,
                                  fdr: str):
    fdr_cutoff = float(f".{fdr}")
    if condition_name == 'max':
        return

    if condition_name == 'highCIN_vs_lowCIN':
        adata = cin_adata
    elif condition_name == 'highCIN_vs_noSTING':
        condition_name = "highCIN_vs_lowCIN"  # FIXME: The sting file uses this name
        adata = sting_adata
    else:
        raise ValueError(f'Unknown condition {condition_name}')

    all_celltypes = set(adata.obs['cell type'])
    full_df = pd.DataFrame()
    for donor_celltype, target_celltype in itertools.combinations(all_celltypes, 2):
        lig_adata = adata[(adata.obs['cell type'] == donor_celltype) &
                              (adata.obs['target'].isin(interactions.ligand))]
        rec_adata = adata[(adata.obs['cell type'] == target_celltype) &
                          (adata.obs['target'].isin(interactions.receptor))]
        # TODO: Expression value?
        ligand_df = pd.DataFrame({
            'ligand': lig_adata.obs.target,
            'MAST_log2FC_ligand': lig_adata.obs[f'MAST_log2FC_{condition_name}'],
            'MAST_fdr_ligand': lig_adata.obs[f'MAST_fdr_{condition_name}'],
            'cell_type_ligand_dc1': lig_adata.obs['cell_type_dc1'],
            'DA_ligand': lig_adata.obs['DA_score'],
        })
        ligand_df['cell_type_ligand'] = donor_celltype

        receptor_df = pd.DataFrame({
            'receptor': rec_adata.obs.target,
            'MAST_log2FC_receptor': rec_adata.obs[f'MAST_log2FC_{condition_name}'],
            'MAST_fdr_receptor': rec_adata.obs[f'MAST_fdr_{condition_name}'],
            'cell_type_receptor_dc1': rec_adata.obs['cell_type_dc1'],
            'DA_receptor': rec_adata.obs['DA_score'],
            # Receptor exclusive features
            'numDEG': (rec_adata.layers['fdr'] < fdr_cutoff).sum(axis=1),
            'numSigI1': (rec_adata.layers['fdr.i1'] < fdr_cutoff).sum(axis=1),
        })
        receptor_df['cell_type_receptor'] = target_celltype

        interactions_ct = pd.merge(interactions, ligand_df, on='ligand')
        interactions_ct = pd.merge(interactions_ct, receptor_df, on='receptor')

        full_df = pd.concat([full_df, interactions_ct])

    return full_df


def _compile_interactions(interactions: pd.DataFrame,
                          cin_adata: ad.AnnData,
                          sting_adata: ad.AnnData,
                          condition_name: str,
                          fdr: str):
    prefix = 'data/compiled/interactions'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, condition_name, fdr, 'csv')
    if os.path.exists(file):
        return pd.read_csv(file)

    interactions_df = _combine_obs_for_interactions(interactions, cin_adata, sting_adata, condition_name, fdr)
    interactions_df.to_csv(file, index=False)
    return interactions_df


def _merge_max_interactions(cin_interactions: pd.DataFrame,
                            sting_interactions: pd.DataFrame,
                            fdr: str):
    prefix = 'data/compiled/interactions'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, 'max', fdr, 'csv')
    if os.path.exists(file):
        return

    both_interactions = pd.merge(cin_interactions, sting_interactions,
                                 on=['cell_type_ligand', 'ligand',
                                     'cell_type_receptor', 'receptor'],
                                 suffixes=('_cin', '_sting'))

    mast_receptor_mask = both_interactions['MAST_log2FC_receptor_cin'].abs() > both_interactions['MAST_log2FC_receptor_sting'].abs()
    both_interactions['MAST_log2FC_receptor'] = np.where(mast_receptor_mask, both_interactions['MAST_log2FC_receptor_cin'], both_interactions['MAST_log2FC_receptor_sting'])
    both_interactions['MAST_fdr_receptor'] = np.where(mast_receptor_mask, both_interactions['MAST_fdr_receptor_cin'], both_interactions['MAST_fdr_receptor_sting'])

    mast_ligand_mask = both_interactions['MAST_log2FC_ligand_cin'].abs() > both_interactions['MAST_log2FC_ligand_sting'].abs()
    both_interactions['MAST_log2FC_ligand'] = np.where(mast_ligand_mask, both_interactions['MAST_log2FC_ligand_cin'], both_interactions['MAST_log2FC_ligand_sting'])
    both_interactions['MAST_fdr_ligand'] = np.where(mast_ligand_mask, both_interactions['MAST_fdr_ligand_cin'], both_interactions['MAST_fdr_ligand_sting'])

    both_interactions['cell_type_ligand_dc1'] = np.where(both_interactions['cell_type_ligand_dc1_cin'] > both_interactions['cell_type_ligand_dc1_sting'],
                                                         both_interactions['cell_type_ligand_dc1_cin'], both_interactions['cell_type_ligand_dc1_sting'])
    both_interactions['cell_type_receptor_dc1'] = np.where(both_interactions['cell_type_receptor_dc1_cin'] > both_interactions['cell_type_receptor_dc1_sting'],
                                                            both_interactions['cell_type_receptor_dc1_cin'], both_interactions['cell_type_receptor_dc1_sting'])

    both_interactions['DA_ligand'] = np.where(both_interactions['DA_ligand_cin'] > both_interactions['DA_ligand_sting'],
                                                  both_interactions['DA_ligand_cin'], both_interactions['DA_ligand_sting'])
    both_interactions['DA_receptor'] = np.where(both_interactions['DA_receptor_cin'] > both_interactions['DA_receptor_sting'],
                                                    both_interactions['DA_receptor_cin'], both_interactions['DA_receptor_sting'])

    both_interactions['numDEG'] = np.where(both_interactions['numDEG_cin'] > both_interactions['numDEG_sting'],
                                             both_interactions['numDEG_cin'], both_interactions['numDEG_sting'])
    both_interactions['numSigI1'] = np.where(both_interactions['numSigI1_cin'] > both_interactions['numSigI1_sting'],
                                                  both_interactions['numSigI1_cin'], both_interactions['numSigI1_sting'])

    both_interactions.to_csv(file, index=False)
    return both_interactions


def _compile_ligand_effects(interactions: pd.DataFrame,
                            cin_adata: ad.AnnData,
                            sting_adata: ad.AnnData,
                            condition_name: str):
    prefix = 'data/compiled/ligand_effects'
    os.makedirs(prefix, exist_ok=True)
    prefix += '/' + condition_name + "_"

    if condition_name == 'max':
        return

    if condition_name == 'highCIN_vs_lowCIN':
        adata = cin_adata
    elif condition_name == 'highCIN_vs_noSTING':
        condition_name = "highCIN_vs_lowCIN"  # FIXME: The sting file uses this name
        adata = sting_adata
    else:
        raise ValueError(f'Unknown condition {condition_name}')

    all_celltypes = set(adata.obs['cell type'])
    genes = [g for g in adata.var.index]
    for celltype in all_celltypes:
        for receptor in set(interactions.receptor):
            filename = prefix + celltype.replace('/', '') + "_" + receptor + ".csv"
            if os.path.exists(filename):
                continue
            rec_adata = adata[(adata.obs['cell type'] == celltype) &
                              (adata.obs['target'] == receptor)]
            if rec_adata.n_obs == 0:
                continue
            df = pd.DataFrame({
                'gene': genes,
                'gene_is_ligand': rec_adata.var.index.isin(interactions.ligand),
                'gene_is_receptor': rec_adata.var.index.isin(interactions.receptor),
                'log2FC': rec_adata.layers['log2FC'].flatten(),
                'fdr': rec_adata.layers['fdr'].flatten(),
                'i1.fdr': rec_adata.layers['fdr.i1'].flatten()
            })
            df['cell_type'] = celltype
            df['target'] = receptor
            df['receptor'] = True
            df['ligand'] = False
            df.to_csv(filename, index=False)
        for ligand in set(interactions.ligand):
            filename = prefix + celltype.replace('/', '') + "_" + ligand + ".csv"
            if os.path.exists(filename):
                continue
            lig_adata = adata[(adata.obs['cell type'] == celltype) &
                              (adata.obs['target'] == ligand)]
            if lig_adata.n_obs == 0:
                continue
            df = pd.DataFrame({
                'gene': genes,
                'gene_is_ligand': lig_adata.var.index.isin(interactions.ligand),
                'gene_is_receptor': lig_adata.var.index.isin(interactions.receptor),
                'log2FC': lig_adata.layers['log2FC'].flatten(),
                'fdr': lig_adata.layers['fdr'].flatten(),
                'i1.fdr': lig_adata.layers['fdr.i1'].flatten()
            })
            df['cell_type'] = celltype
            df['target'] = ligand
            df['receptor'] = False
            df['ligand'] = True
            df.to_csv(filename, index=False)


def _merge_max_ligand_effects(cin_condition: str,
                              sting_condition: str):
    prefix = 'data/compiled/ligand_effects/'
    os.makedirs(prefix, exist_ok=True)
    cin_prefix = prefix + cin_condition + "_"
    sting_prefix = prefix + sting_condition + "_"
    prefix += 'max' + "_"

    for cin_file in glob(cin_prefix + "*.csv"):
        sting_file = sting_prefix + cin_file.replace(cin_prefix, '')
        file = prefix + cin_file.replace(cin_prefix, '')
        if not os.path.exists(sting_file) or os.path.exists(file):
            continue

        cin_interactions = pd.read_csv(cin_file)
        sting_interactions = pd.read_csv(sting_file)
        both_ligand_effects = pd.merge(cin_interactions, sting_interactions,
                                       on=['gene', 'cell_type', 'target', 'ligand', 'receptor', 'gene_is_ligand', 'gene_is_receptor'],
                                       suffixes=('_cin', '_sting'))

        logfc_mask = both_ligand_effects['log2FC_cin'].abs() > both_ligand_effects['log2FC_sting'].abs()
        both_ligand_effects['log2FC'] = np.where(logfc_mask, both_ligand_effects['log2FC_cin'], both_ligand_effects['log2FC_sting'])
        both_ligand_effects['fdr'] = np.where(logfc_mask, both_ligand_effects['fdr_cin'], both_ligand_effects['fdr_sting'])
        both_ligand_effects['i1.fdr'] = np.where(both_ligand_effects['i1.fdr_cin'] < both_ligand_effects['i1.fdr_sting'], both_ligand_effects['i1.fdr_cin'], both_ligand_effects['i1.fdr_sting'])

        both_ligand_effects.to_csv(file, index=False)


def _compile_circos(interactions: pd.DataFrame,
                    cin_adata: ad.AnnData,
                    sting_adata: ad.AnnData,
                    condition_name: str,
                    fdr: str):
    prefix = 'data/compiled/circos'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, condition_name, fdr, 'csv')
    fdr_cutoff = float(f".{fdr}")
    if os.path.exists(file):
        return pd.read_csv(file)

    if condition_name == 'max':
        return

    if condition_name == 'highCIN_vs_lowCIN':
        adata = cin_adata
    elif condition_name == 'highCIN_vs_noSTING':
        condition_name = "highCIN_vs_lowCIN"  # FIXME: The sting file uses this name
        adata = sting_adata
    else:
        raise ValueError(f'Unknown condition {condition_name}')

    df = pd.DataFrame({
        'cell_type': adata.obs['cell type'],
        'target': adata.obs['target'],
        'ligand': adata.obs['ligand'],
        'receptor': adata.obs['receptor'],
        'numDEG': (adata.layers['fdr'] < fdr_cutoff).sum(axis=1),
        'numSigI1': (adata.layers['fdr.i1'] < fdr_cutoff).sum(axis=1),
        'MAST_log2FC': adata.obs[f'MAST_log2FC_{condition_name}'],
        'MAST_fdr': adata.obs[f'MAST_fdr_{condition_name}'],
        'cell_type_dc1': adata.obs['cell_type_dc1'],
        'DA_score': adata.obs['DA_score'],
    })
    df.to_csv(file, index=False)
    return df


def _merge_max_circos(cin_circos: pd.DataFrame,
                      sting_circos: pd.DataFrame,
                      fdr: str):
    prefix = 'data/compiled/circos'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, 'max', fdr, 'csv')
    fdr_cutoff = float(f".{fdr}")
    if os.path.exists(file):
        return

    both_circos = pd.merge(cin_circos, sting_circos,
                            on=['cell_type', 'target', 'ligand', 'receptor'],
                            suffixes=('_cin', '_sting'))

    logfc_mask = both_circos['MAST_log2FC_cin'].abs() > both_circos['MAST_log2FC_sting'].abs()
    both_circos['MAST_log2FC'] = np.where(logfc_mask, both_circos['MAST_log2FC_cin'], both_circos['MAST_log2FC_sting'])
    both_circos['MAST_fdr'] = np.where(logfc_mask, both_circos['MAST_fdr_cin'], both_circos['MAST_fdr_sting'])

    both_circos['numDEG'] = np.where(both_circos['numDEG_cin'] > both_circos['numDEG_sting'], both_circos['numDEG_cin'], both_circos['numDEG_sting'])
    both_circos['numSigI1'] = np.where(both_circos['numSigI1_cin'] > both_circos['numSigI1_sting'], both_circos['numSigI1_cin'], both_circos['numSigI1_sting'])

    both_circos['cell_type_dc1'] = np.where(both_circos['cell_type_dc1_cin'] > both_circos['cell_type_dc1_sting'], both_circos['cell_type_dc1_cin'], both_circos['cell_type_dc1_sting'])
    both_circos['DA_score'] = np.where(both_circos['DA_score_cin'] > both_circos['DA_score_sting'], both_circos['DA_score_cin'], both_circos['DA_score_sting'])

    both_circos.to_csv(file, index=False)
    return both_circos


def compile_data(
        interactions_file: str = 'data/allgenes/all_interactions.tsv',
        cin_adata_file: str = 'data/highCIN_vs_lowCIN/saves/degboth.h5ad',
        sting_adata_file: str = 'data/highCIN_vs_noSTING/saves/degboth.h5ad'
):
    """
    Optimize the builtin files for fast loading
    """

    if os.path.exists('data/compiled/is_compiled'):
        return

    interactions = read_ligand_receptor_file(interactions_file)
    cin_adata = read_ct_data(cin_adata_file)
    sting_adata = read_ct_data(sting_adata_file)

    cin_condition = 'highCIN_vs_lowCIN'
    sting_condition = 'highCIN_vs_noSTING'

    for fdr in ['05', '25']:
        cin_interactions = _compile_interactions(interactions, cin_adata, sting_adata, cin_condition, fdr)
        sting_interactions = _compile_interactions(interactions, cin_adata, sting_adata, sting_condition, fdr)
        _merge_max_interactions(cin_interactions, sting_interactions, fdr)
        del cin_interactions
        del sting_interactions
        cin_circos = _compile_circos(interactions, cin_adata, sting_adata, cin_condition, fdr)
        sting_circos = _compile_circos(interactions, cin_adata, sting_adata, sting_condition, fdr)
        _merge_max_circos(cin_circos, sting_circos, fdr)
        del cin_circos
        del sting_circos
        _compile_ligand_effects(interactions, cin_adata, sting_adata, cin_condition)
        _compile_ligand_effects(interactions, cin_adata, sting_adata, sting_condition)
        _merge_max_ligand_effects(cin_condition, sting_condition)
    # touch file to indicate that we've compiled the data
    with open('data/compiled/is_compiled', 'w') as f:
        f.write('compiled')


def read_ct_data(path: str) -> ad.AnnData:
    ct = sc.read(path, backed='r')  # Low memory mode
    return ct


# def prune_ct_data(path: str, inplace=False):
#     """
#     Strip out all unused data from the data file to improve performance.
#     """
#     if inplace:
#         pruned_path = path
#     else:
#         pruned_path = path.replace('.h5ad', '_pruned.h5ad')
#         if os.path.exists(pruned_path):
#             return pruned_path
#
#     ct = sc.read(path, backed='r')  # Low memory mode
#     new_ct = ad.AnnData(np.zeros_like(ct.X, dtype=np.float32))
#     # Select only the data we need
#     new_ct.obs['cell type'] = ct.obs['cell type']
#     new_ct.obs['target'] = ct.obs['target']
#     new_ct.obs['ligand'] = ct.obs['ligand'].astype(int)
#     new_ct.obs['receptor'] = ct.obs['receptor'].astype(int)
#     #new_ct.obs['DA_score'] = ct.obs['DA_score']
#     new_ct.var.index = ct.var.index
#     new_ct.layers['fdr'] = ct.layers['fdr']
#     new_ct.layers['log2FC'] = ct.layers['log2FC']
#     new_ct.layers['p.bonf'] = ct.layers['p.bonf']
#     new_ct.layers['fdr.i1'] = ct.layers['fdr.i1']
#     del ct
#     gene_is_ligand = dict()
#     for gene in new_ct.var.index:
#         gene_is_ligand[gene] = new_ct[new_ct.obs['target'] == gene].obs.head(1)['ligand'].get(0, False)
#     new_ct.uns['gene_is_ligand'] = gene_is_ligand
#     if inplace:
#         os.remove(path)
#     new_ct.write_h5ad(pruned_path)
#     return pruned_path


def read_interactions_from_data_info(data_dict: dict) -> pd.DataFrame:
    inter_file = data_dict['interactions']
    condition_name = data_dict['condition']
    return read_interactions(inter_file, condition_name)


def read_obs_from_data_info(data_dict: dict) -> pd.DataFrame:
    obs_file = data_dict['obs']
    return read_obs(obs_file)


def read_interactions(path: str, condition_name='highCIN_vs_lowCIN') -> pd.DataFrame:
    df = pd.read_table(path)
    df["MAST_fdr_ligand"] = df[f'MAST_fdr_{condition_name}_ligand']
    df["MAST_log2FC_ligand"] = df[f'MAST_log2FC_{condition_name}_ligand']
    df["MAST_fdr_receptor"] = df[f'MAST_fdr_{condition_name}_receptor']
    df["MAST_log2FC_receptor"] = df[f'MAST_log2FC_{condition_name}_receptor']
    df["numSigI1_fdr_receptor"] = df[f'numSigI1_fdr25_receptor']
    df["numDEG_fdr_receptor"] = df[f'numDEG_fdr25_receptor']
    return df


def read_obs(path: str) -> pd.DataFrame:
    df = pd.read_table(path, index_col=0)
    conditions = list(set([x.replace('numSigI1', '').replace('_fdr25', '').replace('_fdr05', '').replace('_fdr', '') for x in
              df.columns if x.startswith('numSigI1')]))
    fdrs = list(set([x.replace('numSigI_', '').split("_")[1] for x in df.columns if x.startswith('numSigI1')]))
    if '_max' in conditions:
        test_conditions = [x for x in conditions if x != '_max']
        for fdr_name in fdrs:
            if f'numSigI1_{fdr_name}_max' not in df.columns:
                df[f'numSigI1_{fdr_name}_max'] = df[[f'numSigI1_{fdr_name}{condition}' for condition in test_conditions]].max(axis=1)
            if f'numDEG_{fdr_name}_max' not in df.columns:
                df[f'numDEG_{fdr_name}_max'] = df[[f'numDEG_{fdr_name}{condition}' for condition in test_conditions]].max(axis=1)
        if 'MAST_fdr_max' not in df.columns:
            df['MAST_log2FC_max'] = df[[f'MAST_log2FC{condition}' for condition in test_conditions]].max(axis=1)
            fdr_mask = ~df[[f'MAST_log2FC{condition}' for condition in test_conditions]].eq(df['MAST_log2FC_max'], axis=0)
            fdr_mask.columns = [f'MAST_fdr{condition}' for condition in test_conditions]
            fdr_mask = fdr_mask.astype(float).replace(1, np.nan)
            df['MAST_fdr_max'] = (df[[f'MAST_fdr{condition}' for condition in test_conditions]] + fdr_mask).max(axis=1)
    return df


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


# def get_diff_abundance(ct_res, celltype, target):
#     return ct_res[(ct_res.obs['cell type'] == celltype) & (ct_res.obs['target'] == target)].obs['DA_score'].iloc[0]


# def get_interaction_fdr(ct_res, celltype, target, gene):
#     subset_ct = ct_res[(ct_res.obs['cell type'] == celltype) & (ct_res.obs['target'] == target)].layers['fdr.i1']
#     return subset_ct[0, ct_res.var.index == gene][0]


# def get_interaction_logfc(ct_res, celltype, target, gene):
#     subset_ct = ct_res[(ct_res.obs['cell type'] == celltype) & (ct_res.obs['target'] == target)].layers['log2FC']
#     return subset_ct[0, ct_res.var.index == gene][0]


# def get_interaction_logfc_fdr(ct_res, celltype, target, gene):
#     subset_ct = ct_res[(ct_res.obs['cell type'] == celltype) & (ct_res.obs['target'] == target)].layers['fdr']
#     return subset_ct[0, ct_res.var.index == gene][0]

