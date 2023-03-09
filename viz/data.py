import itertools
import os
from collections import namedtuple
from typing import Dict

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd


def _filename(prefix, condition, fdr, extension):
    return os.path.join(prefix, f'{condition}_{fdr}.{extension}')


def _normalize_condition_name(condition: str) -> str:
    if condition == 'cin':
        return 'highCIN_vs_lowCIN'
    elif condition == 'sting' or condition == 'max':
        return 'highCIN_vs_noSTING'
    elif condition == 'highCIN_vs_lowCIN' or condition == 'highCIN_vs_noSTING':
        return condition
    else:
        raise ValueError(f'Unknown condition {condition}')


def _normalize_fdr_name(fdr: str) -> str:
    if fdr is None:
        return None
    return fdr.split('.')[-1].replace("fdr", "")


def read_circos_file(condition: str, fdr: str) -> pd.DataFrame:
    prefix = 'data/compiled/circos'
    df = pd.read_csv(_filename(prefix, _normalize_condition_name(condition), _normalize_fdr_name(fdr), 'csv'))
    return df


def read_ligand_effect_for(condition: str, target: str) -> pd.DataFrame:
    prefix = 'data/compiled/ligand_effects/'
    condition = _normalize_condition_name(condition)
    prefix += condition + '_' + target + '.csv'
    if not os.path.exists(prefix):
        return None
    df = pd.read_csv(prefix)
    return df


def read_interactions_file(condition: str, fdr: str) -> pd.DataFrame:
    prefix = 'data/compiled/interactions'
    df = pd.read_csv(_filename(prefix, _normalize_condition_name(condition), _normalize_fdr_name(fdr), 'csv'))
    return df


def read_ligand_receptor_file(filename: str = 'data/allgenes/all_interactions.tsv') -> pd.DataFrame:
    df = pd.read_csv(filename, sep='\t')[['ligand', 'receptor']].drop_duplicates()
    return df


def calculate_expression(celltype, target, exp_adata: ad.AnnData):
    filter = (exp_adata.obs['Cell Type'] == celltype)
    expressing_cells = (exp_adata[filter, target].layers['X'].flatten() > 0).sum()
    return expressing_cells / filter.sum()


def calculate_expressions(adata: ad.AnnData, exp_adata: ad.AnnData):
    expressions = []
    for i, row in adata.obs.iterrows():
        celltype = row['cell type']
        target = row['target']
        expressions.append(calculate_expression(celltype, target, exp_adata))
    return expressions


def calculate_aggregated_stats(adata: ad.AnnData, fdr: str) -> pd.DataFrame:
    fdr = float("0." + fdr)

    data = list()
    for i, row in adata.obs.iterrows():
        celltype = row['cell type']
        target = row['target']
        logfc = row['MAST_log2FC_highCIN_vs_lowCIN']
        logfc_fdr = row['MAST_fdr_highCIN_vs_lowCIN']
        selected = adata[(adata.obs['cell type'] == celltype) & (adata.obs['target'] == target)]
        numDEG = int((selected.layers['fdr'].flatten() < fdr).sum())
        numSigI1 = int((selected.layers['fdr'].flatten() < fdr).sum())
        data.append({
            'cell_type': celltype,
            'target': target,
            'numDEG': numDEG,
            'numSigI1': numSigI1,
            'MAST_log2FC': logfc,
            'MAST_fdr': logfc_fdr
        })

    return pd.DataFrame(data)


def _combine_obs_for_interactions(interactions: pd.DataFrame,
                                  exp_adata: ad.AnnData,
                                  cin_adata: ad.AnnData,
                                  sting_adata: ad.AnnData,
                                  condition_name: str,
                                  fdr: str):
    if condition_name == 'highCIN_vs_noSTING':
        condition_name = 'max'

    cin_agg = None
    sting_agg = None
    main_agg = None
    main_adata = None
    if condition_name == 'highCIN_vs_noSTING' or condition_name == 'max':
        sting_agg = calculate_aggregated_stats(sting_adata, fdr)
        main_agg = sting_agg
        main_adata = sting_adata

    if condition_name == 'highCIN_vs_lowCIN' or condition_name == 'max':
        cin_agg = calculate_aggregated_stats(cin_adata, fdr)
        main_agg = cin_agg
        main_adata = cin_adata

    assert main_adata is not None

    if condition_name == 'max':
        main_agg = pd.merge(cin_agg, sting_agg, on=['cell_type', 'target'], how='inner', suffixes=('_cin', '_sting'))
        main_agg['numSigI1'] = np.where(main_agg['numSigI1_cin'] > main_agg['numSigI1_sting'],
                                        main_agg['numSigI1_cin'], main_agg['numSigI1_sting'])
        main_agg['numDEG'] = np.where(main_agg['numDEG_cin'] > main_agg['numDEG_sting'],
                                        main_agg['numDEG_cin'], main_agg['numDEG_sting'])
        logfc_selector = main_agg['MAST_log2FC_cin'] > main_agg['MAST_log2FC_sting']
        main_agg['MAST_log2FC'] = np.where(logfc_selector, main_agg['MAST_log2FC_cin'], main_agg['MAST_log2FC_sting'])
        main_agg['MAST_fdr'] = np.where(logfc_selector, main_agg['MAST_fdr_cin'], main_agg['MAST_fdr_sting'])
    # We are only using this for receptors
    main_agg['cell_type_receptor'] = main_agg['cell_type']
    main_agg['receptor'] = main_agg['target']

    all_celltypes = set(main_adata.obs['cell type'])
    full_df = pd.DataFrame()
    for donor_celltype, target_celltype in itertools.product(all_celltypes, repeat=2):
        lig_adata = main_adata[(main_adata.obs['cell type'] == donor_celltype) &
                               (main_adata.obs['target'].isin(interactions.ligand))]
        rec_adata = main_adata[(main_adata.obs['cell type'] == target_celltype) &
                               (main_adata.obs['target'].isin(interactions.receptor))]

        ligand_df = pd.DataFrame({
            'ligand': lig_adata.obs.target,
            'MAST_log2FC_ligand': lig_adata.obs[f'MAST_log2FC_highCIN_vs_lowCIN'],
            'MAST_fdr_ligand': lig_adata.obs["MAST_fdr_highCIN_vs_lowCIN"],
            'cell_type_ligand_dc1': lig_adata.obs['cell_type_dc1'],
            'DA_ligand': lig_adata.obs['DA_score'],
            'expression_ligand': calculate_expressions(lig_adata, exp_adata),
        })
        ligand_df['cell_type_ligand'] = donor_celltype

        agg_df = main_agg[main_agg['cell_type_receptor'] == target_celltype]

        receptor_df = pd.DataFrame({
            'receptor': rec_adata.obs.target,
            'MAST_log2FC_receptor': rec_adata.obs[f'MAST_log2FC_highCIN_vs_lowCIN'],
            'MAST_fdr_receptor': rec_adata.obs["MAST_fdr_highCIN_vs_lowCIN"],
            'cell_type_receptor_dc1': rec_adata.obs['cell_type_dc1'],
            'DA_receptor': rec_adata.obs['DA_score'],
            'expression_receptor': calculate_expressions(rec_adata, exp_adata)
        })
        receptor_df['cell_type_receptor'] = target_celltype
        receptor_df = pd.merge(receptor_df, agg_df, on=['cell_type_receptor', 'receptor'], how='inner', suffixes=('', '_old'))

        interactions_ct = pd.merge(interactions, ligand_df, on='ligand')
        interactions_ct = pd.merge(interactions_ct, receptor_df, on='receptor')

        full_df = pd.concat([full_df, interactions_ct])

    return full_df


def _compile_interactions(interactions: pd.DataFrame,
                          exp_adata: ad.AnnData,
                          cin_adata: ad.AnnData,
                          sting_adata: ad.AnnData,
                          condition_name: str,
                          fdr: str):
    prefix = 'data/compiled/interactions'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, condition_name, fdr, 'csv')

    if os.path.exists(file):
        return

    interactions_df = _combine_obs_for_interactions(interactions, exp_adata, cin_adata, sting_adata, condition_name, fdr)
    interactions_df.to_csv(file, index=False)


def _compile_ligand_effects(interactions: pd.DataFrame,
                            exp_adata: ad.AnnData,
                            adata: ad.AnnData,
                            condition_name: str,
                            fdr: str):
    prefix = 'data/compiled/ligand_effects'
    os.makedirs(prefix, exist_ok=True)
    prefix += '/' + condition_name + "_"

    if condition_name == 'highCIN_vs_noSTING':
        raise ValueError("highCIN_vs_noSTING is not supported")

    all_celltypes = set(adata.obs['cell type'])
    genes = [g for g in adata.var.index]
    for receptor in set(interactions.receptor):
        filename = prefix + receptor + ".csv"
        if os.path.exists(filename):
            continue
        df = pd.DataFrame()
        for celltype in all_celltypes:
            rec_adata = adata[(adata.obs['cell type'] == celltype) &
                              (adata.obs['target'] == receptor)]
            if rec_adata.n_obs == 0:
                continue
            rec_expression = calculate_expression(celltype, receptor, exp_adata)
            cell_df = pd.DataFrame({
                'gene': genes,
                'gene_is_ligand': rec_adata.var.index.isin(interactions.ligand),
                'gene_is_receptor': rec_adata.var.index.isin(interactions.receptor),
                'log2FC': rec_adata.layers['log2FC'].flatten(),
                'fdr': rec_adata.layers['fdr'].flatten(),
                'i1.fdr': rec_adata.layers['fdr.i1'].flatten()
            })
            cell_df['cell_type'] = celltype
            cell_df['target'] = receptor
            cell_df['target_expression'] = rec_expression
            cell_df['receptor'] = True
            cell_df['ligand'] = receptor in interactions.ligand
            cell_df['MAST_log2FC'] = rec_adata.obs[f'MAST_log2FC_highCIN_vs_lowCIN'].values[0]
            cell_df['MAST_fdr'] = rec_adata.obs["MAST_fdr_highCIN_vs_lowCIN"].values[0]
            df = pd.concat([df, cell_df])
        if df.shape[0] == 0:
            continue
        else:
            df.to_csv(filename, index=False)

    for ligand in set(interactions.ligand):
        filename = prefix + ligand + ".csv"
        if os.path.exists(filename):
            continue
        df = pd.DataFrame()
        for celltype in all_celltypes:
            lig_adata = adata[(adata.obs['cell type'] == celltype) &
                              (adata.obs['target'] == ligand)]
            if lig_adata.n_obs == 0:
                continue
            lig_expression = calculate_expression(celltype, ligand, exp_adata)
            cell_df = pd.DataFrame({
                'gene': genes,
                'gene_is_ligand': lig_adata.var.index.isin(interactions.ligand),
                'gene_is_receptor': lig_adata.var.index.isin(interactions.receptor),
                'log2FC': lig_adata.layers['log2FC'].flatten(),
                'fdr': lig_adata.layers['fdr'].flatten(),
                'i1.fdr': lig_adata.layers['fdr.i1'].flatten()
            })
            cell_df['cell_type'] = celltype
            cell_df['target'] = ligand
            cell_df['target_expression'] = lig_expression
            cell_df['receptor'] = ligand in interactions.receptor
            cell_df['ligand'] = True
            cell_df['MAST_log2FC'] = lig_adata.obs[f'MAST_log2FC_highCIN_vs_lowCIN'].values[0]
            cell_df['MAST_fdr'] = lig_adata.obs["MAST_fdr_highCIN_vs_lowCIN"].values[0]
            df = pd.concat([df, cell_df])
        if df.shape[0] == 0:
            continue
        else:
            df.to_csv(filename, index=False)


def _compile_circos(interactions: pd.DataFrame,
                    exp_adata: ad.AnnData,
                    cin_adata: ad.AnnData,
                    sting_adata: ad.AnnData,
                    condition_name: str,
                    fdr: str):
    prefix = 'data/compiled/circos'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, condition_name, fdr, 'csv')

    if os.path.exists(file):
        return

    if condition_name == 'highCIN_vs_noSTING':
        condition_name = "max"

    cin_agg = None
    sting_agg = None
    main_agg = None
    main_adata = None
    if condition_name == 'highCIN_vs_noSTING' or condition_name == 'max':
        sting_agg = calculate_aggregated_stats(sting_adata, fdr)
        main_agg = sting_agg
        main_adata = sting_adata

    if condition_name == 'highCIN_vs_lowCIN' or condition_name == 'max':
        cin_agg = calculate_aggregated_stats(cin_adata, fdr)
        main_agg = cin_agg
        main_adata = cin_adata

    assert main_adata is not None

    if condition_name == 'max':
        main_agg = pd.merge(cin_agg, sting_agg, on=['cell_type', 'target'], how='inner', suffixes=('_cin', '_sting'))
        main_agg['numSigI1'] = np.where(main_agg['numSigI1_cin'] > main_agg['numSigI1_sting'],
                                        main_agg['numSigI1_cin'], main_agg['numSigI1_sting'])
        main_agg['numDEG'] = np.where(main_agg['numDEG_cin'] > main_agg['numDEG_sting'],
                                        main_agg['numDEG_cin'], main_agg['numDEG_sting'])
        logfc_selector = main_agg['MAST_log2FC_cin'] > main_agg['MAST_log2FC_sting']
        main_agg['MAST_log2FC'] = np.where(logfc_selector, main_agg['MAST_log2FC_cin'], main_agg['MAST_log2FC_sting'])
        main_agg['MAST_fdr'] = np.where(logfc_selector, main_agg['MAST_fdr_cin'], main_agg['MAST_fdr_sting'])

    df = pd.DataFrame({
        'cell_type': main_adata.obs['cell type'],
        'target': main_adata.obs['target'],
        'ligand': main_adata.obs['ligand'],
        'receptor': main_adata.obs['receptor'],
        'cell_type_dc1': main_adata.obs['cell_type_dc1'],
        'DA_score': main_adata.obs['DA_score'],
    })

    df = pd.merge(df, main_agg, on=['cell_type', 'target'], how='right', suffixes=('_old', ''))

    df.to_csv(file, index=False)
    return df


def read_ct_data(path: str) -> ad.AnnData:
    ct = sc.read(path, backed='r')  # Low memory mode
    return ct


def compile_data(
        interactions_file: str = 'data/allgenes/all_interactions.tsv',
        expression_adata_file: str = 'data/highCIN_vs_lowCIN/saves/adata.h5ad',
        cin_adata_file: str = 'data/allgenes/saves/deg_newfdr.h5ad',
        cin_sting_adata_file: str = 'data/CIN_and_STING/CIN_and_STING_deg.h5ad'
):
    """
    Optimize the builtin files for fast loading
    """

    if os.path.exists('data/compiled/is_compiled'):
        return

    interactions = read_ligand_receptor_file(interactions_file)
    cin_adata = read_ct_data(cin_adata_file)
    cin_sting_adata = read_ct_data(cin_sting_adata_file)
    exp_adata = read_ct_data(expression_adata_file)

    cin_condition = 'highCIN_vs_lowCIN'
    cin_sting_condition = 'highCIN_vs_noSTING'

    # Compile a bunch of FDRs
    for fdr in range(1, 26):
        fdr = str(fdr).zfill(2)
        print("Compiling FDR", fdr)
        _compile_interactions(interactions, exp_adata, cin_adata, cin_sting_adata, cin_condition, fdr)
        _compile_interactions(interactions, exp_adata, cin_adata, cin_sting_adata, cin_sting_condition, fdr)

        _compile_circos(interactions, exp_adata, cin_adata, cin_sting_adata, cin_condition, fdr)
        _compile_circos(interactions, exp_adata, cin_adata, cin_sting_adata, cin_sting_condition, fdr)

        _compile_ligand_effects(interactions, exp_adata, cin_adata, cin_condition, fdr)
        #_compile_ligand_effects(interactions, exp_adata, cin_sting_adata, cin_sting_condition, fdr)

    # touch file to indicate that we've compiled the data
    with open('data/compiled/is_compiled', 'w') as f:
        f.write('compiled')


if __name__ == '__main__':
    compile_data(
        interactions_file="/workdir/mt269/ContactTracing/CIN_TME/allgenes/all_interactions.tsv",
        expression_adata_file="/workdir/mt269/ContactTracing/CIN_TME/highCIN_vs_lowCIN/saves/adata.h5ad",
        cin_adata_file="/workdir/mt269/ContactTracing/CIN_TME/allgenes/saves/deg_newfdr.h5ad",
        cin_sting_adata_file="/workdir/mt269/ContactTracing/CIN_TME/highCIN_vs_noSTING/deg_fixed.h5ad"
    )
