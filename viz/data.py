import functools
import itertools
import os
from collections import namedtuple
from typing import Dict, List

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd


# Exclude these cell types
EXCLUDED_CELL_TYPES = ['Osteoclasts']


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
        return condition  # Assume it's already normalized


def _normalize_fdr_name(fdr: str) -> str:
    return fdr.split('.')[-1].replace("fdr", "")


def read_circos_file(condition: str, fdr: str) -> pd.DataFrame:
    prefix = 'data/compiled/circos'
    df = pd.read_feather(_filename(prefix, _normalize_condition_name(condition), _normalize_fdr_name(fdr), 'feather'))
    df = df[~df['cell_type'].isin(EXCLUDED_CELL_TYPES)]
    return df


def read_ligand_effect_for(condition: str, target: str) -> pd.DataFrame:
    prefix = 'data/compiled/ligand_effects/'
    condition = _normalize_condition_name(condition)
    prefix += condition + '_' + target + '.feather'
    if not os.path.exists(prefix):
        return None
    df = pd.read_feather(prefix)
    df = df[~df['cell_type'].isin(EXCLUDED_CELL_TYPES)]
    return df


def read_interactions_file(condition: str, fdr: str) -> pd.DataFrame:
    prefix = 'data/compiled/interactions'
    df = pd.read_feather(_filename(prefix, _normalize_condition_name(condition), _normalize_fdr_name(fdr), 'feather'))
    df = df[(~df['cell_type_receptor'].isin(EXCLUDED_CELL_TYPES)) & (~df['cell_type_ligand'].isin(EXCLUDED_CELL_TYPES))]
    return df


def read_ligand_receptor_file(filename: str = 'data/compiled/interactions.feather') -> pd.DataFrame:
    func = pd.read_feather if filename.endswith('.feather') else functools.partial(pd.read_csv, sep='\t')
    df = func(filename)[['ligand', 'receptor']].drop_duplicates()
    return df


def calculate_expression(celltype, target, exp_adata: ad.AnnData):
    filter = (exp_adata.obs['Cell Type'] == celltype)
    expressing_cells = (exp_adata[filter, target].layers['X'].flatten() > 0).sum()
    return np.nan_to_num((expressing_cells / filter.sum()).item())


def calculate_expressions(adata: ad.AnnData, exp_adata: ad.AnnData):
    expressions = []
    for i, row in adata.obs.iterrows():
        celltype = row['cell type']
        target = row['target']
        expressions.append(calculate_expression(celltype, target, exp_adata))
    return expressions


def calculate_aggregated_stats(adata: ad.AnnData,
                               fdr: str) -> pd.DataFrame:
    fdr = float("0." + fdr)

    data = list()
    for i, row in adata.obs.iterrows():
        celltype = row['cell type']
        target = row['target']
        logfc = row['MAST_log2FC_highCIN_vs_lowCIN']
        logfc_fdr = row['MAST_fdr_highCIN_vs_lowCIN']
        selected = adata[(adata.obs['cell type'] == celltype) & (adata.obs['target'] == target)]
        numDEG = int((selected.layers['fdr'].flatten() < fdr).sum())
        numSigI1 = int((selected.layers['fdr.i1'].flatten() < fdr).sum())
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
        logfc_selector = main_agg['MAST_log2FC_cin'].abs() > main_agg['MAST_log2FC_sting'].abs()
        main_agg['MAST_log2FC'] = np.where(logfc_selector, main_agg['MAST_log2FC_cin'], main_agg['MAST_log2FC_sting'])
        main_agg['MAST_fdr'] = np.where(logfc_selector, main_agg['MAST_fdr_cin'], main_agg['MAST_fdr_sting'])

        # Adjust values if directions are different
        main_agg['numSigI1'] = np.where((main_agg['numSigI1_cin'] == 0) | (main_agg['numSigI1_sting'] == 0), 0, main_agg['numSigI1'])
        main_agg['numDEG'] = np.where((main_agg['numDEG_cin'] == 0) | (main_agg['numDEG_sting'] == 0), 0, main_agg['numDEG'])
        main_agg['MAST_log2FC'] = np.where(np.sign(main_agg['MAST_log2FC_cin']) != np.sign(main_agg['MAST_log2FC_sting']), 0, main_agg['MAST_log2FC'])
        main_agg['MAST_fdr'] = np.where(np.sign(main_agg['MAST_log2FC_cin']) != np.sign(main_agg['MAST_log2FC_sting']), 1, main_agg['MAST_fdr'])

    all_celltypes = set(main_adata.obs['cell type'])
    full_df = pd.DataFrame()
    for donor_celltype, target_celltype in itertools.product(all_celltypes, repeat=2):
        lig_adata = main_adata[(main_adata.obs['cell type'] == donor_celltype) &
                               (main_adata.obs['target'].isin(interactions.ligand))]
        rec_adata = main_adata[(main_adata.obs['cell type'] == target_celltype) &
                               (main_adata.obs['target'].isin(interactions.receptor))]

        ligand_df = pd.DataFrame({
            'ligand': lig_adata.obs.target,
            'cell_type_ligand_dc1': lig_adata.obs['cell_type_dc1'],
            'DA_ligand': lig_adata.obs['DA_score'],
            'expression_ligand': calculate_expressions(lig_adata, exp_adata),
        })
        ligand_df['cell_type_ligand'] = donor_celltype
        agg_df = main_agg[main_agg['cell_type'] == donor_celltype]
        # Don't need cell type in the df anymore
        agg_df = agg_df.drop(columns=['cell_type', 'numSigI1', 'numDEG'])
        ligand_df = pd.merge(ligand_df, agg_df, left_on='ligand', right_on='target', how='inner', suffixes=('', '_ignore'))

        receptor_df = pd.DataFrame({
            'receptor': rec_adata.obs.target,
            'cell_type_receptor_dc1': rec_adata.obs['cell_type_dc1'],
            'DA_receptor': rec_adata.obs['DA_score'],
            'expression_receptor': calculate_expressions(rec_adata, exp_adata)
        })
        receptor_df['cell_type_receptor'] = target_celltype
        agg_df = main_agg[main_agg['cell_type'] == target_celltype]
        # Don't need cell type in the df anymore
        agg_df = agg_df.drop(columns=['cell_type'])
        receptor_df = pd.merge(receptor_df, agg_df, left_on='receptor', right_on='target', how='inner', suffixes=('', '_ignore'))

        interactions_ct = pd.merge(interactions, ligand_df, on='ligand', suffixes=('', '_ligand'))
        interactions_ct = pd.merge(interactions_ct, receptor_df, on='receptor', suffixes=('_ligand', '_receptor'))

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
    file = _filename(prefix, condition_name, fdr, 'feather')

    if os.path.exists(file):
        return

    interactions_df = _combine_obs_for_interactions(interactions, exp_adata, cin_adata, sting_adata, condition_name, fdr)
    interactions_df.reset_index().to_feather(file, compression='zstd')


def _custom_aggregated_stats(adata: ad.AnnData,
                             deg_adata: ad.AnnData,
                             target_stats: pd.DataFrame,
                             fdr: str) -> pd.DataFrame:
    fdr = float("0." + fdr)
    data = list()
    for i, row in adata.obs.iterrows():
        celltype = row['cell type']
        target = row['target']
        expression = row['fracExp']
        selected_info = target_stats[(target_stats['cell type'] == celltype) & (target_stats['target'] == target)]
        if selected_info.shape[0] == 0:
            logfc = 0
            logfc_fdr = 1
        else:
            logfc = selected_info['log2FC'].values[0]
            logfc_fdr = selected_info['fdr'].values[0]
        selected_deg = deg_adata[(deg_adata.obs['cell type'] == celltype) & (deg_adata.obs['receptor'] == target)]
        if selected_deg.shape[0] == 0:
            numSigI1 = 0
        else:
            numSigI1 = int((selected_deg.layers['fdr.i1'].flatten() < fdr).sum())
        selected_adata = adata[(adata.obs['cell type'] == celltype) & (adata.obs['target'] == target)]
        if selected_adata.shape[0] == 0:
            numDEG = 0
        else:
            numDEG = int((selected_adata.layers['fdr'].flatten() < fdr).sum())
        data.append({
            'cell_type': celltype,
            'target': target,
            'expression': expression,
            'numDEG': numDEG,
            'numSigI1': numSigI1,
            'MAST_log2FC': logfc,
            'MAST_fdr': logfc_fdr
        })

    return pd.DataFrame(data)


def _compile_custom_interactions(adata: ad.AnnData,
                                 target_stats: pd.DataFrame,
                                 interactions: pd.DataFrame,
                                 agg_stats: pd.DataFrame,
                                 condition_name: str,
                                 fdr: str):
    prefix = 'data/compiled/interactions'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, 'custom', fdr, 'feather')

    all_celltypes = set(adata.obs['cell type'])
    full_df = pd.DataFrame()
    for donor_celltype, target_celltype in itertools.product(all_celltypes, repeat=2):
        lig_data = target_stats[(target_stats['cell type'] == donor_celltype) &
                                (target_stats['target'].isin(interactions.ligand))]
        # Only include select columns
        lig_data = lig_data[['target', 'cell type', 'fracExp', 'cell_type_dc1']]
        lig_data.rename(columns={'target': 'ligand',
                                 'cell type': 'cell_type_ligand',
                                 'fracExp': 'expression_ligand'}, inplace=True)
        rec_data = target_stats[(target_stats['cell type'] == target_celltype) &
                                (target_stats['target'].isin(interactions.receptor))]
        # Only include select columns
        rec_data = rec_data[['target', 'cell type', 'fracExp', 'fdr']]
        rec_data.rename(columns={'target': 'receptor',
                                 'cell type': 'cell_type_receptor',
                                 'fracExp': 'expression_receptor'}, inplace=True)
        agg_df = agg_stats[agg_stats['cell_type'] == donor_celltype]
        # Don't need cell type in the df anymore
        agg_df = agg_df.drop(columns=['cell_type', 'numSigI1', 'numDEG', 'expression'])
        lig_data = pd.merge(lig_data, agg_df, left_on='ligand', right_on='target', how='inner', suffixes=('', '_ignore'))
        agg_df = agg_stats[agg_stats['cell_type'] == target_celltype]
        # Don't need cell type in the df anymore
        agg_df = agg_df.drop(columns=['cell_type'])
        rec_data = pd.merge(rec_data, agg_df, left_on='receptor', right_on='target', how='inner', suffixes=('', '_ignore'))

        interactions_ct = pd.merge(interactions, lig_data, on='ligand', suffixes=('', '_ligand'))
        interactions_ct = pd.merge(interactions_ct, rec_data, on='receptor', suffixes=('_ligand', '_receptor'))

        if interactions_ct.shape[0] == 0:
            continue
        full_df = pd.concat([full_df, interactions_ct])
    full_df['DA_receptor'] = 0
    full_df['DA_ligand'] = 0
    full_df.reset_index().to_feather(file, compression='zstd')


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
    genes = {g for g in adata.var.index} | {g for g in interactions.receptor} | {g for g in interactions.ligand}

    ct_gene_expression = {'max': dict()}
    for celltype in all_celltypes:
        ct_gene_expression[celltype] = dict()
        for gene in genes:
            ct_gene_expression[celltype][gene] = calculate_expression(celltype, gene, exp_adata)
    # Max Expressions
    for gene in genes:
        ct_gene_expression['max'][gene] = max([ct_gene_expression[celltype].get(gene, 0) for celltype in all_celltypes])

    for receptor in set(interactions.receptor):
        filename = prefix + receptor + ".feather"
        if os.path.exists(filename):
            continue
        df = pd.DataFrame()
        for celltype in all_celltypes:
            rec_adata = adata[(adata.obs['cell type'] == celltype) &
                              (adata.obs['target'] == receptor)]
            if rec_adata.n_obs == 0:
                continue
            rec_expression = ct_gene_expression[celltype].get(receptor, 0)

            cell_df = pd.DataFrame({
                'gene': adata.var.index,
                'gene_is_ligand': rec_adata.var.index.isin(interactions.ligand),
                'gene_is_receptor': rec_adata.var.index.isin(interactions.receptor),
                'gene_expression': [ct_gene_expression[celltype].get(g, 0) for g in adata.var.index],
                'log2FC': rec_adata.layers['log2FC'].flatten(),
                'fdr': rec_adata.layers['fdr'].flatten(),
                'i1.fdr': rec_adata.layers['fdr.i1'].flatten()
            })

            # The max expression of genes in the cell type's microenvironment
            # If the gene is a receptor, this represents the max expression of its ligands from any cell type
            # If the gene is not a receptor, it is 0
            gene_max_ligand_expression = []
            for gene, is_receptor in zip(cell_df['gene'], cell_df['gene_is_receptor']):
                if is_receptor:
                    ligands = interactions[interactions.receptor == gene].ligand.tolist()
                    gene_max_ligand_expression.append(max([ct_gene_expression['max'].get(l, 0) for l in ligands]))
                else:
                    gene_max_ligand_expression.append(0)
            cell_df['gene_max_ligand_expression'] = gene_max_ligand_expression

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
            df.reset_index().to_feather(filename, compression='zstd')

    for ligand in set(interactions.ligand):
        filename = prefix + ligand + ".feather"
        if os.path.exists(filename):
            continue
        df = pd.DataFrame()
        for celltype in all_celltypes:
            lig_adata = adata[(adata.obs['cell type'] == celltype) &
                              (adata.obs['target'] == ligand)]
            if lig_adata.n_obs == 0:
                continue
            lig_expression = ct_gene_expression[celltype].get(ligand, 0)
            cell_df = pd.DataFrame({
                'gene': adata.var.index,
                'gene_is_ligand': lig_adata.var.index.isin(interactions.ligand),
                'gene_is_receptor': lig_adata.var.index.isin(interactions.receptor),
                'gene_expression': [ct_gene_expression[celltype].get(g, 0) for g in adata.var.index],
                'log2FC': lig_adata.layers['log2FC'].flatten(),
                'fdr': lig_adata.layers['fdr'].flatten(),
                'i1.fdr': lig_adata.layers['fdr.i1'].flatten()
            })

            # The max expression of genes in the cell type's microenvironment
            # If the gene is a receptor, this represents the max expression of its ligands from any cell type
            # If the gene is not a receptor, it is 0
            gene_max_ligand_expression = []
            for gene, is_receptor in zip(cell_df['gene'], cell_df['gene_is_receptor']):
                if is_receptor:
                    ligands = interactions[interactions.receptor == gene].ligand.tolist()
                    gene_max_ligand_expression.append(max([ct_gene_expression['max'].get(l, 0) for l in ligands]))
                else:
                    gene_max_ligand_expression.append(0)
            cell_df['gene_max_ligand_expression'] = gene_max_ligand_expression

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
            df.reset_index().to_feather(filename, compression='zstd')


def _compile_custom_ligand_effects(interactions: pd.DataFrame,
                                   adata: ad.AnnData,
                                   deg_adata: ad.AnnData,
                                   target_stats: pd.DataFrame,
                                   condition_name: str):
    prefix = 'data/compiled/ligand_effects'
    os.makedirs(prefix, exist_ok=True)
    prefix += '/' + 'custom' + "_"

    all_celltypes = set(adata.obs['cell type'])

    for receptor in set(interactions.receptor):
        filename = prefix + receptor + ".feather"
        df = pd.DataFrame()
        for celltype in all_celltypes:
            rec_adata = adata[(adata.obs['cell type'] == celltype) &
                                (adata.obs['target'] == receptor)]
            rec_deg_adata = deg_adata[(deg_adata.obs['cell type'] == celltype) &
                                        (deg_adata.obs['receptor'] == receptor), rec_adata.var.index.values]
            if rec_adata.n_obs == 0 or rec_adata.n_vars == 0:
                continue

            cell_df = pd.DataFrame({
                'gene': rec_adata.var.index.values,
                'log2FC': rec_adata.layers['lfc'].flatten(),
                'fdr': rec_adata.layers['fdr'].flatten(),
            })
            cell_df['gene_expression'] = 1
            cell_df['gene_is_ligand'] = cell_df['gene'].isin(interactions.ligand)
            cell_df['gene_is_receptor'] = cell_df['gene'].isin(interactions.receptor)
            cell_df['i1.fdr'] = rec_deg_adata.layers['fdr.i1'].flatten()
            cell_df['gene_max_ligand_expression'] = 1
            cell_df['cell_type'] = celltype
            cell_df['target'] = receptor
            cell_df['target_expression'] = 1
            cell_df['receptor'] = True
            cell_df['ligand'] = receptor in interactions.ligand
            filtered_target_stats = target_stats[(target_stats['cell type'] == celltype) &
                                                    (target_stats['target'] == receptor)]
            cell_df['MAST_log2FC'] = filtered_target_stats['log2FC'].values[0] if filtered_target_stats.shape[0] > 0 else 0
            cell_df['MAST_fdr'] = filtered_target_stats['fdr'].values[0] if filtered_target_stats.shape[0] > 0 else 1
            df = pd.concat([df, cell_df])
        if df.shape[0] == 0:
            continue
        else:
            df.reset_index().to_feather(filename, compression='zstd')

    for ligand in set(interactions.ligand):
        if ligand in interactions.receptor:
            continue

        filename = prefix + ligand + ".feather"
        df = pd.DataFrame()
        for celltype in all_celltypes:
            lig_adata = adata[(adata.obs['cell type'] == celltype) &
                              (adata.obs['target'] == ligand)]
            if lig_adata.n_obs == 0 or lig_adata.n_vars == 0:
                continue
            cell_df = pd.DataFrame({
                'gene': lig_adata.var.index.values,
                'log2FC': lig_adata.layers['lfc'].flatten(),
                'fdr': lig_adata.layers['fdr'].flatten(),
            })
            cell_df['gene_expression'] = 1
            cell_df['gene_is_ligand'] = cell_df['gene'].isin(interactions.ligand)
            cell_df['gene_is_receptor'] = cell_df['gene'].isin(interactions.receptor)
            cell_df['i1.fdr'] = 1
            cell_df['gene_max_ligand_expression'] = 1
            cell_df['cell_type'] = celltype
            cell_df['target'] = ligand
            cell_df['target_expression'] = 1
            cell_df['receptor'] = ligand in interactions.receptor
            cell_df['ligand'] = True
            filtered_target_stats = target_stats[(target_stats['cell type'] == celltype) &
                                                    (target_stats['target'] == ligand)]
            cell_df['MAST_log2FC'] = filtered_target_stats['log2FC'].values[0] if filtered_target_stats.shape[0] > 0 else 0
            cell_df['MAST_fdr'] = filtered_target_stats['fdr'].values[0] if filtered_target_stats.shape[0] > 0 else 1
            df = pd.concat([df, cell_df])
        if df.shape[0] == 0:
            continue
        else:
            df.reset_index().to_feather(filename, compression='zstd')


def _compile_circos(interactions: pd.DataFrame,
                    exp_adata: ad.AnnData,
                    cin_adata: ad.AnnData,
                    sting_adata: ad.AnnData,
                    condition_name: str,
                    fdr: str):
    prefix = 'data/compiled/circos'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, condition_name, fdr, 'feather')

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
        logfc_selector = main_agg['MAST_log2FC_cin'].abs() > main_agg['MAST_log2FC_sting'].abs()
        main_agg['MAST_log2FC'] = np.where(logfc_selector, main_agg['MAST_log2FC_cin'], main_agg['MAST_log2FC_sting'])
        main_agg['MAST_fdr'] = np.where(logfc_selector, main_agg['MAST_fdr_cin'], main_agg['MAST_fdr_sting'])

        # Adjust values if directions are different
        main_agg['numSigI1'] = np.where((main_agg['numSigI1_cin'] == 0) | (main_agg['numSigI1_sting'] == 0), 0, main_agg['numSigI1'])
        main_agg['numDEG'] = np.where((main_agg['numDEG_cin'] == 0) | (main_agg['numDEG_sting'] == 0), 0, main_agg['numDEG'])
        main_agg['MAST_log2FC'] = np.where(np.sign(main_agg['MAST_log2FC_cin']) != np.sign(main_agg['MAST_log2FC_sting']), 0, main_agg['MAST_log2FC'])
        main_agg['MAST_fdr'] = np.where(np.sign(main_agg['MAST_log2FC_cin']) != np.sign(main_agg['MAST_log2FC_sting']), 1, main_agg['MAST_fdr'])

    df = pd.DataFrame({
        'cell_type': main_adata.obs['cell type'],
        'target': main_adata.obs['target'],
        'ligand': main_adata.obs['ligand'],
        'receptor': main_adata.obs['receptor'],
        'cell_type_dc1': main_adata.obs['cell_type_dc1'],
        'DA_score': main_adata.obs['DA_score'],
    })

    df = pd.merge(df, main_agg, on=['cell_type', 'target'], how='right', suffixes=('_old', ''))
    # Replace numSigI1 for ligands that are not receptors to 0
    df['numSigI1'] = np.where(~df['receptor'], df['numSigI1'], 0)

    df.reset_index().to_feather(file, compression='zstd')
    return df


def _compile_custom_circos(interactions: pd.DataFrame,
                           target_stats: pd.DataFrame,
                           agg_stats: pd.DataFrame,
                           condition_name: str,
                           fdr: str):
    prefix = 'data/compiled/circos'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, 'custom', fdr, 'feather')

    df = pd.DataFrame({
        'cell_type': target_stats['cell type'],
        'target': target_stats['target'],
        'ligand': target_stats['ligand'],
        'receptor': target_stats['receptor'],
        'cell_type_dc1': target_stats['cell_type_dc1']
    })
    df['DA_score'] = 0

    df = pd.merge(df, agg_stats, on=['cell_type', 'target'], how='right', suffixes=('_old', ''))
    df['receptor'] = df['target'].isin(interactions.receptor)
    df['ligand'] = df['target'].isin(interactions.ligand)
    df['MAST_log2FC'].fillna(0, inplace=True)
    df['MAST_fdr'].fillna(1, inplace=True)
    df['cell_type_dc1'].fillna(0, inplace=True)
    df['DA_score'].fillna(0, inplace=True)
    # Replace numSigI1 for ligands that are not receptors to 0
    df['numSigI1'] = np.where(~df['receptor'], df['numSigI1'], 0)

    df.reset_index().to_feather(file, compression='zstd')
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
    if not os.path.exists('data/compiled/interactions.feather'):
        # Compile to a feather file (no need to compress since this is a pretty small file)
        os.makedirs('data/compiled', exist_ok=True)
        interactions.reset_index().to_feather('data/compiled/interactions.feather')

    cin_adata = read_ct_data(cin_adata_file)
    cin_sting_adata = read_ct_data(cin_sting_adata_file)
    exp_adata = read_ct_data(expression_adata_file)

    cin_condition = 'highCIN_vs_lowCIN'
    cin_sting_condition = 'highCIN_vs_noSTING'

    # Compile a bunch of FDRs
    for fdr in reversed(range(1, 26)):
        fdr = str(fdr).zfill(2)
        print("Compiling FDR", fdr, flush=True)
        _compile_interactions(interactions, exp_adata, cin_adata, cin_sting_adata, cin_condition, fdr)
        _compile_interactions(interactions, exp_adata, cin_adata, cin_sting_adata, cin_sting_condition, fdr)

        _compile_circos(interactions, exp_adata, cin_adata, cin_sting_adata, cin_condition, fdr)
        _compile_circos(interactions, exp_adata, cin_adata, cin_sting_adata, cin_sting_condition, fdr)

        _compile_ligand_effects(interactions, exp_adata, cin_adata, cin_condition, fdr)
        #_compile_ligand_effects(interactions, exp_adata, cin_sting_adata, cin_sting_condition, fdr)

    # touch file to indicate that we've compiled the data
    with open('data/compiled/is_compiled', 'w') as f:
        f.write('compiled')


def compile_custom_dataset(adata: ad.AnnData,
                           deg_adata: ad.AnnData,
                           target_stats: pd.DataFrame,
                           interactions: pd.DataFrame,
                           condition_name: str):
    """
    Compile a custom dataset. Intended for importing non-CIN-TME data.
    """
    # Compile to a feather file (no need to compress since this is a pretty small file)
    os.makedirs('data/compiled', exist_ok=True)
    interactions.reset_index().to_feather('data/compiled/interactions.feather')

    # Compile a bunch of FDRs
    for fdr in reversed(range(1, 26)):
        fdr = str(fdr).zfill(2)
        print("Compiling FDR", fdr, flush=True)

        agg_stats = _custom_aggregated_stats(adata, deg_adata, target_stats, fdr)

        _compile_custom_interactions(adata, target_stats, interactions, agg_stats, condition_name, fdr)

        _compile_custom_circos(interactions, target_stats, agg_stats, condition_name, fdr)

    _compile_custom_ligand_effects(interactions, adata, deg_adata, target_stats, condition_name)

    # touch file to indicate that we've compiled the data
    with open('data/compiled/is_compiled', 'w') as f:
        f.write('compiled')
    with open('data/compiled/custom_data_cell_type_list', 'w') as f:
        f.write('\n'.join(set(adata.obs['cell type'])))
    with open('data/compiled/is_custom', 'w') as f:
        f.write('custom')


def using_custom_data() -> bool:
    return os.path.exists('data/compiled/is_custom')


def get_custom_celltypes() -> List[str]:
    with open('data/compiled/custom_data_cell_type_list', 'r') as f:
        return f.read().splitlines()


if __name__ == '__main__':
    compile_data(
        interactions_file="/workdir/mt269/ContactTracing/CIN_TME/allgenes/all_interactions.tsv",
        expression_adata_file="/workdir/mt269/ContactTracing/CIN_TME/highCIN_vs_lowCIN/saves/adata.h5ad",
        cin_adata_file="/workdir/mt269/ContactTracing/CIN_TME/allgenes/saves/deg_newfdr.h5ad",
        cin_sting_adata_file="/workdir/mt269/ContactTracing/CIN_TME/highCIN_vs_noSTING/deg_fixed.h5ad"
    )
