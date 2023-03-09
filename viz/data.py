import itertools
import os

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd


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


def _combine_obs_for_interactions(interactions: pd.DataFrame,
                                  exp_adata: ad.AnnData,
                                  adata: ad.AnnData,
                                  condition_name: str,
                                  fdr: str):
    fdr = f"fdr{fdr}"

    if condition_name == 'highCIN_vs_noSTING':
        condition_name = 'max'

    all_celltypes = set(adata.obs['cell type'])
    full_df = pd.DataFrame()
    for donor_celltype, target_celltype in itertools.product(all_celltypes, repeat=2):
        lig_adata = adata[(adata.obs['cell type'] == donor_celltype) &
                              (adata.obs['target'].isin(interactions.ligand))]
        rec_adata = adata[(adata.obs['cell type'] == target_celltype) &
                          (adata.obs['target'].isin(interactions.receptor))]

        ligand_df = pd.DataFrame({
            'ligand': lig_adata.obs.target,
            'MAST_log2FC_ligand': lig_adata.obs[f'MAST_log2FC_{condition_name}'],
            'MAST_fdr_ligand': lig_adata.obs["MAST_fdr" if condition_name == 'max' else f'MAST_fdr_{condition_name}'],
            'cell_type_ligand_dc1': lig_adata.obs['cell_type_dc1'],
            'DA_ligand': lig_adata.obs['DA_score'],
            'expression_ligand': calculate_expressions(lig_adata, exp_adata),
        })
        ligand_df['cell_type_ligand'] = donor_celltype

        if condition_name == 'max':
            numdeg = np.where(
                rec_adata.obs[f'numDEG_{fdr}_highCIN_vs_lowCIN'] > rec_adata.obs[f'numDEG_{fdr}_highCIN_vs_noSTING'],
                rec_adata.obs[f'numDEG_{fdr}_highCIN_vs_lowCIN'], rec_adata.obs[f'numDEG_{fdr}_highCIN_vs_noSTING'])
        else:
            numdeg = rec_adata.obs[f'numDEG_{fdr}']

        receptor_df = pd.DataFrame({
            'receptor': rec_adata.obs.target,
            'MAST_log2FC_receptor': rec_adata.obs[f'MAST_log2FC_{condition_name}'],
            'MAST_fdr_receptor': rec_adata.obs["MAST_fdr" if condition_name == 'max' else f'MAST_fdr_{condition_name}'],
            'cell_type_receptor_dc1': rec_adata.obs['cell_type_dc1'],
            'DA_receptor': rec_adata.obs['DA_score'],
            'expression_receptor': calculate_expressions(rec_adata, exp_adata),
            # Receptor exclusive features
            'numDEG': numdeg,
            'numSigI1': rec_adata.obs[f'numSigI1_{fdr}_max' if condition_name == 'max' else f'numSigI1_{fdr}'],
        })
        receptor_df['cell_type_receptor'] = target_celltype

        interactions_ct = pd.merge(interactions, ligand_df, on='ligand')
        interactions_ct = pd.merge(interactions_ct, receptor_df, on='receptor')

        full_df = pd.concat([full_df, interactions_ct])

    return full_df


def _compile_interactions(interactions: pd.DataFrame,
                          exp_adata: ad.AnnData,
                          adata: ad.AnnData,
                          condition_name: str,
                          fdr: str):
    prefix = 'data/compiled/interactions'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, condition_name, fdr, 'csv')

    if os.path.exists(file):
        return

    interactions_df = _combine_obs_for_interactions(interactions, exp_adata, adata, condition_name, fdr)
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
            cell_df['MAST_log2FC'] = rec_adata.obs[f'MAST_log2FC_{condition_name}'].values[0]
            cell_df['MAST_fdr'] = rec_adata.obs["MAST_fdr" if condition_name == 'max' else f'MAST_fdr_{condition_name}'].values[0]
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
            cell_df['MAST_log2FC'] = lig_adata.obs[f'MAST_log2FC_{condition_name}'].values[0]
            cell_df['MAST_fdr'] = lig_adata.obs["MAST_fdr" if condition_name == 'max' else f'MAST_fdr_{condition_name}'].values[0]
            df = pd.concat([df, cell_df])
        if df.shape[0] == 0:
            continue
        else:
            df.to_csv(filename, index=False)


def _compile_circos(interactions: pd.DataFrame,
                    exp_adata: ad.AnnData,
                    adata: ad.AnnData,
                    condition_name: str,
                    fdr: str):
    prefix = 'data/compiled/circos'
    os.makedirs(prefix, exist_ok=True)
    file = _filename(prefix, condition_name, fdr, 'csv')

    if os.path.exists(file):
        return

    fdr = f"fdr{fdr}"

    if condition_name == 'highCIN_vs_noSTING':
        condition_name = "max"

    if condition_name == 'max':
        numdeg = np.where(adata.obs[f'numDEG_{fdr}_highCIN_vs_lowCIN'] > adata.obs[f'numDEG_{fdr}_highCIN_vs_noSTING'],
                          adata.obs[f'numDEG_{fdr}_highCIN_vs_lowCIN'], adata.obs[f'numDEG_{fdr}_highCIN_vs_noSTING'])
    else:
        numdeg = adata.obs[f'numDEG_{fdr}']

    df = pd.DataFrame({
        'cell_type': adata.obs['cell type'],
        'target': adata.obs['target'],
        'ligand': adata.obs['ligand'],
        'receptor': adata.obs['receptor'],
        'numDEG': numdeg,
        'numSigI1': adata.obs[f'numSigI1_{fdr}_max' if condition_name == 'max' else f'numSigI1_{fdr}'],
        'MAST_log2FC': adata.obs[f'MAST_log2FC_{condition_name}'],
        'MAST_fdr': adata.obs["MAST_fdr" if condition_name == 'max' else f'MAST_fdr_{condition_name}'],
        'cell_type_dc1': adata.obs['cell_type_dc1'],
        'DA_score': adata.obs['DA_score'],
    })
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

    for fdr in ['05', '25']:
        print("Compiling FDR", fdr)
        _compile_interactions(interactions, exp_adata, cin_adata, cin_condition, fdr)
        _compile_interactions(interactions, exp_adata, cin_sting_adata, cin_sting_condition, fdr)

        _compile_circos(interactions, exp_adata, cin_adata, cin_condition, fdr)
        _compile_circos(interactions, exp_adata, cin_sting_adata, cin_sting_condition, fdr)

        _compile_ligand_effects(interactions, exp_adata, cin_adata, cin_condition, fdr)
        #_compile_ligand_effects(interactions, exp_adata, cin_sting_adata, cin_sting_condition, fdr)

    # touch file to indicate that we've compiled the data
    with open('data/compiled/is_compiled', 'w') as f:
        f.write('compiled')


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


if __name__ == '__main__':
    compile_data()
