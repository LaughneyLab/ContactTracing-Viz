import numpy as np
import pandas as pd

from viz.data import read_circos_file, read_interactions_file


def main():
    my_data_outer = read_circos_file('highCIN_vs_noSTING', 'fdr25')
    my_data_inter = read_interactions_file('highCIN_vs_noSTING', 'fdr25')

    min_numsigi1 = 10
    logfc_pvalue_cutoff = 0.05
    min_numdeg = 0
    min_logfc = 0.12

    # Filters from web.py
    my_data_outer = my_data_outer[~my_data_outer['cell_type_dc1'].isna()]

    #outer_data['cell_type_dc1'] = outer_data['cell_type_dc1'].fillna(0)
    # Ligand / Receptor filter
    my_data_outer = my_data_outer[(my_data_outer['receptor'] & (my_data_outer['numSigI1_cin'] > 0) & (my_data_outer['numSigI1_sting'] > 0))
                            | (my_data_outer['ligand'] & (my_data_outer['MAST_fdr_cin'] < logfc_pvalue_cutoff) &
                               (my_data_outer['MAST_fdr_sting'] < logfc_pvalue_cutoff) &
                               (np.sign(my_data_outer['MAST_log2FC_cin']) == np.sign(my_data_outer['MAST_log2FC_sting'])))]

    # Only select interactions present in outer_data
    inter_ligand_index = pd.MultiIndex.from_frame(my_data_inter[['ligand', 'cell_type_ligand']])
    outer_ligand_index = pd.MultiIndex.from_frame(my_data_outer[['target', 'cell_type']])
    my_data_inter = my_data_inter.loc[inter_ligand_index.isin(outer_ligand_index)]
    inter_receptor_index = pd.MultiIndex.from_frame(my_data_inter[['receptor', 'cell_type_receptor']])
    outer_receptor_index = pd.MultiIndex.from_frame(my_data_outer[['target', 'cell_type']])
    my_data_inter = my_data_inter.loc[inter_receptor_index.isin(outer_receptor_index)]
    del inter_ligand_index, outer_ligand_index

    # Filter outer data to just receptors and ligands with interactions
    my_data_outer = my_data_outer[my_data_outer['receptor'] | my_data_outer['target'].isin(my_data_inter['ligand'])]

    my_data_inter = my_data_inter[(my_data_inter['numSigI1'] >= min_numsigi1) & (my_data_inter['numDEG'] >= min_numdeg)
                            & (my_data_inter['MAST_fdr_ligand'] < logfc_pvalue_cutoff) & (my_data_inter['MAST_log2FC_ligand'].abs() > min_logfc)]

    links = pd.read_csv('old_data/minInt10_fdr25_logFC0.12/links_tabular.tsv', sep='\t')
    links['cell_type_receptor'] = links['cell type_receptor']
    links['cell_type_ligand'] = links['cell type_ligand']
    links['numSigI1'] = links['numSigI1_fdr25_max_receptor']

    nodes = pd.read_csv('old_data/minInt10_fdr25_logFC0.12/circle_plot_tabular.tsv', sep='\t')
    nodes['cell_type'] = nodes['cell type']
    nodes['numSigI1'] = nodes['numSigI1_fdr25_max']

    links_data = pd.merge(links, my_data_inter, 'outer', ['cell_type_receptor', 'receptor', 'cell_type_ligand', 'ligand'], suffixes=('_old', '_new'))
    # Require at least one NA value
    #links_data = links_data[links_data.isna().any(axis=1)]
    links_data['old_missing'] = links_data['numSigI1_fdr25_max_receptor'].isna()
    links_data['new_missing'] = links_data['numSigI1_new'].isna()

    nodes_data = pd.merge(nodes, my_data_outer, 'outer', ['cell_type', 'target'], suffixes=('_old', '_new'))
    #nodes_data = nodes_data[nodes_data.isna().any(axis=1)]
    nodes_data['old_missing'] = nodes_data['numSigI1_fdr25_max'].isna()
    nodes_data['new_missing'] = nodes_data['numSigI1_new'].isna()

    links_data.to_csv('old_data/links_intersection.csv', index=False)
    nodes_data.to_csv('old_data/nodes_intersection.csv', index=False)


if __name__ == '__main__':
    main()
