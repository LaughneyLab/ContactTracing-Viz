import pandas as pd

from viz.data import read_interactions, read_ct_data, prune_ct_data
from viz.figures import pseudotime_interaction_propagation_graph, bipartite_graph, polar_receptor_figure, \
    cascading_effects_figure
from viz.util import mouse_colors


def make_pairwise_fig():
    interactions = read_interactions("data/mouse_ranked_interactions.tsv")
    relevant_interactions = pd.read_excel("links_tabular_minInt10_fdr25_logFC0p12_EE.xlsx", sheet_name="tumor-ligands-shaping-TME")
    relevant_interactions['cell_type_receptor'] = relevant_interactions['cell type_receptor']
    relevant_interactions['cell_type_ligand'] = relevant_interactions['cell type_ligand']

    # Select only interactions that are in the relevant interactions
    interactions = interactions.merge(relevant_interactions,
                                      on=['cell_type_ligand', 'cell_type_receptor', "ligand", "receptor"],
                                      how="inner", suffixes=("", "_y"))

    fig = bipartite_graph(
        df=interactions,
        cell1="Tumor cells",
        cell2="Macrophages/mMDSC",
        cell3=None,
        numInteractions=0,
        min_logfc_bipartite=0,
        logfc_fdr_bipartite_cutoff=1
    )

    fig.write_image("pairwise.png")
    fig.write_image("pairwise.svg")
    fig.write_html("pairwise.html")


def make_ligand_spread_fig():
    # TODO: Radial plots for each receptor/cell type, with ligands as spokes
    # Then draw arrows connecting each plot to the next (positive log2FC)
    adata = read_ct_data("data/mouse_deg_fixed_pruned.h5ad")
    interactions = read_interactions("data/mouse_ranked_interactions.tsv")
    relevant_interactions = pd.read_excel("links_tabular_minInt10_fdr25_logFC0p12_EE.xlsx",
                                          sheet_name="tumor-ligands-shaping-TME")  # sheet_name="links_tabular")
    relevant_interactions['cell_type_receptor'] = relevant_interactions['cell type_receptor']
    relevant_interactions['cell_type_ligand'] = relevant_interactions['cell type_ligand']

    # Select only interactions that are in the relevant interactions
    interactions = interactions.merge(relevant_interactions,
                                      on=['cell_type_ligand', 'cell_type_receptor', "ligand", "receptor"],
                                      how="inner", suffixes=("", "_y"))

    fig = pseudotime_interaction_propagation_graph(
        ct=adata,
        orig_df=interactions,
        seed_cell="Tumor cells",
        seed_ligands=["Ccl2", "S100a8", "Serpine2", "S100a9", "Cxcl1", "Apoe", "Col7a1", "Il11", "Timp1"],
        iterations=25,
        interaction_fdr_cutoff=1,
        min_logfc=0,
        logfc_fdr_cutoff=1,
        min_expression=0,
        layout="timeline"
    )

    fig.write_image("ligand_spread.png")
    fig.write_image("ligand_spread.svg")
    fig.write_html("ligand_spread.html")


def test_radial_plot():
    adata = read_ct_data("data/mouse_deg_fixed_pruned.h5ad")
    interactions = read_interactions("data/mouse_ranked_interactions.tsv")
    relevant_interactions = pd.read_excel("links_tabular_minInt10_fdr25_logFC0p12_EE.xlsx",
                                          sheet_name="links_tabular")
    relevant_interactions['cell_type_receptor'] = relevant_interactions['cell type_receptor']
    relevant_interactions['cell_type_ligand'] = relevant_interactions['cell type_ligand']

    # Select only interactions that are in the relevant interactions
    #interactions = interactions.merge(relevant_interactions,
    #                                  on=['cell_type_ligand', 'cell_type_receptor', "ligand", "receptor"],
    #                                  how="inner", suffixes=("", "_y"))
    fig, _ = polar_receptor_figure(adata, interactions, "PMN/gMDSC", "Trem2", min_logfc=0.05, max_fdr=.25, main_node_color=mouse_colors["Macrophages/mMDSC"])
    fig.show()


def radial_networks():
    adata = read_ct_data("data/mouse_deg_fixed_pruned.h5ad")
    interactions = read_interactions("data/mouse_ranked_interactions.tsv")
    relevant_interactions = pd.read_excel("links_tabular_minInt10_fdr25_logFC0p12_EE.xlsx",
                                          sheet_name="links_tabular")
    relevant_interactions['cell_type_receptor'] = relevant_interactions['cell type_receptor']
    relevant_interactions['cell_type_ligand'] = relevant_interactions['cell type_ligand']

    # Select only interactions that are in the relevant interactions
    #interactions = interactions.merge(relevant_interactions,
    #                                 on=['cell_type_ligand', 'cell_type_receptor', "ligand", "receptor"],
    #                                 how="inner", suffixes=("", "_y"))

    # fig = cascading_effects_figure(adata, interactions,
    #                               ["Ccl2", "S100a8", "Serpine2", "S100a9", "Cxcl1", "Apoe", "Col7a1", "Il11", "Timp1"],
    #                               'Tumor cells', min_logfc=0, max_fdr=.25, iterations=1)
    #fig.show()

    cascading_effects_figure(adata, interactions,
                             ["Apoe"], 'Tumor cells',
                             min_logfc=0, max_fdr=.25, numSigI1=10,
                             iterations=2, celltype_filters=['Tumor cells', 'Macrophages/mMDSC']).show()


def main():
    #prune_ct_data("data/mouse_deg_fixed.h5ad")
    #make_pairwise_fig()
    #test_radial_plot()
    #make_ligand_spread_fig()
    radial_networks()


if __name__ == '__main__':
    main()
