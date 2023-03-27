import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

from viz.web import make_tooltip


def conditions_def(text):
    return make_tooltip(text, [
        html.I("ContactTracing"),
        " compares the transcriptional profiles of genes between two conditions. In our study, we compared "
        "CINlow vs CINhigh tumors and CINhigh vs STING-KO tumors. Presented are the results for the CINlow vs CINhigh "
        "comparison and the intersection between the CINlow vs CINhigh and CINhigh vs STING-KO comparisons."
    ])


def ligand_log2fc_def(text):
    return make_tooltip(text, [
        "As part of the tests ",
        html.I("ContactTracing"),
        " performs, differential availability of a ligand between conditions is required to identify condition-specific "
        "effects from the TME."
    ])


def diffusion_component_def(text):
    return make_tooltip(text, "The Diffusion component metric depicted is an embedding of each gene's differential expression, "
                              "it is used for ordering genes along the Circos plot.")


def differential_abundance_def(text):
    return make_tooltip(text, "The differential abundance metric depicted represents the relative correlation of each genes expression "
                              "across each cell type's differential abundance between conditions.")


def deg_test_def(text):
    return make_tooltip(text,
                        [
                            html.I("ContactTracing"),
                            " uses MAST to identify whether a downstream gene's expression is significantly differentially"
                            " expressed according to the activation of a given receptor."
                        ])


def interaction_test_def(text):
    return make_tooltip(text,
                        [
                            html.I("ContactTracing"),
                            " uses MAST to test for whether a given gene's expression is significantly dependent on "
                            "the condition-specific activation of a given receptor."
                        ])


def interaction_effects_def(text):
    return make_tooltip(text,
                        "An interaction effect is a gene which has expression that is effected by "
                        "condition-specific activation of a given receptor.")


def circos_help():
    return [
        html.H5(["The ",
                 html.I('ContactTracing'),
                 " Circos plot aggregates all relevant interaction information into a single figure."]),
        html.Hr(),
        html.P(html.H6(html.Strong("Filter Options"))),
        html.Ul([
            html.Li(html.P([
                "Interaction Effect FDR Cutoff: ",
                html.I("ContactTracing"),
                " identifies the downstream genes that have induced expression shifts from the activation of a receptor"
                " within a given cell of interest; this modulates the threshold for classifying induced expression"
                " change as significant after a Benjamini-Hochberg FDR correction. The default value of 0.25 reflects"
                " what was chosen for evaluating the intersection between CIN and STING-dependent effects. A cutoff"
                " of 0.05 might be more sensible when evaluating CIN-dependent effects on its own."
            ])),
            html.Li(html.P([
                "Ligand Log2FC FDR Cutoff: ",
                html.I("ContactTracing"),
                " requires differential availability of a ligand across conditions to identify condition-specific"
                " effects from the Tumor Microenvironment. This cutoff represents the threshold for identifying whether"
                " a ligand's differential expression between conditions significantly differs from 0 after a"
                " Benjamini-Hochberg FDR correction."
            ])),
            html.Li(html.P([
                "Minimum Interaction Effect: As ",
                html.I("ContactTracing"),
                " ContactTracing identifies the condition-specific effects of each receptor's activation by evaluating "
                "all possible genes, this value indicates a filter for the minimum number of condition-specific "
                "activations of genes by each given receptor using the selected Interaction Effect FDR Cutoff. The "
                "interaction effect, in part, determines whether a receptor has interactions depicted in the "
                "\"ribbons\" of the Circos plot."
            ])),
            html.Li(html.P([
                "Minimum Affected Downstream Genes: The number of differentially expressed genes in response to"
                " receptor activation. This value is similar to the Minimum Interaction Effect parameter, but rather"
                " than requiring condition-specific responses, this measures the strength of the Log2FC induced by a"
                " receptor's activation. Note that the FDR-adjusted p-value threshold is the same as the Minimum"
                " Interaction Effect parameter. "
            ])),
            html.Li(html.P([
                "Chord Minimum Ligand abs(Log2FC): Whereas the Ligand Log2FC FDR Cutoff determines whether genes get "
                "included in the plot, this parameter only affects whether a ribbon is drawn to indicate strong "
                "interaction effects within the Circos figure. Increasing this value leads to fewer ribbons drawn "
                "between ligands and receptors, while the opposite is true if decreased."
            ])),
            html.Li(html.P([
                "Genes of Interest: If specified, the chords drawn corresponding to a ligand or receptor listed will "
                "be emphasized (below is an example). Multiple genes can be included by separating names with commas. ",
                html.Em("IMPORTANT: This field is CASE SENSITIVE.")
            ])),
            html.Li(html.P([
                "Interaction Set: We have evaluated both CIN-dependent and CIN/STING-dependent interaction effects. "
                "This toggle lets users instantly see the plot under both conditions. Note that the CIN/STING "
                "effects represent the maximum value between shared interactions across both conditions and require "
                "interaction effects to have the same directionality."
            ]))
        ]),
        html.P(html.H6(html.Strong("How to Interpret the Circos Plot"))),
        html.P("The final figure is very information-dense, so it can be challenging to interpret initially. "
               "Thus, we will explore the diagram layer by layer. "),
        html.P("First, we have sections of the rings grouped by cell type annotation. NOTE: If a cell type name is too "
               "large, it is excluded. However, it is possible to identify the cell type by hovering over each section "
               "of the plot in the browser."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help1.png'), style={'width': '40%', 'align': 'center'},
                          alt="Cell Types", className='mx-auto'), className='text-center'),
        html.P("The next layer represents the value of the first Diffusion Component, as calculated by Palantir. This "
               "value represents the euclidean space embedding of the differential expression score (the Log2FC "
               "between conditions multiplied by the negative Log10-transformed p-value) of a given gene for a cell "
               "type. The diffusion component allows for a consistent ordering of genes along the Circos rings."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help2.png'), style={'width': '40%', 'align': 'center'},
                          alt="Diffusion Components", className='mx-auto'), className='text-center'),
        html.P("Following Diffusion Component 1, the next ring represents the gene-wise differential abundance of a "
               "gene between conditions. The differential abundance is calculated as the Pearson correlation "
               "coefficient of each gene's imputed expression across cells against the corresponding cell type's "
               "differential abundance across conditions as calculated by MILO. This value is then normalized to "
               "range from -1 to 1."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help3.png'), style={'width': '40%', 'align': 'center'},
                          alt="Differential Abundance", className='mx-auto'), className='text-center'),
        html.P("The final ring represents a histogram scaled to the number of significant interaction effects for each "
               "receptor."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help4.png'), style={'width': '40%', 'align': 'center'},
                          alt="Interaction Effects", className='mx-auto'), className='text-center'),
        html.P("Lastly, all interactions that meet the selected requirements and have the same directionality are "
               "depicted as ribbons connecting each receptor to its associated ligands across cell types. Depending on "
               "the number of ribbons, it may be difficult to read some gene labels. Therefore it is possible to hover "
               "over each ribbon to get precise details. The color of these ribbons reflects a ligand's Log2FC between "
               "conditions, and the thickness of the ribbon is relative to the number of significant interaction "
               "effects for the receptor."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help5.png'), style={'width': '40%', 'align': 'center'},
                          alt="Ribbons", className='mx-auto'), className='text-center'),
        html.P("Additionally, the ribbons can be filtered to emphasize specific ligands or receptors of interest using "
               "the aforementioned \"Genes of Interest\" field."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help6.png'), style={'width': '40%', 'align': 'center'},
                          alt="Ribbons Highlighted", className='mx-auto'), className='text-center')
    ]


def ligand_effects_help():
    return [
        html.H5(["The ligand cascade effects plot visualizes all potential downstream effects from TME-specific"
                 " identified by ",
                 html.I("ContactTracing"),
                 "."]),
        html.Hr(),
        html.P(html.H6(html.Strong("Filter Options"))),
        html.Ul([
            html.Li(html.P([
                "Emitting Cell Type: The cell type for which specified ligands are emitted."
            ])),
            html.Li(html.P([
                "Emitted Ligands: The ligand(s) to follow in the figure; must be emitted from the selected cell type."
                " Multiple genes can be included by separating names with commas. IMPORTANT: This field is CASE"
                " SENSITIVE."
            ])),
            html.Li(html.P([
                "Network Building Iterations: The number of levels of downstream interactions to include. Note that"
                " typically, the number of downstream genes per level can increase exponentially. So starting with a"
                " low value and increasing as needed is recommended since larger values require additional time to"
                " render."
            ])),
            html.Li(html.P([
                "Interaction Effect FDR Cutoff: ",
                html.I("ContactTracing"),
                " identifies the downstream genes that have induced expression shifts from the activation of a "
                " receptor within a given cell of interest; this modulates the threshold for classifying induced"
                " expression change as significant after a Benjamini-Hochberg FDR correction. The default value of 0.25"
                " reflects what was chosen for evaluating the intersection between CIN and STING-dependent effects. A"
                " cutoff of 0.05 might be more sensible when evaluating CIN-dependent effects on its own."
            ])),
            html.Li(html.P([
                "Log2FC FDR Cutoff: ",
                html.I("ContactTracing"),
                " requires differential availability of a ligand across conditions to identify condition-specific"
                " effects from the Tumor Microenvironment. Additionally, as ",
                html.I("ContactTracing"),
                " is able to determine the induced Log2FC of a receptor's activation on a downstream gene, this option"
                " will also be applied to induced Log2FC tests. The cutoff represents the threshold for identifying whether"
                " a gene's differential expression between conditions significantly differs from 0 after a"
                " Benjamini-Hochberg FDR correction."
            ])),
            html.Li(html.P([
                "Minimum Expression: Allows for Ligands/Receptors to be filtered according to a minimum expression"
                " value. This value, however, does not represent counts. Instead, \"expression\" refers to the fraction"
                " of the cells from a given cell type with non-zero counts of a particular gene. Therefore, expression"
                " values range from 0-1, where zero means that no cells of a cell type express the gene, and one means"
                " that all cells of a cell type express the gene."
            ])),
            html.Li(html.P([
                "Minimum abs(Log2FC): While the Log2FC FDR Cutoff option filters interactions according"
                " to whether a gene's differential expression between conditions significantly differs from zero,"
                " this filter additionally allows for further refinement by requiring a minimum absolute Log2FC. The"
                " default value reflects the cutoff used for the Circos plot."
            ])),
            html.Li(html.P([
                "Interaction Set: While we have evaluated both CIN-dependent and CIN/STING-dependent interaction"
                " effects, the CIN/STING intersection data set does not have statistical power for examining"
                " individually induced downstream effects of interactions. Therefore, only CIN-dependent interaction"
                " effects are enabled."
            ]))
        ]),
        html.P(html.H6(html.Strong("How to Interpret the Plot"))),
        html.P("This figure illustrates the potential cascading effects of receptor activation in the tumor"
               " microenvironment. It does so by generating a directed network of gene interactions by following the"
               " steps listed below:"),
        html.Ol([
            html.Li([html.P([
                    "Starting from a given cell type and a set of ligands, identify all the corresponding receptors that"
                    " pass the provided filters and have a significant condition-specific interaction effect on downstream"
                    " signaling, then draw an arrow connecting the ligands to the receptors. Finally, add 1 to the number"
                    " of network building iterations."
                ]),
                html.Div(html.Img(src=dash.get_asset_url("ligand_effects_help1.png"), style={'width': '30%', 'align': 'center'}, alt='Step 1', className='mx-auto'), className='text-center')
            ]),
            html.Li([html.P([
                    "Identify significantly differentially regulated ligands for each added receptor in response to the"
                    " receptor's activation. Additionally, identify further receptors with a significant interaction effect"
                    " and a ligand available in the microenvironment. Finally, add 1 to the number of network building"
                    " iterations."
                ]),
                html.Div(html.Img(src=dash.get_asset_url("ligand_effects_help2.png"), style={'width': '25%', 'align': 'center'}, alt='Step 2', className='mx-auto'), className='text-center')
            ]),
            html.Li([html.P([
                    "Repeat steps 1 and 2 until the number of network building iterations equals what was specified by"
                    " the user."
                ]),
                html.Div(html.Img(src=dash.get_asset_url("ligand_effects_help3.png"), style={'width': '25%', 'align': 'center'}, alt='Step 3', className='mx-auto'), className='text-center')
            ])
        ]),
        html.P(["As depicted above, each of these iterations is represented by a column of ligands and receptors. In"
                " addition, these nodes have arrows connecting to their corresponding receptors and ligands. If the"
                " arrow is connecting from a receptor to a ligand, the color and arrow thickness represents the Log2FC"
                " induced by receptor activation. However, if the arrow is connecting from a ligand to a receptor, the"
                " color is equivalent to a Log2FC of 0 as ligands are not expected to regulate receptor expression"
                " differentially."]),
        html.P(["Additionally, the color of nodes represents the cell type in which each gene is expressed; the shape"
                " of the node represents whether the gene is a ligand, receptor, or both; and the size of each node is"
                " relative to the total induced Log2FC from receptor activation to ligands."])
    ]


def interactions_help():
    return [
        html.H5(["The pairwise interactions plot highlights TME-specific interactions identified by ",
                 html.I("ContactTracing"),
                 " between cells of interest."]),
        html.Hr(),
        html.P(html.H6(html.Strong('Filter Options'))),
        html.Ul([
            html.Li(html.P([
                "Interaction Effect FDR Cutoff: ",
                html.I("ContactTracing"),
                " identifies the downstream genes that have induced expression shifts from the activation of a "
                "receptor within a given cell of interest; this modulates the threshold for classifying induced "
                "expression change as significant after a Benjamini-Hochberg FDR correction. The default value of "
                "0.25 reflects what was chosen for evaluating the intersection between CIN and STING-dependent "
                "effects. A cutoff of 0.05 might be more sensible when evaluating CIN-dependent effects on its own."
            ])),
            html.Li(html.P([
                "Ligand Log2FC FDR Cutoff: ",
                html.I("ContactTracing"),
                " requires differential availability of a ligand across conditions to identify condition-specific "
                "effects from the Tumor Microenvironment. This cutoff represents the threshold for identifying "
                "whether a ligand's differential expression between conditions significantly differs from 0 after a "
                "Benjamini-Hochberg FDR correction."
            ])),
            html.Li(html.P([
                "Minimum Interaction Effect: As ",
                html.I("ContactTracing"),
                " identifies the condition-specific effects of each receptor's activation by evaluating all possible "
                "genes, this value indicates a filter for the minimum number of condition-specific activations of "
                "genes by each given receptor using the selected Interaction Effect FDR Cutoff."
            ])),
            html.Li(html.P([
                "Minimum Expression: Allows for Ligands/Receptors to be filtered according to a minimum expression "
                "value. This value, however, does not represent counts. Instead, \"expression\" refers to the "
                "fraction of the cells from a given cell type with non-zero counts of a particular gene. Therefore, "
                "expression values range from 0-1, where zero means that no cells of a cell type express the gene, "
                "and one means that all cells of a cell type express the gene."
            ])),
            html.Li(html.P([
                "Minimum Ligand abs(Log2FC): While the Ligand Log2FC FDR Cutoff option filters interactions according "
                "to whether a ligand's differential expression between conditions significantly differs from zero, "
                "this filter additionally allows for further refinement by requiring a minimum absolute Log2FC. The "
                "default value reflects the cutoff used for the Circos plot."
            ])),
            html.Li(html.P([
                "Interaction Directionality: By default, the figure generated only depicts pairwise interactions "
                "between cell types from left to right. In other words, for every cell type depicted, interactions "
                "shown depict ligands emitted from the cell type directly to the left of a given cell type. "
                "Bidirectional interactions can also be enabled, which can visualize the two-way cross-talk of cells."
            ])),
            html.Li(html.P([
                "First Cell Type: The first cell type of interest to depict interactions between."
            ])),
            html.Li(html.P([
                "Second Cell Type: The second cell type of interest to depict interactions between."
            ])),
            html.Li(html.P([
                "Third Cell Type: The third cell type to include in the figure; this is optional. If not specified, "
                "only two cell types will be shown.",
            ])),
            html.Li(html.P([
                "Interaction Set: We have evaluated both CIN-dependent and CIN/STING-dependent interaction effects. "
                "This toggle lets users instantly see the plot under both conditions. Note that the CIN/STING "
                "effects represent the maximum value between shared interactions across both conditions and require "
                "interaction effects to have the same directionality."
            ]))
        ]),
        html.P(html.H6(html.Strong('How to Interpret the Plot'))),
        html.P("Selected cell types will be represented as columns in the figure according to the order selected "
               "(ex. first cell type is the leftmost column), with selected genes listed within it. These genes may "
               "be either a ligand (circle), receptor (square), or a gene that is both a ligand and receptor (diamond "
               "with dot). Additionally, these genes are colored according to the proportion of cells of the given "
               "cell type expressing each gene. Note: when many genes are selected, it is possible for labels to "
               "overlap, in which case the user can toggle the label of each gene by clicking on the representative "
               "node. Columns of genes are illustrated below:"),
        html.Div(html.Img(src=dash.get_asset_url('pairwise_help1.png'), style={'width': '20%', 'align': 'center'}, alt='Cell Type Columns', className='mx-auto'), className='text-center'),
        html.P("To illustrate interactions between ligands and receptors, arrows are drawn between pairs that have "
               "interactions that meet the user-defined filters. The colors of these arrows indicate the Log2FC "
               "strength and directionality of the ligand between conditions (red corresponds to up-regulation, and "
               "blue corresponds to down-regulation). Additionally, each arrow has a thickness relative to the number "
               "of significant interaction effects in the receptor. Below is an example of the arrows in a "
               "unidirectional interactions plot. However, users can allow interactions to start at either cell type "
               "by selecting the \"bidirectional\" Interaction Directionality setting."),
        html.Div(html.Img(src=dash.get_asset_url('pairwise_help2.png'), style={'width': '20%', 'align': 'center'}, alt='Arrows', className='mx-auto'), className='text-center')
    ]