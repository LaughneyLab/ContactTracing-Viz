import os.path

import dash
from dash import dcc, callback, Output, Input, State, html
import dash_bootstrap_components as dbc

from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output, \
    make_custom_slider, make_fdr_slider

if __name__ != '__main__':
    dash.register_page(__name__,
                       path='/interactions',
                       name='Interactions',
                       order=2)


def build_interface() -> list:
    # Ignore if this file is run directly
    if __name__ == '__main__':
        return []

    from viz.figures import DEFAULT_INTERACTIONS_ARGS, INTERACTIONS_SAVE_LOCATION

    import pickle
    default_plots = None
    try:
        with open(INTERACTIONS_SAVE_LOCATION, 'rb') as f:
            default_plots = pickle.load(f)
    except:
        pass

    controls = control_panel(
        [
            control_panel_element("Interaction Effect FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  make_fdr_slider('bipartite_inter_fdr',
                                                  DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'])),
            control_panel_element("Ligand Log2FC FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  make_fdr_slider('bipartite_logfc_fdr',
                                                  DEFAULT_INTERACTIONS_ARGS['bipartite_logfc_fdr'])),
        ], [
            control_panel_element("Minimum Interaction Effect", "The minimum number of significant target gene interactions.",
                                  make_custom_slider(
                                      id="min_numsigi1_bipartite",
                                      min=0,
                                      max=100,  # Fill in
                                      step=1,
                                      value=DEFAULT_INTERACTIONS_ARGS['min_numsigi1_bipartite']
                                  )),
            control_panel_element("Minimum Expression", "The minimum fraction of cells expressing genes to be considered.",
                                  make_custom_slider(
                                      id="min_expression_bipartite",
                                      min=0,
                                      max=1,
                                      step=0.01,
                                      value=DEFAULT_INTERACTIONS_ARGS['min_expression_bipartite']
                                  ))
        ], [
            control_panel_element("Minimum Ligand abs(Log2FC)", "The minimum log2FC for ligands between conditions.",
                                  make_custom_slider(
                                      id="min_logfc_bipartite",
                                      min=0,
                                      max=1,  # Fill in
                                      step=0.01,
                                      value=DEFAULT_INTERACTIONS_ARGS['min_logfc_bipartite']
                                  )),
            control_panel_element("Interaction Directionality", "Whether to consider interactions between all cell types, or only allow interactions to flow from left to right.",
                                  dbc.RadioItems(
                                      id='bidirectional_bipartite',
                                      options=[{'label': 'Bidirectional Interactions', 'value': True},
                                               {'label': 'Unidirectional Interactions', 'value': False}],
                                      value=DEFAULT_INTERACTIONS_ARGS['bidirectional_bipartite'],
                                      persistence=False
                                  )),
        ], [
            control_panel_element("First Cell Type", "The first cell type to examine interactions between.",
                                  dbc.Select(
                                      id='first_celltype',
                                      options=[{'label': ct, 'value': ct} for ct in
                                               ['Tumor cells',
                                                 'Macrophages/mMDSC',
                                                 'PMN/gMDSC',
                                                 'T cells',
                                                 'B cells',
                                                 'NK cells',
                                                 'cDC',
                                                 'pDC',
                                                 'Fibroblast cells',
                                                 'Endothelial cells',
                                                 'Osteoclasts',
                                                 'Mast cells']],
                                      value=DEFAULT_INTERACTIONS_ARGS['first_celltype']
                                  )),
            control_panel_element("Second Cell Type", "The second cell type to examine interactions between (need not be unique).",
                                  dbc.Select(
                                      id='second_celltype',
                                    options=[{'label': ct, 'value': ct} for ct in
                                             ['Tumor cells',
                                             'Macrophages/mMDSC',
                                             'PMN/gMDSC',
                                             'T cells',
                                             'B cells',
                                             'NK cells',
                                             'cDC',
                                             'pDC',
                                             'Fibroblast cells',
                                             'Endothelial cells',
                                             'Osteoclasts',
                                             'Mast cells']],
                                    value=DEFAULT_INTERACTIONS_ARGS['second_celltype']
                                  )),
            control_panel_element("Third Cell Type", "If specified, include interactions across a third cell type.",
                                  dbc.Select(
                                      id='third_celltype',
                                      options=[{'label': ct, 'value': ct} for ct in
                                               ['(None)',
                                                'Tumor cells',
                                                'Macrophages/mMDSC',
                                                'PMN/gMDSC',
                                                'T cells',
                                                'B cells',
                                                'NK cells',
                                                'cDC',
                                                'pDC',
                                                'Fibroblast cells',
                                                'Endothelial cells',
                                                'Osteoclasts',
                                                'Mast cells']],
                                      value=DEFAULT_INTERACTIONS_ARGS['third_celltype']
                                  )),
        ], [
            control_panel_element("Interaction Set", "Biological condition to compare.",
                                  dbc.RadioItems(
                                      id='inter_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'CIN & STING Max Effect', 'value': 'max'}],
                                      value=DEFAULT_INTERACTIONS_ARGS['inter_set'],
                                      persistence=False
                                  )),
            control_panel_element("Plot", "You must click submit to update the plot.",
                                  html.Div(
                                      dbc.Button(
                                          "Submit",
                                          id="submit-button-bipartite",
                                          size="lg",
                                          color="primary",
                                          className='me-1',
                                          n_clicks=0
                                      ),
                                  className='text-center d-grid gap-2')),
        ]
    )

    results = figure_output(
        title="Cell Type Interactions",
        footer="Circle = Ligand, Square = Receptor, Diamond = Ligand and Receptor",
        element=dcc.Graph(id="celltype-interaction-graph",
                          figure=default_plots[1] if default_plots is not None else {},
                          config={
                              'displaylogo': False,
                              'showTips': True,
                              'toImageButtonOptions': {
                                  'format': 'svg',
                                  'filename': 'interactions_image',
                                  'height': 1275,
                                  'width': 1275*.8,
                                  'scale': 1
                              },
                              'watermark': False
                          }),
        help_info=make_help_info()
    )

    return [
        controls,
        results,
        dcc.Store(id='cin_bipartite_plot', data=default_plots[0] if default_plots is not None else {}),
        dcc.Store(id='max_bipartite_plot', data=default_plots[1] if default_plots is not None else {}),
    ]


def make_help_info():
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
                "Interaction Set: We have evaluated both CIN-dependent and CIN & STING-dependent interaction effects. "
                "This toggle lets users instantly see the plot under both conditions. Note that the CIN & STING "
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


@callback(
    Output('celltype-interaction-graph', 'figure'),
    Input('inter_set', 'value'),
    Input('cin_bipartite_plot', 'data'),
    Input('max_bipartite_plot', 'data'),
    prevent_initial_call=True
)
def update_bipartite_plot(inter_set, cin_bipartite_plot, max_bipartite_plot):
    if inter_set == 'cin':
        return cin_bipartite_plot
    else:
        return max_bipartite_plot


@callback(
    Output('cin_bipartite_plot', 'data'),
    Output('max_bipartite_plot', 'data'),
    Input('submit-button-bipartite', 'n_clicks'),
    State('first_celltype', 'value'),
    State('second_celltype', 'value'),
    State('third_celltype', 'value'),
    State('min_logfc_bipartite', 'data'),
    State('min_expression_bipartite', 'data'),
    State('min_numsigi1_bipartite', 'data'),
    State('bipartite_inter_fdr', 'data'),
    State('bipartite_logfc_fdr', 'data'),
    State('bidirectional_bipartite', 'value'),
    interval=500,
    cache_args_to_ignore=['submit-button-bipartite', 'n_clicks'],
    background=True,  # Run in background
    prevent_initial_call=True,
    running=[  # Disable the button while the callback is running
        (Output('submit-button-bipartite', 'disabled'), True, False),
        (Output('progress-bar', 'style'), {'visibility': 'visible'}, {'visibility': 'hidden'})
    ],
    progress=[  # Show a progress bar while the callback is running
        Output('progress-bar', "value"),
        Output('progress-bar', "max")
    ]
)
def make_graph(set_progress, n_clicks,
               first_ct, second_ct, third_ct,
               min_logfc, min_expression, min_numSigI1,
               inter_fdr, logfc_fdr, bidirectional):
    if n_clicks == 0:
        from dash.exceptions import PreventUpdate
        raise PreventUpdate

    from viz.data import read_interactions_file
    from viz.figures import bipartite_graph, DEFAULT_INTERACTIONS_ARGS, INTERACTIONS_SAVE_LOCATION

    set_progress((0, 2))

    # Check if arguments match default, if so return the pre-computed default
    if (first_ct == DEFAULT_INTERACTIONS_ARGS['first_celltype'] and
            second_ct == DEFAULT_INTERACTIONS_ARGS['second_celltype'] and
            third_ct == DEFAULT_INTERACTIONS_ARGS['third_celltype'] and
            min_logfc == DEFAULT_INTERACTIONS_ARGS['min_logfc_bipartite'] and
            min_expression == DEFAULT_INTERACTIONS_ARGS['min_expression_bipartite'] and
            min_numSigI1 == DEFAULT_INTERACTIONS_ARGS['min_numsigi1_bipartite'] and
            inter_fdr == DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'] and
            logfc_fdr == DEFAULT_INTERACTIONS_ARGS['bipartite_logfc_fdr'] and
            bidirectional == DEFAULT_INTERACTIONS_ARGS['bidirectional_bipartite']):
        import pickle
        try:
            with open(INTERACTIONS_SAVE_LOCATION, 'rb') as f:
                return pickle.load(f)
        except:
            pass

    # Do some work
    cin_interactions = read_interactions_file('cin', inter_fdr)
    sting_interactions = read_interactions_file('sting', inter_fdr)

    cin_fig = bipartite_graph(
        df=cin_interactions,
        cell1=first_ct,
        cell2=second_ct,
        cell3=third_ct,
        numInteractions=min_numSigI1,
        min_logfc_bipartite=min_logfc,
        min_expression_bipartite=min_expression,
        logfc_fdr_bipartite_cutoff=float(logfc_fdr.replace('fdr', '.')),
        bidirectional=bidirectional
    )

    set_progress((1, 2))

    sting_fig = bipartite_graph(
        df=sting_interactions,
        cell1=first_ct,
        cell2=second_ct,
        cell3=third_ct,
        numInteractions=min_numSigI1,
        min_logfc_bipartite=min_logfc,
        min_expression_bipartite=min_expression,
        logfc_fdr_bipartite_cutoff=float(logfc_fdr.replace('fdr', '.')),
        bidirectional=bidirectional
    )

    set_progress((2, 2))

    return [cin_fig, sting_fig]


# @callback(
#     Output('min_logfc_bipartite', 'max'),
#     Output('min_numsigi1_bipartite', 'max'),
#     Input('bipartite_inter_fdr', 'value'),
#     Input('bipartite_logfc_fdr', 'value'),
#     interval=10,
#     background=False,  # Run in background,
#     running=[  # Disable the button while the callback is running
#         (Output('submit-button-bipartite', 'disabled'), True, False),
#         (Output('spinner-holder', 'children'), [
#             dbc.Spinner(color='primary',
#                         size='md',
#                         fullscreen=True,
#                         type='grow')],
#          []),
#     ]
# )
# def initialize_options(inter_fdr, logfc_fdr):
#     from viz.data import read_interactions_file
#
#     fdr = inter_fdr
#
#     # Use the max set
#     interactions = read_interactions_file('highCIN_vs_noSTING', fdr)
#
#     max_logfc = max(max(interactions.MAST_log2FC_ligand.abs()), max(interactions.MAST_log2FC_receptor.abs()))
#     max_inter = max(interactions.numSigI1)
#
    return max_logfc, max_inter


layout = [
    interactive_panel(wrap_icon('fa-arrows-left-right-to-line', 'Cell Type to Cell Type Interactions'), "Identify crosstalk between cell types.",
                      *build_interface()
                      )
]


# If run as a script, compile the default plot
if __name__ == '__main__':
    from viz.figures import bipartite_graph, DEFAULT_INTERACTIONS_ARGS, INTERACTIONS_SAVE_LOCATION
    if not os.path.exists(INTERACTIONS_SAVE_LOCATION):
        from viz.data import read_interactions_file
        import pickle
        cin_interactions = read_interactions_file('cin', DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'])
        sting_interactions = read_interactions_file('max', DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'])

        cin_fig = bipartite_graph(
            df=cin_interactions,
            cell1=DEFAULT_INTERACTIONS_ARGS['first_celltype'],
            cell2=DEFAULT_INTERACTIONS_ARGS['second_celltype'],
            cell3=DEFAULT_INTERACTIONS_ARGS['third_celltype'],
            numInteractions=DEFAULT_INTERACTIONS_ARGS['min_numsigi1_bipartite'],
            min_logfc_bipartite=DEFAULT_INTERACTIONS_ARGS['min_logfc_bipartite'],
            min_expression_bipartite=DEFAULT_INTERACTIONS_ARGS['min_expression_bipartite'],
            logfc_fdr_bipartite_cutoff=float(DEFAULT_INTERACTIONS_ARGS['bipartite_logfc_fdr'].replace('fdr', '.')),
            bidirectional=DEFAULT_INTERACTIONS_ARGS['bidirectional_bipartite']
        )

        sting_fig = bipartite_graph(
            df=sting_interactions,
            cell1=DEFAULT_INTERACTIONS_ARGS['first_celltype'],
            cell2=DEFAULT_INTERACTIONS_ARGS['second_celltype'],
            cell3=DEFAULT_INTERACTIONS_ARGS['third_celltype'],
            numInteractions=DEFAULT_INTERACTIONS_ARGS['min_numsigi1_bipartite'],
            min_logfc_bipartite=DEFAULT_INTERACTIONS_ARGS['min_logfc_bipartite'],
            min_expression_bipartite=DEFAULT_INTERACTIONS_ARGS['min_expression_bipartite'],
            logfc_fdr_bipartite_cutoff=float(DEFAULT_INTERACTIONS_ARGS['bipartite_logfc_fdr'].replace('fdr', '.')),
            bidirectional=DEFAULT_INTERACTIONS_ARGS['bidirectional_bipartite']
        )

        with open(INTERACTIONS_SAVE_LOCATION, 'wb') as f:
            pickle.dump([cin_fig, sting_fig], f)
