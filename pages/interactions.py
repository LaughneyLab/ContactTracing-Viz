import os.path

import dash
from dash import dcc, callback, Output, Input, State
import dash_bootstrap_components as dbc

from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output, \
    make_custom_slider, make_fdr_slider

if __name__ != '__main__':
    dash.register_page(__name__,
                       path='/interactions',
                       name='Interactions',
                       order=2)


def build_interface() -> list:
    from viz.figures import DEFAULT_INTERACTIONS_ARGS, INTERACTIONS_SAVE_LOCATION

    import pickle
    default_plot = {}
    try:
        with open(INTERACTIONS_SAVE_LOCATION, 'rb') as f:
            default_plot = pickle.load(f)
    except:
        pass

    controls = control_panel(
        [
            control_panel_element("Interaction Set", "Biological condition to compare.",
                                  dbc.RadioItems(
                                      id='inter_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'CIN & STING Max Effect', 'value': 'max'}],
                                      value=DEFAULT_INTERACTIONS_ARGS['inter_set']
                                  ))
        ], [
            control_panel_element("Minimum numSigI1", "The minimum number of significant target gene interactions.",
                                  make_custom_slider(
                                      id="min_numsigi1_bipartite",
                                      min=0,
                                      max=10,  # Fill in
                                      step=1,
                                      value=DEFAULT_INTERACTIONS_ARGS['min_numsigi1_bipartite']
                                  ))
        ], [
            control_panel_element("Minimum abs(log2FC)", "The minimum log2FC for either ligands or receptors between conditions.",
                                  make_custom_slider(
                                      id="min_logfc_bipartite",
                                      min=0,
                                      max=1,  # Fill in
                                      step=0.01,
                                      value=DEFAULT_INTERACTIONS_ARGS['min_logfc_bipartite']
                                  ))
        ], [
            control_panel_element("Minimum Expression", "The minimum fraction of cells expressing genes to be considered.",
                                  make_custom_slider(
                                      id="min_expression_bipartite",
                                      min=0,
                                      max=1,
                                      step=0.01,
                                      value=DEFAULT_INTERACTIONS_ARGS['min_expression_bipartite']
                                  ))
        ], [
            control_panel_element("Interaction Effect FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  make_fdr_slider('bipartite_inter_fdr', DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'])),
            control_panel_element("log2FC FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  make_fdr_slider('bipartite_logfc_fdr', DEFAULT_INTERACTIONS_ARGS['bipartite_logfc_fdr'])),
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
            control_panel_element("Plot", "",
                                  dbc.Button(
                                      "Submit",
                                      id="submit-button-bipartite",
                                      size="lg",
                                      color="primary",
                                      className='me-1',
                                      n_clicks=0
                                  ))
        ]
    )

    results = figure_output(
        title="Cell Type Interactions",
        footer="Circle = Ligand, Square = Receptor, Diamond = Ligand and Receptor",
        element=dcc.Graph(id="celltype-interaction-graph",
                          figure=default_plot,
                          config={
                              'displaylogo': False,
                              'showTips': True,
                              'toImageButtonOptions': {
                                  'format': 'png',
                                  'filename': 'exported_image',
                                  'height': 800,
                                  'width': 1200,
                                  'scale': 6
                              },
                              'watermark': False
                          })
    )

    return [
        controls,
        results
    ]


@callback(
    Output('celltype-interaction-graph', 'figure'),
    Input('submit-button-bipartite', 'n_clicks'),
    State('inter_set', 'value'),
    State('first_celltype', 'value'),
    State('second_celltype', 'value'),
    State('third_celltype', 'value'),
    State('min_logfc_bipartite', 'data'),
    State('min_expression_bipartite', 'data'),
    State('min_numsigi1_bipartite', 'data'),
    State('bipartite_inter_fdr', 'data'),
    State('bipartite_logfc_fdr', 'data'),
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
               inter_set,
               first_ct, second_ct, third_ct,
               min_logfc, min_expression, min_numSigI1,
               inter_fdr, logfc_fdr):
    if n_clicks == 0:
        from dash.exceptions import PreventUpdate
        raise PreventUpdate

    from viz.data import read_interactions_file
    from viz.figures import bipartite_graph, DEFAULT_INTERACTIONS_ARGS, INTERACTIONS_SAVE_LOCATION

    set_progress((50, 100))

    # Check if arguments match default, if so return the pre-computed default
    if (inter_set == DEFAULT_INTERACTIONS_ARGS['inter_set'] and
            first_ct == DEFAULT_INTERACTIONS_ARGS['first_celltype'] and
            second_ct == DEFAULT_INTERACTIONS_ARGS['second_celltype'] and
            third_ct == DEFAULT_INTERACTIONS_ARGS['third_celltype'] and
            min_logfc == DEFAULT_INTERACTIONS_ARGS['min_logfc_bipartite'] and
            min_expression == DEFAULT_INTERACTIONS_ARGS['min_expression_bipartite'] and
            min_numSigI1 == DEFAULT_INTERACTIONS_ARGS['min_numsigi1_bipartite'] and
            inter_fdr == DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'] and
            logfc_fdr == DEFAULT_INTERACTIONS_ARGS['bipartite_logfc_fdr']):
        import pickle
        try:
            with open(INTERACTIONS_SAVE_LOCATION, 'rb') as f:
                return pickle.load(f)
        except:
            pass

    # Do some work
    interactions = read_interactions_file(inter_set, inter_fdr)

    fig = bipartite_graph(
        df=interactions,
        cell1=first_ct,
        cell2=second_ct,
        cell3=third_ct,
        numInteractions=min_numSigI1,
        min_logfc_bipartite=min_logfc,
        min_expression_bipartite=min_expression,
        logfc_fdr_bipartite_cutoff=float(logfc_fdr.replace('fdr', '.'))
    )

    return fig


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
        interactions = read_interactions_file(DEFAULT_INTERACTIONS_ARGS['inter_set'], DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'])

        fig = bipartite_graph(
            df=interactions,
            cell1=DEFAULT_INTERACTIONS_ARGS['first_celltype'],
            cell2=DEFAULT_INTERACTIONS_ARGS['second_celltype'],
            cell3=DEFAULT_INTERACTIONS_ARGS['third_celltype'],
            numInteractions=DEFAULT_INTERACTIONS_ARGS['min_numsigi1_bipartite'],
            min_logfc_bipartite=DEFAULT_INTERACTIONS_ARGS['min_logfc_bipartite'],
            min_expression_bipartite=DEFAULT_INTERACTIONS_ARGS['min_expression_bipartite'],
            logfc_fdr_bipartite_cutoff=float(DEFAULT_INTERACTIONS_ARGS['bipartite_logfc_fdr'].replace('fdr', '.'))
        )

        with open(INTERACTIONS_SAVE_LOCATION, 'wb') as f:
            pickle.dump(fig, f)
