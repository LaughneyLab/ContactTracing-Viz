import os.path

import dash
from dash import dcc, callback, Output, Input, State, html, clientside_callback
import dash_bootstrap_components as dbc
from plotly import graph_objects as go

from viz.docs import interactions_help, interaction_effects_def, ligand_log2fc_def, conditions_def
from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output, \
    make_custom_slider, make_fdr_slider

if __name__ != '__main__':
    dash.register_page(__name__,
                       path='/interactions',
                       name='Pairwise Interactions',
                       order=2)


def build_interface() -> list:
    # Ignore if this file is run directly
    if __name__ == '__main__':
        return []

    from viz.figures import DEFAULT_INTERACTIONS_ARGS, INTERACTIONS_SAVE_LOCATION
    from viz.data import using_custom_data, get_custom_celltypes

    import pickle
    default_plots = None
    try:
        if not using_custom_data():
            with open(INTERACTIONS_SAVE_LOCATION, 'rb') as f:
                default_plots = pickle.load(f)
    except:
        pass

    custom_celltypes = get_custom_celltypes() if using_custom_data() else []

    controls = control_panel("submit-button-bipartite",
         [
             control_panel_element("Interaction Set",
                                   [
                                       "The ",
                                       *conditions_def("biological condition"),
                                       " to compare."
                                   ],
                                   dbc.RadioItems(
                                       id='inter_set',
                                       options=[{'label': 'Custom', 'value': 'custom'}],
                                       value='custom',
                                       persistence=False
                                   ) if using_custom_data() else dbc.RadioItems(
                                       id='inter_set',
                                       options=[{'label': 'CIN-Dependent Interactions', 'value': 'cin'},
                                                {'label': 'CIN & STING-Dependent Interactions', 'value': 'max'}],
                                       value=DEFAULT_INTERACTIONS_ARGS['inter_set'],
                                       persistence=False
                                   )),
             control_panel_element("Plot", "You must click submit to update the plot.",
                                   html.Div(
                                       dbc.Row([dbc.Col(
                                           dbc.Button(
                                               wrap_icon("fas fa-play", "Submit"),
                                               id="submit-button-bipartite",
                                               size="lg",
                                               color="primary",
                                               className='me-1',
                                               n_clicks=0),
                                           width='auto', align='left'),
                                           dbc.Col(
                                               dbc.Button(wrap_icon('fas fa-rotate-left', 'Reset'),
                                                          id='reset_default_interactions_button',
                                                          color='secondary', size='lg',
                                                          className='float-end dark', n_clicks=0),
                                               align='right', width='auto')
                                       ]),
                                       className='text-center d-grid d-md-block gap-1 col-gap-2'))
         ], [
             control_panel_element("First Cell Type",
                                   "The first cell type to examine interactions between.",
                                   dbc.Select(
                                       id='first_celltype',
                                       options=[{'label': ct, 'value': ct} for ct in custom_celltypes],
                                       value=custom_celltypes[0]
                                   ) if using_custom_data() else dbc.Select(
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
                                                 'Mast cells']],
                                       value=DEFAULT_INTERACTIONS_ARGS['first_celltype']
                                   )),
             control_panel_element("Second Cell Type",
                                   "The second cell type to examine interactions between (need not be unique).",
                                   dbc.Select(
                                       id='second_celltype',
                                       options=[{'label': ct, 'value': ct} for ct in custom_celltypes],
                                       value=custom_celltypes[2]
                                   ) if using_custom_data() else dbc.Select(
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
                                                 'Mast cells']],
                                       value=DEFAULT_INTERACTIONS_ARGS['second_celltype']
                                   )),
             control_panel_element("Third Cell Type",
                                   "If specified, include interactions across a third cell type.",
                                   dbc.Select(
                                       id='third_celltype',
                                       options=[{'label': ct, 'value': ct} for ct in (['(None)']+custom_celltypes)],
                                       value='(None)'
                                   ) if using_custom_data() else dbc.Select(
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
                                                 'Mast cells']],
                                       value='(None)'
                                   )),
         ], [
            control_panel_element("Interaction Effect FDR Cutoff",
                                  [
                                      "FDR-adjusted requirements for ",
                                      *interaction_effects_def("interaction effects"),
                                      "."
                                  ],
                                  make_fdr_slider('bipartite_inter_fdr',
                                                  DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'])),
            control_panel_element("Ligand Log2FC FDR Cutoff",
                                  [
                                      "FDR-adjusted requirements for ",
                                      *ligand_log2fc_def("ligand differential expression"),
                                      "."
                                  ],
                                  make_fdr_slider('bipartite_logfc_fdr',
                                                  DEFAULT_INTERACTIONS_ARGS['bipartite_logfc_fdr'])),
        ], [
            control_panel_element("Minimum Interaction Effect",
                                  [
                                      "The minimum number of significant target gene ",
                                      *interaction_effects_def("interaction effects"),
                                      "."
                                  ],
                                  make_custom_slider(
                                      id="min_numsigi1_bipartite",
                                      min=0,
                                      max=100,  # Fill in
                                      step=1,
                                      value=DEFAULT_INTERACTIONS_ARGS['min_numsigi1_bipartite']
                                  )),
            control_panel_element("Minimum Expression",
                                  "The minimum fraction of cells expressing genes to be considered.",
                                  make_custom_slider(
                                      id="min_expression_bipartite",
                                      min=0,
                                      max=1,
                                      step=0.01,
                                      value=DEFAULT_INTERACTIONS_ARGS['min_expression_bipartite']
                                  ))
        ], [
            control_panel_element("Minimum Ligand abs(Log2FC)",
                                  [
                                      'The minimum ',
                                       *ligand_log2fc_def('ligand log2FC'),
                                       ' required for a chord to be drawn.'
                                  ],
                                  make_custom_slider(
                                      id="min_logfc_bipartite",
                                      min=0,
                                      max=1,  # Fill in
                                      step=0.01,
                                      value=DEFAULT_INTERACTIONS_ARGS['min_logfc_bipartite']
                                  )),
            control_panel_element("Interaction Directionality",
                                  "Whether to consider ligand/receptor pairs between all cell types, or only display these pairs from cell types from left to right.",
                                  dbc.RadioItems(
                                      id='bidirectional_bipartite',
                                      options=[{'label': 'Bidirectional Interactions', 'value': True},
                                               {'label': 'Unidirectional Interactions', 'value': False}],
                                      value=DEFAULT_INTERACTIONS_ARGS['bidirectional_bipartite'],
                                      persistence=False
                                  )),
        ],
    )

    results = figure_output(
        title="Pairwise Cell-Cell Interactions",
        subtitle=dbc.Alert(["The pairwise cell type interactions plot highlights the condition-specific signals sent by a donor "
                "cell type (emitting ligands) to a target cell type (expressing receptors). This figure gives a "
                "focused view of cell-cell interactions of interest."], color='dark'),
        footer="Circle = Ligand, Square = Receptor, Diamond = Ligand and Receptor",
        element=dcc.Graph(id="celltype-interaction-graph",
                          figure=default_plots[1] if default_plots is not None else go.Figure(data=[go.Scatter(x=[], y=[])]),
                          config={
                              'displaylogo': False,
                              'showTips': True,
                              'toImageButtonOptions': {
                                  'format': 'svg',
                                  'filename': 'interactions_image',
                                  'height': 1020 * 1.3,
                                  'width': 1020,
                                  'scale': 1
                              },
                              'watermark': False
                          }),
        help_info=interactions_help(),
        download_btn_id="download-bipartite-interactions-btn",
    )

    return [
        results,
        html.Br(),
        html.Hr(),
        html.Br(),
        controls,
        dcc.Store(id='bipartite-placeholder', storage_type='memory', data=None),
        dcc.Store(id='cin_bipartite_plot', data=default_plots[0] if default_plots is not None else {}),
        dcc.Store(id='max_bipartite_plot', data=default_plots[1] if default_plots is not None else {}),
        dcc.Download(id="download-bipartite-interactions"),
        # DeferScript(src=dash.get_asset_url('plotly_hooks.js'))
    ]


# When the graph is updated, inject custom hooks
clientside_callback("""
function(_) {
    plotlyInjection(1.3);
    return null;
}
""",
Output('bipartite-placeholder', 'data'),
Input("celltype-interaction-graph", 'children'),
)


@callback(
    Output('first_celltype', 'value', allow_duplicate=True),
    Output('second_celltype', 'value', allow_duplicate=True),
    Output('third_celltype', 'value', allow_duplicate=True),
    Output('min_logfc_bipartite', 'data', allow_duplicate=True),
    Output('min_expression_bipartite', 'data', allow_duplicate=True),
    Output('min_numsigi1_bipartite', 'data', allow_duplicate=True),
    Output('bipartite_inter_fdr', 'data', allow_duplicate=True),
    Output('bipartite_logfc_fdr', 'data', allow_duplicate=True),
    Output('bidirectional_bipartite', 'value', allow_duplicate=True),
    Output('submit-button-bipartite', 'n_clicks'),
    Input('reset_default_interactions_button', 'n_clicks'),
    State('submit-button-bipartite', 'n_clicks'),
    prevent_initial_call=True
)
def reset_to_defaults(n_clicks, submission_button_clicks):
    from dash.exceptions import PreventUpdate
    from viz.figures import DEFAULT_INTERACTIONS_ARGS
    if n_clicks > 0:
        return DEFAULT_INTERACTIONS_ARGS['first_celltype'], \
                DEFAULT_INTERACTIONS_ARGS['second_celltype'], \
                DEFAULT_INTERACTIONS_ARGS['third_celltype'], \
                DEFAULT_INTERACTIONS_ARGS['min_logfc_bipartite'], \
                DEFAULT_INTERACTIONS_ARGS['min_expression_bipartite'], \
                DEFAULT_INTERACTIONS_ARGS['min_numsigi1_bipartite'], \
                DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'], \
                DEFAULT_INTERACTIONS_ARGS['bipartite_logfc_fdr'], \
                DEFAULT_INTERACTIONS_ARGS['bidirectional_bipartite'], \
                submission_button_clicks+1  # Trigger the plotting callback to run on default values
    raise PreventUpdate


@callback(
    Output('celltype-interaction-graph', 'figure'),
    Input('inter_set', 'value'),
    Input('cin_bipartite_plot', 'data'),
    Input('max_bipartite_plot', 'data'),
    prevent_initial_call=True
)
def update_bipartite_plot(inter_set, cin_bipartite_plot, max_bipartite_plot):
    if inter_set == 'cin' or inter_set == 'custom':
        return cin_bipartite_plot
    else:
        return max_bipartite_plot


@callback(
    Output('download-bipartite-interactions', 'data'),
    Input('download-bipartite-interactions-btn', 'n_clicks'),
    State('inter_set', 'value'),
    State('first_celltype', 'value'),
    State('second_celltype', 'value'),
    State('third_celltype', 'value'),
    State('min_logfc_bipartite', 'data'),
    State('min_expression_bipartite', 'data'),
    State('min_numsigi1_bipartite', 'data'),
    State('bipartite_inter_fdr', 'data'),
    State('bipartite_logfc_fdr', 'data'),
    State('bidirectional_bipartite', 'value'),
    prevent_initial_call=True,
)
def download_data(n_clicks,
                  inter_set,
                  cell1, cell2, cell3,
                  min_logfc, min_expression, min_numSigI1,
                  inter_fdr, logfc_fdr, bidirectional
                  ):
    if n_clicks <= 0:
        from dash.exceptions import PreventUpdate
        raise PreventUpdate

    from viz.data import read_interactions_file, using_custom_data
    if using_custom_data():
        df = read_interactions_file('custom', inter_fdr)
    else:
        df = read_interactions_file('cin' if inter_set == 'cin' else 'sting', inter_fdr)
    logfc_fdr_cutoff = float(logfc_fdr.replace('fdr', '.'))

    if cell3 is not None and 'none' in str(cell3).lower():
        cell3 = None

    if cell3 is None:
        selected = df[((df.cell_type_ligand == cell1) & (df.cell_type_receptor == cell2)) |
                      (bidirectional and (df.cell_type_ligand == cell2) & (df.cell_type_receptor == cell1))]
    else:
        selected = df[((df.cell_type_ligand == cell1) & (df.cell_type_receptor == cell2)) |
                      (bidirectional and (df.cell_type_ligand == cell2) & (df.cell_type_receptor == cell1)) |
                      ((df.cell_type_ligand == cell2) & (df.cell_type_receptor == cell3)) |
                      (bidirectional and (df.cell_type_ligand == cell3) & (df.cell_type_receptor == cell2)) |
                      (bidirectional & (df.cell_type_ligand == cell1) & (df.cell_type_receptor == cell3)) |
                      (bidirectional & (df.cell_type_ligand == cell3) & (df.cell_type_receptor == cell1))]

    import numpy as np
    selected = selected[(selected.numSigI1 >= min_numSigI1) &
                        (selected.expression_ligand >= min_expression) &
                        (selected.expression_receptor >= min_expression) &
                        (selected["MAST_fdr_ligand"] < logfc_fdr_cutoff) &
                        (selected["MAST_fdr_receptor"] < logfc_fdr_cutoff) &
                        (np.abs(selected["MAST_log2FC_ligand"]) > min_logfc)
                        ].sort_values(
        by="numSigI1", ascending=False)

    # Only keep the columns we want to download
    selected = selected[["cell_type_ligand", "cell_type_receptor", "ligand", "receptor", "numSigI1", "numDEG",
                            "expression_ligand", "expression_receptor", "MAST_log2FC_ligand", "MAST_log2FC_receptor",
                            "MAST_fdr_ligand", "MAST_fdr_receptor"]].reset_index(drop=True)

    return dcc.send_data_frame(selected.to_csv, f"{inter_set}_interactions.csv")


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
    cache_args_to_ignore=[0],
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
    from viz.data import read_interactions_file, using_custom_data
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
            if not using_custom_data():
                with open(INTERACTIONS_SAVE_LOCATION, 'rb') as f:
                    # return [(Serverside(p) if p is not None else None) for p in pickle.load(f)]  FIXME
                    return pickle.load(f)
        except:
            pass

    # Do some work
    if using_custom_data():
        cin_interactions = read_interactions_file('custom', inter_fdr)
        sting_interactions = None
    else:
        cin_interactions = read_interactions_file('cin', inter_fdr)
        sting_interactions = read_interactions_file('sting', inter_fdr)

    cin_fig = bipartite_graph(
        df=cin_interactions,
        cin_only=True,
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

    if using_custom_data():
        sting_fig = None
    else:
        sting_fig = bipartite_graph(
            df=sting_interactions,
            cin_only=False,
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

    # return [Serverside(cin_fig), Serverside(sting_fig)]  FIXME
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
    interactive_panel(html.Div(),  # wrap_icon('fa-arrows-left-right-to-line', 'Cell Type to Cell Type Interactions'),
                      html.Div(),  # "Identify crosstalk between cell types.",
                      *build_interface()
                      )
]


# If run as a script, compile the default plot
if __name__ == '__main__':
    from viz.figures import bipartite_graph, DEFAULT_INTERACTIONS_ARGS, INTERACTIONS_SAVE_LOCATION
    from viz.data import read_interactions_file
    import pickle
    cin_interactions = read_interactions_file('cin', DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'])
    sting_interactions = read_interactions_file('max', DEFAULT_INTERACTIONS_ARGS['bipartite_inter_fdr'])

    cin_fig = bipartite_graph(
        df=cin_interactions,
        cin_only=True,
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
        cin_only=False,
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
