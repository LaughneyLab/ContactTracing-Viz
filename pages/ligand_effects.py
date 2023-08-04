import os.path

import dash
import dash_bootstrap_components as dbc
from dash_extensions import DeferScript
from dash_extensions.enrich import callback, Output, Input, dcc, State, html, Serverside, clientside_callback
from plotly import graph_objects as go

from viz.docs import ligand_effects_help, interaction_test_def, conditions_def, deg_test_def
from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output, \
    make_custom_slider, make_fdr_slider

if __name__ != '__main__':
    dash.register_page(__name__,
                       path='/ligand-effects',
                       name='Ligand Effects',
                       order=3)


def build_interface() -> list:
    # Ignore if this file is run directly
    if __name__ == '__main__':
        return []

    from viz.figures import DEFAULT_LIGAND_EFFECT_ARGS, LIGAND_EFFECT_SAVE_LOCATION
    from viz.data import using_custom_data, get_custom_celltypes

    import pickle
    default_plots = None
    try:
        if not using_custom_data():
            with open(LIGAND_EFFECT_SAVE_LOCATION, 'rb') as f:
                default_plots = pickle.load(f)
    except:
        pass

    custom_celltypes = get_custom_celltypes() if using_custom_data() else []

    controls = control_panel("submit_button",
         [
             control_panel_element("Interaction Set",
                                   [
                                       "The ",
                                       *conditions_def("biological condition"),
                                       " to compare."
                                   ],
                                   dbc.RadioItems(
                                       id='effect_set',
                                       options=[{'label': 'Custom', 'value': 'custom'}],
                                       value='custom',
                                       persistence=False
                                   ) if using_custom_data() else dbc.RadioItems(
                                       id='effect_set',
                                       options=[{'label': 'CIN-Dependent Interactions', 'value': 'cin'},
                                                {'label': 'CIN/STING-Dependent Interactions', 'value': 'max',
                                                 'disabled': True}
                                                ],
                                       value=DEFAULT_LIGAND_EFFECT_ARGS['effect_set'],
                                       persistence=False
                                   )),
             control_panel_element("Plot", "You must click submit to update the plot.",
                                   html.Div(
                                       dbc.Row([dbc.Col(
                                           dbc.Button(
                                               wrap_icon("fas fa-play", "Submit"),
                                               id='submit_button',
                                               size='lg',
                                               color="primary",
                                               className='me-1',
                                               n_clicks=0),
                                           width='auto', align='left'),
                                           dbc.Col(
                                               dbc.Button(wrap_icon('fas fa-rotate-left', 'Reset'),
                                                          id='reset_default_ligand_effects_button',
                                                          color='secondary',
                                                          size='lg', className='float-end dark',
                                                          n_clicks=0),
                                               align='right', width='auto')
                                       ]),
                                       className='text-center d-grid d-md-block gap-1 col-gap-2'))
         ],
        [
            control_panel_element('Emitting Cell Type', 'Select the initial cell type.',
                                  dbc.Select(
                                      id='cell_type',
                                      options=[{'label': ct, 'value': ct} for ct in custom_celltypes],
                                      value=custom_celltypes[0]
                                  ) if using_custom_data() else dbc.Select(
                                    id='cell_type',
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
                                    value=DEFAULT_LIGAND_EFFECT_ARGS['cell_type']
                                  )),
            control_panel_element("Emitted Ligands",
                                  'Comma-separated list of ligands to emit from the selected cell type.',
                                  dbc.Input(
                                      id='ligands',
                                      autofocus=True,
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['ligands'],
                                      placeholder='Example: Ccl2,Apoe'
                                  )),
            control_panel_element("Network Building Iterations",
                                  "This controls how many downstream interactions may be detected.",
                                  dbc.Input(
                                      id='iterations',
                                      type='number',
                                      debounce=True,
                                      max=5,
                                      min=1,
                                      step=1,
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['iterations'],
                                      persistence=True, persistence_type='session'
                                  )),
        ], [
            control_panel_element("Interaction Effect FDR Cutoff",
                                  [
                                      "The maximum ",
                                      *interaction_test_def("interaction test"),
                                      " FDR to consider."
                                  ],
                                  make_fdr_slider("interaction_fdr", DEFAULT_LIGAND_EFFECT_ARGS['interaction_fdr'])),
            control_panel_element("Log2FC FDR Cutoff",
                                  [
                                      "The FDR-adjusted p-value cutoff for determining if a ",
                                      *deg_test_def("gene's receptor-induced log2FC"),
                                      " value is non-zero."
                                  ],
                                  make_fdr_slider("logfc_fdr", DEFAULT_LIGAND_EFFECT_ARGS['logfc_fdr']))
        ], [
            control_panel_element("Minimum Expression", "The minimum expression of target genes to be considered.",
                                  make_custom_slider(
                                      id='min_expression',
                                      max=1,
                                      min=0,
                                      step=0.01,
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['min_expression']
                                  )),
            control_panel_element("Minimum abs(Log2FC)",
                                  [
                                      "The minimum induced log2FC for gene ",
                                      *interaction_test_def("interactions"),
                                      "."
                                  ],
                                  make_custom_slider(
                                      id='min_logfc',
                                      max=2,
                                      min=0,
                                      step=0.01,
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['min_logfc']
                                  )),
        ]
    )

    results = figure_output(
        title='CIN-dependent Ligand Cascading Effects',
        footer="Circle = Ligand, Square = Receptor, Diamond = Ligand and Receptor",
        element=dcc.Graph(id='network_graph',
                          figure=default_plots[0] if default_plots else go.Figure(data=[go.Scatter(x=[], y=[])]),
                          config={
                             'displaylogo': False,
                             'showTips': True,
                             'toImageButtonOptions': {
                                'format': 'svg',  # one of png, svg, jpeg, webp
                                'filename': 'ligand_effects_image',
                                 'height': 1025 * 1.2,
                                 'width': 1025,
                                'scale': 1  # Multiply title/legend/axis/canvas sizes by this factor
                             },
                             'watermark': False
                          }),
        help_info=ligand_effects_help(),
        download_btn_id="download_ligand_effects_network_btn"
    )

    return [
        results,
        html.Br(),
        html.Hr(),
        html.Br(),
        controls,
        dcc.Store(id='cin_network_plot', storage_type='memory', data=default_plots[0] if default_plots is not None else {}),
        dcc.Store(id='sting_network_plot', storage_type='memory', data=default_plots[1] if default_plots is not None else {}),
        dcc.Download(id="download_ligand_effects_network"),
        # DeferScript(src=dash.get_asset_url('plotly_hooks.js'))
    ]


# When the graph is updated, inject custom hooks
clientside_callback("""
function(ignore) {
    plotlyInjection(1.2);
}
""", Input("network_graph", 'children'))


@callback(
    Output('cell_type', 'value', allow_duplicate=True),
    Output('ligands', 'value', allow_duplicate=True),
    Output('interaction_fdr', 'data', allow_duplicate=True),
    Output('min_logfc', 'data', allow_duplicate=True),
    Output('min_expression', 'data', allow_duplicate=True),
    Output('logfc_fdr', 'data', allow_duplicate=True),
    Output('iterations', 'value', allow_duplicate=True),
    Output('submit_button', 'n_clicks'),
    Input('reset_default_ligand_effects_button', 'n_clicks'),
    State('submit_button', 'n_clicks'),
    prevent_initial_call=True
)
def reset_to_defaults(n_clicks, submission_button_clicks):
    from dash.exceptions import PreventUpdate
    from viz.figures import DEFAULT_LIGAND_EFFECT_ARGS
    if n_clicks > 0:
        return [DEFAULT_LIGAND_EFFECT_ARGS['cell_type'],
            DEFAULT_LIGAND_EFFECT_ARGS['ligands'],
            DEFAULT_LIGAND_EFFECT_ARGS['interaction_fdr'],
            DEFAULT_LIGAND_EFFECT_ARGS['min_logfc'],
            DEFAULT_LIGAND_EFFECT_ARGS['min_expression'],
            DEFAULT_LIGAND_EFFECT_ARGS['logfc_fdr'],
            DEFAULT_LIGAND_EFFECT_ARGS['iterations'],
            submission_button_clicks+1]  # Trigger the plotting callback to run on default values
    raise PreventUpdate


@callback(
    Output('download_ligand_effects_network', 'data'),
    Input('download_ligand_effects_network_btn', 'n_clicks'),
    State('cell_type', 'value'),
    State('ligands', 'value'),
    State('interaction_fdr', 'data'),
    State('min_logfc', 'data'),
    State('min_expression', 'data'),
    State('logfc_fdr', 'data'),
    State('iterations', 'value'),
    prevent_initial_call=True
)
def download_data(n_clicks,
                  cell_type, ligands,
                  interaction_fdr, min_logfc, min_expression,
                  logfc_fdr, iterations
                  ):
    if n_clicks <= 0:
        from dash.exceptions import PreventUpdate
        raise PreventUpdate

    from viz.data import using_custom_data
    from viz.figures import pseudotime_interaction_propagation_graph, DEFAULT_LIGAND_EFFECT_ARGS
    from networkx import to_pandas_edgelist
    graph = pseudotime_interaction_propagation_graph(
        effect_set='custom' if using_custom_data() else DEFAULT_LIGAND_EFFECT_ARGS['effect_set'],  # Keep this fixed
        seed_cell=cell_type,
        seed_ligands=ligands,
        iterations=int(iterations),
        interaction_fdr_cutoff=float("." + interaction_fdr[3:]),
        min_logfc=float(min_logfc),
        min_expression=float(min_expression),
        logfc_fdr_cutoff=float("." + logfc_fdr[3:]),
        return_only_raw_data=True
    )
    df = to_pandas_edgelist(graph)

    return dcc.send_data_frame(df.to_csv, "ligand_effects_network.csv")


@callback(
    Output('network_graph', 'figure'),
    Input('effect_set', 'value'),
    Input('cin_network_plot', 'data'),
    Input('sting_network_plot', 'data'),
    prevent_initial_call=True
)
def update_network_figure(effect_set, cin_network_plot, sting_network_plot):
    if effect_set == 'cin' or effect_set == 'custom':
        return cin_network_plot
    else:
        return sting_network_plot


@callback(
    Output('cin_network_plot', 'data'),
    Output('sting_network_plot', 'data'),
    Input('submit_button', 'n_clicks'),
    State('cell_type', 'value'),
    State('ligands', 'value'),
    State('interaction_fdr', 'data'),
    State('min_logfc', 'data'),
    State('min_expression', 'data'),
    State('logfc_fdr', 'data'),
    State('iterations', 'value'),
    interval=500,
    cache_args_to_ignore=[0],
    background=True,  # Run in background
    prevent_initial_call=True,
    running=[  # Disable the button while the callback is running
        (Output('submit_button', 'disabled'), True, False),
        (Output('progress-bar', 'style'), {'visibility': 'visible'}, {'visibility': 'hidden'})
    ],
    progress=[  # Show a progress bar while the callback is running
        Output('progress-bar', "value"),
        Output('progress-bar', "max")
    ]
)
def make_graph(set_progress, n_clicks,
               cell_type, ligands,
               interaction_fdr, min_logfc, min_expression,
               logfc_fdr, iterations):
    from viz.figures import pseudotime_interaction_propagation_graph, DEFAULT_LIGAND_EFFECT_ARGS, LIGAND_EFFECT_SAVE_LOCATION
    from viz.data import using_custom_data

    set_progress((0, iterations))

    # Check if arguments match default, if so return the pre-computed default
    if (cell_type == DEFAULT_LIGAND_EFFECT_ARGS['cell_type'] and
        ligands == DEFAULT_LIGAND_EFFECT_ARGS['ligands'] and
        interaction_fdr == DEFAULT_LIGAND_EFFECT_ARGS['interaction_fdr'] and
        min_logfc == DEFAULT_LIGAND_EFFECT_ARGS['min_logfc'] and
        min_expression == DEFAULT_LIGAND_EFFECT_ARGS['min_expression'] and
        logfc_fdr == DEFAULT_LIGAND_EFFECT_ARGS['logfc_fdr'] and
        iterations == DEFAULT_LIGAND_EFFECT_ARGS['iterations']):
        import pickle
        try:
            if not using_custom_data():
                with open(LIGAND_EFFECT_SAVE_LOCATION, 'rb') as f:
                    # return [(Serverside(p) if p is not None else None) for p in pickle.load(f)]  FIXME
                    return pickle.load(f)
        except:
            pass

    fig = pseudotime_interaction_propagation_graph(
        effect_set='custom' if using_custom_data() else DEFAULT_LIGAND_EFFECT_ARGS['effect_set'],  # Keep this fixed
        seed_cell=cell_type,
        seed_ligands=ligands,
        iterations=int(iterations),
        interaction_fdr_cutoff=float("." + interaction_fdr[3:]),
        min_logfc=float(min_logfc),
        min_expression=float(min_expression),
        logfc_fdr_cutoff=float("." + logfc_fdr[3:]),
        set_progress_callback=set_progress
    )

    # return [Serverside(fig), None] FIXME
    return [fig, None]


# @callback(
#     Output('min_logfc', 'max'),
#     Input('effect_set', 'value'),
#     interval=10,
#     background=False,  # Run in background,
#     running=[  # Disable the button while the callback is running
#         (Output('submit_button', 'disabled'), True, False),
#         (Output('spinner-holder', 'children'), [
#             dbc.Spinner(color='primary',
#                         size='md',
#                         fullscreen=True,
#                         type='grow')],
#          []),
#     ]
# )
# def initialize_options(effect_set):
#     max_logfc = 2
#     return max_logfc


layout = [
    interactive_panel(html.Div(),  # wrap_icon('fa-maximize', 'Downstream Cascading Effects'),
                      html.Div(),  # "Predict the downstream effects of ligands across a microenvironment.",
                      *build_interface()
                      )
]


# If run as a script, compile the default plot
if __name__ == '__main__':
    from viz.figures import pseudotime_interaction_propagation_graph, DEFAULT_LIGAND_EFFECT_ARGS, LIGAND_EFFECT_SAVE_LOCATION
    import pickle

    fig = pseudotime_interaction_propagation_graph(
        effect_set=DEFAULT_LIGAND_EFFECT_ARGS['effect_set'],  # Calculate both sets if we decide to implement a selector
        seed_cell=DEFAULT_LIGAND_EFFECT_ARGS['cell_type'],
        seed_ligands=DEFAULT_LIGAND_EFFECT_ARGS['ligands'],
        iterations=int(DEFAULT_LIGAND_EFFECT_ARGS['iterations']),
        interaction_fdr_cutoff=float(DEFAULT_LIGAND_EFFECT_ARGS['interaction_fdr']),
        min_logfc=float(DEFAULT_LIGAND_EFFECT_ARGS['min_logfc']),
        min_expression=float(DEFAULT_LIGAND_EFFECT_ARGS['min_expression']),
        logfc_fdr_cutoff=float(DEFAULT_LIGAND_EFFECT_ARGS['logfc_fdr']),
        set_progress_callback=None
    )

    with open(LIGAND_EFFECT_SAVE_LOCATION, 'wb') as f:
        pickle.dump([fig, None], f)
