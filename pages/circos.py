import dash
from dash_extensions import DeferScript
from dash_extensions.enrich import html, callback, Output, Input, State, dcc, clientside_callback, Serverside
import dash_bootstrap_components as dbc

from viz.docs import circos_help, interaction_effects_def, diffusion_component_def, differential_abundance_def, \
    deg_test_def, ligand_log2fc_def, conditions_def
from viz.web import interactive_panel, wrap_icon, control_panel, figure_output, control_panel_element, \
    make_custom_slider, make_fdr_slider

if __name__ != '__main__':
    dash.register_page(__name__,
                       path='/circos',
                       name='Circos',
                       order=1)


def build_interface() -> list:
    # Ignore if this file is run directly
    if __name__ == '__main__':
        return []

    from viz.figures import DEFAULT_CIRCOS_ARGS, CIRCOS_SAVE_LOCATION

    import pickle
    default_plots = None
    try:
        with open(CIRCOS_SAVE_LOCATION, 'rb') as f:
            default_plots = pickle.load(f)
    except:
        pass

    controls = control_panel("submit-button-circos",
         [
             control_panel_element("Interaction Set",
                                   [
                                       "The ",
                                       *conditions_def("biological condition"),
                                       " to compare."
                                   ],
                                   dbc.RadioItems(
                                       id='circos_set',
                                       options=[{'label': 'CIN-Dependent Interactions', 'value': 'cin'},
                                                {'label': 'CIN/STING-Dependent Interactions',
                                                 'value': 'sting'}],
                                       value=DEFAULT_CIRCOS_ARGS['circos_set'],
                                       persistence=False
                                   )),
            control_panel_element("Genes of Interest",
                                  'Comma-separated list of genes to highlight in the plot.',
                                  dbc.Input(
                                      id='genes',
                                      autofocus=False,
                                      value=DEFAULT_CIRCOS_ARGS['genes'],
                                      placeholder='Example: Ccl2,Apoe',
                                  )),
             control_panel_element("Plot", "You must click submit to update the plot.",
                                   html.Div(
                                       dbc.Row([dbc.Col(
                                           dbc.Button(
                                               wrap_icon("fas fa-play", "Submit"),
                                               id="submit-button-circos",
                                               size="lg",
                                               color="primary",
                                               className='me-1',
                                               n_clicks=0),
                                           width='auto', align='left'),
                                           dbc.Col(
                                               dbc.Button(wrap_icon('fas fa-rotate-left', 'Reset'),
                                                          id='reset_default_circos_button',
                                                          color='secondary', size='lg',
                                                          className='float-end dark', n_clicks=0),
                                               align='right', width='auto')
                                       ]),
                                       className='text-center d-grid gap-2'))
         ], [
            control_panel_element("Interaction Effect FDR Cutoff",
                                  [
                                      "FDR-adjusted requirements for ",
                                      *interaction_effects_def("interaction effects"),
                                      "."
                                  ],
                                  make_fdr_slider('inter_circos_fdr', DEFAULT_CIRCOS_ARGS['inter_circos_fdr'])),
            control_panel_element("Ligand Log2FC FDR Cutoff",
                                  [
                                      "FDR-adjusted requirements for ",
                                      *ligand_log2fc_def("ligand differential expression"),
                                      "."
                                  ],
                                  make_fdr_slider('logfc_circos_fdr', DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'])),
        ], [
             control_panel_element('Minimum Interaction Effect',
                                   [
                                       'Minimum number of significant ',
                                       *interaction_effects_def("interaction effects"),
                                       ' for a receptor to be included.'
                                   ],
                                   make_custom_slider(
                                       id='circos_min_numsigi1',
                                       min=0,
                                       max=100,  # Fill in later
                                       value=DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
                                       step=1
                                   )),
             control_panel_element('Minimum Affected Downstream Genes',
                                   [
                                       'Minimum number of ',
                                       *deg_test_def('differentially expressed genes conditioned on a target gene'),
                                       '.'
                                   ],
                                   make_custom_slider(
                                       id='circos_min_numdeg',
                                       min=0,
                                       max=100,  # Fill in later
                                       value=DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                                       step=1
                                   ))
        ], [
            control_panel_element('Chord Minimum Ligand abs(Log2FC)',
                                  [
                                      'The minimum ',
                                       *ligand_log2fc_def('ligand log2FC'),
                                       ' required for a chord to be drawn.'
                                  ],
                                  make_custom_slider(
                                      id='circos_min_ligand_logfc',
                                      min=0,
                                      max=2,  # Fill in
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'],
                                      step=0.01
                                  ))
        ]
    )

    results = figure_output(
        title="Interactive Circos Plot",
        footer=[
            "Layers from outside ring to center: Cell Type, ",
            *diffusion_component_def("Diffusion Component 1"),
            ", ",
            *differential_abundance_def("Differential Abundance"),
            ", Number of significant ",
            *interaction_effects_def("interaction effects"),
            ", ligand/receptor interactions. Tip: You can copy the information you are hovering over by hitting Ctrl+C and clear highlighting/zoom by double-clicking."
        ],
        element=html.Div([
            circos_toolbar(),
            html.Div(
                id="circos-graph-holder",
                children=default_plots[1] if default_plots is not None else None
            )
        ]),
        help_info=circos_help()
    )

    return [
        results,
        html.Br(),
        html.Hr(),
        html.Br(),
        controls,
        dcc.Store(id='cin_circos_plot', storage_type='memory', data=[default_plots[0] if default_plots is not None else None]),
        dcc.Store(id='sting_circos_plot', storage_type='memory', data=[default_plots[1] if default_plots is not None else None]),
        # DeferScript(src=dash.get_asset_url('circos_hooks.js'))
    ]


def circos_toolbar():
    reset_button = dbc.Button(wrap_icon('fas fa-undo', 'Reset'), outline=True, color='dark', id='circos-reset-button')
    zoom_in_button = dbc.Button(wrap_icon('fas fa-search-plus', 'Zoom In'), outline=True, color='dark', id='circos-zoom-in-button')
    zoom_out_button = dbc.Button(wrap_icon('fas fa-search-minus', 'Zoom Out'), outline=True, color='dark', id='circos-zoom-out-button')
    download_button = dbc.Button(wrap_icon('fas fa-download', 'Download'), outline=True, color='dark', id='circos-download-button')
    return html.Div(
        dbc.ButtonGroup(
            [reset_button, zoom_in_button, zoom_out_button, download_button],
            size='sm'
        ), id='circos-toolbar'
    )


# When the graph is updated, inject custom hooks
clientside_callback("""
function(_) {
    circosInjection();
}
""", Input('circos-graph-holder', 'children'))

clientside_callback("""
function(n_clicks) {
    if (n_clicks > 0) {
        downloadCircosSvg();
    }
}
""", Input('circos-download-button', 'n_clicks'), prevent_initial_call=True)

clientside_callback("""
function(n_clicks) {
    if (n_clicks > 0) {
        resetCircosTransform();
    }
}
""", Input('circos-reset-button', 'n_clicks'), prevent_initial_call=True)

clientside_callback("""
function(n_clicks) {
    if (n_clicks > 0) {
        zoomInCircos();
    }
}
""", Input('circos-zoom-in-button', 'n_clicks'), prevent_initial_call=True)

clientside_callback("""
function(n_clicks) {
    if (n_clicks > 0) {
        zoomOutCircos();
    }
}
""", Input('circos-zoom-out-button', 'n_clicks'), prevent_initial_call=True)

clientside_callback("""
function(genes, circos_set) {
    highlightGenes(genes);
}""", Input('genes', 'value'), Input('circos-graph-holder', 'children'), prevent_initial_call=True)


@callback(
    Output('inter_circos_fdr', 'data', allow_duplicate=True),
    Output('logfc_circos_fdr', 'data', allow_duplicate=True),
    Output('circos_min_numsigi1', 'data', allow_duplicate=True),
    Output('circos_min_numdeg', 'data', allow_duplicate=True),
    Output('circos_min_ligand_logfc', 'data', allow_duplicate=True),
    Output('genes', 'value', allow_duplicate=True),
    Output('submit-button-circos', 'n_clicks'),
    Input('reset_default_circos_button', 'n_clicks'),
    State('submit-button-circos', 'n_clicks'),
    prevent_initial_call=True
)
def reset_to_defaults(n_clicks, submission_button_clicks):
    from dash.exceptions import PreventUpdate
    from viz.figures import DEFAULT_CIRCOS_ARGS
    if n_clicks > 0:
        return DEFAULT_CIRCOS_ARGS['inter_circos_fdr'], \
                DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'], \
                DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'], \
                DEFAULT_CIRCOS_ARGS['circos_min_numdeg'], \
                DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'], \
                DEFAULT_CIRCOS_ARGS['genes'], \
                submission_button_clicks+1  # Trigger the plotting callback to run on default values
    raise PreventUpdate


@callback(
    Output('circos-graph-holder', 'children'),
    Input('circos_set', 'value'),
    Input('cin_circos_plot', 'data'),
    Input('sting_circos_plot', 'data')
)
def update_circos_plot(circos_set, cin_circos_plot, sting_circos_plot):
    if circos_set == 'cin':
        return cin_circos_plot
    else:
        return sting_circos_plot


@callback(
    Output('cin_circos_plot', 'data'),
    Output('sting_circos_plot', 'data'),
    Input('submit-button-circos', 'n_clicks'),
    State('inter_circos_fdr', 'data'),
    State('logfc_circos_fdr', 'data'),
    State('circos_min_numsigi1', 'data'),
    State('circos_min_numdeg', 'data'),
    State('circos_min_ligand_logfc', 'data'),
    # State('genes', 'value'),
    background=True,
    prevent_initial_call=True,
    interval=500,
    cache_args_to_ignore=[0],
    running=[
        (Output('submit-button-circos', 'disabled'), True, False),
        (Output('progress-bar', 'style'), {'visibility': 'visible'}, {'visibility': 'hidden'}),
    ],
    progress=[
        Output('progress-bar', 'value'),
        Output('progress-bar', 'max'),
    ]
)
def make_circos_plot(set_progress, n_clicks,
                     inter_circos_fdr, logfc_circos_fdr,
                     min_numsigi1, min_numdeg, min_chord_ligand_logfc #, gene_list
                     ):
    #if gene_list is not None and len(gene_list.strip()) == 0:
    #    gene_list = None
    gene_list = None  # TODO: Remove, gene list is now implemented client-side

    set_progress((0, 14))
    from viz.data import read_circos_file, read_interactions_file
    from viz.web import make_circos_figure
    from viz.figures import DEFAULT_CIRCOS_ARGS, CIRCOS_SAVE_LOCATION

    # Check if arguments match default, if so return the pre-computed default
    if (inter_circos_fdr == DEFAULT_CIRCOS_ARGS['inter_circos_fdr'] and
        logfc_circos_fdr == DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'] and
        min_numsigi1 == DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'] and
        min_numdeg == DEFAULT_CIRCOS_ARGS['circos_min_numdeg'] and
        min_chord_ligand_logfc == DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'] and
        gene_list == DEFAULT_CIRCOS_ARGS['genes']):
        import pickle
        try:
            with open(CIRCOS_SAVE_LOCATION, 'rb') as f:
                # return [(Serverside(p) if p is not None else None) for p in pickle.load(f)]  FIXME
                return pickle.load(f)
        except:
            pass

    #data = read_circos_file(circos_set, inter_circos_fdr)
    data = read_circos_file('sting', inter_circos_fdr)  # Keep outer data consistent
    cin_inter_data = read_interactions_file('cin', inter_circos_fdr)
    sting_inter_data = read_interactions_file('sting', inter_circos_fdr)

    cin_circos = make_circos_figure(set_progress, 0,
                               data,
                               cin_inter_data,
                               True,
                               logfc_circos_fdr,
                               min_numsigi1,
                               min_numdeg,
                               min_chord_ligand_logfc,
                               gene_list)

    sting_circos = make_circos_figure(set_progress, 1,
                               data,
                               sting_inter_data,
                               False,
                               logfc_circos_fdr,
                               min_numsigi1,
                               min_numdeg,
                               min_chord_ligand_logfc,
                               gene_list)

    set_progress((14, 14))

    # return [Serverside(cin_circos), Serverside(sting_circos)]  FIXME
    return [cin_circos, sting_circos]


# @callback(
#     Output('circos_min_numsigi1', 'max'),
#     Output('circos_min_numdeg', 'max'),
#     Output('circos_min_ligand_logfc', 'max'),
#     Input('inter_circos_fdr', 'data'),
#     Input('logfc_circos_fdr', 'data'),
#     interval=10,
#     background=False,
#     running=[
#         (Output('submit-button-circos', 'disabled'), True, False),
#         (Output('spinner-holder', 'children'), [
#             dbc.Spinner(color='primary', size='md', fullscreen=True, type='grow')
#         ], [])
#     ]
# )
# def initialize_options(inter_circos_fdr, logfc_circos_fdr):
#     from viz.data import read_circos_file
#
#     #fdr = inter_circos_fdr
#
#     # Read the maximum values from the circos file to determine range of sliders
#     #max_obs = read_circos_file('highCIN_vs_noSTING', 'fdr25')
#
#     max_deg = 250 #  max_obs['numDEG'].max()
#     max_numsigi1 = 250 #  max_obs['numSigI1'].max()
#     max_logfc = 2 #  max_obs['MAST_log2FC'].abs().max()
#
#     return max_numsigi1, max_deg, max_logfc


layout = [
    interactive_panel(wrap_icon('fa-circle-dot', 'Circos Plot of Interactions'),
                      "Visualize all interactions across an experiment.",
                      *build_interface())
]


# If run as a script, compile the default plot
if __name__ == '__main__':
    from viz.figures import DEFAULT_CIRCOS_ARGS, CIRCOS_SAVE_LOCATION
    from viz.data import read_circos_file, read_interactions_file
    from viz.web import make_circos_figure
    import pickle

    data = read_circos_file('sting', DEFAULT_CIRCOS_ARGS['inter_circos_fdr'])
    cin_inter_data = read_interactions_file('cin', DEFAULT_CIRCOS_ARGS['inter_circos_fdr'])
    sting_inter_data = read_interactions_file('sting', DEFAULT_CIRCOS_ARGS['inter_circos_fdr'])
    cin_fig = make_circos_figure(None, 0,
                             data,
                             cin_inter_data,
                             True,
                             DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'],
                             DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
                             DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                             DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'],
                             DEFAULT_CIRCOS_ARGS['genes'])
    sting_fig = make_circos_figure(None, 1,
                                 data,
                                 sting_inter_data,
                                 False,
                                 DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'],
                                 DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
                                 DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                                 DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'],
                                 DEFAULT_CIRCOS_ARGS['genes'])

    # Dump the figure to a pickle file
    with open(CIRCOS_SAVE_LOCATION, 'wb') as f:
        pickle.dump([cin_fig, sting_fig], f)
