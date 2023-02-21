import math

import dash
import networkx as nx
import numpy as np
from dash import html, dcc, callback, Output, Input, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
from matplotlib import cm, colors

from viz.data import read_interactions
from viz.figures import bipartite_graph
from viz.util import multipartite_layout, get_quiver_arrows
from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output

dash.register_page(__name__,
                   path='/interactions',
                   name='Interactions',
                   order=1)


def build_interface() -> list:
    controls = control_panel(
        [
            control_panel_element("Minimum numSigI1", "The minimum number of significant target gene interactions.",
                                  dcc.Slider(
                                      id="min_numsigi1_bipartite",
                                      min=0,
                                      max=10,  # Fill in
                                      step=1,
                                      value=1,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className="form-range"
                                  )),
            control_panel_element("Minimum Expression", "The minimum expression of each protein.",
                                  dcc.Slider(
                                      id="min_expression_bipartite",
                                      min=0,
                                      max=1,
                                      step=0.1,
                                      value=0,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className="form-range"
                                  ))
        ], [
            control_panel_element("Minimum abs(log2FC)", "The minimum log2FC for either ligands or receptors between conditions.",
                                  dcc.Slider(
                                      id="min_logfc_bipartite",
                                      min=0,
                                      max=1,  # Fill in
                                      step=0.1,
                                      value=0,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className="form-range"
                                  )),
            control_panel_element("log2FC FDR Cutoff", "The FDR-corrected p-value cutoff for whether the log2FC is not 0.",
                                  dcc.Slider(
                                      id="logfc_fdr_bipartite",
                                      min=0,
                                      max=1,
                                      step=0.01,
                                      value=0.05,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=True, persistence_type='session',
                                      className="form-range"
                                  ))
        ], [
            control_panel_element("First Cell Type", "The first cell type to examine interactions between.",
                                  dbc.Select(
                                      id='first_celltype',
                                      options=[],  # Filled in by callback
                                      persistence=False
                                  )),
            control_panel_element("Second Cell Type", "The second cell type to examine interactions between (need not be unique).",
                                  dbc.Select(
                                      id='second_celltype',
                                      options=[],  # Filled in by callback
                                      persistence=False
                                  )),
            control_panel_element("Third Cell Type", "If specified, include interactions across a third cell type.",
                                  dbc.Select(
                                      id='third_celltype',
                                      options=['(None)'],  # Filled in by callback
                                      value='(None)',
                                      persistence=False
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
                          animate=True,
                          #figure={},  # Filled in by callback
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
    Input('data-session', 'data'),
    Input('submit-button-bipartite', 'n_clicks'),
    State('first_celltype', 'value'),
    State('second_celltype', 'value'),
    State('third_celltype', 'value'),
    State('min_expression_bipartite', 'value'),
    State('min_logfc_bipartite', 'value'),
    State('logfc_fdr_bipartite', 'value'),
    State('min_numsigi1_bipartite', 'value'),
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
def make_graph(set_progress, data, n_clicks, first_ct, second_ct, third_ct, min_expression_bipartite, min_logfc_bipartite, logfc_fdr_bipartite, min_numSigI1):
    if data is None or n_clicks == 0:
        raise PreventUpdate
    set_progress((0, 100))

    # Do some work

    file = data['path']
    interaction_file = data['tsv']
    #adata = read_ct_data(file)
    interactions = read_interactions(interaction_file)
    fig = bipartite_graph(
        df=interactions,
        cell1=first_ct,
        cell2=second_ct,
        cell3=third_ct,
        numInteractions=min_numSigI1,
        min_logfc_bipartite=min_logfc_bipartite,
        logfc_fdr_bipartite_cutoff=logfc_fdr_bipartite
    )

    return fig


@callback(
    Output('first_celltype', 'options'),
    Output('first_celltype', 'value'),
    Output('second_celltype', 'options'),
    Output('second_celltype', 'value'),
    Output('third_celltype', 'options'),
    Output('third_celltype', 'value'),
    Output('min_expression_bipartite', 'max'),
    Output('min_expression_bipartite', 'value'),
    Output('min_logfc_bipartite', 'max'),
    Output('min_logfc_bipartite', 'value'),
    Output('min_numsigi1_bipartite', 'max'),
    Output('min_numsigi1_bipartite', 'value'),
    Input('data-session', 'data'),
    background=True,  # Run in background,
    running=[  # Disable the button while the callback is running
        (Output('submit-button-bipartite', 'disabled'), True, False),
        (Output('spinner-holder', 'children'), [
            dbc.Spinner(color='primary',
                        size='md',
                        fullscreen=True,
                        type='grow')],
         []),
    ]
)
def initialize_options(data):
    from viz.data import read_ct_data, read_interactions
    # print(data)
    if data is None:
        return [], '', \
               [], '', \
               [{'label': '(None)', 'value': '(None)'}], '(None)', \
               1, 0, \
               1, 0, \
               10, 0
    file = data['path']
    interaction_file = data['tsv']
    #adata = read_ct_data(file)
    interactions = read_interactions(interaction_file)
    max_exp = max(max(interactions.expression_ligand), max(interactions.expression_receptor))
    max_logfc = max(max(interactions.MAST_log2FC_ligand.abs()), max(interactions.MAST_log2FC_receptor.abs()))
    max_inter = max(interactions.numSigI1_fdr_receptor)
    celltypes = list(sorted(set(interactions.cell_type_ligand) | set(interactions.cell_type_receptor)))
    return [{'label': ct, 'value': ct} for ct in celltypes], celltypes[0], \
           [{'label': ct, 'value': ct} for ct in celltypes], celltypes[1], \
           [{'label': '(None)', 'value': '(None)'}] + [{'label': ct, 'value': ct} for ct in celltypes], '(None)', \
           max_exp, 0, \
           max_logfc, 0, \
           max_inter, 0


layout = [
    interactive_panel(wrap_icon('fa-arrows-left-right-to-line', 'Cell Type to Cell Type Interactions'),
                      *build_interface()
                      )
]
