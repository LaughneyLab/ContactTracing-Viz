import math

import anndata as ad
import dash
import dash_bootstrap_components as dbc
from dash import html, callback, Output, Input, dcc, State
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import networkx as nx
from dash.exceptions import PreventUpdate
from matplotlib import colors, cm

from viz.data import get_diff_abundance, get_interaction_fdr, get_effective_logfc, get_downstream_ligands, read_ct_data, \
    read_interactions
from viz.figures import pseudotime_interaction_propagation_graph
from viz.util import enhance_plotly_export, celltype_to_colors, get_quiver_arrows
from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output

dash.register_page(__name__,
                   path='/ligand-effects',
                   name='Ligand Effects',
                   order=1)


def build_interface() -> list:
    controls = control_panel(
        [
            control_panel_element("Network Layout", 'abc',
                                  dbc.Select(
                                      id='network-layout',
                                      options=[
                                          {'label': 'Planar Layout', 'value': 'planar'},
                                          {'label': 'Spring Layout', 'value': 'spring'},
                                          {'label': 'Circular Layout', 'value': 'circular'},
                                          {'label': 'Timeline Layout', 'value': 'timeline'}
                                      ],
                                      persistence=True, persistence_type='session', value='timeline'
                                  )),
            control_panel_element('Emitting Cell Type', 'abc',
                                  dbc.Select(
                                    id='cell-type',
                                    options=[],  # Filled in by callback
                                    persistence=False
                                  )),
            control_panel_element("Emitted Ligands", 'abc',
                                  dbc.Input(
                                      id='ligands',
                                      autofocus=True,
                                      placeholder='Example: Ccl2,Apoe',
                                      persistence=False
                                  ))
        ], [
            control_panel_element("Minimum Expression", "abc",
                                  dcc.Slider(
                                      id='min_expression',
                                      min=0,
                                      max=1,
                                      step=0.01,
                                      value=0,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className='form-range'
                                  )),
            control_panel_element("Interaction FDR Cutoff", "abc",
                                  dcc.Slider(
                                      id='interaction_fdr',
                                      max=1,
                                      min=0,
                                      step=0.01,
                                      value=0.05,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=True, persistence_type='session',
                                      className='form-range'
                                  ))
        ], [
            control_panel_element("Minimum abs(Log2FC)", "abc",
                                  dcc.Slider(
                                      id='min_logfc',
                                      max=1,
                                      min=0,
                                      step=0.01,
                                      value=0,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className='form-range'
                                  )),
            control_panel_element("log2FC FDR Cutoff", "abc",
                                  dcc.Slider(
                                      id='logfc_fdr',
                                      max=1,
                                      min=0,
                                      step=0.01,
                                      value=0.05,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=True, persistence_type='session',
                                      className='form-range'
                                  ))
        ], [
            control_panel_element("Network Building Iterations", "abc",
                                  dbc.Input(
                                      id='iterations',
                                      type='number',
                                      debounce=True,
                                      max=100,
                                      min=1,
                                      step=1,
                                      value=10,
                                      persistence=True, persistence_type='session'
                                  )),
            control_panel_element("Plot", "abc",
                                  dbc.Button(
                                      "Submit",
                                      id='submit-button',
                                      size='lg',
                                      color="primary",
                                      className='me-1',
                                      n_clicks=0
                                  ))
        ]
    )

    results = figure_output(
        title='Ligand Network Figure',
        footer="This is a figure",
        element=dcc.Graph(id='network-graph',
                          animate=True,
                          #figure={},  # Filled in by callback
                          config={
                             'displaylogo': False,
                             'showTips': True,
                             'toImageButtonOptions': {
                                'format': 'png',  # one of png, svg, jpeg, webp
                                'filename': 'exported_image',
                                'height': 800,
                                'width': 1200,
                                'scale': 6  # Multiply title/legend/axis/canvas sizes by this factor
                             },
                             'watermark': False
                          })
    )

    return [
        controls,
        results
    ]


@callback(
    Output('network-graph', 'figure'),
    Input('data-session', 'data'),
    Input('submit-button', 'n_clicks'),
    State('network-layout', 'value'),
    State('cell-type', 'value'),
    State('ligands', 'value'),
    State('min_expression', 'value'),
    State('interaction_fdr', 'value'),
    State('min_logfc', 'value'),
    State('logfc_fdr', 'value'),
    State('iterations', 'value'),
    background=True,  # Run in background
    prevent_initial_call=True,
    running=[  # Disable the button while the callback is running
        (Output('submit-button', 'disabled'), True, False),
        (Output('progress-bar', 'style'), {'visibility': 'visible'}, {'visibility': 'hidden'})
    ],
    progress=[  # Show a progress bar while the callback is running
        Output('progress-bar', "value"),
        Output('progress-bar', "max")
    ]
)
def make_graph(set_progress, data, n_clicks, network_layout, cell_type, ligands, min_expression, interaction_fdr, min_logfc, logfc_fdr, iterations):
    if data is None or n_clicks == 0 or ligands is None:
        raise PreventUpdate
    set_progress((0, 100))

    # Do some work

    file = data['path']
    interaction_file = data['tsv']
    adata = read_ct_data(file)
    interactions = read_interactions(interaction_file)
    fig = pseudotime_interaction_propagation_graph(
        ct=adata,
        orig_df=interactions,
        seed_cell=cell_type,
        seed_ligands=ligands,
        iterations=int(iterations),
        interaction_fdr_cutoff=float(interaction_fdr),
        min_logfc=float(min_logfc),
        logfc_fdr_cutoff=float(logfc_fdr),
        min_expression=int(min_expression),
        layout=network_layout,
        set_progress_callback=set_progress
    )

    return fig


@callback(
    Output('cell-type', 'options'),
    Output('cell-type', 'value'),
    Output('min_expression', 'max'),
    Output('min_expression', 'value'),
    Output('min_logfc', 'max'),
    Output('min_logfc', 'value'),
    Input('data-session', 'data'),
    background=True,  # Run in background,
    running=[  # Disable the button while the callback is running
        (Output('submit-button', 'disabled'), True, False),
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
    print(data)
    if data is None:
        return [], '', \
               1, 0, \
               1, 0
    file = data['path']
    interaction_file = data['tsv']
 #   adata = read_ct_data(file)
    interactions = read_interactions(interaction_file)
    max_exp = max(max(interactions.expression_ligand), max(interactions.expression_receptor))
    max_logfc = max(max(interactions.MAST_log2FC_ligand.abs()), max(interactions.MAST_log2FC_receptor.abs()))
    celltypes = list(sorted(set(interactions.cell_type_ligand) | set(interactions.cell_type_receptor)))
  #  celltypes = list(adata.obs['cell type'].unique())
    return [{'label': ct, 'value': ct} for ct in celltypes], celltypes[0], \
           max_exp, 0, \
           max_logfc, 0


layout = [
    interactive_panel(wrap_icon('fa-maximize', 'Cell Type to Cell Type Interactions'),
                      *build_interface()
                      )
]
