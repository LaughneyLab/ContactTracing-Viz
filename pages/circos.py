from typing import Tuple, List

import dash
import numpy as np
import pandas as pd
from dash import html, callback, Output, Input, State, dcc
import dash_bio as dashbio
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate

from viz.data import read_interactions
from viz.util import celltype_to_colors, ColorTransformer
from viz.web import interactive_panel, wrap_icon, control_panel, figure_output, control_panel_element, get_obs_columns

dash.register_page(__name__,
                   path='/circos',
                   name='Circos',
                   order=2)


def build_interface() -> list:
    controls = control_panel(
        [
            control_panel_element("FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  dbc.Select(
                                      id='circos_fdr',
                                      options=[],  # Fill in later
                                      persistence=False
                                  )),
            control_panel_element("Interaction Set", "Biological condition to compare.",
                                  dbc.Select(
                                      id='circos_interaction_set',
                                      options=[],  # Fill in later
                                      persistence=False
                                  ))
        ], [
            control_panel_element('Minimum numSigI1', 'Minimum number of significant interactions for a receptor to be included.',
                                  dcc.Slider(
                                      id='circos_min_numsigi1',
                                      min=0,
                                      max=10,  # Fill in
                                      step=1,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className='form-range'
                                  )),
            control_panel_element('Minimum numDEG', 'Minimum number of differentially expressed genes conditioned on a target gene.',
                                  dcc.Slider(
                                      id='circos_min_numdeg',
                                      min=0,
                                      max=10,  # Fill in
                                      step=1,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className='form-range'
                                  ))
        ], [
            control_panel_element('Chord Minimum Ligand abs(log2FC)', 'The minimum ligand log2FC required for a chord to be drawn.',
                                  dcc.Slider(
                                      id='circos_min_ligand_logfc',
                                      min=0,
                                      max=1,  # Fill in
                                      step=0.1,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className='form-range'
                                  )),
            control_panel_element('Chord Minimum Receptor numSigI1', 'The minimum number of significant interactions for a chord to be drawn.',
                                  dcc.Slider(
                                      id='circos_min_receptor_numsigi1',
                                      min=0,
                                      max=10,  # Fill in
                                      step=1,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className='form-range'
                                  ))
        ], [
            control_panel_element("Plot", "",
                                  dbc.Button(
                                      "Submit",
                                      id="submit-button-circos",
                                      size="lg",
                                      color="primary",
                                      className='me-1',
                                      n_clicks=0
                                  ))
        ]
    )

    results = figure_output(
        title="Interactive Circos Plot",
        footer="Layers from outside ring to center: Cell Type, Diffusion Component Value, Differential Abundance, Number of significant interactions, Strongest ligand/receptor interactions",
        element=html.Div(
            id="circos-graph-holder",
            children=[]
        )
    )

    return [
        controls,
        results
    ]


@callback(
    Output('circos-graph-holder', 'children'),
    Input('data-session', 'data'),
    Input('submit-button-circos', 'n_clicks'),
    State('circos_fdr', 'value'),
    State('circos_interaction_set', 'value'),
    State('circos_min_numsigi1', 'value'),
    State('circos_min_numdeg', 'value'),
    State('circos_min_ligand_logfc', 'value'),
    State('circos_min_receptor_numsigi1', 'value'),
    background=True,
    prevent_initial_call=True,
    running=[
        (Output('submit-button-circos', 'disabled'), True, False),
        (Output('progress-bar', 'style'), {'visibility': 'visible'}, {'visibility': 'hidden'}),
    ],
    progress=[
        Output('progress-bar', 'value'),
        Output('progress-bar', 'max'),
    ]
)
def make_circos_plot(set_progress, data, n_clicks,
                     fdr_cutoff, interaction_set, min_numsigi1, min_numdeg, min_chord_ligand_logfc, min_chord_numsigi1):
    if data is None or n_clicks == 0:
        raise PreventUpdate
    set_progress((0, 7))
    from viz.data import read_interactions, read_obs
    from viz.web import make_circos_figure
    file = data['path']
    interaction_file = data['tsv']
    obs_file = data['obs']

    obs = read_obs(obs_file)
    interactions = read_interactions(interaction_file)
    # Get distinct pairs of ligands and receptors
    lr_pairs = interactions[['ligand', 'receptor']].drop_duplicates()
    # Convert to list of tuples
    lr_pairs = [tuple(x) for x in lr_pairs.values]

    return [make_circos_figure(set_progress, obs, lr_pairs, interaction_set, fdr_cutoff, min_numsigi1, min_numdeg, min_chord_ligand_logfc, min_chord_numsigi1)]


@callback(
    Output('circos_fdr', 'options'),
    Output('circos_fdr', 'value'),
    Output('circos_interaction_set', 'options'),
    Output('circos_interaction_set', 'value'),
    Output('circos_min_numsigi1', 'max'),
    Output('circos_min_numsigi1', 'value'),
    Output('circos_min_numdeg', 'max'),
    Output('circos_min_numdeg', 'value'),
    Output('circos_min_ligand_logfc', 'max'),
    Output('circos_min_ligand_logfc', 'value'),
    Output('circos_min_receptor_numsigi1', 'max'),
    Output('circos_min_receptor_numsigi1', 'value'),
    Input('data-session', 'data'),
    background=True,
    running=[
        (Output('submit-button-circos', 'disabled'), True, False),
        (Output('spinner-holder', 'children'), [
            dbc.Spinner(color='primary', size='md', fullscreen=True, type='grow')
        ], [])
    ]
)
def initialize_options(data):
    if data is None:
        return [{'label': '0.05', 'value': 'fdr05'}], 'fdr05', \
               [{'label': '', 'value': ''}], '', \
                10, 0, \
                10, 0, \
                10, 0, \
                10, 0

    import pandas as pd
    from viz.data import read_obs
    from viz.web import get_obs_columns
    file = data['path']
    interaction_file = data['tsv']
    obs_file = data['obs']
    obs = read_obs(obs_file)

    fdr_values = [
        {'label': '0.05', 'value': 'fdr05'}
    ]
    if len(obs.columns[obs.columns.str.contains('fdr25')]) > 0:
        fdr_values = [{
            'label': '0.25',
            'value': 'fdr25'
        }] + fdr_values
    default_fdr = fdr_values[0]['value']

    conditions = list(set([x.replace('numSigI1', '').replace('_fdr25', '').replace('_fdr05', '').replace('_fdr', '') for x in obs.columns if x.startswith('numSigI1')]))
    condition_names = [c.replace('_', ' ').strip().replace("max", "Intersection Max") for c in conditions]
    condition = conditions[0] if '_max' not in conditions else '_max'

    # columns
    cell_type, target, receptor, ligand, numDEG, numSigI1, MAST_log2FC, MAST_fdr, cell_type_dc1, DA_score = get_obs_columns(default_fdr, condition)

    return fdr_values, default_fdr, \
        [{'label': name, 'value': cond} for name, cond in zip(condition_names, conditions)], condition, \
        obs[obs[receptor]][numSigI1].max(), 0, \
        obs[numDEG].max(), 0, \
        obs[obs[ligand]][MAST_log2FC].abs().max(), obs[obs[ligand]][MAST_log2FC].abs().max()*.1, \
        obs[obs[receptor]][numSigI1].max(), obs[obs[receptor]][numSigI1].max()*.1


layout = [
    interactive_panel(wrap_icon('fa-circle-dot', 'Circos Plot of Interactions'), *build_interface())
]

