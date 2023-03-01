import os

import dash
from dash import html, callback, Output, Input, State, dcc
import dash_bootstrap_components as dbc

from viz.web import interactive_panel, wrap_icon, control_panel, figure_output, control_panel_element

if __name__ != '__main__':
    dash.register_page(__name__,
                       path='/circos',
                       name='Circos',
                       order=1)


def build_interface() -> list:
    from viz.figures import DEFAULT_CIRCOS_ARGS, CIRCOS_SAVE_LOCATION

    import pickle
    default_plot = None
    try:
        with open(CIRCOS_SAVE_LOCATION, 'rb') as f:
            default_plot = pickle.load(f)
    except:
        pass

    controls = control_panel(
        [
            control_panel_element("Interaction Effect FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  dbc.Select(
                                      id='inter_circos_fdr',
                                      options=[{'label': '0.05', 'value': 'fdr05'},
                                               {'label': '0.25', 'value': 'fdr25'}],
                                      value=DEFAULT_CIRCOS_ARGS['inter_circos_fdr']
                                  )),
            control_panel_element("log2FC FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  dbc.Select(
                                      id='logfc_circos_fdr',
                                      options=[{'label': '0.05', 'value': 'fdr05'},
                                               {'label': '0.25', 'value': 'fdr25'}],
                                      value=DEFAULT_CIRCOS_ARGS['logfc_circos_fdr']
                                  ))
        ], [
            control_panel_element("Outer Interaction Set", "Biological condition to compare.",
                                  dbc.Select(
                                      id='circos_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'STING-Dependent Effect', 'value': 'sting'}],
                                      value=DEFAULT_CIRCOS_ARGS['circos_set']
                                  ))
        ], [
            control_panel_element('Minimum numSigI1', 'Minimum number of significant interactions for a receptor to be included.',
                                  dcc.Slider(
                                      id='circos_min_numsigi1',
                                      min=0,
                                      max=10,  # Fill in later
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
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
                                      max=10,  # Fill in later
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                                      step=1,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className='form-range'
                                  ))
        ], [
            control_panel_element('Chord Ligand abs(log2FC) Cutoff', 'The minimum ligand log2FC required for a chord to be drawn.',
                                  dcc.Slider(
                                      id='circos_min_ligand_logfc',
                                      min=0,
                                      max=1,  # Fill in
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'],
                                      step=0.01,
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
            children=default_plot
        )
    )

    return [
        controls,
        results
    ]


@callback(
    Output('circos-graph-holder', 'children'),
    Input('submit-button-circos', 'n_clicks'),
    State('inter_circos_fdr', 'value'),
    State('logfc_circos_fdr', 'value'),
    State('circos_set', 'value'),
    State('circos_min_numsigi1', 'value'),
    State('circos_min_numdeg', 'value'),
    State('circos_min_ligand_logfc', 'value'),
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
def make_circos_plot(set_progress, n_clicks,
                     inter_circos_fdr, logfc_circos_fdr, circos_set,
                     min_numsigi1, min_numdeg, min_chord_ligand_logfc):
    if n_clicks == 0:
        from dash.exceptions import PreventUpdate
        raise PreventUpdate

    set_progress((0, 7))
    from viz.data import read_circos_file, read_ligand_receptor_file, read_interactions_file
    from viz.web import make_circos_figure
    from viz.figures import DEFAULT_CIRCOS_ARGS, CIRCOS_SAVE_LOCATION

    # Check if arguments match default, if so return the pre-computed default
    if (inter_circos_fdr == DEFAULT_CIRCOS_ARGS['inter_circos_fdr'] and
        logfc_circos_fdr == DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'] and
        circos_set == DEFAULT_CIRCOS_ARGS['circos_set'] and
        min_numsigi1 == DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'] and
        min_numdeg == DEFAULT_CIRCOS_ARGS['circos_min_numdeg'] and
        min_chord_ligand_logfc == DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc']):
        import pickle
        try:
            with open(CIRCOS_SAVE_LOCATION, 'rb') as f:
                return [pickle.load(f)]
        except:
            pass

    data = read_circos_file(circos_set, inter_circos_fdr)
    inter_data = read_interactions_file(circos_set, inter_circos_fdr)

    # Get distinct pairs of ligands and receptors
    lr_pairs = read_ligand_receptor_file()
    # Convert to list of tuples
    lr_pairs = [tuple(x) for x in lr_pairs.values]

    return [make_circos_figure(set_progress,
                               lr_pairs,
                               data,
                               inter_data,
                               logfc_circos_fdr,
                               min_numsigi1,
                               min_numdeg,
                               min_chord_ligand_logfc)]


@callback(
    Output('circos_min_numsigi1', 'max'),
    Output('circos_min_numdeg', 'max'),
    Output('circos_min_ligand_logfc', 'max'),
    Input('inter_circos_fdr', 'value'),
    Input('logfc_circos_fdr', 'value'),
    background=True,
    running=[
        (Output('submit-button-circos', 'disabled'), True, False),
        (Output('spinner-holder', 'children'), [
            dbc.Spinner(color='primary', size='md', fullscreen=True, type='grow')
        ], [])
    ]
)
def initialize_options(inter_circos_fdr, logfc_circos_fdr):
    from viz.data import read_circos_file

    fdr = inter_circos_fdr

    # Read the maximum values from the circos file to determine range of sliders
    max_obs = read_circos_file('highCIN_vs_noSTING', fdr)

    max_deg = max_obs['numDEG'].max()
    max_numsigi1 = max_obs['numSigI1'].max()
    max_logfc = max_obs['MAST_log2FC'].abs().max()

    return max_numsigi1, max_deg, max_logfc


layout = [
    interactive_panel(wrap_icon('fa-circle-dot', 'Circos Plot of Interactions'), *build_interface())
]


# If run as a script, compile the default plot
if __name__ == '__main__':
    from viz.figures import DEFAULT_CIRCOS_ARGS, CIRCOS_SAVE_LOCATION

    if not os.path.exists(CIRCOS_SAVE_LOCATION):
        from viz.data import read_circos_file, read_ligand_receptor_file, read_interactions_file
        from viz.web import make_circos_figure
        import pickle

        data = read_circos_file(DEFAULT_CIRCOS_ARGS['circos_set'], DEFAULT_CIRCOS_ARGS['inter_circos_fdr'])
        inter_data = read_interactions_file(DEFAULT_CIRCOS_ARGS['circos_set'], DEFAULT_CIRCOS_ARGS['inter_circos_fdr'])

        # Get distinct pairs of ligands and receptors
        lr_pairs = read_ligand_receptor_file()
        # Convert to list of tuples
        lr_pairs = [tuple(x) for x in lr_pairs.values]

        fig = make_circos_figure(None,
                                 lr_pairs,
                                 data,
                                 inter_data,
                                 DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'],
                                 DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
                                 DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                                 DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'])

        # Dump the figure to a pickle file
        with open(CIRCOS_SAVE_LOCATION, 'wb') as f:
            pickle.dump(fig, f)
