import dash
from dash import html, callback, Output, Input, State, dcc
import dash_bootstrap_components as dbc

from viz.web import interactive_panel, wrap_icon, control_panel, figure_output, control_panel_element

dash.register_page(__name__,
                   path='/circos',
                   name='Circos',
                   order=1)


def build_interface() -> list:
    controls = control_panel(
        [
            control_panel_element("Interaction Effect FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  dbc.Select(
                                      id='inter_circos_fdr',
                                      options=[{'label': '0.05', 'value': 'fdr05'},
                                               {'label': '0.25', 'value': 'fdr25'}],
                                      value='fdr25'  # Default from fig 4
                                  )),
            control_panel_element("log2FC FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  dbc.Select(
                                      id='logfc_circos_fdr',
                                      options=[{'label': '0.05', 'value': 'fdr05'},
                                               {'label': '0.25', 'value': 'fdr25'}],
                                      value='fdr05'  # Default from fig 4
                                  ))
        ], [
            control_panel_element("Outer Interaction Set", "Biological condition to compare.",
                                  dbc.Select(
                                      id='circos_outer_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'STING-Dependent Effect', 'value': 'sting'}],
                                      value='cin'  # Default from fig 4
                                  )),
            control_panel_element("Inner Interaction Set", "Biological condition to compare.",
                                  dbc.Select(
                                      id='circos_inner_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'STING-Dependent Effect', 'value': 'sting'}],
                                      value='sting'  # Default from fig 4
                                  ))
        ], [
            control_panel_element('Minimum numSigI1', 'Minimum number of significant interactions for a receptor to be included.',
                                  dcc.Slider(
                                      id='circos_min_numsigi1',
                                      min=0,
                                      max=10,  # Fill in later
                                      value=1,  # Default from fig 4
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
                                      value=0,    # Default from fig 4
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
                                      value=0.12,  # Default from fig 4
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
            children=[]
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
    State('circos_outer_set', 'value'),
    State('circos_inner_set', 'value'),
    State('circos_min_numsigi1', 'value'),
    State('circos_min_numdeg', 'value'),
    State('circos_min_ligand_logfc', 'value'),
    background=True,
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
                     inter_circos_fdr, logfc_circos_fdr, outer_set, inner_set,
                     min_numsigi1, min_numdeg, min_chord_ligand_logfc):
    set_progress((0, 7))
    from viz.data import read_circos_file, read_ligand_receptor_file
    from viz.web import make_circos_figure

    outer_data = read_circos_file(outer_set, inter_circos_fdr)
    if inner_set == outer_set:
        max_data = None
        is_max = False
    else:
        max_data = read_circos_file('max', inter_circos_fdr)
        is_max = True


    # Get distinct pairs of ligands and receptors
    lr_pairs = read_ligand_receptor_file()
    # Convert to list of tuples
    lr_pairs = [tuple(x) for x in lr_pairs.values]

    return [make_circos_figure(set_progress,
                               lr_pairs,
                               outer_data,
                               max_data,
                               is_max,
                               inter_circos_fdr,
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
    max_obs = read_circos_file('max', fdr)

    max_deg = max_obs['numDEG'].max()
    max_numsigi1 = max_obs['numSigI1'].max()
    max_logfc = max_obs['MAST_log2FC'].abs().max()

    return max_numsigi1, max_deg, max_logfc


layout = [
    interactive_panel(wrap_icon('fa-circle-dot', 'Circos Plot of Interactions'), *build_interface())
]

