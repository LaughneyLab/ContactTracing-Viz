import dash
import dash_bootstrap_components as dbc
from dash import html, callback, Output, Input, dcc, State

from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output

dash.register_page(__name__,
                   path='/ligand-effects',
                   name='Ligand Effects',
                   order=3)


def build_interface() -> list:
    controls = control_panel(
        [
            control_panel_element("Interaction Set", "Biological condition to compare.",
                                  dbc.Select(
                                      id='effect_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'STING-Dependent Effect', 'value': 'sting'},
                                               {'label': 'CIN & STING Max Effect', 'value': 'max'}],
                                      value='cin'
                                  )),
            control_panel_element("Network Layout", 'Select how you would like to structure nodes.',
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
            control_panel_element('Emitting Cell Type', 'Select the initial cell type.',
                                  dbc.Select(
                                    id='cell-type',
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
                                    value='Tumor cells'
                                  )),
            control_panel_element("Emitted Ligands", 'List of ligands to emit from the selected cell type.',
                                  dbc.Input(
                                      id='ligands',
                                      autofocus=True,
                                      value='Ccl2,Apoe',
                                      placeholder='Example: Ccl2,Apoe'
                                  ))
        ], [
            control_panel_element("Interaction FDR Cutoff", "The maximum interaction test FDR to consider.",
                                  dcc.Slider(
                                      id='interaction_fdr',
                                      max=1,
                                      min=0,
                                      step=0.01,
                                      value=0.05,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=True,
                                      persistence_type='session',
                                      className='form-range'
                                  ))
        ], [
            control_panel_element("Minimum abs(Log2FC)", "The minimum induced log2FC for targets.",
                                  dcc.Slider(
                                      id='min_logfc',
                                      max=1,
                                      min=0,
                                      step=0.01,
                                      value=0.1,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className='form-range'
                                  )),
            control_panel_element("log2FC FDR Cutoff", "The FDR-adjusted cutoff for determining if a log2FC value is non-zero.",
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
            control_panel_element("Network Building Iterations", "This controls how many downstream interactions may be detected.",
                                  dbc.Input(
                                      id='iterations',
                                      type='number',
                                      debounce=True,
                                      max=5,
                                      min=1,
                                      step=1,
                                      value=3,
                                      persistence=True, persistence_type='session'
                                  )),
            control_panel_element("Plot", "",
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
        footer="Circle = Ligand, Square = Receptor, Diamond = Ligand and Receptor",
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
    Input('submit-button', 'n_clicks'),
    State('effect_set', 'value'),
    State('network-layout', 'value'),
    State('cell-type', 'value'),
    State('ligands', 'value'),
    State('interaction_fdr', 'value'),
    State('min_logfc', 'value'),
    State('logfc_fdr', 'value'),
    State('iterations', 'value'),
    background=True,  # Run in background
    running=[  # Disable the button while the callback is running
        (Output('submit-button', 'disabled'), True, False),
        (Output('progress-bar', 'style'), {'visibility': 'visible'}, {'visibility': 'hidden'})
    ],
    progress=[  # Show a progress bar while the callback is running
        Output('progress-bar', "value"),
        Output('progress-bar', "max")
    ]
)
def make_graph(set_progress, n_clicks,
               effect_set, network_layout,
               cell_type, ligands,
               interaction_fdr, min_logfc,
               logfc_fdr, iterations):
    from viz.figures import pseudotime_interaction_propagation_graph

    set_progress((0, 100))

    fig = pseudotime_interaction_propagation_graph(
        effect_set=effect_set,
        seed_cell=cell_type,
        seed_ligands=ligands,
        iterations=int(iterations),
        interaction_fdr_cutoff=float(interaction_fdr),
        min_logfc=float(min_logfc),
        logfc_fdr_cutoff=float(logfc_fdr),
        layout=network_layout,
        set_progress_callback=set_progress
    )

    return fig


@callback(
    Output('min_logfc', 'max'),
    Input('effect_set', 'value'),
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
def initialize_options(effect_set):
    # TODO Should we scan to find true max logfc?
    max_logfc = 2
    return max_logfc


layout = [
    interactive_panel(wrap_icon('fa-maximize', 'Downstream Ligand Effects'),
                      *build_interface()
                      )
]
