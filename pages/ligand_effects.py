import os.path

import dash
import dash_bootstrap_components as dbc
from dash import callback, Output, Input, dcc, State

from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output

if __name__ != '__main__':
    dash.register_page(__name__,
                       path='/ligand-effects',
                       name='Ligand Effects',
                       order=3)


def build_interface() -> list:
    from viz.figures import DEFAULT_LIGAND_EFFECT_ARGS, LIGAND_EFFECT_SAVE_LOCATION

    import pickle
    default_plot = None
    try:
        with open(LIGAND_EFFECT_SAVE_LOCATION, 'rb') as f:
            default_plot = pickle.load(f)
    except:
        pass

    controls = control_panel(
        [
            control_panel_element("Interaction Set", "Biological condition to compare.",
                                  dbc.Select(
                                      id='effect_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'STING-Dependent Effect', 'value': 'sting'},
                                               {'label': 'CIN & STING Max Effect', 'value': 'max'}],
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['effect_set']
                                  )),
            control_panel_element("Network Layout", 'Select how you would like to structure nodes.',
                                  dbc.Select(
                                      id='network_layout',
                                      options=[
                                          {'label': 'Planar Layout', 'value': 'planar'},
                                          {'label': 'Spring Layout', 'value': 'spring'},
                                          {'label': 'Circular Layout', 'value': 'circular'},
                                          {'label': 'Timeline Layout', 'value': 'timeline'}
                                      ],
                                      persistence=True, persistence_type='session',
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['network_layout']
                                  )),
            control_panel_element('Emitting Cell Type', 'Select the initial cell type.',
                                  dbc.Select(
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
            control_panel_element("Emitted Ligands", 'Comma-separated list of ligands to emit from the selected cell type.',
                                  dbc.Input(
                                      id='ligands',
                                      autofocus=True,
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['ligands'],
                                      placeholder='Example: Ccl2,Apoe'
                                  ))
        ], [
            control_panel_element("Interaction FDR Cutoff", "The maximum interaction test FDR to consider.",
                                  dcc.Slider(
                                      id='interaction_fdr',
                                      max=1,
                                      min=0,
                                      step=0.01,
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['interaction_fdr'],
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
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['min_logfc'],
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
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['logfc_fdr'],
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
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['iterations'],
                                      persistence=True, persistence_type='session'
                                  )),
            control_panel_element("Plot", "",
                                  dbc.Button(
                                      "Submit",
                                      id='submit_button',
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
        element=dcc.Graph(id='network_graph',
                          animate=True,
                          figure=default_plot,
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
    Output('network_graph', 'figure'),
    Input('submit_button', 'n_clicks'),
    State('effect_set', 'value'),
    State('network_layout', 'value'),
    State('cell_type', 'value'),
    State('ligands', 'value'),
    State('interaction_fdr', 'value'),
    State('min_logfc', 'value'),
    State('logfc_fdr', 'value'),
    State('iterations', 'value'),
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
               effect_set, network_layout,
               cell_type, ligands,
               interaction_fdr, min_logfc,
               logfc_fdr, iterations):
    if n_clicks == 0:
        from dash.exceptions import PreventUpdate
        raise PreventUpdate

    from viz.figures import pseudotime_interaction_propagation_graph, DEFAULT_LIGAND_EFFECT_ARGS, LIGAND_EFFECT_SAVE_LOCATION

    set_progress((0, 100))

    # Check if arguments match default, if so return the pre-computed default
    if (effect_set == DEFAULT_LIGAND_EFFECT_ARGS['effect_set'] and
        network_layout == DEFAULT_LIGAND_EFFECT_ARGS['network_layout'] and
        cell_type == DEFAULT_LIGAND_EFFECT_ARGS['cell_type'] and
        ligands == DEFAULT_LIGAND_EFFECT_ARGS['ligands'] and
        interaction_fdr == DEFAULT_LIGAND_EFFECT_ARGS['interaction_fdr'] and
        min_logfc == DEFAULT_LIGAND_EFFECT_ARGS['min_logfc'] and
        logfc_fdr == DEFAULT_LIGAND_EFFECT_ARGS['logfc_fdr'] and
        iterations == DEFAULT_LIGAND_EFFECT_ARGS['iterations']):
        import pickle
        try:
            with open(LIGAND_EFFECT_SAVE_LOCATION, 'rb') as f:
                return pickle.load(f)
        except:
            pass

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
        (Output('submit_button', 'disabled'), True, False),
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


# If run as a script, compile the default plot
if __name__ == '__main__':
    from viz.figures import pseudotime_interaction_propagation_graph, DEFAULT_LIGAND_EFFECT_ARGS, LIGAND_EFFECT_SAVE_LOCATION

    if not os.path.exists(LIGAND_EFFECT_SAVE_LOCATION):
        import pickle

        fig = pseudotime_interaction_propagation_graph(
            effect_set=DEFAULT_LIGAND_EFFECT_ARGS['effect_set'],
            seed_cell=DEFAULT_LIGAND_EFFECT_ARGS['cell_type'],
            seed_ligands=DEFAULT_LIGAND_EFFECT_ARGS['ligands'],
            iterations=int(DEFAULT_LIGAND_EFFECT_ARGS['iterations']),
            interaction_fdr_cutoff=float(DEFAULT_LIGAND_EFFECT_ARGS['interaction_fdr']),
            min_logfc=float(DEFAULT_LIGAND_EFFECT_ARGS['min_logfc']),
            logfc_fdr_cutoff=float(DEFAULT_LIGAND_EFFECT_ARGS['logfc_fdr']),
            layout=DEFAULT_LIGAND_EFFECT_ARGS['network_layout'],
            set_progress_callback=None
        )

        with open(LIGAND_EFFECT_SAVE_LOCATION, 'wb') as f:
            pickle.dump(fig, f)
