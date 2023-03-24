import os.path

import dash
import dash_bootstrap_components as dbc
from dash import callback, Output, Input, dcc, State, html

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

    import pickle
    default_plots = None
    try:
        with open(LIGAND_EFFECT_SAVE_LOCATION, 'rb') as f:
            default_plots = pickle.load(f)
    except:
        pass

    controls = control_panel(
        [
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
                                  )),
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
        ], [
            control_panel_element("Interaction Effect FDR Cutoff", "The maximum interaction test FDR to consider.",
                                  make_fdr_slider("interaction_fdr", DEFAULT_LIGAND_EFFECT_ARGS['interaction_fdr'])),
            control_panel_element("Log2FC FDR Cutoff",
                                  "The FDR-adjusted p-value cutoff for determining if a gene's log2FC value is non-zero.",
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
            control_panel_element("Minimum abs(Log2FC)", "The minimum induced log2FC for gene interactions.",
                                  make_custom_slider(
                                      id='min_logfc',
                                      max=2,
                                      min=0,
                                      step=0.01,
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['min_logfc']
                                  )),
        ], [
            control_panel_element("Interaction Set", "Biological condition to compare.",
                                  dbc.RadioItems(
                                      id='effect_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'CIN & STING Max Effect', 'value': 'max', 'disabled': True}
                                               ],
                                      value=DEFAULT_LIGAND_EFFECT_ARGS['effect_set'],
                                      persistence=False
                                  )),
            control_panel_element("Plot", "You must click submit to update the plot.",
                                  html.Div(
                                      dbc.Button(
                                          "Submit",
                                          id='submit_button',
                                          size='lg',
                                          color="primary",
                                          className='me-1',
                                          n_clicks=0
                                      ),
                                      className='text-center d-grid gap-2')),
        ]
    )

    results = figure_output(
        title='Ligand Network Figure',
        footer="Circle = Ligand, Square = Receptor, Diamond = Ligand and Receptor",
        element=dcc.Graph(id='network_graph',
                          figure=default_plots[0] if default_plots else {},
                          config={
                             'displaylogo': False,
                             'showTips': True,
                             'toImageButtonOptions': {
                                'format': 'svg',  # one of png, svg, jpeg, webp
                                'filename': 'ligand_effects_image',
                                'height': 1800,
                                'width': 1200,
                                'scale': 1  # Multiply title/legend/axis/canvas sizes by this factor
                             },
                             'watermark': False
                          }),
        help_info=make_help_info()
    )

    return [
        controls,
        results,
        dcc.Store(id='cin_network_plot', storage_type='memory', data=default_plots[0] if default_plots is not None else {}),
        dcc.Store(id='sting_network_plot', storage_type='memory', data=default_plots[1] if default_plots is not None else {})
    ]


def make_help_info():
    return [
        html.H5(["The ligand cascade effects plot visualizes all potential downstream effects from TME-specific"
                 " identified by ",
                 html.I("ContactTracing"),
                 "."]),
        html.Hr(),
        html.P(html.H6(html.Strong("Filter Options"))),
        html.Ul([
            html.Li(html.P([
                "Emitting Cell Type: The cell type for which specified ligands are emitted."
            ])),
            html.Li(html.P([
                "Emitted Ligands: The ligand(s) to follow in the figure; must be emitted from the selected cell type."
                " Multiple genes can be included by separating names with commas. IMPORTANT: This field is CASE"
                " SENSITIVE."
            ])),
            html.Li(html.P([
                "Network Building Iterations: The number of levels of downstream interactions to include. Note that"
                " typically, the number of downstream genes per level can increase exponentially. So starting with a"
                " low value and increasing as needed is recommended since larger values require additional time to"
                " render."
            ])),
            html.Li(html.P([
                "Interaction Effect FDR Cutoff: ",
                html.I("ContactTracing"),
                " identifies the downstream genes that have induced expression shifts from the activation of a "
                " receptor within a given cell of interest; this modulates the threshold for classifying induced"
                " expression change as significant after a Benjamini-Hochberg FDR correction. The default value of 0.25"
                " reflects what was chosen for evaluating the intersection between CIN and STING-dependent effects. A"
                " cutoff of 0.05 might be more sensible when evaluating CIN-dependent effects on its own."
            ])),
            html.Li(html.P([
                "Log2FC FDR Cutoff: ",
                html.I("ContactTracing"),
                " requires differential availability of a ligand across conditions to identify condition-specific"
                " effects from the Tumor Microenvironment. Additionally, as ",
                html.I("ContactTracing"),
                " is able to determine the induced Log2FC of a receptor's activation on a downstream gene, this option"
                " will also be applied to induced Log2FC tests. The cutoff represents the threshold for identifying whether"
                " a gene's differential expression between conditions significantly differs from 0 after a"
                " Benjamini-Hochberg FDR correction."
            ])),
            html.Li(html.P([
                "Minimum Expression: Allows for Ligands/Receptors to be filtered according to a minimum expression"
                " value. This value, however, does not represent counts. Instead, \"expression\" refers to the fraction"
                " of the cells from a given cell type with non-zero counts of a particular gene. Therefore, expression"
                " values range from 0-1, where zero means that no cells of a cell type express the gene, and one means"
                " that all cells of a cell type express the gene."
            ])),
            html.Li(html.P([
                "Minimum abs(Log2FC): While the Log2FC FDR Cutoff option filters interactions according"
                " to whether a gene's differential expression between conditions significantly differs from zero,"
                " this filter additionally allows for further refinement by requiring a minimum absolute Log2FC. The"
                " default value reflects the cutoff used for the Circos plot."
            ])),
            html.Li(html.P([
                "Interaction Set: While we have evaluated both CIN-dependent and CIN & STING-dependent interaction"
                " effects, the CIN & STING intersection data set does not have statistical power for examining"
                " individually induced downstream effects of interactions. Therefore, only CIN-dependent interaction"
                " effects are enabled."
            ]))
        ]),
        html.P(html.H6(html.Strong("How to Interpret the Plot"))),
        html.P("This figure illustrates the potential cascading effects of receptor activation in the tumor"
               " microenvironment. It does so by generating a directed network of gene interactions by following the"
               " steps listed below:"),
        html.Ol([
            html.Li([html.P([
                    "Starting from a given cell type and a set of ligands, identify all the corresponding receptors that"
                    " pass the provided filters and have a significant condition-specific interaction effect on downstream"
                    " signaling, then draw an arrow connecting the ligands to the receptors. Finally, add 1 to the number"
                    " of network building iterations."
                ]),
                html.Div(html.Img(src=dash.get_asset_url("ligand_effects_help1.png"), style={'width': '30%', 'align': 'center'}, alt='Step 1', className='mx-auto'), className='text-center')
            ]),
            html.Li([html.P([
                    "Identify significantly differentially regulated ligands for each added receptor in response to the"
                    " receptor's activation. Additionally, identify further receptors with a significant interaction effect"
                    " and a ligand available in the microenvironment. Finally, add 1 to the number of network building"
                    " iterations."
                ]),
                html.Div(html.Img(src=dash.get_asset_url("ligand_effects_help2.png"), style={'width': '25%', 'align': 'center'}, alt='Step 2', className='mx-auto'), className='text-center')
            ]),
            html.Li([html.P([
                    "Repeat steps 1 and 2 until the number of network building iterations equals what was specified by"
                    " the user."
                ]),
                html.Div(html.Img(src=dash.get_asset_url("ligand_effects_help3.png"), style={'width': '25%', 'align': 'center'}, alt='Step 3', className='mx-auto'), className='text-center')
            ])
        ]),
        html.P(["As depicted above, each of these iterations is represented by a column of ligands and receptors. In"
                " addition, these nodes have arrows connecting to their corresponding receptors and ligands. If the"
                " arrow is connecting from a receptor to a ligand, the color and arrow thickness represents the Log2FC"
                " induced by receptor activation. However, if the arrow is connecting from a ligand to a receptor, the"
                " color is equivalent to a Log2FC of 0 as ligands are not expected to regulate receptor expression"
                " differentially."]),
        html.P(["Additionally, the color of nodes represents the cell type in which each gene is expressed; the shape"
                " of the node represents whether the gene is a ligand, receptor, or both; and the size of each node is"
                " relative to the total induced Log2FC from receptor activation to ligands."])
    ]


@callback(
    Output('network_graph', 'figure'),
    Input('effect_set', 'value'),
    Input('cin_network_plot', 'data'),
    Input('sting_network_plot', 'data'),
    prevent_initial_call=True
)
def update_network_figure(effect_set, cin_network_plot, sting_network_plot):
    if effect_set == 'cin':
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
    cache_args_to_ignore=['submit_button', 'n_clicks'],
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
    if n_clicks == 0:
        from dash.exceptions import PreventUpdate
        raise PreventUpdate

    from viz.figures import pseudotime_interaction_propagation_graph, DEFAULT_LIGAND_EFFECT_ARGS, LIGAND_EFFECT_SAVE_LOCATION

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
            with open(LIGAND_EFFECT_SAVE_LOCATION, 'rb') as f:
                return pickle.load(f)
        except:
            pass

    fig = pseudotime_interaction_propagation_graph(
        effect_set=DEFAULT_LIGAND_EFFECT_ARGS['effect_set'],  # Keep this fixed
        seed_cell=cell_type,
        seed_ligands=ligands,
        iterations=int(iterations),
        interaction_fdr_cutoff=float("." + interaction_fdr[3:]),
        min_logfc=float(min_logfc),
        min_expression=float(min_expression),
        logfc_fdr_cutoff=float("." + logfc_fdr[3:]),
        set_progress_callback=set_progress
    )

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
    interactive_panel(wrap_icon('fa-maximize', 'Downstream Ligand Effects'), "Predict the downstream effects of ligands across a microenvironment.",
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
