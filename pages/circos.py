import os

import dash
from dash import html, callback, Output, Input, State, dcc
import dash_bootstrap_components as dbc

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

    controls = control_panel(
        [
            control_panel_element("Interaction Effect FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  make_fdr_slider('inter_circos_fdr', DEFAULT_CIRCOS_ARGS['inter_circos_fdr'])),
            control_panel_element("Ligand Log2FC FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  make_fdr_slider('logfc_circos_fdr', DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'])),
        ], [
             control_panel_element('Minimum Interaction Effect', 'Minimum number of significant interactions for a receptor to be included.',
                                  make_custom_slider(
                                      id='circos_min_numsigi1',
                                      min=0,
                                      max=100,  # Fill in later
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
                                      step=1
                                  )),
            control_panel_element('Minimum Affected Downstream Genes', 'Minimum number of differentially expressed genes conditioned on a target gene.',
                                  make_custom_slider(
                                      id='circos_min_numdeg',
                                      min=0,
                                      max=100,  # Fill in later
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                                      step=1
                                  ))
        ], [
            control_panel_element('Chord Minimum Ligand abs(Log2FC)', 'The minimum ligand log2FC required for a chord to be drawn.',
                                  make_custom_slider(
                                      id='circos_min_ligand_logfc',
                                      min=0,
                                      max=2,  # Fill in
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'],
                                      step=0.01
                                  )),
            control_panel_element("Genes of Interest",
                                  'Comma-separated list of genes to highlight interactions for in the plot.',
                                  dbc.Input(
                                      id='genes',
                                      autofocus=True,
                                      value=DEFAULT_CIRCOS_ARGS['genes'],
                                      placeholder='Example: Ccl2,Apoe'
                                  ))
        ], [
            control_panel_element("Interaction Set", "Biological condition to compare.",
                                  dbc.RadioItems(
                                      id='circos_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'CIN & STING Max Effect', 'value': 'sting'}],
                                      value=DEFAULT_CIRCOS_ARGS['circos_set'],
                                      persistence=False
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
                                                         id='reset_default_circos_button', color='secondary', size='lg', className='float-end dark', n_clicks=0),
                                          align='right', width='auto')
                                      ]),
                                      className='text-center d-grid gap-2'))
        ]
    )

    results = figure_output(
        title="Interactive Circos Plot",
        footer="Layers from outside ring to center: Cell Type, Diffusion Component Value, Differential Abundance, Number of significant interactions, Strongest ligand/receptor interactions",
        element=html.Div(
            id="circos-graph-holder",
            children=default_plots[1] if default_plots is not None else None
        ),
        help_info=make_help_info()
    )

    return [
        controls,
        results,
        dcc.Store(id='cin_circos_plot', storage_type='memory', data=[default_plots[0] if default_plots is not None else None]),
        dcc.Store(id='sting_circos_plot', storage_type='memory', data=[default_plots[1] if default_plots is not None else None])
    ]


def make_help_info():
    return [
        html.H5(["The ",
                 html.I('ContactTracing'),
                 " Circos plot aggregates all relevant interaction information into a single figure."]),
        html.Hr(),
        html.P(html.H6(html.Strong("Filter Options"))),
        html.Ul([
            html.Li(html.P([
                "Interaction Effect FDR Cutoff: ",
                html.I("ContactTracing"),
                " identifies the downstream genes that have induced expression shifts from the activation of a receptor"
                " within a given cell of interest; this modulates the threshold for classifying induced expression"
                " change as significant after a Benjamini-Hochberg FDR correction. The default value of 0.25 reflects"
                " what was chosen for evaluating the intersection between CIN and STING-dependent effects. A cutoff"
                " of 0.05 might be more sensible when evaluating CIN-dependent effects on its own."
            ])),
            html.Li(html.P([
                "Ligand Log2FC FDR Cutoff: ",
                html.I("ContactTracing"),
                " requires differential availability of a ligand across conditions to identify condition-specific"
                " effects from the Tumor Microenvironment. This cutoff represents the threshold for identifying whether"
                " a ligand's differential expression between conditions significantly differs from 0 after a"
                " Benjamini-Hochberg FDR correction."
            ])),
            html.Li(html.P([
                "Minimum Interaction Effect: As ",
                html.I("ContactTracing"),
                " ContactTracing identifies the condition-specific effects of each receptor's activation by evaluating "
                "all possible genes, this value indicates a filter for the minimum number of condition-specific "
                "activations of genes by each given receptor using the selected Interaction Effect FDR Cutoff. The "
                "interaction effect, in part, determines whether a receptor has interactions depicted in the "
                "\"ribbons\" of the Circos plot."
            ])),
            html.Li(html.P([
                "Minimum Affected Downstream Genes: The number of differentially expressed genes in response to"
                " receptor activation. This value is similar to the Minimum Interaction Effect parameter, but rather"
                " than requiring condition-specific responses, this measures the strength of the Log2FC induced by a"
                " receptor's activation. Note that the FDR-adjusted p-value threshold is the same as the Minimum"
                " Interaction Effect parameter. "
            ])),
            html.Li(html.P([
                "Chord Minimum Ligand abs(Log2FC): Whereas the Ligand Log2FC FDR Cutoff determines whether genes get "
                "included in the plot, this parameter only affects whether a ribbon is drawn to indicate strong "
                "interaction effects within the Circos figure. Increasing this value leads to fewer ribbons drawn "
                "between ligands and receptors, while the opposite is true if decreased."
            ])),
            html.Li(html.P([
                "Genes of Interest: If specified, the chords drawn corresponding to a ligand or receptor listed will "
                "be emphasized (below is an example). Multiple genes can be included by separating names with commas. ",
                html.Em("IMPORTANT: This field is CASE SENSITIVE.")
            ])),
            html.Li(html.P([
                "Interaction Set: We have evaluated both CIN-dependent and CIN & STING-dependent interaction effects. "
                "This toggle lets users instantly see the plot under both conditions. Note that the CIN & STING "
                "effects represent the maximum value between shared interactions across both conditions and require "
                "interaction effects to have the same directionality."
            ]))
        ]),
        html.P(html.H6(html.Strong("How to Interpret the Circos Plot"))),
        html.P("The final figure is very information-dense, so it can be challenging to interpret initially. "
               "Thus, we will explore the diagram layer by layer. "),
        html.P("First, we have sections of the rings grouped by cell type annotation. NOTE: If a cell type name is too "
               "large, it is excluded. However, it is possible to identify the cell type by hovering over each section "
               "of the plot in the browser."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help1.png'), style={'width': '40%', 'align': 'center'}, alt="Cell Types", className='mx-auto'), className='text-center'),
        html.P("The next layer represents the value of the first Diffusion Component, as calculated by Palantir. This "
               "value represents the euclidean space embedding of the differential expression score (the Log2FC "
               "between conditions multiplied by the negative Log10-transformed p-value) of a given gene for a cell "
               "type. The diffusion component allows for a consistent ordering of genes along the Circos rings."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help2.png'), style={'width': '40%', 'align': 'center'}, alt="Diffusion Components", className='mx-auto'), className='text-center'),
        html.P("Following Diffusion Component 1, the next ring represents the gene-wise differential abundance of a "
               "gene between conditions. The differential abundance is calculated as the Pearson correlation "
               "coefficient of each gene's imputed expression across cells against the corresponding cell type's "
               "differential abundance across conditions as calculated by MILO. This value is then normalized to "
               "range from -1 to 1."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help3.png'), style={'width': '40%', 'align': 'center'}, alt="Differential Abundance", className='mx-auto'), className='text-center'),
        html.P("The final ring represents a histogram scaled to the number of significant interaction effects for each "
               "receptor."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help4.png'), style={'width': '40%', 'align': 'center'}, alt="Interaction Effects", className='mx-auto'), className='text-center'),
        html.P("Lastly, all interactions that meet the selected requirements and have the same directionality are "
               "depicted as ribbons connecting each receptor to its associated ligands across cell types. Depending on "
               "the number of ribbons, it may be difficult to read some gene labels. Therefore it is possible to hover "
               "over each ribbon to get precise details. The color of these ribbons reflects a ligand's Log2FC between "
               "conditions, and the thickness of the ribbon is relative to the number of significant interaction "
               "effects for the receptor."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help5.png'), style={'width': '40%', 'align': 'center'}, alt="Ribbons", className='mx-auto'), className='text-center'),
        html.P("Additionally, the ribbons can be filtered to emphasize specific ligands or receptors of interest using "
               "the aforementioned \"Genes of Interest\" field."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help6.png'), style={'width': '40%', 'align': 'center'}, alt="Ribbons Highlighted", className='mx-auto'), className='text-center')
    ]


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
    Input('sting_circos_plot', 'data'),
    prevent_initial_call=True
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
    State('genes', 'value'),
    background=True,
    prevent_initial_call=True,
    interval=500,
    cache_args_to_ignore=['submit-button-circos', 'n_clicks'],
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
                     min_numsigi1, min_numdeg, min_chord_ligand_logfc, gene_list):
    if n_clicks == 0:
        from dash.exceptions import PreventUpdate
        raise PreventUpdate

    if gene_list is not None and len(gene_list.strip()) == 0:
        gene_list = None

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
                               logfc_circos_fdr,
                               min_numsigi1,
                               min_numdeg,
                               min_chord_ligand_logfc,
                               gene_list)

    sting_circos = make_circos_figure(set_progress, 1,
                               data,
                               sting_inter_data,
                               logfc_circos_fdr,
                               min_numsigi1,
                               min_numdeg,
                               min_chord_ligand_logfc,
                               gene_list)

    set_progress((14, 14))

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
                             DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'],
                             DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
                             DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                             DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'],
                             DEFAULT_CIRCOS_ARGS['genes'])
    sting_fig = make_circos_figure(None, 1,
                                 data,
                                 sting_inter_data,
                                 DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'],
                                 DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
                                 DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                                 DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'],
                                 DEFAULT_CIRCOS_ARGS['genes'])

    # Dump the figure to a pickle file
    with open(CIRCOS_SAVE_LOCATION, 'wb') as f:
        pickle.dump([cin_fig, sting_fig], f)
