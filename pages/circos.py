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
            control_panel_element("log2FC FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  make_fdr_slider('logfc_circos_fdr', DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'])),
        ], [
             control_panel_element('Minimum numSigI1', 'Minimum number of significant interactions for a receptor to be included.',
                                  make_custom_slider(
                                      id='circos_min_numsigi1',
                                      min=0,
                                      max=250,  # Fill in later
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
                                      step=1
                                  )),
            control_panel_element('Minimum numDEG', 'Minimum number of differentially expressed genes conditioned on a target gene.',
                                  make_custom_slider(
                                      id='circos_min_numdeg',
                                      min=0,
                                      max=250,  # Fill in later
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                                      step=1
                                  ))
        ], [
            control_panel_element('Chord Ligand abs(log2FC) Cutoff', 'The minimum ligand log2FC required for a chord to be drawn.',
                                  make_custom_slider(
                                      id='circos_min_ligand_logfc',
                                      min=0,
                                      max=2,  # Fill in
                                      value=DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'],
                                      step=0.01
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
            control_panel_element("Plot", "",
                                  html.Div(
                                      dbc.Button(
                                          "Submit",
                                          id="submit-button-circos",
                                          size="lg",
                                          color="primary",
                                          className='me-1',
                                          n_clicks=0
                                      ),
                                      className='text-center d-grid gap-2'))
        ]
    )

    results = figure_output(
        title="Interactive Circos Plot",
        footer="Layers from outside ring to center: Cell Type, Diffusion Component Value, Differential Abundance, Number of significant interactions, Strongest ligand/receptor interactions",
        element=html.Div(
            id="circos-graph-holder",
            children=default_plots[1] if default_plots is not None else None
        )
    )

    return [
        controls,
        results,
        dcc.Store(id='cin_circos_plot', storage_type='memory', data=[default_plots[0] if default_plots is not None else None]),
        dcc.Store(id='sting_circos_plot', storage_type='memory', data=[default_plots[1] if default_plots is not None else None])
    ]


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
                     min_numsigi1, min_numdeg, min_chord_ligand_logfc):
    if n_clicks == 0:
        from dash.exceptions import PreventUpdate
        raise PreventUpdate

    set_progress((0, 14))
    from viz.data import read_circos_file, read_interactions_file
    from viz.web import make_circos_figure
    from viz.figures import DEFAULT_CIRCOS_ARGS, CIRCOS_SAVE_LOCATION

    # Check if arguments match default, if so return the pre-computed default
    if (inter_circos_fdr == DEFAULT_CIRCOS_ARGS['inter_circos_fdr'] and
        logfc_circos_fdr == DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'] and
        min_numsigi1 == DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'] and
        min_numdeg == DEFAULT_CIRCOS_ARGS['circos_min_numdeg'] and
        min_chord_ligand_logfc == DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc']):
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
                               min_chord_ligand_logfc)

    sting_circos = make_circos_figure(set_progress, 1,
                               data,
                               sting_inter_data,
                               logfc_circos_fdr,
                               min_numsigi1,
                               min_numdeg,
                               min_chord_ligand_logfc)

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

    if not os.path.exists(CIRCOS_SAVE_LOCATION):
        from viz.data import read_circos_file, read_ligand_receptor_file, read_interactions_file
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
                                 DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'])
        sting_fig = make_circos_figure(None, 1,
                                     data,
                                     sting_inter_data,
                                     DEFAULT_CIRCOS_ARGS['logfc_circos_fdr'],
                                     DEFAULT_CIRCOS_ARGS['circos_min_numsigi1'],
                                     DEFAULT_CIRCOS_ARGS['circos_min_numdeg'],
                                     DEFAULT_CIRCOS_ARGS['circos_min_ligand_logfc'])

        # Dump the figure to a pickle file
        with open(CIRCOS_SAVE_LOCATION, 'wb') as f:
            pickle.dump([cin_fig, sting_fig], f)
