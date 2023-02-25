import dash
from dash import html, dcc, callback, Output, Input, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate

from viz.figures import bipartite_graph
from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output

dash.register_page(__name__,
                   path='/interactions',
                   name='Interactions',
                   order=2)


def build_interface() -> list:
    controls = control_panel(
        [
            control_panel_element("Interaction Set", "Biological condition to compare.",
                                  dbc.Select(
                                      id='inter_set',
                                      options=[{'label': 'CIN-Dependent Effect', 'value': 'cin'},
                                               {'label': 'STING-Dependent Effect', 'value': 'sting'},
                                               {'label': 'CIN & STING Max Effect', 'value': 'max'}],
                                      value='cin'
                                  )),
            control_panel_element("Minimum numSigI1", "The minimum number of significant target gene interactions.",
                                  dcc.Slider(
                                      id="min_numsigi1_bipartite",
                                      min=0,
                                      max=10,  # Fill in
                                      step=1,
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
                                      step=0.01,
                                      value=0,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className="form-range"
                                  )),
            control_panel_element("Interaction Effect FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  dbc.Select(
                                      id='bipartite_inter_fdr',
                                      options=[{'label': '0.05', 'value': 'fdr05'},
                                               {'label': '0.25', 'value': 'fdr25'}],
                                      value='fdr05'
                                  )),
            control_panel_element("log2FC FDR Cutoff", "FDR-adjusted requirements for interaction effects.",
                                  dbc.Select(
                                      id='bipartite_logfc_fdr',
                                      options=[{'label': '0.05', 'value': 'fdr05'},
                                               {'label': '0.25', 'value': 'fdr25'}],
                                      value='fdr05'
                                  ))
        ], [
            control_panel_element("First Cell Type", "The first cell type to examine interactions between.",
                                  dbc.Select(
                                      id='first_celltype',
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
            control_panel_element("Second Cell Type", "The second cell type to examine interactions between (need not be unique).",
                                  dbc.Select(
                                      id='second_celltype',
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
                                    value='Macrophages/mMDSC'
                                  )),
            control_panel_element("Third Cell Type", "If specified, include interactions across a third cell type.",
                                  dbc.Select(
                                      id='third_celltype',
                                      options=[{'label': ct, 'value': ct} for ct in
                                               ['(None)',
                                                'Tumor cells',
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
                                      value='(None)',
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
    Input('submit-button-bipartite', 'n_clicks'),
    State('inter_set', 'value'),
    State('first_celltype', 'value'),
    State('second_celltype', 'value'),
    State('third_celltype', 'value'),
    State('min_logfc_bipartite', 'value'),
    State('min_numsigi1_bipartite', 'value'),
    State('bipartite_inter_fdr', 'value'),
    State('bipartite_logfc_fdr', 'value'),
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
def make_graph(set_progress, n_clicks,
               inter_set,
               first_ct, second_ct, third_ct,
               min_logfc, min_numSigI1,
               inter_fdr, logfc_fdr):
    from viz.data import read_interactions_file

    set_progress((0, 100))

    # Do some work
    interactions = read_interactions_file(inter_set, inter_fdr)

    fig = bipartite_graph(
        df=interactions,
        cell1=first_ct,
        cell2=second_ct,
        cell3=third_ct,
        numInteractions=min_numSigI1,
        min_logfc_bipartite=min_logfc,
        logfc_fdr_bipartite_cutoff=float(logfc_fdr.replace('fdr', '.'))
    )

    return fig


@callback(
    Output('min_logfc_bipartite', 'max'),
    Output('min_numsigi1_bipartite', 'max'),
    Input('bipartite_inter_fdr', 'value'),
    Input('bipartite_logfc_fdr', 'value'),
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
def initialize_options(inter_fdr, logfc_fdr):
    from viz.data import read_interactions_file

    fdr = inter_fdr

    # Use the max set
    interactions = read_interactions_file('max', fdr)

    max_logfc = max(max(interactions.MAST_log2FC_ligand.abs()), max(interactions.MAST_log2FC_receptor.abs()))
    max_inter = max(interactions.numSigI1)

    return max_logfc, max_inter


layout = [
    interactive_panel(wrap_icon('fa-arrows-left-right-to-line', 'Cell Type to Cell Type Interactions'),
                      *build_interface()
                      )
]
