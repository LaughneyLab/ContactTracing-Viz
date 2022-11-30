import math

import dash
import networkx as nx
import numpy as np
from dash import html, dcc, callback, Output, Input, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
from matplotlib import cm, colors

from viz.data import read_interactions
from viz.util import multipartite_layout, get_quiver_arrows
from viz.web import interactive_panel, wrap_icon, control_panel, control_panel_element, figure_output

dash.register_page(__name__,
                   path='/interactions',
                   name='Interactions',
                   order=1)


def build_interface() -> list:
    controls = control_panel(
        [
            control_panel_element("Minimum numSigI1", "abc",
                                  dcc.Slider(
                                      id="min_numsigi1_bipartite",
                                      min=0,
                                      max=10,  # Fill in
                                      step=1,
                                      value=1,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className="form-range"
                                  )),
            control_panel_element("Minimum Expression", "abc",
                                  dcc.Slider(
                                      id="min_expression_bipartite",
                                      min=0,
                                      max=1,
                                      step=0.1,
                                      value=0,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className="form-range"
                                  ))
        ], [
            control_panel_element("Minimum abs(log2FC)", "abc",
                                  dcc.Slider(
                                      id="min_logfc_bipartite",
                                      min=0,
                                      max=1,  # Fill in
                                      step=0.1,
                                      value=0,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=False,
                                      className="form-range"
                                  )),
            control_panel_element("log2FC FDR Cutoff", "abc",
                                  dcc.Slider(
                                      id="logfc_fdr_bipartite",
                                      min=0,
                                      max=1,
                                      step=0.01,
                                      value=0.05,
                                      marks=None,
                                      tooltip={'placement': 'bottom'},
                                      persistence=True, persistence_type='session',
                                      className="form-range"
                                  ))
        ], [
            control_panel_element("First Cell Type", "abc",
                                  dbc.Select(
                                      id='first_celltype',
                                      options=[],  # Filled in by callback
                                      persistence=False
                                  )),
            control_panel_element("Second Cell Type", "abc",
                                  dbc.Select(
                                      id='second_celltype',
                                      options=[],  # Filled in by callback
                                      persistence=False
                                  )),
            control_panel_element("Third Cell Type", "abc",
                                  dbc.Select(
                                      id='third_celltype',
                                      options=['(None)'],  # Filled in by callback
                                      value='(None)',
                                      persistence=False
                                  )),
        ], [
            control_panel_element("Plot", "abc",
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
        footer="This is a figure",
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
    Input('data-session', 'data'),
    Input('submit-button-bipartite', 'n_clicks'),
    State('first_celltype', 'value'),
    State('second_celltype', 'value'),
    State('third_celltype', 'value'),
    State('min_expression_bipartite', 'value'),
    State('min_logfc_bipartite', 'value'),
    State('logfc_fdr_bipartite', 'value'),
    State('min_numsigi1_bipartite', 'value'),
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
def make_graph(set_progress, data, n_clicks, first_ct, second_ct, third_ct, min_expression_bipartite, min_logfc_bipartite, logfc_fdr_bipartite, min_numSigI1):
    if data is None or n_clicks == 0:
        raise PreventUpdate
    set_progress((0, 100))

    # Do some work

    file = data['path']
    interaction_file = data['tsv']
    #adata = read_ct_data(file)
    interactions = read_interactions(interaction_file)
    fig = bipartite_graph(
        df=interactions,
        cell1=first_ct,
        cell2=second_ct,
        cell3=third_ct,
        numInteractions=min_numSigI1,
        min_logfc_bipartite=min_logfc_bipartite,
        logfc_fdr_bipartite_cutoff=logfc_fdr_bipartite
    )

    return fig


@callback(
    Output('first_celltype', 'options'),
    Output('first_celltype', 'value'),
    Output('second_celltype', 'options'),
    Output('second_celltype', 'value'),
    Output('third_celltype', 'options'),
    Output('third_celltype', 'value'),
    Output('min_expression_bipartite', 'max'),
    Output('min_expression_bipartite', 'value'),
    Output('min_logfc_bipartite', 'max'),
    Output('min_logfc_bipartite', 'value'),
    Output('min_numsigi1_bipartite', 'max'),
    Output('min_numsigi1_bipartite', 'value'),
    Input('data-session', 'data'),
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
def initialize_options(data):
    from viz.data import read_ct_data, read_interactions
    print(data)
    if data is None:
        return [], '', \
               [], '', \
               [{'label': '(None)', 'value': '(None)'}], '(None)', \
               1, 0, \
               1, 0, \
               10, 0
    file = data['path']
    interaction_file = data['tsv']
    #adata = read_ct_data(file)
    interactions = read_interactions(interaction_file)
    max_exp = max(max(interactions.expression_ligand), max(interactions.expression_receptor))
    max_logfc = max(max(interactions.MAST_log2FC_ligand.abs()), max(interactions.MAST_log2FC_receptor.abs()))
    max_inter = max(interactions.numSigI1_fdr05_receptor)
    celltypes = list(sorted(set(interactions.cell_type_ligand) | set(interactions.cell_type_receptor)))
    return [{'label': ct, 'value': ct} for ct in celltypes], celltypes[0], \
           [{'label': ct, 'value': ct} for ct in celltypes], celltypes[1], \
           [{'label': '(None)', 'value': '(None)'}] + [{'label': ct, 'value': ct} for ct in celltypes], '(None)', \
           max_exp, 0, \
           max_logfc, 0, \
           max_inter, 0


def bipartite_graph(df,
                    cell1,
                    cell2,
                    cell3=None,
                    numInteractions=50,
                    min_logfc_bipartite=0,
                    logfc_fdr_bipartite_cutoff=0.05,
                    min_expression_bipartite=0,
                    allow_overlap=True):
    scaleratio = .9
    height = 1000

    if 'none' in str(cell3).lower():
        selected = df[((df.cell_type_ligand == cell1) & (df.cell_type_receptor == cell2)) |
                      ((df.cell_type_ligand == cell2) & (df.cell_type_receptor == cell1))]
    else:
        selected = df[((df.cell_type_ligand == cell1) & (df.cell_type_receptor == cell2)) |
                      ((df.cell_type_ligand == cell2) & (df.cell_type_receptor == cell1)) |
                      ((df.cell_type_ligand == cell2) & (df.cell_type_receptor == cell3)) |
                      ((df.cell_type_ligand == cell3) & (df.cell_type_receptor == cell2)) |
                      (allow_overlap & (df.cell_type_ligand == cell1) & (df.cell_type_receptor == cell3)) |
                      (allow_overlap & (df.cell_type_ligand == cell3) & (df.cell_type_receptor == cell1))]

    selected = selected[(selected.numSigI1_fdr05_receptor >= numInteractions) &
                        (selected["MAST_fdr_ligand"] < logfc_fdr_bipartite_cutoff) &
                        (selected["MAST_fdr_receptor"] < logfc_fdr_bipartite_cutoff) &
                        (np.abs(selected["MAST_log2FC_ligand"]) > min_logfc_bipartite) &
                        (np.abs(selected["MAST_log2FC_receptor"]) > min_logfc_bipartite) &
                        (np.abs(selected["expression_receptor"]) > min_expression_bipartite) &
                        (np.abs(selected["expression_ligand"]) > min_expression_bipartite)].sort_values(
        by="numDEG_fdr05_receptor", ascending=False)

    G = nx.DiGraph()
    for i, row in selected.iterrows():
        lig_ct = row['cell_type_ligand']
        rec_ct = row['cell_type_receptor']
        lig = row['ligand']
        rec = row['receptor']

        lig_node = lig_ct + "_" + lig
        rec_node = rec_ct + "_" + rec

        if not G.has_node(lig_node):
            G.add_node(lig_node, ligand=lig in df.ligand.values, receptor=lig in df.receptor.values, name=lig,
                       celltype=lig_ct,
                       expression=row['expression_ligand'],
                       logfc_fdr_bipartite=row["MAST_fdr_ligand"],
                       logfc=row["MAST_log2FC_receptor"])
        if not G.has_node(rec_node):
            G.add_node(rec_node, ligand=rec in df.ligand.values, receptor=rec in df.receptor.values, name=rec,
                       celltype=rec_ct,
                       expression=row['expression_receptor'],
                       logfc_fdr_bipartite=row["MAST_fdr_receptor"],
                       logfc=row["MAST_log2FC_receptor"])

        G.add_edge(lig_node, rec_node,
                   numInteractions=row['numSigI1_fdr05_receptor'],
                   ligand_fdr=row["MAST_fdr_ligand"],
                   ligand_logfc=row["MAST_log2FC_ligand"],
                   receptor_fdr=row["MAST_fdr_receptor"],
                   receptor_logfc=row["MAST_log2FC_receptor"],
                   numDEG=row['numDEG_fdr05_receptor'])

    ct_ordering = [cell1, cell2]
    if cell3:
        ct_ordering.append(cell3)
    pos = multipartite_layout(G, subset_key='celltype', scale=1, space_mult_x=30 * scaleratio, space_mult_y=60,
                              ordering=ct_ordering)

    # start_x = []
    # start_y = []
    # end_x = []
    # end_y = []
    # edge_x = []
    # edge_y = []
    max_inter = None
    min_inter = None
    min_logfc_bipartite = None
    max_logfc = None
    for edge in G.edges():
        inters = G[edge[0]][edge[1]]['numInteractions']
        logfc = G[edge[0]][edge[1]]['ligand_logfc']
        if max_inter is None:
            min_inter = max_inter = inters
            min_logfc_bipartite = max_logfc = logfc
        else:
            if max_inter < inters:
                max_inter = inters
            if min_inter > inters:
                min_inter = inters
            if logfc < min_logfc_bipartite:
                min_logfc_bipartite = logfc
            if logfc > max_logfc:
                max_logfc = logfc

    cmap = cm.get_cmap('bwr')
    color_step = 255
    norm = colors.Normalize(vmin=0, vmax=color_step)
    colorscale = []
    for i in range(0, color_step):
        colorscale.append(colors.rgb2hex(colors.colorConverter.to_rgb(cmap(norm(i)))))

    # colorscale = [x.hex for x in list(Color(start_color).range_to(Color(end_color), color_step))]

    min_color_threshold = -0.2  # min_logfc_bipartite
    max_color_threshold = 0.2  # max_logfc

    def get_color(value, minimum, maximum):
        step_size = (max_color_threshold - min_color_threshold) / (color_step - 1)
        index = (value - min_color_threshold) // step_size
        return colorscale[min(max(int(index), 0), color_step - 1)]

    # Can't hover over lines directly, so add ivisible points
    # See https://stackoverflow.com/a/46039229/5179044
    interpolation_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers',
        hoverinfo='text',
        marker=dict(
            opacity=0,
            cmin=min_color_threshold,
            cmax=max_color_threshold,
            color=[],
            colorbar=dict(
                title="Ligand LogFC",
                xanchor='left',
                yanchor='top',
                thickness=25,
                titleside='right',
                x=1.15,
                y=1
            ),
            colorscale=[[i / (len(colorscale) - 1), col] for (i, col) in enumerate(colorscale)]
        )
    )

    edge_traces = []

    EDGE_PADDING = .08
    MIN_WIDTH = 2
    MAX_WIDTH = 6
    INTERPOLATION_POINTS = 20
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        num_inter = G[edge[0]][edge[1]]['numInteractions']
        logfc = G[edge[0]][edge[1]]['ligand_logfc']

        dif_x = x1 - x0
        dif_y = y1 - y0
        # total_len = math.hypot(dif_x, dif_y)
        angle = math.atan2(dif_y, dif_x / scaleratio)
        deltaX = EDGE_PADDING * math.cos(angle)
        deltaY = EDGE_PADDING * math.sin(angle)

        x0 += deltaX * .6
        y0 += deltaY * .6

        x1 -= deltaX * .4
        y1 -= deltaY * .4

        # edge_x += [x0, x1, None]
        # edge_y += [y0, y1, None]
        # start_x.append(x0)
        # start_y.append(y0)
        # end_x.append(x1)
        # end_y.append(y1)

        edge_x = [x0, x1, None]
        edge_y = [y0, y1, None]
        arrow_x, arrow_y = get_quiver_arrows([x0], [y0], [x1], [y1], scaleratio)

        color = get_color(logfc, min_logfc_bipartite, max_logfc)
        arrow_xs = edge_x + arrow_x
        arrow_ys = edge_y + arrow_y
        arrow_width = MIN_WIDTH + (MAX_WIDTH - MIN_WIDTH) * ((num_inter - min_inter + 1) / (max_inter - min_inter + 1))
        edge_trace = go.Scatter(
            x=arrow_xs, y=arrow_ys,
            fill="toself",
            fillcolor=color,
            line=dict(
                width=arrow_width,
                color=color
            ),
            hoverinfo='none',
            mode='lines'
        )
        hover_text = f"Ligand: {G.nodes[edge[0]]['name']}<br>Receptor: {G.nodes[edge[1]]['name']}<br>Interactions: {num_inter}<br>Ligand logFC: {logfc}"
        edge_trace.text = [hover_text] * (len(edge_x) + len(arrow_x))
        edge_traces.append(edge_trace)

        step_x = dif_x / INTERPOLATION_POINTS
        step_y = dif_y / INTERPOLATION_POINTS

        for point in [(x0 + i * step_x, y0 + i * step_y) for i in range(0, INTERPOLATION_POINTS)]:
            interpolation_trace['x'] += (point[0],)
            interpolation_trace['y'] += (point[1],)
            interpolation_trace['text'] += (hover_text,)
            interpolation_trace['marker']['color'] += (color,)

    # arrow_x, arrow_y = get_quiver_arrows(start_x, start_y, end_x, end_y, scaleratio)

    # edge_trace = go.Scatter(
    #    x=edge_x + arrow_x, y=edge_y + arrow_y,
    #    line=dict(
    #        width=1.5,
    #        color='#888'
    #    ),
    #    hoverinfo='none',
    #    mode='lines'
    # )

    # edge_trace = create_quiver(start_x, start_y, end_x, end_y)

    node_x = []
    node_y = []
    symbols = []
    layer_anchors = dict()
    for node in G.nodes():
        layer_name = G.nodes[node]['celltype']
        ligand = G.nodes[node]['ligand']
        receptor = G.nodes[node]['receptor']
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

        assert ligand or receptor

        if ligand and receptor:
            symbols.append('diamond-dot')
        elif ligand:
            symbols.append("circle")
        elif receptor:
            symbols.append('square')

        if layer_name not in layer_anchors:
            layer_anchors[layer_name] = (x, y, 1)
        else:
            curr_anchor = layer_anchors[layer_name]
            if y > curr_anchor[1]:
                layer_anchors[layer_name] = (x, y, 1 + curr_anchor[2])

    annotations = []
    celltype_background_points = []
    yrange = max(node_y) - min(node_y) if len(node_y) > 0 else 0
    if len(layer_anchors) == 0:
        xrange = 0
    else:
        xrange = (max(node_x) - min(node_x) if len(node_x) > 0 else 0) / len(layer_anchors)
    for layer_name, (anchor_x, anchor_y, layer_size) in layer_anchors.items():
        annotations.append(dict(
            x=anchor_x,
            y=max(node_y) + (yrange / 12) if len(node_y) > 0 else 0,
            xref="x",
            yref="y",
            text=layer_name,
            showarrow=False,
            font=dict(
                size=16,
                color="white"
            ),
            align="center",
            bordercolor="black",
            borderwidth=2,
            borderpad=4,
            bgcolor="grey"
        ))

        celltype_background_points += [
            (anchor_x - .4 * xrange, min(node_y) - .05 * yrange),
            (anchor_x - .4 * xrange, max(node_y) + .05 * yrange),
            (anchor_x + .4 * xrange, max(node_y) + .05 * yrange),
            (anchor_x + .4 * xrange, min(node_y) - .05 * yrange),
            (anchor_x - .4 * xrange, min(node_y) - .05 * yrange),
            (None, None)
        ]

    celltype_bg_trace = go.Scatter(
        x=[c[0] for c in celltype_background_points],
        y=[c[1] for c in celltype_background_points],
        fill="toself",
        fillcolor='#f8f3dd',
        line=dict(
            width=2,
            color='grey'
        ),
        hoverinfo='none',
        mode='lines'
    )

    first_layer = cell1
    last_layer = cell2 if len(layer_anchors) == 2 else cell3
    middle_layer = cell2 if len(layer_anchors) > 2 else None

    for node in G.nodes():
        layer_name = G.nodes[node]['celltype']
        layer_size = layer_anchors[layer_name][2]
        x, y = pos[node]

        xshift = 0
        yshift = 0
        # if layer_name == last_layer:
        #    xshift += 80
        # elif layer_name == first_layer:
        #    xshift -= 80
        # elif layer_name == middle_layer:
        #    yshift += 20
        yshift += 20

        align = "center"
        if layer_name == last_layer:
            align = "left"
        elif layer_name == first_layer:
            align = "right"

        clicktoshow = "onoff"
        # if (layer_size > 40 and layer_name != middle_layer) or (layer_size > 20 and layer_name == middle_layer):

        annotations.append(dict(
            x=x,
            xshift=xshift,
            y=y,
            yshift=yshift,
            xref="x",
            yref="y",
            clicktoshow=clicktoshow,
            text=G.nodes[node]['name'],
            showarrow=False,
            font=dict(
                size=14,
                color="black"
            ),
            align=align
        ))

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        marker_symbol=symbols,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='gray',
            reversescale=True,
            color=[],
            size=15,
            line=dict(width=1, color='black'),
            colorbar=dict(
                thickness=25,
                title="Expression",
                xanchor='left',
                titleside='right'
            )
        )
    )

    node_expressions = []
    node_text = []
    for node in G.nodes:
        expression = G.nodes[node]['expression']
        name = G.nodes[node]['name']
        logfc_fdr_bipartite = G.nodes[node]['logfc_fdr_bipartite']
        logfc = G.nodes[node]['logfc']
        node_expressions.append(expression)
        ligand = G.nodes[node]['ligand']
        receptor = G.nodes[node]['receptor']
        if ligand and receptor:
            prot_type = "Ligand/Receptor"
        elif ligand:
            prot_type = "Ligand"
        else:
            prot_type = "Receptor"

        node_text.append(
            f'Name: {name}<br>Type: {prot_type}<br>Expression: {expression}<br>logFC: {logfc}<br>logFC FDR: {logfc_fdr_bipartite}')

    node_trace.marker.color = node_expressions
    node_trace.text = node_text

    if len(node_expressions) == 0:
        annotations.append(dict(
            x=0,
            y=0,
            xref="x",
            yref="y",
            clicktoshow=False,
            text="No Interactions (try loosening the filters?)",
            showarrow=False,
            font=dict(
                size=24,
                color="red"
            ),
            align="center"
        ))

    fig = go.Figure(data=[celltype_bg_trace] + edge_traces + [node_trace, interpolation_trace],
                    layout=go.Layout(
                        title='&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Detected Interactions',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=40, l=5, r=5, t=40),
                        autosize=False,
                        width=height * scaleratio,
                        height=height,
                        plot_bgcolor='white',
                        annotations=annotations,
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, scaleratio=scaleratio,
                                   scaleanchor="x"))
                    )
#    return enhance_plotly_export(fig, height, scaleratio)
    return fig


layout = [
    interactive_panel(wrap_icon('fa-arrows-left-right-to-line', 'Cell Type to Cell Type Interactions'),
                      *build_interface()
                      )
]
