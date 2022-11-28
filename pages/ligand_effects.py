import math

import anndata as ad
import dash
import dash_bootstrap_components as dbc
from dash import html, callback, Output, Input
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly
import networkx as nx
from matplotlib import colors, cm

from viz.data import get_diff_abundance, get_interaction_fdr, get_effective_logfc, get_downstream_ligands, read_ct_data
from viz.util import enhance_plotly_export, celltype_to_colors, get_quiver_arrows
from viz.web import interactive_panel, wrap_icon

dash.register_page(__name__,
                   path='/ligand-effects',
                   name='Ligand Effects',
                   order=1)


def build_interface() -> list:
    controls = [  # Each CardGroup is a row
        dbc.CardGroup([
            dbc.Card(
                dbc.CardHeader("Network Layout"),
                dbc.CardBody([
                    dbc.Select(
                        id='network-layout',
                        options=[
                            {'label': 'Planar Layout', 'value': 'planar'},
                            {'label': 'Spring Layout', 'value': 'spring'},
                            {'label': 'Circular Layout', 'value': 'circular'},
                        ],
                        persistence=True, persistence_type='session', value='planar'
                    ),
                    html.P("abc", className='card-text'),
                ]),
                outline=True, color='light'
            ),
            dbc.Card(
                dbc.CardHeader("Emitting Cell Type"),
                dbc.CardBody([
                    dbc.Select(
                        id='cell-type',
                        options=[],  # Filled in by callback
                        persistence=False
                    ),
                    html.P("abc", className='card-text'),
                ]),
                outline=True, color='light'
            ),
            dbc.Card(
                dbc.CardHeader("Emitted Ligands"),
                dbc.CardBody([
                    dbc.Input(
                        id='ligands',
                        autofocus=True,
                        placeholder='Example: Ccl2,Apoe',
                        persistence=False
                    ),
                    html.P("abc", className='card-text'),
                ]),
                outline=True, color='light'
            ),
        ]),
        dbc.CardGroup([
            dbc.Card(
                dbc.CardHeader("Minimum Expression"),
                dbc.CardBody([
                    dbc.Slider()
                    html.P("abc", className='card-text'),
                ]),
                outline=True, color='light'
            ),
            dbc.Card(
                dbc.CardHeader("Interaction FDR cutoff"),
                dbc.CardBody([
                    html.P("abc", className='card-text'),
                ]),
                outline=True, color='light'
            ),
        ]),
        dbc.CardGroup([
            dbc.Card(
                dbc.CardHeader("Minimum abs(LogFC)"),
                dbc.CardBody([
                    html.P("abc", className='card-text'),
                ]),
                outline=True, color='light'
            ),
            dbc.Card(
                dbc.CardHeader("log2FC FDR cutoff"),
                dbc.CardBody([
                    html.P("abc", className='card-text'),
                ]),
                outline=True, color='light'
            ),
        ]),
        dbc.CardGroup([
            dbc.Card(
                dbc.CardHeader("Network Building Iteratiosn"),
                dbc.CardBody([
                    html.P("abc", className='card-text'),
                ]),
                outline=True, color='light'
            ),
            dbc.Card(
                dbc.CardHeader("Plot"),
                dbc.CardBody([
                    html.P("abc", className='card-text'),
                ]),
                outline=True, color='light'
            )
        ])
    ]

    results = []

    return [
        html.Div(controls),
        html.Div(results)
    ]


@callback(
    Output('cell-type', 'options'),
    Output('cell-type', 'value'),
    Input('data-session', 'data')
)
def initialize_options(data):
    if data is None:
        return [], ''
    file = data['filename']
    adata = read_ct_data(file)
    celltypes = list(adata.obs['celltype'].unique())
    return [{'label': ct, 'value': ct} for ct in celltypes], celltypes[0]

def get_gene_type(df, gene):
    # Get gene type label
    is_ligand = (df.ligand == gene).sum() > 0
    is_receptor = (df.receptor == gene).sum() > 0
    if is_ligand and is_receptor:
        return "Ligand/Receptor"
    elif is_ligand:
        return "Ligand"
    elif is_receptor:
        return "Receptor"
    else:
        return "Gene"


def make_plot_from_graph(G: nx.DiGraph, orig_df, layout="Planar Layout") -> go.Figure:
    scaleratio = 1.5
    height = 800

    if layout == 'Spring Layout':
        pos = nx.spring_layout(G, weight='weight', iterations=150, seed=1234567890, threshold=1e-6)
    elif layout == "Planar Layout":
        pos = nx.planar_layout(G)
    elif layout == "Circular Layout":
        pos = nx.circular_layout(G)
    else:
        raise AssertionError("Can't recognize layout " + layout)

    celltypes = list(sorted(set(orig_df.cell_type_ligand) | set(orig_df.cell_type_receptor)))

    # possible_symbols = ['square', 'circle', 'diamond', 'pentagon', 'hexagon2', 'x', 'hexagram', 'star-triangle-up', 'circle-x', 'octagon']
    possible_colors = [colors.to_hex(c) for c in cm.tab10.colors]

    # celltype2symbols = dict()
    celltype2colors = celltype_to_colors(celltypes)
    # for i, celltype in enumerate(celltypes):
    #    celltype2symbols[celltype] = possible_symbols[i]
    genetype2symbol = {
        'Receptor': 'square',
        'Ligand': 'triangle-up',  # 'arrow-up',
        'Ligand/Receptor': 'pentagon',
        'Gene': 'circle'
    }

    cmap = cm.get_cmap('bwr')
    color_step = 255
    norm = colors.Normalize(vmin=0, vmax=color_step)
    colorscale = []
    for i in range(0, color_step):
        colorscale.append(colors.rgb2hex(colors.colorConverter.to_rgb(cmap(norm(i)))))

    # colorscale = [x.hex for x in list(Color(start_color).range_to(Color(end_color), color_step))]

    min_logfc = None
    max_logfc = None
    max_abs_logfc = None
    min_abs_logfc = None
    for edge in G.edges():
        logfc = G[edge[0]][edge[1]]['logfc']
        if max_logfc is None:
            min_logfc = max_logfc = logfc
            min_abs_logfc = max_abs_logfc = abs(logfc)
        else:
            if logfc < min_logfc:
                min_logfc = logfc
            if logfc > max_logfc:
                max_logfc = logfc
            if abs(logfc) > max_abs_logfc:
                max_abs_logfc = abs(logfc)
            if abs(logfc) < min_abs_logfc:
                min_abs_logfc = abs(logfc)

    min_color_threshold = -0.2  # min_logfc
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
        showlegend=False,
        marker=dict(
            opacity=0,
            cmin=min_color_threshold,
            cmax=max_color_threshold,
            color=[],
            colorbar=dict(
                title="Induced LogFC",
                xanchor='left',
                yanchor='top',
                thickness=25,
                titleside='right',
                x=1,
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
        logfc = G[edge[0]][edge[1]]['logfc']
        ligand = G[edge[0]][edge[1]]['ligand']
        receptor = G[edge[0]][edge[1]]['receptor']

        dif_x = x1 - x0
        dif_y = y1 - y0
        # total_len = math.hypot(dif_x, dif_y)
        angle = math.atan2(dif_y, dif_x / scaleratio)
        deltaX = EDGE_PADDING * math.cos(angle)
        deltaY = EDGE_PADDING * math.sin(angle)

        x0 += deltaX * .1
        y0 += deltaY * .1

        x1 -= deltaX * .3
        y1 -= deltaY * .3

        # edge_x += [x0, x1, None]
        # edge_y += [y0, y1, None]
        # start_x.append(x0)
        # start_y.append(y0)
        # end_x.append(x1)
        # end_y.append(y1)

        edge_x = [x0, x1, None]
        edge_y = [y0, y1, None]
        arrow_x, arrow_y = get_quiver_arrows([x0], [y0], [x1], [y1], scaleratio)

        color = get_color(logfc, min_logfc, max_logfc)
        arrow_xs = edge_x + arrow_x
        arrow_ys = edge_y + arrow_y
        arrow_width = MIN_WIDTH + (MAX_WIDTH - MIN_WIDTH) * math.sqrt(
            (abs(logfc) - min_abs_logfc) / (max_abs_logfc - min_abs_logfc))
        edge_trace = go.Scatter(
            x=arrow_xs, y=arrow_ys,
            fill="toself",
            fillcolor=color,
            showlegend=False,
            line=dict(
                width=arrow_width,
                color=color
            ),
            hoverinfo='none',
            mode='lines'
        )
        hover_text = f"Ligand: {ligand}<br>Receptor: {receptor}<br>logFC Effect: {logfc}"
        edge_trace.text = [hover_text] * (len(edge_x) + len(arrow_x))
        edge_traces.append(edge_trace)

        step_x = dif_x / INTERPOLATION_POINTS
        step_y = dif_y / INTERPOLATION_POINTS

        for point in [(x0 + i * step_x, y0 + i * step_y) for i in range(0, INTERPOLATION_POINTS)]:
            interpolation_trace['x'] += (point[0],)
            interpolation_trace['y'] += (point[1],)
            interpolation_trace['text'] += (hover_text,)
            interpolation_trace['marker']['color'] += (color,)

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    annotations = []
    yrange = max(node_y) - min(node_y) if len(node_y) > 0 else 0
    xrange = max(node_x) - min(node_x) if len(node_x) > 0 else 0

    node_colors = []
    genetype2nodes = dict()
    for node in G.nodes():
        step = G.nodes[node]['t']
        celltype = G.nodes[node]['celltype']
        gene = G.nodes[node]['gene']
        gene_type = G.nodes[node]['gene_type']

        if step == 0:  # Initial gene
            gene_type += '-initial'

        if gene_type not in genetype2nodes:
            genetype2nodes[gene_type] = ([], [], [], [])
        node_colors.append(celltype2colors[celltype])
        x, y = pos[node]

        xshift = 0
        yshift = 0
        yshift += 20

        align = "center"

        clicktoshow = "onoff"

        annotations.append(dict(
            x=x,
            xshift=xshift,
            y=y,
            yshift=yshift,
            xref="x",
            yref="y",
            clicktoshow=clicktoshow,
            text=G.nodes[node]['gene'],
            visible=False,
            showarrow=False,
            font=dict(
                size=14,
                color="black"
            ),
            align=align
        ))

        genetype2nodes[gene_type][0].append(x)
        genetype2nodes[gene_type][1].append(y)
        genetype2nodes[gene_type][2].append(f'Name: {gene}<br>Step: {step}')
        genetype2nodes[gene_type][3].append(celltype)

    nodes = []
    for gene_type in genetype2nodes.keys():
        gene_xs, gene_ys, texts, celltypes = genetype2nodes[gene_type]
        is_initial = gene_type.endswith('-initial')
        if is_initial:
            gene_type = gene_type.replace('-initial', '')
            symbol = genetype2symbol[gene_type] + '-dot'
        else:
            symbol = genetype2symbol[gene_type]
        # if symbol == 'arrow-up':
        # Need to adjust y-values to center the triangle
        # gene_ys = [y*.95 for y in gene_ys]

        node_trace = go.Scatter(
            name=gene_type,
            x=gene_xs, y=gene_ys,
            marker_symbol=symbol,
            mode='markers',
            hoverinfo='text',
            text=texts,
            showlegend=False,
            marker=dict(
                color=[celltype2colors[celltype] for celltype in celltypes],
                size=15 * (1.4 if 'triangle-up' in symbol else 1),  # Triangles visually appear smaller than expected
                line=dict(width=1, color='black')
            )
        )
        nodes.append(node_trace)

    legends = []
    # Dummy traces for legends
    present_celltypes = {G.nodes[node]['celltype'] for node in G.nodes}
    present_genetypes = {G.nodes[node]['gene_type'] for node in G.nodes}
    for celltype in present_celltypes:
        legends.append(go.Scatter(
            name=celltype,
            x=[None], y=[None],
            marker_symbol='circle',
            mode='markers',
            showlegend=True,
            # visible = 'legendonly',
            legendgroup='celltypes',
            legendgrouptitle=dict(
                text='Cell Type'
            ),
            marker=dict(
                color=celltype2colors[celltype],
                size=15,
                line=dict(width=1, color='black')
            )
        ))
    for gene_type in present_genetypes:
        symbol = genetype2symbol[gene_type]
        legends.append(go.Scatter(
            name=gene_type,
            x=[None], y=[None],
            marker_symbol=symbol,
            mode='markers',
            showlegend=True,
            # visible = 'legendonly',
            legendgroup='genetypes',
            legendgrouptitle=dict(
                text='Gene'
            ),
            marker=dict(
                color='white',
                size=15,
                line=dict(width=1, color='black')
            )
        ))

    if len(nodes) == 0:
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

    fig = go.Figure(data=edge_traces + nodes + [interpolation_trace] + legends,
                    layout=go.Layout(
                        title='&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Detected Interactions',
                        titlefont_size=16,
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=40, l=5, r=5, t=40),
                        autosize=False,
                        width=height * scaleratio,
                        height=height,
                        plot_bgcolor='white',
                        annotations=annotations,
                        legend=dict(
                            yanchor="top",
                            y=1,
                            xanchor="left",
                            x=1.15
                        ),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, scaleratio=scaleratio,
                                   scaleanchor="x"))
                    )
    return enhance_plotly_export(fig, height, scaleratio)


def pseudotime_interaction_propagation_graph(ct: ad.AnnData,
                                             orig_df: pd.DataFrame,
                                             seed_cell: str,
                                             seed_ligands: list,
                                             iterations: int = 4,
                                             interaction_fdr_cutoff=0.05,
                                             min_logfc=0.01,
                                             logfc_fdr_cutoff=0.05,
                                             min_expression=0,
                                             layout="Planar Layout"):
    if isinstance(seed_ligands, str):
        ligands = []
        for ligand in seed_ligands.split(","):
            ligands.append(ligand.strip())
        seed_ligands = ligands

    gene2ligand = ct.uns['gene_is_ligand']

    # Basic filters to build the network from
    df = orig_df[(orig_df.MAST_fdr_ligand < logfc_fdr_cutoff) &
                 (orig_df.MAST_fdr_receptor < logfc_fdr_cutoff) &
                 (np.abs(orig_df.MAST_log2FC_ligand) > min_logfc) &
                 (np.abs(orig_df.MAST_log2FC_receptor) > min_logfc) &
                 (np.abs(orig_df.expression_receptor) > min_expression) &
                 (np.abs(orig_df.expression_ligand) > min_expression)]

    Gs = []
    curr_G = nx.DiGraph()
    frontier = set()
    for t in range(iterations):
        if t == 0:
            for seed_ligand in seed_ligands:
                gene_type = get_gene_type(orig_df, seed_ligand)
                seed_node = f"{seed_cell}_{seed_ligand}_ligand"
                curr_G.add_node(seed_node, t=t, gene=seed_ligand, celltype=seed_cell, gene_type=gene_type,
                                diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand))
                frontier.add(seed_node)
        else:
            new_frontier = set()
            for node in frontier:
                celltype = curr_G.nodes[node]['celltype']
                gene_type = curr_G.nodes[node]['gene_type']

                if 'Ligand' in gene_type:
                    ligand = curr_G.nodes[node]['gene']
                    receivers = df[(df.ligand == ligand) & (df.cell_type_ligand == celltype)]
                    for i, row in receivers.iterrows():
                        receptor = row['receptor']
                        receptor_celltype = row['cell_type_receptor']
                        next_node = f"{receptor_celltype}_{receptor}_receptor"

                        if get_interaction_fdr(ct, celltype, ligand, receptor) < interaction_fdr_cutoff:
                            continue

                        if not curr_G.has_node(next_node):
                            curr_G.add_node(next_node, t=t,
                                            gene=receptor,
                                            celltype=receptor_celltype,
                                            gene_type=get_gene_type(orig_df, receptor),
                                            diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand))
                        if not curr_G.has_edge(node, next_node) and not curr_G.has_edge(next_node, node):
                            logfc = get_effective_logfc(row['MAST_log2FC_ligand'], row["MAST_log2FC_receptor"])
                            curr_G.add_edge(node, next_node,
                                            t=t, ligand=node,
                                            receptor=next_node,
                                            logfc=logfc,
                                            weight=abs(logfc),
                                            numDEG=row['numDEG_fdr05_receptor'],
                                            numInteractions=row['numSigI1_fdr05_receptor'])
                            new_frontier.add(next_node)

                if 'Receptor' in gene_type:
                    receptor = curr_G.nodes[node]['gene']
                    down_ligands = get_downstream_ligands(ct, receptor, celltype, min_logfc,
                                                          logfc_fdr_cutoff)
                    senders = df[(df.receptor == receptor) & (df.cell_type_receptor == celltype)]
                    for i, row in senders.iterrows():
                        ligand = row['ligand']
                        if ligand not in down_ligands:
                            continue
                        ligand_celltype = row['cell_type_ligand']
                        next_node = f"{ligand_celltype}_{ligand}_ligand"

                        if get_interaction_fdr(ct, celltype, receptor, ligand) < interaction_fdr_cutoff:
                            continue

                        if not curr_G.has_node(next_node):
                            curr_G.add_node(next_node, t=t,
                                            gene=ligand,
                                            celltype=ligand_celltype,
                                            gene_type=get_gene_type(orig_df, ligand),
                                            diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand))
                        if not curr_G.has_edge(node, next_node) and not curr_G.has_edge(next_node, node):
                            logfc = get_effective_logfc(row["MAST_log2FC_receptor"], row['MAST_log2FC_ligand'])
                            curr_G.add_edge(node, next_node,
                                            t=t, ligand=node,
                                            receptor=next_node,
                                            logfc=logfc,
                                            weight=abs(logfc),
                                            numDEG=row['numDEG_fdr05_receptor'],
                                            numInteractions=row['numSigI1_fdr05_receptor'])
                            new_frontier.add(next_node)

            frontier = new_frontier
        Gs.append(curr_G.copy())

    return make_plot_from_graph(Gs[-1], df, layout=layout)


layout = [
    interactive_panel(wrap_icon('fa-maximize', 'Cell Type to Cell Type Interactions'),
                      *build_interface()
                      )
]
