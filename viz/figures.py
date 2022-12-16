import bisect
import itertools
import math
from dataclasses import dataclass, field
from io import BytesIO
from typing import List, Tuple, Optional, Set

import anndata as ad
import networkx as nx
import numpy as np
from matplotlib import cm, colors
import pandas as pd
from PIL import Image
import plotly.graph_objects as go

from viz.data import get_interaction_fdr, get_diff_abundance, get_downstream_ligands, get_interaction_logfc, \
    get_interaction_logfc_fdr
from viz.util import multipartite_layout, get_quiver_arrows, celltype_to_colors, timeline_layout


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

    selected = selected[(selected.numSigI1_fdr_receptor >= numInteractions) &
                        (selected["MAST_fdr_ligand"] < logfc_fdr_bipartite_cutoff) &
                        (selected["MAST_fdr_receptor"] < logfc_fdr_bipartite_cutoff) &
                        (np.abs(selected["MAST_log2FC_ligand"]) > min_logfc_bipartite) &
                        (np.abs(selected["MAST_log2FC_receptor"]) > min_logfc_bipartite) &
                        (np.abs(selected["expression_receptor"]) > min_expression_bipartite) &
                        (np.abs(selected["expression_ligand"]) > min_expression_bipartite)].sort_values(
        by="numDEG_fdr_receptor", ascending=False)

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
                   numInteractions=row['numSigI1_fdr_receptor'],
                   ligand_fdr=row["MAST_fdr_ligand"],
                   ligand_logfc=row["MAST_log2FC_ligand"],
                   receptor_fdr=row["MAST_fdr_receptor"],
                   receptor_logfc=row["MAST_log2FC_receptor"],
                   numDEG=row['numDEG_fdr_receptor'])

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


def get_all_edge_attr_for_node(G, node, attr):
    successor_attrs = [G[node][s][attr] for s in G.successors(node)]
    predecessor_attrs = [G[p][node][attr] for p in G.predecessors(node)]
    return successor_attrs + predecessor_attrs


def make_plot_from_graph(G: nx.DiGraph, orig_df, layout="planar") -> go.Figure:
    scaleratio = 1.5
    height = 800

    layout = layout.lower()
    if 'spring' in layout:
        pos = nx.spring_layout(G, weight='weight', iterations=150, seed=1234567890, threshold=1e-6)
    elif "planar" in layout:
        pos = nx.planar_layout(G)
    elif "circular" in layout:
        pos = nx.circular_layout(G)
    elif "multipartite" in layout or 'timeline' in layout:
        pos = timeline_layout(G, step_attr='t', scaleratio=scaleratio)
    else:
        raise AssertionError("Can't recognize layout " + layout)

    celltypes = list(sorted(set(orig_df.cell_type_ligand) | set(orig_df.cell_type_receptor)))

    # possible_symbols = ['square', 'circle', 'diamond', 'pentagon', 'hexagon2', 'x', 'hexagram', 'star-triangle-up', 'circle-x', 'octagon']
    # possible_colors = [colors.to_hex(c) for c in cm.tab10.colors]

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
        logfc = G[edge[0]][edge[1]]['response_logfc']
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

    # Can't hover over lines directly, so add invisible points
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
                title="Ligand LogFC",
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
        logfc = G[edge[0]][edge[1]]['response_logfc']
        ligand = G[edge[0]][edge[1]]['ligand']
        receptor = G[edge[0]][edge[1]]['receptor']

        dif_x = x1 - x0
        dif_y = y1 - y0
        # total_len = math.hypot(dif_x, dif_y)
        angle = math.atan2(dif_y, dif_x / scaleratio)
        deltaX = EDGE_PADDING * math.cos(angle)
        deltaY = EDGE_PADDING * math.sin(angle)

        x0 += deltaX * .5
        y0 += deltaY * .5

        x1 -= deltaX * .5
        y1 -= deltaY * .5

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
            (abs(logfc) - min_abs_logfc) / (max_abs_logfc - min_abs_logfc + 1e-8))
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
        hover_text = f"Ligand: {ligand}<br>Receptor: {receptor}<br>Ligand LogFC Effect: {logfc}"
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

        ligand_logfc = max(get_all_edge_attr_for_node(G, node, "ligand_logfc") + [0], key=abs)
        receptor_logfc = max(get_all_edge_attr_for_node(G, node, "receptor_logfc") + [0], key=abs)

        if step == 0:  # Initial gene
            gene_type += '-initial'

        if gene_type not in genetype2nodes:
            genetype2nodes[gene_type] = ([], [], [], [], [], [])
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
            visible=True,
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
        genetype2nodes[gene_type][4].append(ligand_logfc)
        genetype2nodes[gene_type][5].append(receptor_logfc)

    # Get the size range
    max_node_logfc = max(nx.get_edge_attributes(G, 'weight').values())
    min_node_logfc = min(nx.get_edge_attributes(G, 'weight').values())

    MIN_NODE_SIZE = 8
    MAX_NODE_SIZE = 27
    NODE_SIZE_RANGE = MAX_NODE_SIZE - MIN_NODE_SIZE

    nodes = []
    for gene_type in genetype2nodes.keys():
        gene_xs, gene_ys, texts, celltypes, ligand_logfc, receptor_logfc = genetype2nodes[gene_type]

        is_initial = gene_type.endswith('-initial')
        if is_initial:
            gene_type = gene_type.replace('-initial', '')
            symbol = genetype2symbol[gene_type] + '-dot'
        else:
            symbol = genetype2symbol[gene_type]
        # if symbol == 'arrow-up':
        # Need to adjust y-values to center the triangle
        # gene_ys = [y*.95 for y in gene_ys]

        node_sizes = [(MIN_NODE_SIZE + (NODE_SIZE_RANGE*(abs(logfc)/max_node_logfc))) * (1.25 if 'triangle-up' in symbol else 1) for logfc in receptor_logfc]  # Triangles visually appear smaller than expected

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
                opacity=1,
                size=node_sizes,
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
    # legend for node size
    legends.append(go.Scatter(
        name="{:.0E}".format(max_node_logfc),
        x=[None], y=[None],
        marker_symbol='circle',
        mode='markers',
        showlegend=True,
        # visible = 'legendonly',
        legendgroup='nodesize',
        legendgrouptitle=dict(
            text='Ligand-Induced abs(LogFC)'
        ),
        marker=dict(
            color='grey',
            size=MAX_NODE_SIZE,
            line=dict(width=1, color='black')
        )
    ))
    legends.append(go.Scatter(
        name="{:.0E}".format(min_node_logfc),
        x=[None], y=[None],
        marker_symbol='circle',
        mode='markers',
        showlegend=True,
        # visible = 'legendonly',
        legendgroup='nodesize',
        legendgrouptitle=dict(
            text='Ligand-Induced abs(LogFC)'
        ),
        marker=dict(
            color='grey',
            size=MIN_NODE_SIZE,
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
    #return enhance_plotly_export(fig, height, scaleratio)
    return fig


def pseudotime_interaction_propagation_graph(ct: ad.AnnData,
                                             orig_df: pd.DataFrame,
                                             seed_cell: str,
                                             seed_ligands: list,
                                             iterations: int = 4,
                                             interaction_fdr_cutoff=0.05,
                                             min_logfc=0.01,
                                             logfc_fdr_cutoff=0.05,
                                             min_expression=0,
                                             layout="Planar Layout",
                                             set_progress_callback=None):
    if isinstance(seed_ligands, str):
        ligands = []
        for ligand in seed_ligands.split(","):
            ligands.append(ligand.strip())
        seed_ligands = ligands

    gene2ligand = ct.uns['gene_is_ligand']

    # Basic filters to build the network from
    df = orig_df[(orig_df.MAST_fdr_ligand <= logfc_fdr_cutoff) &
                 (orig_df.MAST_fdr_receptor <= logfc_fdr_cutoff) &
                 (np.abs(orig_df.MAST_log2FC_ligand) >= min_logfc) &
                 (np.abs(orig_df.MAST_log2FC_receptor) >= min_logfc) &
                 (np.abs(orig_df.expression_receptor) >= min_expression) &
                 (np.abs(orig_df.expression_ligand) >= min_expression)]

    Gs = []
    curr_G = nx.DiGraph()
    frontier = set()
    for t in range(iterations):
        if t == 0:
            for seed_ligand in seed_ligands:
                gene_type = get_gene_type(orig_df, seed_ligand)
                seed_node = f"{seed_cell}_{seed_ligand}_ligand"
                curr_G.add_node(seed_node, t=t, gene=seed_ligand, celltype=seed_cell, gene_type=gene_type,
                                #diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand)
                                )
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

                        if get_interaction_fdr(ct, celltype, ligand, receptor) > interaction_fdr_cutoff:
                            continue

                        if not curr_G.has_node(next_node):
                            curr_G.add_node(next_node, t=t,
                                            gene=receptor,
                                            celltype=receptor_celltype,
                                            gene_type=get_gene_type(orig_df, receptor),
                                            #diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand)
                                            )
                        if not curr_G.has_edge(node, next_node):
                            curr_G.add_edge(node, next_node,
                                            t=t, ligand=node,
                                            receptor=next_node,
                                            response_logfc=row['MAST_log2FC_receptor'],
                                            ligand_logfc=row['MAST_log2FC_ligand'],
                                            receptor_logfc=row['MAST_log2FC_receptor'],
                                            weight=max(abs(row['MAST_log2FC_ligand']), abs(row['MAST_log2FC_receptor'])),
                                            numDEG=row['numDEG_fdr_receptor'],
                                            numInteractions=row['numSigI1_fdr_receptor'])
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

                        if get_interaction_fdr(ct, celltype, receptor, ligand) > interaction_fdr_cutoff:
                            continue

                        if not curr_G.has_node(next_node):
                            curr_G.add_node(next_node, t=t,
                                            gene=ligand,
                                            celltype=ligand_celltype,
                                            gene_type=get_gene_type(orig_df, ligand),
                                            #diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand)
                                            )
                        if not curr_G.has_edge(node, next_node):
                            curr_G.add_edge(node, next_node,
                                            t=t, ligand=node,
                                            receptor=next_node,
                                            response_logfc=row['MAST_log2FC_ligand'],
                                            ligand_logfc=row['MAST_log2FC_ligand'],
                                            receptor_logfc=row['MAST_log2FC_receptor'],
                                            weight=max(abs(row['MAST_log2FC_ligand']), abs(row['MAST_log2FC_receptor'])),
                                            numDEG=row['numDEG_fdr_receptor'],
                                            numInteractions=row['numSigI1_fdr_receptor'])
                            new_frontier.add(next_node)

            frontier = new_frontier
        Gs.append(curr_G.copy())
        if set_progress_callback:
            set_progress_callback(((t+1)//iterations * 100, 100))
        if len(frontier) == 0:
            break

    return make_plot_from_graph(Gs[-1], df, layout=layout)


def polar2cartesian(r, theta):
    theta = np.deg2rad(theta)
    return r * np.cos(theta), r * np.sin(theta)


@dataclass(unsafe_hash=True)
class PlotRing:
    name: str
    color: str
    height: float
    circumference: float
    children: List['PlotRing'] = field(default_factory=list)
    theta: float = 0  # computed later
    start_height: float = 0  # computed later

    def add_child(self, child):
        # Add child in sorted order
        bisect.insort(self.children, child, key=lambda x: (x.height, x.circumference))

    def get_all_children(self):
        return [self] + list(itertools.chain.from_iterable([c.get_all_children() for c in self.children]))

    def get_max_height(self):
        return self.height + max([c.get_max_height() for c in self.children], default=0)

    def get_number_of_levels(self):
        return max([c.get_number_of_levels() for c in self.children], default=0) + 1

    def iter_rings(self):
        # Breadth-first iteration to get list of lists
        frontier = [self]
        for level in range(self.get_number_of_levels()):
            yield frontier
            frontier = list(itertools.chain.from_iterable([c.children for c in frontier]))
        if len(frontier) > 0:
            raise Exception("Something went wrong")

    def iter_hierarchy(self):
        # Depth-first iteration of (parent, children) tuples
        yield self, self.children
        for child in self.children:
            yield from child.iter_hierarchy()


def plot_polar_data(rings: PlotRing, center_height: float) -> go.Figure:
    # Normalize the data
    # First ring is always the full circle
    rings.circumference = 360
    # Normalize the heights to be bounded by 0 and 1
    max_height = rings.get_max_height()
    # Reset center node height
    rings.height = center_height*max_height
    for ring_list in rings.iter_rings():
        for r in ring_list:
            r.height /= max_height
    # Calculate circumference of each ring group
    for parent, children in rings.iter_hierarchy():
        total_circumference = sum([c.circumference for c in children])
        curr_theta = parent.theta
        last_child = parent
        for child in children:
            child.circumference = parent.circumference * (child.circumference / total_circumference)
            child.theta = curr_theta
            child.start_height = parent.start_height + parent.height
            curr_theta += child.circumference/2 + last_child.circumference/2
            last_child = child

    # Build the figure
    annotations = []
    traces = []
    for level, ring_list in enumerate(rings.iter_rings()):
        if level == 0:  # Handle first circle
            assert len(ring_list) == 1
            node = ring_list[0]
            annotations.append(dict(
                x=0, y=0,
                xref='x', yref='y',
                clicktoshow='onoff',
                text=node.name,
                visible=True,
                showarrow=False,
                font=dict(
                    size=25,
                    color="black"
                ),
                align='center'
            ))
            traces.append(go.Barpolar(
                r=[node.height],
                theta=[0],
                width=[360],
                thetaunit='degrees',
                opacity=1,
                base=0,
                text=node.name,
                marker=dict(
                    color=[node.color],
                    line=dict(
                        color='black',
                        width=1
                    )
                )
            ))
        else:
            thetas = []
            circumferences = []
            bases = []
            colors = []
            heights = []
            names = []
            for node in ring_list:
                base = node.start_height
                theta = node.theta
                height = node.height
                circumference = node.circumference
                color = node.color
                name = node.name

                # Label
                annotation_distance = 2*height/3
                annotation_x, annotation_y = polar2cartesian(annotation_distance, theta)
                annotations.append(dict(
                    x=annotation_x, y=annotation_y,
                    clicktoshow='onoff',
                    text=name,
                    visible=level < 3,
                    showarrow=False,
                    font=dict(
                        size=20,
                        color="black"
                    ),
                    align='center'
                ))

                # Ring info
                thetas.append(theta)
                circumferences.append(circumference)
                bases.append(base)
                colors.append(color)
                heights.append(height)
                names.append(name)

            traces.append(go.Barpolar(
                r0=0,
                r=heights,
                theta0=0,
                theta=thetas,
                width=circumferences,
                thetaunit='degrees',
                opacity=1,
                base=bases,
                text=names,
                marker=dict(
                    color=colors,
                    line=dict(
                        color='black',
                        width=1
                    )
                )
            ))
        
    fig = go.Figure(data=traces, layout=go.Layout(
        template=None,
        showlegend=False,
        hovermode='closest',
        margin=dict(l=0, r=0, t=0, b=0),
        autosize=False,
        barmode='overlay',
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        annotations=annotations,
        polar=dict(
            bgcolor='rgba(0,0,0,0)',
            radialaxis=dict(
                showticklabels=False,
                ticks='',
                range=[0, 1.5],
                showgrid=False,
                showline=False,
                angle=0,
                tick0=0,
                tickangle=90
            ),
            angularaxis=dict(
                ticks='',
                showgrid=False,
                rotation=0,
                tick0=0,
                tickangle=90,
                direction='clockwise'
            ),
            xaxis=dict(
                showgrid=False,
                zeroline=False,
                showticklabels=False,
                range=[-1, 1]
            ),
            yaxis=dict(
                showgrid=False,
                zeroline=False,
                showticklabels=False,
                range=[-1, 1],
                scaleanchor="x"
            )
        )
    ))
    return fig


def radial_ligand_effect_figure(
        ct: ad.AnnData, interactions: pd.DataFrame,
        ligand: str, celltype: str,
        min_logfc: float = 0, max_fdr: float = 0.05,
        max_levels: int = 3, celltype_filter: Set[str] = None
) -> go.Figure:
    all_celltypes = set(ct.obs['cell type'].unique())
    if celltype_filter is None:
        celltype_filter = all_celltypes
    celltype2color = celltype_to_colors(all_celltypes)

    interactions = interactions[(interactions.MAST_log2FC_ligand > min_logfc) &
                                (interactions.MAST_fdr_ligand <= max_fdr) &
                                (interactions.celltype_receptor.isin(celltype_filter)) &
                                (interactions.celltype_ligand.isin(celltype_filter))]

    rings = PlotRing(ligand, celltype2color[celltype], 0, 360)

    def build_downstream_rings(i: int,
                               curr_ring: PlotRing,
                               start: str,
                               start_celltype: str,
                               is_ligand: bool):
        if i >= max_levels:
            return

        if is_ligand:
            # Get and map receptors
            filtered_interactions = interactions[(interactions.ligand == start)]
            for i, row in filtered_interactions.iterrows():
                receptor = row['receptor']
                receptor_celltype = row['cell type_receptor']
                next_ring = PlotRing(
                    f"{receptor}+", celltype2color[receptor_celltype],
                    row['numDEG_fdr25_receptor'], row['numSigI1_fdr25_receptor']
                )
                curr_ring.add_child(next_ring)
                build_downstream_rings(i+1, next_ring, receptor, receptor_celltype, False)
        else:
            # Get and map ligands
            downstream_ligands = get_downstream_ligands(ct, start, start_celltype, min_logfc, max_fdr)
            for lig in downstream_ligands:
                interaction_logfc = get_interaction_logfc(ct, start_celltype, start, lig)
                interaction_fdr = get_interaction_fdr(ct, start_celltype, start, lig)
                interaction_logfc_fdr = get_interaction_logfc_fdr(ct, start_celltype, start, lig)
                if np.isnan(interaction_fdr) or interaction_fdr > max_fdr or interaction_logfc < min_logfc:
                    continue
                next_ring = PlotRing(
                    lig, celltype2color[start_celltype],
                    interaction_logfc, interaction_logfc
                )
                curr_ring.add_child(next_ring)
                build_downstream_rings(i+1, next_ring, lig, start_celltype, True)

    build_downstream_rings(0, rings, ligand, celltype, True)

    # TODO: TEST
    return plot_polar_data(rings, .2)

def polar_receptor_figure(ct: ad.AnnData, interactions: pd.DataFrame,
                          celltype: str, receptor: Optional[str] = None,
                          ligands: Optional[List[str]] = None,
                          min_logfc: float = 0, max_fdr: float = 0.05,
                          main_node_color: str = "white", main_node_height: float = .33) -> Tuple[go.Figure, List[str]]:
    assert min_logfc >= 0

    #selected_inter = interactions[(interactions.receptor == receptor) & (interactions.cell_type_receptor == celltype)]
    #if selected_inter.shape[0] == 0:
    #    return None
    if ligands is None:
        downstream_ligands = get_downstream_ligands(ct, receptor, celltype, min_logfc, max_fdr)
    else:
        downstream_ligands = set(ligands)
    # Require all ligands to have passed previous filters
    #downstream_ligands = [ligand for ligand in downstream_ligands if ligand in selected_inter.ligand.values]
    downstream_ligand_info = dict(
        ligand=list(),
        induced_logfc=list(),
        interaction_effect=list(),
        logfc_fdr=list()
    )
    for ligand in reversed(downstream_ligands):
        if receptor is None:
            interaction_effect = interactions[(interactions.ligand == ligand) & (interactions.cell_type_ligand == celltype)].MAST_log2FC_ligand.max()
            interaction_fdr = 0
            interaction_logfc_fdr = interactions[(interactions.ligand == ligand) & (interactions.cell_type_ligand == celltype) & (interactions.MAST_log2FC_ligand == interaction_effect)].MAST_fdr_ligand.min()
        else:
            interaction_effect = get_interaction_logfc(ct, celltype, receptor, ligand)
            interaction_fdr = get_interaction_fdr(ct, celltype, receptor, ligand)
            interaction_logfc_fdr = get_interaction_logfc_fdr(ct, celltype, receptor, ligand)
        if interaction_fdr > max_fdr or interaction_effect < min_logfc or np.isnan(interaction_fdr):
            continue
        downstream_ligand_info['ligand'].append(ligand)
        downstream_ligand_info['induced_logfc'].append(interaction_effect)
        downstream_ligand_info['interaction_effect'].append(max(-1*np.log(interaction_fdr + 1e-5), 1e-5))
        downstream_ligand_info['logfc_fdr'].append(max(-1*np.log(interaction_logfc_fdr + 1e-5), 1e-5))
    downstream_ligand_info = pd.DataFrame(downstream_ligand_info)

    origin = polar2cartesian(0, 0)

    annotations = []

    traces = []
    if downstream_ligand_info.shape[0] > 0:
        max_logfc = downstream_ligand_info.induced_logfc.max()
        downstream_ligand_info['induced_logfc'] = (downstream_ligand_info['induced_logfc']+.05) / (max_logfc+0.05)
        # Transform interaction effect
        downstream_ligand_info['interaction_effect'] = downstream_ligand_info.interaction_effect / downstream_ligand_info.interaction_effect.sum()
        downstream_ligand_info['logfc_fdr'] = downstream_ligand_info.logfc_fdr / downstream_ligand_info.logfc_fdr.sum()

        downstream_ligand_info['interaction_rotation'] = downstream_ligand_info.interaction_effect * 360
        downstream_ligand_info['interaction_rotation'] = downstream_ligand_info.interaction_rotation - ((downstream_ligand_info.interaction_rotation.sum() - 360)/downstream_ligand_info.shape[0])/2
        downstream_ligand_info['start_angles'] = np.array([0] + list(downstream_ligand_info.interaction_rotation.cumsum())[:-1])
        downstream_ligand_info['end_angles'] = downstream_ligand_info.interaction_rotation.cumsum()
        downstream_ligand_info['midpoint_angles'] = (downstream_ligand_info.start_angles + downstream_ligand_info.end_angles) / 2

    if receptor is not None:
        annotations.append(dict(
           x=origin[0],
           y=origin[1],
           xref="x",
           yref="y",
           clicktoshow="onoff",
           text=receptor,
           visible=True,
           showarrow=False,
           font=dict(
               size=30,
               color="black"
           ),
           align="center"
        ))

    theta_steps = []
    for i, row in downstream_ligand_info.iterrows():
        if i == 0:
            theta_step = 0
        else:
            last_theta = theta_steps[-1]
            last_theta_step = downstream_ligand_info.interaction_rotation.iloc[i-1]/2
            curr_theta_step = downstream_ligand_info.interaction_rotation.iloc[i]/2
            theta_step = last_theta + last_theta_step + curr_theta_step
        theta_steps.append(theta_step)

        ligand = row.ligand
        height = row.induced_logfc
        regularized_distance = 2*height/3   # * 0.5 + main_node_height
        x, y = polar2cartesian(regularized_distance, theta_step)

        annotations.append(dict(
            x=x,
            y=y,
            clicktoshow="onoff",
            text=ligand,
            visible=True,
            showarrow=False,
            font=dict(
                size=30,
                color="black"
            ),
            align="center",
           # textangle=start_angle - 90
        ))

        dict(
            x=origin[0],
            y=origin[1],
            xref="x",
            yref="y",
            clicktoshow="onoff",
            text=f"{receptor}+",
            visible=True,
            showarrow=False,
            font=dict(
                size=14,
                color="black"
            ),
            align="center"
        )

    # Main figure
    if downstream_ligand_info.shape[0] > 0:
        traces.append(go.Barpolar(
            r0=0,
            r=list(downstream_ligand_info.induced_logfc),
            theta0=0,
            #theta=[0]*downstream_ligand_info.shape[0],
            theta=theta_steps,
            width=list(downstream_ligand_info.interaction_rotation),
            thetaunit="degrees",
            opacity=1,
            base=main_node_height,
            text=list(downstream_ligand_info.ligand),
            marker=dict(
                color=(['grey']*downstream_ligand_info.shape[0]),
                line=dict(
                    color='black',
                    width=2
                )
            )
        ))

    traces.append(go.Barpolar(
        r=[main_node_height],
        theta=[0],
        width=[360],
        thetaunit="degrees",
        opacity=1,
        base=0,
        text=[f"{celltype} {receptor}"],
        marker=dict(
            color=[main_node_color],
            line=dict(
                color='black',
                width=2
            )
        )
    ))

    fig = go.Figure(
        traces,
        layout=go.Layout(
            #title=f'&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;{receptor}+ {celltype} Effects',
            #titlefont_size=16,
            template=None,
            showlegend=False,
            hovermode='closest',
            margin=dict(b=0, l=0, r=0, t=0),
            autosize=False,
            barmode='overlay',
            #width=height * scaleratio,
            #height=height,
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            annotations=annotations,
            polar=dict(
                bgcolor='rgba(0,0,0,0)',
                radialaxis=dict(
                    showticklabels=False,
                    ticks='',
                    range=[0, 1.5],
                    showgrid=False,
                    showline=False,
                    angle=0,
                    tick0=0,
                    tickangle=90,
                ),
                angularaxis=dict(
                    showticklabels=False,
                    ticks='',
                    showgrid=False,
                    rotation=0,
                    tick0=0,
                    tickangle=90,
                )
            ),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-1, 1]),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False,
                       scaleanchor="x", range=[-1, 1])
        )
    )

    return fig, downstream_ligand_info['ligand'].tolist()


def cascading_effects_figure(ct: ad.AnnData, interactions: pd.DataFrame,
                             starting_ligands: List[str], celltype: str,
                             min_logfc: float = 0, max_fdr: float = 0.05,
                             numSigI1: int = 0, celltype_filters=None, iterations: int = 15):
    if celltype_filters is None:
        celltype_filters = []

    interactions = interactions[(interactions.MAST_log2FC_ligand > min_logfc) & (interactions.MAST_fdr_ligand <= max_fdr) & (interactions.numSigI1_fdr_receptor >= numSigI1)]

    G = nx.DiGraph()
    celltype2colors = celltype_to_colors(interactions.cell_type_receptor.unique().tolist())
    frontier = set()
    for ligand in starting_ligands:
        G.add_node(f"{ligand}_START", ligands=[ligand], fig=None, receptor=None, celltype=celltype, step=0)
        frontier.add(f"{ligand}_START")
    nodes_ignored = set()
    from tqdm import tqdm
    for i in tqdm(range(iterations)):
        new_frontier = set()
        for node in frontier:
            # We are assuming that the interactions dataframe is already filtered
            ligands = G.nodes[node]['ligands']
            cell = G.nodes[node]['celltype']
            selected_interactions = interactions[interactions.ligand.isin(ligands) & (interactions.cell_type_ligand == cell)]
            for _, row in selected_interactions.iterrows():
                new_node = f"{row.cell_type_receptor}_{row.receptor}"
                if row.cell_type_receptor not in celltype_filters or G.has_edge(node, new_node) or new_node in nodes_ignored:
                    continue

                if not G.has_node(new_node):
                    node_fig, ligands = polar_receptor_figure(ct, interactions, row.cell_type_receptor,
                                                              row.receptor, None, min_logfc, max_fdr,
                                                              main_node_color=celltype2colors[row.cell_type_receptor])

                    if len(ligands) == 0:
                        nodes_ignored.add(new_node)
                        continue

                    G.add_node(new_node,
                               receptor=row.receptor,
                               fig=node_fig,
                               ligands=ligands,
                               celltype=row.cell_type_receptor, step=i+1)
                    new_frontier.add(new_node)

                G.add_edge(node, new_node,
                           logfc=row.MAST_log2FC_ligand,
                           weight=abs(row.MAST_log2FC_ligand),
                           ligand=row.ligand,
                           receptor=row.receptor)

        frontier = new_frontier
        if len(frontier) == 0:
            break

    layout = timeline_layout(G, step_attr='step', scaleratio=1.5)
    #layout = nx.spring_layout(G)
    #shells = [[]]*(max([G.nodes[node]['step'] for node in G.nodes])+1)
    #for node in G.nodes:
    #    shells[G.nodes[node]['step']].append(node)
    #layout = nx.shell_layout(G, shells)

    cmap = cm.get_cmap('viridis')
    color_step = 255
    norm = colors.Normalize(vmin=0, vmax=color_step)
    colorscale = []
    for i in range(0, color_step):
        colorscale.append(colors.rgb2hex(colors.colorConverter.to_rgb(cmap(norm(i)))))
    # Reverse colorscale to have dark colors = high values
    colorscale = colorscale[::-1]

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

    min_color_threshold = 0  # min_logfc
    max_color_threshold = max_logfc  # max_logfc

    def get_color(value, minimum, maximum):
        step_size = (max_color_threshold - min_color_threshold) / (color_step - 1)
        index = (value - min_color_threshold) // step_size
        return colorscale[min(max(int(index), 0), color_step - 1)]


    # Can't hover over lines directly, so add invisible points
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
                title="Ligand LogFC",
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
    EDGE_PADDING = .05
    MIN_WIDTH = 2
    MAX_WIDTH = 6
    INTERPOLATION_POINTS = 20

    for edge in G.edges():
        x0, y0 = layout[edge[0]]
        x1, y1 = layout[edge[1]]
        logfc = G[edge[0]][edge[1]]['logfc']
        ligand = G[edge[0]][edge[1]]['ligand']
        receptor = G[edge[0]][edge[1]]['receptor']

        dif_x = x1 - x0
        dif_y = y1 - y0
        # total_len = math.hypot(dif_x, dif_y)
        angle = math.atan2(dif_y, dif_x)
        deltaX = EDGE_PADDING * math.cos(angle)
        deltaY = EDGE_PADDING * math.sin(angle)

        x0 += deltaX * .25
        y0 += deltaY * .25

        x1 -= deltaX * .75
        y1 -= deltaY * .75

        # edge_x += [x0, x1, None]
        # edge_y += [y0, y1, None]
        # start_x.append(x0)
        # start_y.append(y0)
        # end_x.append(x1)
        # end_y.append(y1)

        edge_x = [x0, x1, None]
        edge_y = [y0, y1, None]
        arrow_x, arrow_y = get_quiver_arrows([x0], [y0], [x1], [y1], 1)

        color = get_color(logfc, min_logfc, max_logfc)
        arrow_xs = edge_x + arrow_x
        arrow_ys = edge_y + arrow_y
        arrow_width = MIN_WIDTH + (MAX_WIDTH - MIN_WIDTH) * math.sqrt(
            (abs(logfc) - min_abs_logfc) / (max_abs_logfc - min_abs_logfc + 1e-8))
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
        hover_text = f"Ligand: {ligand}<br>Receptor: {receptor}<br>Ligand LogFC Effect: {logfc}"
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
    initial_x = []
    initial_y = []
    images = []
    annotations = []
    for node in G.nodes():
        x, y = layout[node]
        if G.nodes[node]['fig'] is None:
            initial_x.append(x)
            initial_y.append(y)
            annotations.append(dict(
                x=x,
                y=y,
                xref='x',
                yref='y',
                clicktoshow='onoff',
                text=G.nodes[node]['ligands'][0],
                visible=True,
                showarrow=False,
                font=dict(
                    size=14,
                    color='black'
                ),
                align='center'
            ))
        else:
            node_x.append(x)
            node_y.append(y)
            images.append(dict(
                source=Image.open(BytesIO(G.nodes[node]['fig'].to_image(format="png", scale=2)), formats=['png']),
                xref="x",
                yref="y",
                xanchor="center",
                yanchor="middle",
                x=x,
                y=y,
                sizex=.66,
                sizey=.66,
                sizing="contain",
                opacity=1,
                #layer="above",
            ))

    legends = []
    for celltype in {G.nodes[node]['celltype'] for node in G.nodes}:
        legends.append(go.Scatter(
            name=celltype,
            x=[None],
            y=[None],
            marker_symbol='circle',
            mode='markers',
            showlegend=True,
            legendgroup='celltypes',
            legendgrouptitle=dict(
                text='Cell Type',
            ),
            marker=dict(
                color=celltype2colors[celltype],
                size=15,
                line=dict(width=1, color='black')
            )
        ))

    initial_trace = go.Scatter(
        name="nodes",
        x=initial_x,
        y=initial_y,
        mode='markers',
        showlegend=False,
        marker=dict(
            color=[celltype2colors[G.nodes[node]['celltype']] for node in G.nodes() if G.nodes[node]['fig'] is None],
            size=50,
            line=dict(width=1, color='black')
        )
    )

    node_trace = go.Scatter(
        name="nodes",
        x=node_x, y=node_y,
        mode='markers',
        showlegend=False,
        marker=dict(
            color='rgba(0,0,0,0)',
            size=100,
            line=dict(width=1, color='rgba(0,0,0,0)')
        )
    )

    fig = go.Figure(data=[interpolation_trace] + edge_traces + [initial_trace, node_trace] + legends,
                    layout=go.Layout(
                        title="&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Receptor Network",
                        titlefont_size=16,
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=30, l=5, r=5, t=30),
                        autosize=False,
                        width=1000,
                        height=800,
                        plot_bgcolor=None,
                        images=images,
                        annotations=annotations,
                        legend=dict(
                            yanchor="top",
                            y=1,
                            xanchor="left",
                            x=1.15
                        ),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, scaleratio=1, scaleanchor="x"))
                    )

    return fig
