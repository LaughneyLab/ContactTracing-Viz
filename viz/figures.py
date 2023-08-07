import bisect
import itertools
import math
from dataclasses import dataclass, field
from functools import total_ordering
from typing import List

import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from viz.data import read_ligand_receptor_file, read_ligand_effect_for
from viz.util import multipartite_layout, get_quiver_arrows, celltype_to_colors, timeline_layout, ColorTransformer

DEFAULT_LIGAND_EFFECT_ARGS = {
    'effect_set': 'cin',
    'cell_type': 'Tumor cells',
    'ligands': 'Ccl2,Apoe',
    'interaction_fdr': 0.05,
    'min_logfc': 0.12,
    'min_expression': 0.0,
    'logfc_fdr': 0.05,
    'iterations': 2
}
LIGAND_EFFECT_SAVE_LOCATION = "compiled_figs/default_ligand_effects.pkl"

DEFAULT_INTERACTIONS_ARGS = {
    'inter_set': 'max',
    'min_numsigi1_bipartite': 10,
    'min_logfc_bipartite': 0.12,
    'min_expression_bipartite': 0.0,
    'bipartite_inter_fdr': 'fdr25',
    'bipartite_logfc_fdr': 'fdr05',
    'first_celltype': 'Tumor cells',
    'second_celltype': 'Macrophages/mMDSC',
    'third_celltype':  '(None)',
    'bidirectional_bipartite': False,
}
INTERACTIONS_SAVE_LOCATION = 'compiled_figs/default_interactions.pkl'

DEFAULT_CIRCOS_ARGS = {
    'inter_circos_fdr': 'fdr25',  # Default from fig 4
    'logfc_circos_fdr': 'fdr05',  # Default from fig 4
    'circos_set': 'sting',  # Default from fig 4
    'circos_min_numsigi1': 10,  # Default from fig 4
    'circos_min_numdeg': 0,  # Default from fig 4
    'circos_min_ligand_logfc': 0.12,  # Default from fig 4
    'genes': None
}
CIRCOS_SAVE_LOCATION = 'compiled_figs/default_circos.pkl'


def bipartite_graph(df,
                    cin_only,
                    cell1,
                    cell2,
                    cell3=None,
                    numInteractions=50,
                    min_logfc_bipartite=0,
                    min_expression_bipartite=0,
                    logfc_fdr_bipartite_cutoff=0.05,
                    bidirectional=True,
                    allow_overlap=True):
    scaleratio = .8
    height = 1275

    if cell3 is not None and 'none' in str(cell3).lower():
        cell3 = None

    if cell3 is None:
        selected = df[((df.cell_type_ligand == cell1) & (df.cell_type_receptor == cell2)) |
                      (bidirectional and (df.cell_type_ligand == cell2) & (df.cell_type_receptor == cell1))]
    else:
        selected = df[((df.cell_type_ligand == cell1) & (df.cell_type_receptor == cell2)) |
                      (bidirectional and (df.cell_type_ligand == cell2) & (df.cell_type_receptor == cell1)) |
                      ((df.cell_type_ligand == cell2) & (df.cell_type_receptor == cell3)) |
                      (bidirectional and (df.cell_type_ligand == cell3) & (df.cell_type_receptor == cell2)) |
                      (bidirectional and allow_overlap & (df.cell_type_ligand == cell1) & (df.cell_type_receptor == cell3)) |
                      (bidirectional and allow_overlap & (df.cell_type_ligand == cell3) & (df.cell_type_receptor == cell1))]

    selected = selected[(selected.numSigI1 >= numInteractions) &
                        (selected.expression_ligand >= min_expression_bipartite) &
                        (selected.expression_receptor >= min_expression_bipartite) &
                        (selected["MAST_fdr_ligand"] < logfc_fdr_bipartite_cutoff) &
                        (selected["MAST_fdr_receptor"] < logfc_fdr_bipartite_cutoff) &
                        (np.abs(selected["MAST_log2FC_ligand"]) > min_logfc_bipartite)
                        ].sort_values(
        by="numDEG", ascending=False)

    G = nx.DiGraph()
    for i, row in selected.iterrows():
        lig_ct = row['cell_type_ligand']
        rec_ct = row['cell_type_receptor']
        lig = row['ligand']
        rec = row['receptor']

        lig_node = lig_ct + "_" + lig
        rec_node = rec_ct + "_" + rec

        if not G.has_node(lig_node):
            column = 0
            if cell3 is None:
                if cell1 == cell2:
                    column = 0  # Ligand from left to right
                elif cell1 == lig_ct:
                    column = 0
                elif cell2 == lig_ct:
                    column = 1
            else:
                if cell1 != cell2 and cell2 != cell3 and cell1 != cell3:
                    # All different
                    if cell1 == lig_ct:
                        column = 0
                    elif cell2 == lig_ct:
                        column = 1
                    elif cell3 == lig_ct:
                        column = 2
                else:
                    if cell1 == lig_ct and cell1 == cell2 and cell2 == cell3:
                        column = 0  # Leftmost
                    elif cell1 == lig_ct and (cell1 == cell2 or cell1 == cell3):
                        column = 0
                    elif cell1 == lig_ct:
                        column = 0
                    elif cell2 == lig_ct:
                        column = 1
                    elif cell3 == lig_ct:
                        column = 2

            G.add_node(lig_node, ligand=lig in df.ligand.values, receptor=lig in df.receptor.values, name=lig,
                       column=column,
                       celltype=lig_ct,
                       logfc_fdr_bipartite=row["MAST_fdr_ligand"],
                       logfc=row["MAST_log2FC_receptor"],
                       expression=row['expression_ligand'])
        if not G.has_node(rec_node):
            column = 0
            if cell3 is None:
                if cell1 == cell2:
                    column = 1  # Receptor from right to left
                elif cell1 == rec_ct:
                    column = 0
                elif cell2 == rec_ct:
                    column = 1
            else:
                if cell1 != cell2 and cell2 != cell3 and cell1 != cell3:
                    # All different
                    if cell1 == rec_ct:
                        column = 0
                    elif cell2 == rec_ct:
                        column = 1
                    elif cell3 == rec_ct:
                        column = 2
                else:
                    if cell1 == rec_ct and cell1 == cell2 and cell2 == cell3:
                        column = 1  # Middle since leftmost is ligand
                    elif cell1 == rec_ct and cell1 == cell2:
                        column = 1  # Middle since leftmost is ligand
                    elif cell1 == rec_ct and cell1 == cell3:
                        column = 2
                    elif cell2 == rec_ct and cell2 == cell3:
                        column = 2
                    elif cell1 == rec_ct:
                        column = 0
                    elif cell2 == rec_ct:
                        column = 1
                    elif cell3 == rec_ct:
                        column = 2

            G.add_node(rec_node, ligand=rec in df.ligand.values, receptor=rec in df.receptor.values, name=rec,
                       column=column,
                       celltype=rec_ct,
                       logfc_fdr_bipartite=row["MAST_fdr_receptor"],
                       logfc=row["MAST_log2FC_receptor"],
                       expression=row['expression_receptor'])

        G.add_edge(lig_node, rec_node,
                   numInteractions=row['numSigI1'],
                   ligand_fdr=row["MAST_fdr_ligand"],
                   ligand_logfc=row["MAST_log2FC_ligand"],
                   receptor_fdr=row["MAST_fdr_receptor"],
                   receptor_logfc=row["MAST_log2FC_receptor"],
                   numDEG=row['numDEG'])

    ct_ordering = [0, 1]
    if cell3:
        ct_ordering.append(2)
    pos = multipartite_layout(G, subset_key='column', scale=1, space_mult_x=30 * scaleratio, space_mult_y=60,
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

    color_step = 255
    colorscale = ColorTransformer(0, color_step, cmap='bwr')

    # colorscale = [x.hex for x in list(Color(start_color).range_to(Color(end_color), color_step))]

    min_color_threshold = -0.2  # min_logfc_bipartite
    max_color_threshold = 0.2  # max_logfc

    def get_color(value, minimum, maximum):
        step_size = (max_color_threshold - min_color_threshold) / (color_step - 1)
        index = (value - min_color_threshold) // step_size
        return colorscale(min(max(int(index), 0), color_step - 1))

    # Can't hover over lines directly, so add invisible points
    # See https://stackoverflow.com/a/46039229/5179044
    interpolation_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        showlegend=False,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            opacity=0,
            cmin=min_color_threshold,
            cmax=max_color_threshold,
            color=[],
            colorbar=dict(
                title="Ligand LogFC",
                xanchor='right',
                yanchor='top',
                thickness=25,
                titleside='right',
                x=1.02,
                y=1,
                lenmode="fraction",
                len=0.25,
                tickmode='array',
                tickvals=[min_color_threshold, max_color_threshold],
                ticktext=["CIN<sup>low</sup>", "CIN<sup>high</sup>" + ("" if cin_only else "/<br>STING<sup>KD</sup>")]
            ),
            colorscale=colorscale.make_plotly_colorscale()
        )
    )

    edge_traces = []

    EDGE_PADDING = .08
    MIN_WIDTH = 2
    MAX_WIDTH = 6
    INTERPOLATION_POINTS = 20
    for edge in G.edges():
        if edge[0] not in pos or edge[1] not in pos:
            continue
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
            showlegend=False,
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
    celltypes = set()
    for node in G.nodes():
        if node not in pos:
            continue
        column = G.nodes[node]['column']
        celltype = G.nodes[node]['celltype']
        celltypes.add(celltype)
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

        if column not in layer_anchors:
            layer_anchors[column] = (x, y, 1, celltype)
        else:
            curr_anchor = layer_anchors[column]
            if y > curr_anchor[1]:
                layer_anchors[column] = (x, y, 1 + curr_anchor[2], celltype)

    annotations = []
    celltype_background_points = []
    yrange = max(node_y) - min(node_y) if len(node_y) > 0 else 0
    if len(layer_anchors) == 0:
        xrange = 0
    else:
        xrange = (max(node_x) - min(node_x) if len(node_x) > 0 else 0) / len(layer_anchors)
    celltype2colors = celltype_to_colors(celltypes)
    for column, (anchor_x, anchor_y, layer_size, celltype) in layer_anchors.items():
        annotations.append(dict(
            x=anchor_x,
            y=max(node_y) + (yrange / 12) if len(node_y) > 0 else 0,
            xref="x",
            yref="y",
            text=celltype,
            showarrow=False,
            font=dict(
                size=16,
                color="white"
            ),
            align="center",
            bordercolor="black",
            borderwidth=2,
            borderpad=4,
            bgcolor=celltype2colors.get(celltype, "grey")
        ))

        celltype_background_points += [
            (anchor_x - .45 * xrange, min(node_y) - .05 * yrange),
            (anchor_x - .45 * xrange, max(node_y) + .05 * yrange),
            (anchor_x + .45 * xrange, max(node_y) + .05 * yrange),
            (anchor_x + .45 * xrange, min(node_y) - .05 * yrange),
            (anchor_x - .45 * xrange, min(node_y) - .05 * yrange),
            (None, None)
        ]

    celltype_bg_trace = go.Scatter(
        x=[c[0] for c in celltype_background_points],
        y=[c[1] for c in celltype_background_points],
        showlegend=False,
        fill="toself",
        fillcolor='#fffff2',
        line=dict(
            width=2,
            color='grey'
        ),
        hoverinfo='none',
        mode='lines'
    )

    first_layer = 0
    last_layer = 1 if len(layer_anchors) == 2 else 2

    for node in G.nodes():
        if node not in pos:
            continue
        column = G.nodes[node]['column']
        layer_size = layer_anchors[column][2]
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
        if column == last_layer:
            align = "left"
        elif column == first_layer:
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
        showlegend=False,
        marker=dict(
            showscale=True,
            colorscale='gray',
            reversescale=True,
            color=[],
            #color='grey',
            size=15,
            cmin=0.0,
            cmax=1.0,
            line=dict(width=1, color='black'),
            colorbar=dict(
                thickness=25,
                title="Expression",
                yanchor='top',
                xanchor='right',
                x=.92,
                y=1,
                titleside='right',
                lenmode="fraction",
                len=0.25,
                tickmode='array',
                tickvals=[0.0, 1.0],
                ticktext=["0%", "100%"]
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
        #node_expressions.append(expression)
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
        node_expressions.append(expression)

    node_trace.marker.color = node_expressions
    node_trace.text = node_text

    if len(node_text) == 0:
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

    legend_traces = [
        go.Scatter(
            name="Ligand",
            x=[None], y=[None],
            marker_symbol="circle",
            mode='markers',
            showlegend=True,
            legendgroup="genetypes",
            legendgrouptitle=dict(
                text="Gene"
            ),
            marker=dict(
                color='white',
                size=15,
                line=dict(width=1, color='black')
            )
        ),
        go.Scatter(
            name="Receptor",
            x=[None], y=[None],
            marker_symbol="square",
            mode='markers',
            showlegend=True,
            legendgroup="genetypes",
            legendgrouptitle=dict(
                text="Gene"
            ),
            marker=dict(
                color='white',
                size=15,
                line=dict(width=1, color='black')
            )
        ),
        go.Scatter(
            name="Ligand/Receptor",
            x=[None], y=[None],
            marker_symbol="diamond-dot",
            mode='markers',
            showlegend=True,
            legendgroup="genetypes",
            legendgrouptitle=dict(
                text="Gene"
            ),
            marker=dict(
                color='white',
                size=15,
                line=dict(width=1, color='black')
            )
        ),
    ]

    fig = go.Figure(data=[celltype_bg_trace] + edge_traces + [node_trace, interpolation_trace] + legend_traces,
                    layout=go.Layout(
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=40, l=5, r=5, t=40),
                        autosize=False,
                        width=height * scaleratio,
                        height=height,
                        plot_bgcolor='rgba(0,0,0,0)',
                        annotations=annotations,
                        legend=dict(
                            yanchor="top",
                            xanchor="right",
                            x=.985,
                            y=.75
                        ),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, scaleratio=scaleratio,
                                   scaleanchor="x"))
                    )
#    return enhance_plotly_export(fig, height, scaleratio)
    return fig


def get_gene_type(df: pd.DataFrame, gene):
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


def make_plot_from_graph(G: nx.DiGraph, celltypes, layout="planar", colormap=None) -> go.Figure:
    max_t = max([G.nodes[n]['t'] for n in G.nodes])
    height = min((max((800, len([n for n in G.nodes if G.nodes[n]['t'] == max_t]) * 40)), 1600))  # Adjust height if too many nodes
    width = 1000
    scaleratio = width / height

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

    color_step = 255
    colorscale = ColorTransformer(0, color_step, cmap='bwr' if not colormap else colormap)

    # colorscale = [x.hex for x in list(Color(start_color).range_to(Color(end_color), color_step))]

    min_logfc = None
    max_logfc = None
    max_abs_logfc = None
    min_abs_logfc = None
    for edge in G.edges():
        logfc = G[edge[0]][edge[1]]['response_logfc']
        if logfc is None:
            logfc = 0
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

    min_color_threshold = -0.5  # min_logfc
    max_color_threshold = 0.5  # max_logfc

    def get_color(value, minimum, maximum):
        if value is None:
            value = 0
        step_size = (max_color_threshold - min_color_threshold) / (color_step - 1)
        index = (value - min_color_threshold) // step_size
        return colorscale(min(max(int(index), 0), color_step - 1))

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
                title="Induced LogFC",
                xanchor='left',
                yanchor='top',
                thickness=25,
                titleside='right',
                x=1.425,
                y=1,
                lenmode='fraction',
                len=0.21
            ),
            colorscale=colorscale.make_plotly_colorscale()
        )
    )
    edge_traces = []
    EDGE_PADDING = .1
    MIN_WIDTH = 3
    MAX_WIDTH = 6
    INTERPOLATION_POINTS = 20

    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        logfc = G[edge[0]][edge[1]]['response_logfc'] or 0
        ligand = G.nodes[G[edge[0]][edge[1]]['ligand']]['gene']
        donor_cell = G.nodes[G[edge[0]][edge[1]]['ligand']]['celltype']
        receptor = G.nodes[G[edge[0]][edge[1]]['receptor']]['gene']
        target_cell = G.nodes[G[edge[0]][edge[1]]['receptor']]['celltype']

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
        hover_text = f"Ligand: {ligand} ({donor_cell})<br>Receptor: {receptor} ({target_cell})<br>"
        if logfc is not None and logfc != 0:
            hover_text += f"Induced LogFC Effect: {logfc:.2f}<br>"
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

    # Get the size range from node weights
    weights = [abs(x or 0) for x in nx.get_node_attributes(G, 'weight').values()]
    max_node_logfc = max(weights)
    min_node_logfc = min(weights)

    node_colors = []
    genetype2nodes = dict()
    for node in G.nodes():
        step = G.nodes[node]['t']
        celltype = G.nodes[node]['celltype']
        gene = G.nodes[node]['gene']
        gene_type = G.nodes[node]['gene_type']

        response_logfc = G.nodes[node].get('weight', 0)
        expression = G.nodes[node].get('expression', 0) or 0
        #ligand_logfc = max(get_all_edge_attr_for_node(G, node, "ligand_logfc") + [0], key=abs)
        #receptor_logfc = max(get_all_edge_attr_for_node(G, node, "receptor_logfc") + [0], key=abs)

        if step == 0:  # Initial gene
            gene_type += '-initial'

        if gene_type not in genetype2nodes:
            genetype2nodes[gene_type] = ([],
                                         [],
                                         [],
                                         [],
                                         [],
                                         #[]
                                         )
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
        genetype2nodes[gene_type][2].append(f'Name: {gene}<br>Cell Type: {celltype}<br>Step: {step}' + (f"<br>Expression: {expression:.2f}" if expression != 0 else ''))
        genetype2nodes[gene_type][3].append(celltype)
        genetype2nodes[gene_type][4].append(response_logfc if (response_logfc is not None and response_logfc != 0) else (min_node_logfc + (max_node_logfc - min_node_logfc) * 0.5))
        #genetype2nodes[gene_type][5].append(receptor_logfc)

    MIN_NODE_SIZE = 14
    MAX_NODE_SIZE = 24
    NODE_SIZE_RANGE = MAX_NODE_SIZE - MIN_NODE_SIZE

    nodes = []
    for gene_type in genetype2nodes.keys():
        gene_xs, gene_ys, texts, celltypes, response_logfc = genetype2nodes[gene_type]

        is_initial = gene_type.endswith('-initial')
        if is_initial:
            gene_type = gene_type.replace('-initial', '')
            symbol = genetype2symbol[gene_type] + '-dot'
        else:
            symbol = genetype2symbol[gene_type]
        # if symbol == 'arrow-up':
        # Need to adjust y-values to center the triangle
        # gene_ys = [y*.95 for y in gene_ys]

        node_sizes = [(MIN_NODE_SIZE + (NODE_SIZE_RANGE*(abs(logfc or (max_node_logfc/2))/max_node_logfc))) * (1.25 if 'triangle-up' in symbol else 1) for logfc in response_logfc]  # Triangles visually appear smaller than expected

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
        name="High",
        x=[None], y=[None],
        marker_symbol='circle',
        mode='markers',
        showlegend=True,
        # visible = 'legendonly',
        legendgroup='nodesize',
        legendgrouptitle=dict(
            text='Induced abs(LogFC)'
        ),
        marker=dict(
            color='grey',
            size=MAX_NODE_SIZE,
            line=dict(width=1, color='black')
        )
    ))
    legends.append(go.Scatter(
        name="Low",
        x=[None], y=[None],
        marker_symbol='circle',
        mode='markers',
        showlegend=True,
        # visible = 'legendonly',
        legendgroup='nodesize',
        legendgrouptitle=dict(
            text='Induced abs(LogFC)'
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
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=10, l=1, r=15, t=30),
                        autosize=False,
                        width=height * scaleratio,
                        height=height,
                        plot_bgcolor='rgba(0,0,0,0)',
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


def pseudotime_interaction_propagation_graph(effect_set: str,
                                             seed_cell: str,
                                             seed_ligands: list,
                                             iterations: int = 4,
                                             interaction_fdr_cutoff=0.05,
                                             min_logfc=0.01,
                                             min_expression=0.1,
                                             logfc_fdr_cutoff=0.05,
                                             set_progress_callback=None,
                                             return_only_raw_data=False):

    iterations = iterations + 1  # +1 because the first iteration are the seeds

    # Basic filters to build the network from
    def _read_filtered_df(celltype, target):
        df = read_ligand_effect_for(effect_set, target)
        if celltype is not None:
            df = df[df['cell_type'] == celltype]
        if df is None:
            return pd.DataFrame()

        # Ligand target filter
        ligand_filter = (
            (df['MAST_log2FC'].abs() >= min_logfc) &  # Require the ligand in question to be differentially expressed
            (df['MAST_fdr'] < logfc_fdr_cutoff) &
            (df['target_expression'] > min_expression) &
            (df['gene_expression'] > min_expression) &
            (df['gene_is_receptor'])
        ) & df['ligand']

        # Receptor target filter
        receptor_filter = (
            ((
                (df['log2FC'].abs() >= min_logfc) &  # Require the downstream gene to be differentially expressed
                (df['fdr'] < logfc_fdr_cutoff) &
                (df['gene_is_ligand'])
            ) | (
                (df['gene_max_ligand_expression'] > min_expression) &  # If a receptor, require the ligand to be available somewhere
                (df['log2FC'].abs() >= min_logfc) &  # Require the downstream gene to be differentially expressed
                (df['fdr'] < logfc_fdr_cutoff) &
                df['gene_is_receptor']
            )) &  # Allow for receptor-receptor interaction effects
            (df['i1.fdr'] < interaction_fdr_cutoff) &
            (df['gene_expression'] > min_expression) &  # Require the downstream gene to be expressed
            (df['target_expression'] > min_expression)
        ) & df['receptor']

        return df[ligand_filter | receptor_filter]

    if isinstance(seed_ligands, str):
        ligands = []
        for ligand in seed_ligands.split(","):
            ligands.append(ligand.strip())
        seed_ligands = ligands

    ligand_receptor_pairs = read_ligand_receptor_file()

    Gs = []
    curr_G = nx.DiGraph()
    frontier = set()
    celltypes = {seed_cell}
    for t in range(iterations):
        if set_progress_callback is not None:
            set_progress_callback((t, iterations))

        # Create dynamic dataframe for the frontier
        df = pd.DataFrame()
        added = set()  # Track genes that we added full sets for
        for node in frontier:
            celltype = curr_G.nodes[node]['celltype']
            gene = curr_G.nodes[node]['gene']
            gene_type = curr_G.nodes[node]['gene_type']
            # Get current node's info
            df = pd.concat([df, _read_filtered_df(celltype, gene)])
            if 'Ligand' in gene_type:
                # If its a ligand, we need all possible receptor info
                receptors = ligand_receptor_pairs[ligand_receptor_pairs['ligand'] == gene]['receptor'].unique()
                for receptor in receptors:
                    if receptor not in added:
                        df = pd.concat([df, _read_filtered_df(None, receptor)])  # Get all possible receptors
                        added.add(receptor)
            if 'Receptor' in gene_type:
                ligands = ligand_receptor_pairs[ligand_receptor_pairs['receptor'] == gene]['ligand'].unique()
                for ligand in ligands:
                    if ligand not in added:
                        df = pd.concat([df, _read_filtered_df(None, ligand)])
                        added.add(ligand)
        df = df.drop_duplicates(ignore_index=True)
        del added

        if t == 0:
            for seed_ligand in seed_ligands:
                gene_type = get_gene_type(ligand_receptor_pairs, seed_ligand)
                seed_node = f"{seed_cell}_{seed_ligand}"
                curr_G.add_node(seed_node, t=t, gene=seed_ligand, celltype=seed_cell, gene_type=gene_type,
                                expression=None,
                                #diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand)
                                )
                frontier.add(seed_node)
        else:
            # Break early if no new nodes or too many nodes
            if df.shape[0] == 0 or len(frontier) > 45:
                break
            new_frontier = set()
            for node in frontier:
                celltype = curr_G.nodes[node]['celltype']
                gene_type = curr_G.nodes[node]['gene_type']

                if 'Ligand' in gene_type:
                    ligand = curr_G.nodes[node]['gene']
                    receivers = ligand_receptor_pairs[ligand_receptor_pairs['ligand'] == ligand]['receptor'].unique()
                    receiver_inters = df[df['target'].isin(receivers)].drop_duplicates(('cell_type', 'target'), ignore_index=True)
                    for i, row in receiver_inters.iterrows():
                        receptor = row['target']
                        receptor_celltype = row['cell_type']
                        expression = row['target_expression']
                        celltypes.add(receptor_celltype)
                        next_node = f"{receptor_celltype}_{receptor}"

                        #if get_interaction_fdr(ct, celltype, ligand, receptor) > interaction_fdr_cutoff:
                        #    continue

                        if not curr_G.has_node(next_node):
                            curr_G.add_node(next_node, t=t,
                                            gene=receptor,
                                            celltype=receptor_celltype,
                                            expression=expression,
                                            gene_type=get_gene_type(ligand_receptor_pairs, receptor),
                                            #diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand)
                                            )

                        # Only add if the edge doesn't exist and that the node is only one time point behind
                        if not curr_G.has_edge(node, next_node) and curr_G.nodes[node]['t']+1 == curr_G.nodes[next_node]['t'] and curr_G.nodes[next_node]['t'] == t:
                            curr_G.add_edge(node, next_node,
                                            t=t, ligand=node,
                                            receptor=next_node,
                                            response_logfc=None,
                                            weight=None  # row['MAST_log2FC']
                                            )
                            new_frontier.add(next_node)

                if 'Receptor' in gene_type:
                    receptor = curr_G.nodes[node]['gene']
                    # Get the downstream ligands
                    down_ligands = df[(df['target'] == receptor) &
                                      (df['cell_type'] == celltype) &
                                      (df['i1.fdr'] < interaction_fdr_cutoff) &
                                      (df['log2FC'].abs() > min_logfc) &
                                      (df['fdr'] < logfc_fdr_cutoff) &
                                      (df['gene_is_ligand'])].drop_duplicates(ignore_index=True)

                    for i, row in down_ligands.iterrows():
                        ligand = row['gene']
                        expression = row['gene_expression']
                        ligand_celltype = celltype
                        celltypes.add(ligand_celltype)
                        next_node = f"{ligand_celltype}_{ligand}"

                        if not curr_G.has_node(next_node):
                            curr_G.add_node(next_node, t=t,
                                            gene=ligand,
                                            celltype=ligand_celltype,
                                            expression=expression,
                                            gene_type=get_gene_type(ligand_receptor_pairs, ligand),
                                            #diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand)
                                            )
                        # Only add if the edge doesn't exist and that the node is only one time point behind
                        if not curr_G.has_edge(node, next_node) and curr_G.nodes[node]['t']+1 == curr_G.nodes[next_node]['t'] and curr_G.nodes[next_node]['t'] == t:
                            curr_G.add_edge(node, next_node,
                                            t=t, ligand=node,
                                            receptor=next_node,
                                            response_logfc=row['log2FC'],
                                            weight=abs(row['log2FC'])
                                            )
                            new_frontier.add(next_node)

                    down_receptors = df[(df['target'] == receptor) &
                                        (df['cell_type'] == celltype) &
                                        (df['gene_max_ligand_expression'] > min_expression) &
                                        (df['i1.fdr'] < interaction_fdr_cutoff) &
                                        (df['gene_is_receptor'])].drop_duplicates(ignore_index=True)
                    # Include downstream receptor if ligand is available
                    for i, row in down_receptors.iterrows():
                        receptor = row['gene']
                        expression = row['gene_expression']
                        receptor_celltype = celltype
                        celltypes.add(receptor_celltype)
                        next_node = f"{receptor_celltype}_{receptor}"

                        if not curr_G.has_node(next_node):
                            curr_G.add_node(next_node, t=t,
                                            gene=receptor,
                                            celltype=receptor_celltype,
                                            expression=expression,
                                            gene_type=get_gene_type(ligand_receptor_pairs, receptor),
                                            #diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand)
                                            )
                        # Only add if the edge doesn't exist and that the node is only one time point behind
                        if not curr_G.has_edge(node, next_node) and curr_G.nodes[node]['t']+1 == curr_G.nodes[next_node]['t'] and curr_G.nodes[next_node]['t'] == t:
                            curr_G.add_edge(node, next_node,
                                            t=t, ligand=node,
                                            receptor=next_node,
                                            response_logfc=row['log2FC'],
                                            weight=abs(row['log2FC'])
                                            )
                            new_frontier.add(next_node)

            frontier = new_frontier

        # Fill in node weights according to total induced logfc
        for node in curr_G.nodes:
            if 'Ligand' not in curr_G.nodes[node]['gene_type']:
                curr_G.nodes[node]['weight'] = 0
                continue
            for pred in curr_G.predecessors(node):
                curr_G.nodes[node]['weight'] = curr_G.nodes[node].get('weight', 0) + (curr_G.edges[pred, node]['weight'] or 0)

        Gs.append(curr_G.copy())
        if len(frontier) == 0:
            break

    if len(Gs) == 0:
        return None

    last_graph = Gs[-1]

    if last_graph.number_of_nodes() == 0 or last_graph.number_of_edges() == 0:
        return None

    if return_only_raw_data:
        return last_graph

    return make_plot_from_graph(last_graph, list(celltypes), layout='timeline', colormap='cividis')


def polar2cartesian(r, theta):
    theta = np.deg2rad(theta)
    return r * np.sin(theta), r * np.cos(theta)


@total_ordering
@dataclass(unsafe_hash=True)
class PlotRing:
    name: str
    color: str
    is_receptor: bool
    height: float
    circumference: float
    children: List['PlotRing'] = field(default_factory=list)
    theta: float = 0  # computed later
    start_height: float = 0  # computed later

    def __lt__(self, other):
        return self.height < other.height or (self.height == other.height and self.circumference < other.circumference)

    def add_child(self, child):
        # Add child in sorted order
        bisect.insort(self.children, child)

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


# def plot_polar_data(rings: PlotRing, legend_info: Dict[str, str], with_text: bool = True) -> go.Figure:
#     # Normalize the data
#     # First ring is always the full circle
#     rings.circumference = 360
#     rings.theta = 180
#     # Reset center node height
#     number_of_levels = rings.get_number_of_levels()
#     rings.height = 1
#     for i, ring in enumerate(rings.iter_rings()):
#         max_ring_height = max([r.height for r in ring], default=1)
#         for r in ring:
#             r.height = (r.height / max_ring_height) * 1/number_of_levels
#     rings.height *= .5
#     # Calculate circumference of each ring group
#     for parent, children in rings.iter_hierarchy():
#         total_child_circumference = sum([c.circumference for c in children])
#         total_group_circumference = parent.circumference
#         curr_theta = parent.theta - total_group_circumference/2
#        #rings_max_height = max([c.height for c in children], default=1)
#         for child in children:
#             child.circumference = total_group_circumference * (child.circumference / total_child_circumference)
#             child.theta = curr_theta + child.circumference/2
#             child.start_height = parent.start_height + parent.height
#             #child.height += child.start_height
#             curr_theta += child.circumference
#
#     max_radius = rings.get_max_height()
#
#     # Build the figure
#     annotations = []
#     traces = []
#     for name, color in legend_info.items():
#         traces.append(go.Scatter(
#             name=name,
#             x=[None],
#             y=[None],
#             marker_symbol='circle',
#             mode='markers',
#             showlegend=True,
#             legendgroup='colors',
#             legendgrouptitle=dict(
#                 text='Ring Color',
#             ),
#             marker=dict(
#                 color=color,
#                 size=15,
#                 line=dict(width=1, color='black')
#             )
#         ))
#
#     font_size = 75 // number_of_levels
#     for level, ring_list in enumerate(rings.iter_rings()):
#         if level == 0:  # Handle first circle
#             assert len(ring_list) == 1
#             node = ring_list[0]
#             if with_text:
#                 annotations.append(dict(
#                     x=0, y=0,
#                     xref='x', yref='y',
#                     clicktoshow='onoff',
#                     text=node.name,
#                     visible=True,
#                     showarrow=False,
#                     font=dict(
#                         size=font_size,
#                         color="black"
#                     ),
#                     align='center'
#                 ))
#             traces.append(go.Barpolar(
#                 r=[node.height],
#                 theta=[0],
#                 width=[360],
#                 thetaunit='degrees',
#                 opacity=1,
#                 base=0,
#                 text=node.name,
#                 showlegend=False,
#                 marker=dict(
#                     color=[node.color],
#                     line=dict(
#                         color='black',
#                         width=2
#                     )
#                 )
#             ))
#         else:
#             thetas = []
#             circumferences = []
#             bases = []
#             colors = []
#             heights = []
#             names = []
#             for node in ring_list:
#                 base = node.start_height
#                 theta = node.theta
#                 height = node.height
#                 circumference = node.circumference
#                 color = node.color
#                 name = node.name
#
#                 # Label
#                 if with_text:
#                     annotation_x, annotation_y = polar2cartesian((3/4)*(height) + base, theta+90)
#                     annotations.append(dict(
#                         x=annotation_x, y=annotation_y,
#                         clicktoshow='onoff',
#                         text=name,
#                         visible=True,
#                         showarrow=False,
#                         font=dict(
#                             size=font_size,
#                             color="black"
#                         ),
#                         align='center'
#                     ))
#
#                 # Ring info
#                 thetas.append(theta)
#                 circumferences.append(circumference)
#                 bases.append(base)
#                 colors.append(color)
#                 heights.append(height)
#                 names.append(name)
#
#             traces.append(go.Barpolar(
#                 r0=0,
#                 r=heights,
#                 theta0=0,
#                 theta=thetas,
#                 width=circumferences,
#                 thetaunit='degrees',
#                 opacity=1,
#                 base=bases,
#                 text=names,
#                 showlegend=False,
#                 marker=dict(
#                     color=colors,
#                     line=dict(
#                         color='black',
#                         width=2
#                     )
#                 )
#             ))
#
#     fig = go.Figure(data=traces, layout=go.Layout(
#         template=None,
#         showlegend=True,
#         hovermode='closest',
#         margin=dict(l=0, r=0, t=0, b=0),
#         autosize=False,
#         barmode='overlay',
#         plot_bgcolor='rgba(0,0,0,0)',
#         paper_bgcolor='rgba(0,0,0,0)',
#         width=1000,
#         height=1000,
#         annotations=annotations,
#         legend=dict(
#             yanchor="top",
#             y=1,
#             xanchor="left",
#             x=1.15
#         ),
#         polar=dict(
#             bgcolor='rgba(0,0,0,0)',
#             radialaxis=dict(
#                 showticklabels=False,
#                 ticks='',
#                 range=[0, max_radius + 0.1],
#                 showgrid=False,
#                 showline=False,
#                 angle=0,
#                 tick0=0,
#                 tickangle=90
#             ),
#             angularaxis=dict(
#                 showticklabels=False,
#                 ticks='',
#                 showgrid=False,
#                 rotation=0,
#                 tick0=0,
#                 tickangle=90,
#                 direction='clockwise'
#         )),
#         xaxis=dict(
#             showgrid=False,
#             zeroline=False,
#             showticklabels=False,
#             range=[-1*max_radius - 0.1, max_radius + 0.1]
#         ),
#         yaxis=dict(
#             showgrid=False,
#             zeroline=False,
#             showticklabels=False,
#             range=[-1*max_radius - 0.1, max_radius + 0.1],
#             scaleanchor="x"
#         )
#     ))
#     return fig


# def radial_ligand_effect_figure(
#         ct: ad.AnnData, interactions: pd.DataFrame,
#         ligand: str, celltype: str,
#         min_logfc: float = 0, max_fdr: float = 0.05,
#         min_numSigI1: int = 1,
#         max_levels: int = 3, with_text: bool = True, initial_celltype_filter: Set[str] = None
# ) -> go.Figure:
#     all_celltypes = set(ct.obs['cell type'].unique())
#     if initial_celltype_filter is None:
#         initial_celltype_filter = all_celltypes
#     celltype2color = celltype_to_colors(all_celltypes)
#
#     #interactions = interactions[(interactions.MAST_log2FC_ligand >= min_logfc) &
#     #                            (interactions.MAST_fdr_ligand <= max_fdr)]
#
#     rings = PlotRing(ligand, celltype2color[celltype], False, 0, 360)
#     legend_data = dict()
#
#     def build_downstream_rings(i: int,
#                                curr_ring: PlotRing,
#                                start: str,
#                                start_celltype: str,
#                                is_ligand: bool):
#         if i >= max_levels:
#             return
#
#         legend_data[start_celltype] = celltype2color[start_celltype]
#
#         filtered_interactions = interactions
#         if i < 2:  # Initial filtering
#             filtered_interactions = filtered_interactions[
#                 (filtered_interactions.cell_type_receptor.isin(initial_celltype_filter)) &
#                 (filtered_interactions.cell_type_ligand.isin(initial_celltype_filter))]
#
#         if is_ligand:
#             # Get and map receptors
#             filtered_interactions = filtered_interactions[(filtered_interactions.ligand == start)].drop_duplicates(subset=['receptor', 'cell_type_receptor'])
#             for _, row in filtered_interactions.iterrows():
#                 receptor = row['receptor']
#                 receptor_celltype = row['cell_type_receptor']
#                 numDEG = row['numDEG_fdr25_receptor']
#                 numSigI1 = row['numSigI1_fdr25_receptor']
#                 if numSigI1 < min_numSigI1:
#                     continue
#                 next_ring = PlotRing(
#                     f"{receptor}+", celltype2color[receptor_celltype], True,
#                     1, 1
#                 )
#                 curr_ring.add_child(next_ring)
#                 build_downstream_rings(i+1, next_ring, receptor, receptor_celltype, False)
#         else:
#             # Get and map ligands
#             downstream_ligands = get_downstream_ligands(ct, start, start_celltype, min_logfc, max_fdr)
#             for lig in downstream_ligands:
#                 interaction_logfc = get_interaction_logfc(ct, start_celltype, start, lig)
#                 interaction_fdr = get_interaction_fdr(ct, start_celltype, start, lig)
#                 interaction_logfc_fdr = get_interaction_logfc_fdr(ct, start_celltype, start, lig)
#                 if np.isnan(interaction_fdr) or interaction_fdr > max_fdr or interaction_logfc < min_logfc:
#                     continue
#                 next_ring = PlotRing(
#                     lig, 'grey', False, # celltype2color[start_celltype],
#                     interaction_logfc, max(-1 * np.log(interaction_fdr + 1e-10), 1e-5)
#                 )
#                 curr_ring.add_child(next_ring)
#                 build_downstream_rings(i+1, next_ring, lig, start_celltype, True)
#
#     build_downstream_rings(0, rings, ligand, celltype, True)
#
#     # Scale receptor rings by number of downstream effects
#     for level in rings.iter_rings():
#         if len(level) < 2:
#             continue
#         if not level[0].is_receptor:
#             continue
#         number_of_downstream = [len(r.get_all_children())+1 for r in level]
#         total_downstream = sum(number_of_downstream)
#         for r, levels in zip(level, number_of_downstream):
#             r.circumference = levels / total_downstream
#
#     legend_data['Ligand'] = 'grey'
#
#     return plot_polar_data(rings, legend_data, with_text)
