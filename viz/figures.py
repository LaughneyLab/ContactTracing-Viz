import math

import anndata as ad
import networkx as nx
import numpy as np
from matplotlib import cm, colors
import pandas as pd
import plotly.graph_objects as go

from viz.data import get_interaction_fdr, get_effective_logfc, get_diff_abundance, get_downstream_ligands
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

                        if get_interaction_fdr(ct, celltype, ligand, receptor) < interaction_fdr_cutoff:
                            continue

                        if not curr_G.has_node(next_node):
                            curr_G.add_node(next_node, t=t,
                                            gene=receptor,
                                            celltype=receptor_celltype,
                                            gene_type=get_gene_type(orig_df, receptor),
                                            #diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand)
                                            )
                        if not curr_G.has_edge(node, next_node) and not curr_G.has_edge(next_node, node):
                            logfc = get_effective_logfc(row['MAST_log2FC_ligand'], row["MAST_log2FC_receptor"])
                            curr_G.add_edge(node, next_node,
                                            t=t, ligand=node,
                                            receptor=next_node,
                                            logfc=logfc,
                                            weight=abs(logfc),
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

                        if get_interaction_fdr(ct, celltype, receptor, ligand) < interaction_fdr_cutoff:
                            continue

                        if not curr_G.has_node(next_node):
                            curr_G.add_node(next_node, t=t,
                                            gene=ligand,
                                            celltype=ligand_celltype,
                                            gene_type=get_gene_type(orig_df, ligand),
                                            #diff_abundance=get_diff_abundance(ct, seed_cell, seed_ligand)
                                            )
                        if not curr_G.has_edge(node, next_node) and not curr_G.has_edge(next_node, node):
                            logfc = get_effective_logfc(row["MAST_log2FC_receptor"], row['MAST_log2FC_ligand'])
                            curr_G.add_edge(node, next_node,
                                            t=t, ligand=node,
                                            receptor=next_node,
                                            logfc=logfc,
                                            weight=abs(logfc),
                                            numDEG=row['numDEG_fdr_receptor'],
                                            numInteractions=row['numSigI1_fdr_receptor'])
                            new_frontier.add(next_node)

            frontier = new_frontier
        Gs.append(curr_G.copy())
        if set_progress_callback:
            set_progress_callback(((t+1)//iterations * 100, 100))

    return make_plot_from_graph(Gs[-1], df, layout=layout)
