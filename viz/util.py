from typing import Set, Sized, Callable, Iterable

import networkx as nx
import numpy as np
from matplotlib import cm, colors
import math
import itertools


class ColorTransformer(Sized, Callable, Iterable):

    def __init__(self, min=0, max=255, cmap='bwr'):
        self.min = min
        self.max = max
        self.cmap = cm.get_cmap(cmap)
        self.norm = colors.Normalize(vmin=min, vmax=max, clip=True)

    def __call__(self, value):
        return colors.rgb2hex(colors.colorConverter.to_rgb(self.cmap(self.norm(value))))

    def __len__(self) -> int:
        color_range = self.max - self.min
        return 255 if color_range <= 1 else int(color_range)

    def __iter__(self):
        for i in range(len(self)):
            yield self(i)


def saturate_color(color, saturation):
    color = colors.to_rgb(color)
    color = [min(c * saturation, 1) for c in color]
    return colors.to_hex(color)


mouse_colors = {
    "Tumor cells": '#a3d106',
    "T cells": '#f47d06',
    "PMN/gMDSC": '#11eafb',
    "B cells": '#f9d319',
    "pDC": '#51119a',
    "cDC": '#c438fb',
    "Fibroblast cells": '#226c04',
    "NK cells": '#ea0d0c',
    "Macrophages/mMDSC": '#0047cf'
}


def celltype_to_colors(celltypes: list):
    default_colors = [colors.to_hex(c) for c in cm.tab20.colors]
    # Color labels from paper

    celltype2color = dict()
    for i, ct in enumerate(celltypes):
        celltype2color[ct] = mouse_colors.get(ct, default_colors[i])
    return celltype2color


def enhance_plotly_export(fig, height, scaleratio):
    old_show = fig.show

    def new_show(*args, **kwargs):
        kwargs['config'] = {
            'toImageButtonOptions': {
                'format': 'png',  # one of png, svg, jpeg, webp
                'filename': 'exported_image',
                'height': height,
                'width': height * scaleratio,
                'scale': 6  # Multiply title/legend/axis/canvas sizes by this factor
            }
        }
        return old_show(*args, **kwargs)

    setattr(fig, 'show', new_show)
    return fig


# Initial placement of nodes in layers
def space_elements(elements, attributes):
    if len(elements) <= 2:
        return elements

    # Visually space out high degree nodes
    sorted_elements = list(sorted(zip(elements, attributes), key=lambda x: x[1], reverse=True))
    sorted_elements = [x[0] for x in sorted_elements]

    # Approximate positions
    positions = []
    finished = False
    for i in range(1, len(sorted_elements)):
        range_obj = range(i) if i % 2 == 1 else range(i - 1, -1, -1)
        for j in range_obj:
            positions.append(((2*j) + 1) / ((2**i)))
            if len(positions) == len(sorted_elements):
                finished = True
                break
        if finished:
            break

    element2position = dict(zip(sorted_elements, positions))
    return list(sorted(elements, key=lambda e: element2position[e]))


# Bias node placement based on the initial placement
def space_elements_from_previous(G, elements, prev_elements):
    elem2score = dict()
    prev2degree = {e: G.degree(e) for e in prev_elements}
    total_degree = sum(prev2degree.values())
    for elem in elements:
        elem2score[elem] = 0
        for prev_elem in enumerate(prev_elements):
            if G.has_edge(prev_elem, elem) or G.has_edge(elem, prev_elem):
                elem2score[elem] += prev2degree[prev_elem] / total_degree

    initial_sorted_elements = list(sorted(elements, key=lambda x: elem2score[x], reverse=False))

    sorted_elements = []
    for elem in prev_elements:
        for elem2 in initial_sorted_elements:
            if (not G.has_edge(elem, elem2) and not G.has_edge(elem2, elem)) or elem2 in sorted_elements:
                continue
            sorted_elements.append(elem2)

    return sorted_elements


def downstream_degree(G, node):
    return len(list(nx.descendants(G, node)))


# Layout nodes to represent a timeline
def timeline_layout(G, step_attr="step", scaleratio=1.0, scale=30):
    return multipartite_layout(G, subset_key=step_attr, scale=1, space_mult_x=scale * scaleratio, space_mult_y=scale*2, ordering=list(range(max(nx.get_node_attributes(G, step_attr).values()))))
    node2step = nx.get_node_attributes(G, step_attr)
    node2descendent_count = dict()
    for node in G.nodes:
        node2descendent_count[node] = len(list(nx.descendants(G, node)))
    all_steps = sorted(list(set(node2step.values())))

    step2ordered_nodes = dict()
    for step in all_steps:
        nodes = [n for n in G.nodes if node2step[n] == step]
        nodes = space_elements(nodes, [node2descendent_count[n] for n in nodes])
        step2ordered_nodes[step] = nodes


# Adapted from networkx to fix a weird crash, and add space multiplier between layers
# Also allow for setting the explict ordering of layers
def multipartite_layout(G, subset_key="subset", align="vertical", scale=1, center=None, space_mult_x=1, space_mult_y=1,
                        ordering=None, weigh_degree=True):
    """Position nodes in layers of straight lines.

    Parameters
    ----------
    G : NetworkX graph or list of nodes
        A position will be assigned to every node in G.

    subset_key : string (default='subset')
        Key of node data to be used as layer subset.

    align : string (default='vertical')
        The alignment of nodes. Vertical or horizontal.

    scale : number (default: 1)
        Scale factor for positions.

    center : array-like or None
        Coordinate pair around which to center the layout.

    Returns
    -------
    pos : dict
        A dictionary of positions keyed by node.

    Examples
    --------
    >>> G = nx.complete_multipartite_graph(28, 16, 10)
    >>> pos = nx.multipartite_layout(G)

    Notes
    -----
    This algorithm currently only works in two dimensions and does not
    try to minimize edge crossings.

    Network does not need to be a complete multipartite graph. As long as nodes
    have subset_key data, they will be placed in the corresponding layers.

    """
    if align not in ("vertical", "horizontal"):
        msg = "align must be either vertical or horizontal."
        raise ValueError(msg)

    G, center = nx.drawing.layout._process_params(G, center=center, dim=2)
    if len(G) == 0:
        return {}

    layers = {}
    for v, data in G.nodes(data=True):
        try:
            layer = data[subset_key]
        except KeyError:
            msg = "all nodes must have subset_key (default='subset') as data"
            raise ValueError(msg)
        layers[layer] = [v] + layers.get(layer, [])

    if ordering and len(ordering) == len(layers):
        layers_list = []
        for layer in ordering:
            layers_list.append((layer, layers[layer]))
        layers = layers_list
        del layers_list
    else:
        # Sort by layer, if possible
        try:
            layers = sorted(layers.items())
        except TypeError:
            layers = list(layers.items())

    pos = None
    nodes = []
    width = len(layers)
    height = max(len(l) for l in layers)
    prev_layer = None
    for i, (_, layer) in enumerate(layers):
        layer_elements = len(layer)
        xs = np.repeat(i * space_mult_x, layer_elements)
        dy = height / layer_elements
        # ys = np.arange(0, height, dy, dtype=float) * space_mult_y
        ys = np.linspace(0, height, xs.shape[0], dtype=float) * space_mult_y

        offset = ((width - 1) * space_mult_x / 2, (height - 1) * space_mult_y / 2)

        if layer_elements == 1:
            ys[0] = (height - 1) * space_mult_y

        layer_pos = np.column_stack([xs, ys]) - offset

        if pos is None:
            pos = layer_pos
        else:
            pos = np.concatenate([pos, layer_pos])
        if weigh_degree:
            if prev_layer is None:
                layer = space_elements(layer, [downstream_degree(G, n) for n in layer])
            else:
                layer = space_elements_from_previous(G, layer, prev_layer)
            prev_layer = layer
        nodes.extend(layer)
    pos = nx.rescale_layout(pos, scale=scale) + center
    if align == "horizontal":
        pos = pos[:, ::-1]  # swap x and y coords
    pos = dict(zip(nodes, pos))
    return pos


# Adapted from plotly quiver plots
# https://github.com/plotly/plotly.py/blob/master/packages/python/plotly/plotly/figure_factory/_quiver.py#L127
def get_quiver_arrows(start_x, start_y, end_x, end_y, scaleratio=1):
    """
    Creates lists of x and y values to plot the arrows
    Gets length of each barb then calculates the length of each side of
    the arrow. Gets angle of barb and applies angle to each side of the
    arrowhead. Next uses arrow_scale to scale the length of arrowhead and
    creates x and y values for arrowhead point1 and point2. Finally x and y
    values for point1, endpoint and point2s for each arrowhead are
    separated by a None and zipped to create lists of x and y values for
    the arrows.
    :rtype: (list, list) arrow_x, arrow_y: list of point1, endpoint, point2
        x_values separated by a None to create the arrowhead and list of
        point1, endpoint, point2 y_values separated by a None to create
        the barb of the arrow.
    """
    # Get the distance in both the x and y axis for the lines
    dif_x = [i - j for i, j in zip(end_x, start_x)]
    dif_y = [i - j for i, j in zip(end_y, start_y)]

    # Get barb lengths(default arrow length = 30% barb length)
    barb_len = [None] * len(start_x)
    for index in range(len(barb_len)):
        barb_len[index] = math.hypot(dif_x[index] / scaleratio, dif_y[index])

    # Make arrow lengths
    arrow_scale = 0.02
    arrow_len = [None] * len(start_x)
    # arrow_len = [i * arrow_scale for i in barb_len] orig dynamic barbs
    # arrow_len = [min(barb_len) * arrow_scale for i in barb_len]
    arrow_len = [.025 for i in barb_len]

    # Get barb angles
    barb_ang = [None] * len(start_x)
    for index in range(len(barb_ang)):
        barb_ang[index] = math.atan2(dif_y[index], dif_x[index] / scaleratio)

    # Set angles to create arrow
    angle = math.pi / 6
    ang1 = [i + angle for i in barb_ang]
    ang2 = [i - angle for i in barb_ang]

    cos_ang1 = [None] * len(ang1)
    for index in range(len(ang1)):
        cos_ang1[index] = math.cos(ang1[index])
    seg1_x = [i * j for i, j in zip(arrow_len, cos_ang1)]

    sin_ang1 = [None] * len(ang1)
    for index in range(len(ang1)):
        sin_ang1[index] = math.sin(ang1[index])
    seg1_y = [i * j for i, j in zip(arrow_len, sin_ang1)]

    cos_ang2 = [None] * len(ang2)
    for index in range(len(ang2)):
        cos_ang2[index] = math.cos(ang2[index])
    seg2_x = [i * j for i, j in zip(arrow_len, cos_ang2)]

    sin_ang2 = [None] * len(ang2)
    for index in range(len(ang2)):
        sin_ang2[index] = math.sin(ang2[index])
    seg2_y = [i * j for i, j in zip(arrow_len, sin_ang2)]

    # Set coordinates to create arrow
    for index in range(len(end_x)):
        point1_x = [i - j * scaleratio for i, j in zip(end_x, seg1_x)]
        point1_y = [i - j for i, j in zip(end_y, seg1_y)]
        point2_x = [i - j * scaleratio for i, j in zip(end_x, seg2_x)]
        point2_y = [i - j for i, j in zip(end_y, seg2_y)]

    # Combine lists to create arrow
    empty = [None] * len(end_x)
    arrow_x = list(itertools.chain.from_iterable(zip(point1_x, end_x, point2_x, empty)))
    arrow_y = list(itertools.chain.from_iterable(zip(point1_y, end_y, point2_y, empty)))
    return arrow_x, arrow_y


def scaled_logistic(x, min, max):
    # Based on: https://stackoverflow.com/a/29863846
    sig = np.exp(-np.logaddexp(0, -x))
    scaled = min + (max - min) * sig
    return float(scaled)
