from typing import List, Tuple, Optional

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
from dash import html, Output, Input
import dash_bio as dashbio

from viz.util import ColorTransformer, celltype_to_colors, saturate_color


def jumbotron(title, main_content, sub_content, *additional_content, dark=False):
    return html.Div(dbc.Container([
        html.H1(title, className='display-3'),
        html.P(main_content, className='lead'),
        html.Hr(className='my-2'),
        html.P(sub_content),
        *[html.P(content) for content in additional_content]
    ], fluid=True, className='py-3'), className="p-3 rounded-3 " + ('bg-dark text-white' if dark else 'bg-light'))


def control_panel_element(title, description, input, footer=None, outline=True) -> dbc.Card:
    elements = []
    if title:
        elements.append(dbc.CardHeader(title))
    elements.append(dbc.CardBody([
        input,
        html.P(description, className="card-text")
    ]))
    if footer:
        elements.append(dbc.CardFooter(footer))

    return dbc.Card(elements,
                    outline=outline,
                    color='light')


def control_panel(*element_rows: List[dbc.Card]) -> html.Div:
    return html.Div(
        [dbc.Row(dbc.CardGroup(cols)) for cols in element_rows],
        className='mb-3'
    )


def figure_output(title, footer, element) -> html.Div:
    return html.Div([
        dbc.Card([
            dbc.CardHeader(title),
            dbc.CardBody([
                html.P([
                    html.Div(id='spinner-holder'),
                    dbc.Progress(id='progress-bar',
                                 striped=True,
                                 animated=True,
                                 value=0,
                                 style={'visibility': 'hidden'}),
                    element
                ], className='card-text')
            ]),
            dbc.CardFooter(footer)
        ], color='light')
    ])


def interactive_panel(title, *content, dark=False):
    return dbc.Card([
        dbc.CardHeader(title),
        dbc.CardBody([
            html.Div(content)
        ])
    ], className='mb-3 ' + ('bg-dark text-white' if dark else 'bg-light'))


def wrap_icon(icon, *content, right=False, high_margin='.75em', low_margin='.18em'):
    elements = [(html.Span(c) if isinstance(c, str) else c) for c in content]
    icon = html.I(className='fas '+icon, style={'margin-left': high_margin if right else low_margin, 'margin-right': low_margin if right else high_margin})
    if right:
        elements.append(icon)
    else:
        elements.insert(0, icon)
    return html.Span(elements)


def make_data_redirect_buttons():
    return dbc.Row([
        dbc.Col(dbc.Button(wrap_icon('fa-arrows-left-right-to-line', 'View Cell Type Interactions'), id='ct-interaction-button', href='/interactions')),
        dbc.Col(dbc.Button(wrap_icon('fa-maximize', 'View Ligand Effects'), id='ligand-effect-button', href='/ligand-effects')),
        dbc.Col(dbc.Button(wrap_icon('fa-circle-dot', 'View Circos'), id='circos-button', href='/circos')),
    ])


def make_circos_figure(set_progress,
                       lr_pairs: List[Tuple[str, str]],
                       outer_data: pd.DataFrame,
                       max_data: Optional[pd.DataFrame],
                       is_max: bool,
                       inter_fdr: str,
                       logfc_fdr: str,
                       min_numsigi1: int,
                       min_numdeg: int,
                       min_logfc: float):
    set_progress((1, 7))

    #inter_pvalue_cutoff = float("0." + inter_fdr.replace('fdr', ''))  # pval from name
    logfc_pvalue_cutoff = float("0." + logfc_fdr.replace('fdr', ''))  # pval from name

    # Build outer ring for cell types
    celltypes = outer_data['cell_type'].unique()
    celltype2color = celltype_to_colors(celltypes)
    celltype2id = dict()
    celltype2targets = dict()

    from viz.data import _get_original_data
    old = _get_original_data()
    differing = pd.merge(old, outer_data, on=['cell_type', 'target'], how='left', suffixes=('', '_new'))

    # Join inner data to max data if necessary
    if is_max:
        outer_data = pd.merge(outer_data, max_data, on=['cell_type', 'target', 'ligand', 'receptor'], how='left', suffixes=('', '_max'))
    del max_data

    # Initial obs filters
    outer_data = outer_data[~outer_data['cell_type_dc1'].isna()]
    #outer_data['cell_type_dc1'] = outer_data['cell_type_dc1'].fillna(0)
    outer_data = outer_data[(outer_data['receptor'] & (outer_data['numSigI1_max' if is_max else 'numSigI1'] >= min_numsigi1) & (outer_data['numDEG_max' if is_max else 'numDEG'] >= min_numdeg))
                            | ((outer_data['ligand'] & ~outer_data['receptor']) & (outer_data['MAST_log2FC_max' if is_max else 'MAST_log2FC'].abs() > min_logfc) & (outer_data['MAST_fdr_max' if is_max else 'MAST_fdr'] < logfc_pvalue_cutoff))]

    outer_data['logfc_direction'] = np.sign(outer_data['MAST_log2FC'])

    if is_max:
        # Require same direction of change
        outer_data['logfc_direction_cin'] = np.sign(outer_data['MAST_log2FC_cin'])
        outer_data['logfc_direction_sting'] = np.sign(outer_data['MAST_log2FC_sting'])
        outer_data = outer_data[(~outer_data['ligand']) | (outer_data['logfc_direction_cin'] == outer_data['logfc_direction_sting'])]

    all_receptors = set(outer_data[outer_data['receptor']].target.unique())
    all_ligands = set(outer_data[outer_data['ligand']].target.unique())
    lr_pairs = [(l, r) for l, r in lr_pairs if (l in all_ligands and r in all_receptors)]
    all_ligands = set(l for l, r in lr_pairs)
    all_receptors = set(r for l, r in lr_pairs)

    # Filter outer_data to just the selected ligands and receptors
    outer_data = outer_data[(outer_data['target'].isin(all_ligands)) | (outer_data['target'].isin(all_receptors))]

    # Filter obs to just the selected ligands and receptors
    set_progress((2, 7))
    layout = []
    for celltype in celltypes:
        color = celltype2color[celltype]
        label = celltype
        id = celltype.replace(" ", "_").replace("-", "_").replace("/", "_")
        celltype2id[celltype] = id
        celltype2targets[celltype] = sorted(outer_data[outer_data['cell_type'] == celltype].target, key=lambda t: outer_data[(outer_data['cell_type'] == celltype) & (outer_data['target'] == t)]['cell_type_dc1'].values[0])
        layout.append({
            'id': id,
            'label': label,
            'color': color,
            'len': len(celltype2targets[celltype]),
        })
    layout = list(sorted(layout, key=lambda x: x['len'], reverse=True))
    set_progress((3, 7))

    # Build next ring for DC1 heatmap
    max_dc1 = outer_data['cell_type_dc1'].max()
    min_dc1 = outer_data['cell_type_dc1'].min()
    colormap = ColorTransformer(min_dc1, max_dc1, 'cividis')
    diffusion_data = []
    for celltype in celltypes:
        id = celltype2id[celltype]
        targets = celltype2targets[celltype]
        for i, t in enumerate(targets):
            dc1 = outer_data[(outer_data['cell_type'] == celltype) & (outer_data['target'] == t)]['cell_type_dc1'].values[0]
            # Normalize to 0-1
            #dc1 = (dc1 - min_dc1) / (max_dc1 - min_dc1)
            diffusion_data.append({
                'block_id': id,
                'start': i,
                'end': (i + 1),
                'value': dc1,
                'value_text': f"DC1={dc1:.2f}",
                'color': colormap(dc1),
                'target': t
            })
    set_progress((4, 7))

    # Next ring for Differential abundance
    max_da = outer_data['DA_score'].max()
    min_da = outer_data['DA_score'].min()
    colormap = ColorTransformer(min_da, max_da, 'seismic')
    da_data = []
    for celltype in celltypes:
        id = celltype2id[celltype]
        targets = celltype2targets[celltype]
        for i, t in enumerate(targets):
            da = outer_data[(outer_data['cell_type'] == celltype) & (outer_data['target'] == t)]['DA_score'].values[0]
            da_data.append({
                'block_id': id,
                'start': i,
                'end': (i + 1),
                'value': da,
                'value_text': f"DA={da:.2f}",
                'color': colormap(da),
                'target': t
            })
    set_progress((5, 7))

    # Next ring for magnitude of CIN-dependent effect
    receptor_info = outer_data[outer_data['receptor']]
    max_receptor_numSigI1 = receptor_info['numSigI1'].max()
    min_receptor_numSigI1 = receptor_info['numSigI1'].min()
    numSigI1_data = []
    for celltype in celltypes:
        id = celltype2id[celltype]
        color = saturate_color(celltype2color[celltype], .75)
        targets = celltype2targets[celltype]
        for i, t in enumerate(targets):
            # Set to 0 if target is not a receptor
            if t not in receptor_info['target'].values:
                numSigI1_data.append({
                    'block_id': id,
                    'start': i,
                    'end': (i + 1),
                    'value': 0,
                    'value_text': "",
                    'color': color,
                    'target': t
                })
            else:
                lig_effect = receptor_info[receptor_info['target'] == t]['numSigI1'].values[0]
                # Normalize to 0 minimum
                lig_effect -= min_receptor_numSigI1
                numSigI1_data.append({
                    'block_id': id,
                    'start': i,
                    'end': (i + 1),
                    'value': lig_effect,
                    'value_text': f"numSigI1={lig_effect:d}",
                    'color': color,
                    'target': t
                })
    set_progress((6, 7))

    # Next ring for chords connecting ligands to receptors
    chord_receptors = outer_data[outer_data['receptor']]
    chord_ligands = outer_data[outer_data['ligand']]
    colormap = ColorTransformer(-0.2, 0.2, 'bwr')
    chord_data = []
    text_data = dict()
    for (l, r) in lr_pairs:
        if l not in chord_ligands['target'].values or r not in chord_receptors['target'].values:
            continue
        for _, ligand_row in chord_ligands[chord_ligands['target'] == l].iterrows():
            for _, receptor_row in chord_receptors[chord_receptors['target'] == r].iterrows():
                source_celltype = ligand_row['cell_type']
                target_celltype = receptor_row['cell_type']
                source_position = celltype2targets[source_celltype].index(l)
                target_position = celltype2targets[target_celltype].index(r)

                # Max thickness of ribbons on either end of the target position
                thickness = max((receptor_row['numSigI1_max' if is_max else 'numSigI1']/max_receptor_numSigI1), 1)

                lig = ligand_row['target']
                rec = receptor_row['target']
                text_data[(source_celltype, lig)] = {
                    'block_id': celltype2id[source_celltype],
                    'position': source_position+.5,
                    'value': lig
                }
                text_data[(target_celltype, rec)] = {
                    'block_id': celltype2id[target_celltype],
                    'position': target_position+.5,
                    'value': rec
                }

                chord_data.append({
                    'color': colormap(ligand_row['MAST_log2FC_max' if is_max else 'MAST_log2FC']),
                    'value_text': f"{lig}/{rec} MAST_log2FC={ligand_row['MAST_log2FC_max' if is_max else 'MAST_log2FC']:.2f}",
                    'logfc': ligand_row['MAST_log2FC_max' if is_max else 'MAST_log2FC'],
                    'source': {
                        'id': celltype2id[source_celltype],
                        'start': source_position - thickness,
                        'end': source_position + thickness,
                    },
                    'target': {
                        'id': celltype2id[target_celltype],
                        'start': target_position - thickness,
                        'end': target_position + thickness,
                    },
                })
    # Sort to place red chords on top
    chord_data = sorted(chord_data, key=lambda c: c['logfc'], reverse=False)

    set_progress((7, 7))
    ring_width = 50
    return dashbio.Circos(
        enableDownloadSVG=True,
        enableZoomPan=True,
        layout=layout,
        selectEvent={
            "0": "hover",
            #"1": "hover",
            #"2": "hover",
            #"3": "hover",
            "4": "hover"
        },
        tracks=[{
            'type': 'TEXT',
            'data': list(text_data.values()),
            'config': {
                'innerRadius': 2.8*ring_width,
                'outerRadius': 3.2*ring_width,
                'style': {
                    'font-size': 6
                }
            }
        }, {
            'type': 'CHORDS',
            'data': chord_data,
            'config': {
                'color': {
                    'name': 'color',
                },
                'opacity': 0.9,
                'radius': 2.75*ring_width,
                'tooltipContent': {
                    'source': 'source',
                    'sourceID': 'id',
                    'target': 'target',
                    'targetID': 'id',
                    'targetEnd': 'value_text'
                }
            }
        }, {
            'type': 'HISTOGRAM',
            'data': numSigI1_data,
            'config': {
                'color': {
                    'name': 'color',
                },
                'innerRadius': 3*ring_width,
                'outerRadius': 4*ring_width,
                'direction': 'in',
                'min': 0,
                'max': max_receptor_numSigI1 - min_receptor_numSigI1,
                'tooltipContent': {
                    'name': 'value_text'
                }
            }
        }, {
            'type': 'HEATMAP',
            'data': da_data,
            'config': {
                'color': {
                    'name': 'color',
                },
                'innerRadius': 4*ring_width,
                'outerRadius': 5*ring_width,
                'tooltipContent': {
                    'name': 'value_text'
                }
            }
        }, {
            'type': 'HEATMAP',
            'data': diffusion_data,
            'config': {
                'color': {
                    'name': 'color',  # Redirect to color property
                },
                'innerRadius': 5*ring_width,
                'outerRadius': 6*ring_width,
                'tooltipContent': {
                    'name': 'value_text'
               }
            }
        }],
        config={
            "labels": {
                "display": True
            },
            "ticks": {
                "display": False
            },
            "innerRadius": 6*ring_width,
            "outerRadius": 7*ring_width,
            'size': 800
            #"cornerRadius": 4,
        }
    )
