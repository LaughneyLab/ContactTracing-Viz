from typing import List, Tuple

import dash_bootstrap_components as dbc
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
        dbc.Col(dbc.Button(wrap_icon('fa-link-slash', 'Unselect Data'), id='unselect-button', n_clicks=0))  # Callback in pages/home
    ])


def get_obs_columns(fdr_name: str, condition_name: str) -> Tuple[str, ...]:
    cell_type = "cell type"
    target = "target"
    receptor = "receptor"
    ligand = "ligand"
    numDEG = f"numDEG_{fdr_name}{condition_name}"
    numSigI1 = f"numSigI1_{fdr_name}{condition_name}"
    MAST_log2FC = f"MAST_log2FC{condition_name}"
    MAST_fdr = f"MAST_fdr{condition_name}"
    cell_type_dc1 = "cell_type_dc1"
    DA_score = "DA_score"
    return cell_type, target, receptor, ligand, numDEG, numSigI1, MAST_log2FC, MAST_fdr, cell_type_dc1, DA_score


def make_circos_figure(set_progress, obs: pd.DataFrame, lr_pairs: List[Tuple[str, str]], condition_name: str,
                       fdr_name: str, min_numsigi1: int, min_numdeg: int, chord_logfc_cutoff: float, chord_numsigi1_cutoff: int):
    set_progress((1, 7))

    cell_type, target, receptor, ligand, numDEG, numSigI1, MAST_log2FC, MAST_fdr, cell_type_dc1, DA_score = get_obs_columns(fdr_name, condition_name)

    pvalue_cutoff = float("0." + fdr_name[3:])  # pval from name

    # Build outer ring for cell types
    celltypes = obs[cell_type].unique()
    celltype2color = celltype_to_colors(celltypes)
    celltype2id = dict()
    celltype2targets = dict()
    # Initial obs filters
    obs = obs[(obs[MAST_fdr] < pvalue_cutoff) & (~obs[cell_type_dc1].isna())]
    all_receptors = obs[(obs[receptor]) & (obs[numSigI1] >= min_numsigi1) & (obs[numDEG] >= min_numdeg)]
    all_ligands = []
    for (l, r) in lr_pairs:
        if r in all_receptors[target].values:
            all_ligands.append(l)
    # Filter obs to just the selected ligands and receptors
    obs = obs[(obs[target].isin(all_ligands)) | (obs[target].isin(all_receptors[target].values))]
    set_progress((2, 7))
    layout = []
    for celltype in celltypes:
        color = celltype2color[celltype]
        label = celltype
        id = celltype.replace(" ", "_").replace("-", "_").replace("/", "_")
        celltype2id[celltype] = id
        celltype2targets[celltype] = sorted(obs[obs[cell_type] == celltype].target, key=lambda t: obs[(obs[cell_type] == celltype) & (obs[target] == t)][cell_type_dc1].values[0])
        layout.append({
            'id': id,
            'label': label,
            'color': color,
            'len': len(celltype2targets[celltype]),
        })
    layout = list(sorted(layout, key=lambda x: x['len'], reverse=True))
    set_progress((3, 7))

    # Build next ring for DC1 heatmap
    max_dc1 = obs[cell_type_dc1].max()
    min_dc1 = obs[cell_type_dc1].min()
    colormap = ColorTransformer(min_dc1, max_dc1, 'cividis')
    diffusion_data = []
    for celltype in celltypes:
        id = celltype2id[celltype]
        targets = celltype2targets[celltype]
        for i, t in enumerate(targets):
            dc1 = obs[(obs[cell_type] == celltype) & (obs[target] == t)][cell_type_dc1].values[0]
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
    max_da = obs[DA_score].max()
    min_da = obs[DA_score].min()
    colormap = ColorTransformer(min_da, max_da, 'seismic')
    da_data = []
    for celltype in celltypes:
        id = celltype2id[celltype]
        targets = celltype2targets[celltype]
        for i, t in enumerate(targets):
            da = obs[(obs[cell_type] == celltype) & (obs[target] == t)][DA_score].values[0]
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
    max_receptor_numSigI1 = all_receptors[numSigI1].max()
    min_receptor_numSigI1 = all_receptors[numSigI1].min()
    numSigI1_data = []
    for celltype in celltypes:
        id = celltype2id[celltype]
        color = saturate_color(celltype2color[celltype], .75)
        targets = celltype2targets[celltype]
        for i, t in enumerate(targets):
            # Set to 0 if target is not a receptor
            if t not in all_receptors[target].values:
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
                lig_effect = all_receptors[all_receptors[target] == t][numSigI1].values[0]
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
    chord_receptors = obs[(obs[receptor]) & (obs[numSigI1] >= chord_numsigi1_cutoff)]
    chord_ligand = obs[(obs[ligand]) & (obs[MAST_log2FC].abs() >= chord_logfc_cutoff) & (obs[MAST_fdr] <= pvalue_cutoff)]
    colormap = ColorTransformer(-0.2, 0.2, 'bwr')
    chord_data = []
    text_data = dict()
    for (l, r) in lr_pairs:
        if l not in chord_ligand[target].values or r not in chord_receptors[target].values:
            continue
        for _, ligand_row in chord_ligand[chord_ligand[target] == l].iterrows():
            for _, receptor_row in chord_receptors[chord_receptors[target] == r].iterrows():
                source_celltype = ligand_row[cell_type]
                target_celltype = receptor_row[cell_type]
                source_position = celltype2targets[source_celltype].index(l)
                target_position = celltype2targets[target_celltype].index(r)

                # Max thickness of ribbons on either end of the target position
                thickness = max(10 * (receptor_row[numSigI1]/obs[numSigI1].max()), 1)

                lig = ligand_row[target]
                rec = receptor_row[target]
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
                    'color': colormap(ligand_row[MAST_log2FC]),
                    'value_text': f"{lig}/{rec} MAST_log2FC={ligand_row[MAST_log2FC]:.2f}",
                    'logfc': ligand_row[MAST_log2FC],
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
