from typing import List, Tuple, Optional
from uuid import uuid4

import dash
import dash_bootstrap_components as dbc
import numpy as np
import plotly.graph_objects as go
import pandas as pd
from dash_extensions.enrich import html, Output, Input, dcc, callback, State
import dash_bio as dashbio
from dash.exceptions import PreventUpdate

from viz.util import ColorTransformer, celltype_to_colors, saturate_color, smooth_step, \
    brighten_color, brighten_and_saturate_color


def make_tooltip(content, tooltip_content):
    # uuid4 is used to generate a unique id for the tooltip
    tooltip_id = str(uuid4())
    content_span = html.Span(content,
                             id=tooltip_id,
                             style={'display': 'inline', 'textDecoration': 'underline dotted', 'cursor': 'help'})
    tooltip = dbc.Tooltip(tooltip_content,
                          target=tooltip_id,
                          placement='auto')

    return content_span, tooltip


def _slider_transform(x):
    if isinstance(x, str):
        if 'fdr' in x:
            return x
        else:
            return f"fdr{x.split('.')[1]}"
    else:
        return "fdr" + f"{x:.2f}".split(".")[1]


def make_fdr_slider(id: str, value) -> html.Div:
    if isinstance(value, str):
        value = float("." + value[3:])
    return make_custom_slider(
        id, 0.01, 0.25, value, 0.01,
        _slider_transform
    )


def make_custom_slider(id: str, min, max, value, step, transform=None) -> html.Div:
    # Generate a slider with a textbox
    # Use the 'data' field of this component to get the value
    slider = dcc.Slider(
        id=id + '-slider',
        min=min,
        max=max,
        value=value,
        step=step,
        marks=None,
        tooltip={'placement': 'bottom'},
        className='form-range ct-slider'
    )
    input = dbc.Input(id=id + '-input',
                      type='number',
                      min=min,
                      max=max,
                      step=step,
                      value=value)
    value_store = dcc.Store(id=id, data=transform(value) if transform else value)
    slider_element = dbc.Container([
        value_store,
        dbc.Row([
            dbc.Col(slider),
            dbc.Col(input, className='slider_text_input')
        ], align='center')
    ], fluid=True)

    @callback(
        Output(value_store, 'data'),
        Output(input, 'value'),
        Output(slider, 'value'),
        Input(input, 'value'),
        Input(slider, 'value'),
        Input(value_store, 'data')
    )
    def update_input(text_input, slider_input, slider_store):
        if dash.callback_context.triggered_id == id + '-input':
            # Input should always be in range
            if text_input is None or (text_input < min or text_input > max):
                raise PreventUpdate
            # Round to the nearest multiple of the step
            rounded_input = round(text_input / step) * step
            # If step is an integer, cast to int
            if step == int(step):
                rounded_input = int(rounded_input)
            return transform(rounded_input) if transform else rounded_input, rounded_input, rounded_input
        elif dash.callback_context.triggered_id == id + '-slider':
            # Slider should always be in range
            return transform(slider_input) if transform else slider_input, slider_input, slider_input
        elif dash.callback_context.triggered_id == id:  # Store updated
            return transform(slider_store) if transform else slider_store, slider_input, slider_input

        return transform(value) if transform else value, value, value

    return slider_element


def jumbotron(title, main_content, sub_content, *additional_content, dark=False):
    return html.Div(dbc.Container([
        html.H1(title, className='display-3 d-none d-md-block'),
        html.P(),
        html.Div(main_content, className='lead'),
        html.Hr(className='my-2'),
        html.P(),
        html.Div(sub_content),
        *[html.Div([html.P(), content]) for content in additional_content],
    ], fluid=True, className='py-3'), className="p-3 rounded-3 " + ('dark text-white' if dark else 'light'))


def control_panel_element(title, description, input, footer=None, outline=True) -> dbc.Card:
    elements = []
    if title:
        elements.append(dbc.CardHeader(title))
    elements.append(dbc.CardBody([
        input,
        html.P(),
        html.Div(description, className="card-text")
    ]))
    if footer:
        elements.append(dbc.CardFooter(footer))

    return dbc.Card(elements,
                    outline=outline,
                    color='light')


def control_panel(submit_btn_id: str, *element_rows: List[dbc.Card]) -> html.Div:
    rows = [dbc.Row(dbc.CardGroup(cols)) for cols in element_rows]
    # Add all rows to a container that will be within an Accordion except for the
    # first row (which contains the submit button!)

    options_accordion = dbc.Accordion(
        dbc.AccordionItem(dbc.Container(rows[1:], fluid=True), title="Figure Options", item_id='options'),
        start_collapsed=False, flush=False, persistence=False
    )

    if submit_btn_id:
        @callback(
            Output(options_accordion, 'active_item'),
            Input(submit_btn_id, 'n_clicks')
        )
        def hide_options(n_clicks):
            if n_clicks and n_clicks > 0:
                return []
            raise PreventUpdate

    return html.Div([
        rows[0],
        options_accordion
    ],
        className='mb-3'
    )


def figure_output(title, footer, element, help_info, download_btn_id, outline=True) -> html.Div:
    if download_btn_id:
        download_button = dbc.Button(wrap_icon('fa-file-csv', '   ', low_margin='', high_margin='', right=False), id=download_btn_id,
                                     color='primary', outline=True, className='float-end', type='button', n_clicks=0)
    else:
        download_button = None
    help_button = dbc.Button(wrap_icon('fa-circle-question', low_margin='', high_margin='', right=True),
                             color='primary', outline=True, className='float-end', type='button', n_clicks=0)
    close_button = dbc.Button("Close", className="ms-auto", n_clicks=0)
    help_modal = dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle(wrap_icon('fa-info-circle', "Plot Information"))),
        dbc.ModalBody(help_info),
        dbc.ModalFooter(close_button)
    ], is_open=False, size='xl', scrollable=True, centered=True)

    @callback(
        Output(help_modal, "is_open"),
        Input(help_button, "n_clicks"),
        Input(close_button, "n_clicks"),
        State(help_modal, "is_open")
    )
    def toggle_help(open_click, close_click, is_open):
        if open_click or close_click:
            return not is_open
        return is_open

    return html.Div([
        help_modal,
        dbc.Card([
            dbc.CardHeader(dbc.Row([
                dbc.Col(html.H5(title, className='card-title')),
                html.Div(className="w-100 d-md-none"),
                dbc.Col(dbc.ButtonGroup(
                    ([download_button] if download_button else []) + [help_button],
                    className='help-button-group'
                ), width=2),
            ])),
            dbc.CardBody([
                html.P(),
                html.Div([
                    html.Div(id='spinner-holder'),
                    dbc.Progress(id='progress-bar',
                                 striped=True,
                                 animated=True,
                                 value=0,
                                 style={'visibility': 'hidden'}),
                    element,
                    html.Br()
                ])
            ]),
            dbc.CardFooter(footer)
        ], color='light', outline=outline)
    ])


def interactive_panel(title, subtitle, *content, dark=False):
    return html.Div(dbc.Container([
        html.H3(title, className='display-7'),
        html.P(),
        html.Div(subtitle, className='lead'),
        html.Hr(className='my-2'),
        html.Br(),
        html.Div(dbc.Container(content), className='card-text')
    ], fluid=True, className='py-3'), className="p-3 rounded-3 " + ('dark text-white' if dark else 'light'))


def wrap_icon(icon, *content, right=False, high_margin='.75em', low_margin='.18em'):
    elements = [(html.Span(c) if isinstance(c, str) else c) for c in content]
    icon = html.I(className='fas '+icon, style={'marginLeft': high_margin if right else low_margin, 'marginRight': low_margin if right else high_margin})
    if right:
        elements.append(icon)
    else:
        elements.insert(0, icon)
    return html.Span(elements)


def make_data_redirect_buttons():
    return html.Div([
                dbc.Button(wrap_icon('fa-circle-dot', 'View Circos'), id='circos-button', href='/circos'),
                html.Span(" "),
                dbc.Button(wrap_icon('fa-arrows-left-right-to-line', 'View Cell Type Interactions'),
                           id='ct-interaction-button', href='/interactions'),
                html.Span(" "),
                dbc.Button(wrap_icon('fa-maximize', 'View Ligand Effects'), id='ligand-effect-button',
                           href='/ligand-effects')
            ], className='text-center d-grid d-md-block gap-1 col-gap-2')


def make_circos_figure(set_progress, progress_offset: int,
                       outer_data: pd.DataFrame,
                       inter_data: pd.DataFrame,
                       cin_only: bool,
                       logfc_fdr: str,
                       min_numsigi1: int,
                       min_numdeg: int,
                       min_logfc: float,
                       highlighted_genes: str):
    if set_progress is not None:
        set_progress((1+(progress_offset*7), 7*(progress_offset+1)))

    if highlighted_genes is not None and highlighted_genes != '':
        highlighted_genes = [g.strip() for g in highlighted_genes.split(',')]
        should_highlight = len(highlighted_genes) > 0
    else:
        highlighted_genes = []
        should_highlight = False

    #inter_pvalue_cutoff = float("0." + inter_fdr.replace('fdr', ''))  # pval from name
    logfc_pvalue_cutoff = float("0." + logfc_fdr.replace('fdr', ''))  # pval from name

    # Build outer ring for cell types
    celltypes = outer_data['cell_type'].unique()
    celltype2color = celltype_to_colors(celltypes)
    celltype2id = dict()
    celltype2targets = dict()

    # Initial obs filters
    outer_data = outer_data[~outer_data['cell_type_dc1'].isna()]

    #outer_data['cell_type_dc1'] = outer_data['cell_type_dc1'].fillna(0)
    # Ligand / Receptor filter
    outer_data = outer_data[(outer_data['receptor'] & (outer_data['numSigI1_cin'] > 0) & (outer_data['numSigI1_sting'] > 0))
                            | (outer_data['ligand'] & (outer_data['MAST_fdr_cin'] < logfc_pvalue_cutoff) & (outer_data['MAST_fdr_sting'] < logfc_pvalue_cutoff) & (np.sign(outer_data['MAST_log2FC_cin']) == np.sign(outer_data['MAST_log2FC_sting'])))]

    # Only select interactions present in outer_data
    inter_ligand_index = pd.MultiIndex.from_frame(inter_data[['ligand', 'cell_type_ligand']])
    outer_ligand_index = pd.MultiIndex.from_frame(outer_data[['target', 'cell_type']])
    inter_data = inter_data.loc[inter_ligand_index.isin(outer_ligand_index)]
    inter_receptor_index = pd.MultiIndex.from_frame(inter_data[['receptor', 'cell_type_receptor']])
    outer_receptor_index = pd.MultiIndex.from_frame(outer_data[['target', 'cell_type']])
    inter_data = inter_data.loc[inter_receptor_index.isin(outer_receptor_index)]
    del inter_ligand_index, outer_ligand_index, inter_receptor_index, outer_receptor_index

    # Filter outer data to just receptors and ligands with interactions
    outer_data = outer_data[outer_data['receptor'] | outer_data['target'].isin(inter_data['ligand'])]

    inter_data = inter_data[(inter_data['numSigI1'] >= min_numsigi1) & (inter_data['numDEG'] >= min_numdeg)
                            & (inter_data['MAST_fdr_ligand'] < logfc_pvalue_cutoff) & (inter_data['MAST_log2FC_ligand'].abs() > min_logfc)]

    # Annotate nodes to highlight
    if should_highlight:
        outer_data['highlight'] = outer_data['target'].isin(highlighted_genes)
        inter_data['highlight'] = inter_data['ligand'].isin(highlighted_genes) | inter_data['receptor'].isin(highlighted_genes)

    def maybe_brighten(colormap, value, *genes):
        if should_highlight:
            if not any([g in genes for g in highlighted_genes]):
                return brighten_and_saturate_color(colormap(value) if colormap is not None else value, 1.75, .5)
            else:
                return brighten_and_saturate_color(colormap(value) if colormap is not None else value, 0.95, 1.05)
        else:
            return colormap(value) if colormap is not None else value

    celltypes = outer_data['cell_type'].unique()

    # Filter obs to just the selected ligands and receptors
    if set_progress is not None:
        set_progress((2+(progress_offset*7), 7*(progress_offset+1)))

    layout = []
    for celltype in celltypes:
        color = celltype2color[celltype]
        label = celltype
        id = celltype.replace(" ", "_").replace("-", "_").replace("/", "_")
        celltype2id[celltype] = id
        celltype2targets[celltype] = outer_data[outer_data['cell_type'] == celltype].sort_values(['cell_type_dc1', 'DA_score']).target.tolist()
        target_count = len(celltype2targets[celltype])
        layout.append({
            'id': id,
            'label': label if target_count > len(label)*2 else "",
            'color': color,
            'len': target_count,
            'celltype': celltype,
            'receiving_celltypes': inter_data[inter_data['cell_type_ligand'] == celltype].cell_type_receptor.unique().tolist(),
            'sending_celltypes': inter_data[inter_data['cell_type_receptor'] == celltype].cell_type_ligand.unique().tolist(),
            'direct_indirect_targets': list(set(inter_data[inter_data['cell_type_ligand'] == celltype].ligand.str.lower().tolist() + inter_data[inter_data['cell_type_receptor'] == celltype].receptor.str.lower().tolist())),
        })
    # Hard coded ordering according to the paper
    ct2order = {
        'Macrophages/mMDSC': 8,
        'Tumor cells': 7,
        'Endothelial cells': 6,
        'T cells': 5,
        'PMN/gMDSC': 4,
        'B cells': 3,
        'pDC': 2,
        'cDC': 1,
        'NK cells': 0
    }
    layout = list(sorted(layout, key=lambda x: (ct2order.get(x['celltype'], -1), x['len']), reverse=True))
    if set_progress is not None:
        set_progress((3+(progress_offset*7), 7*(progress_offset+1)))

    def get_data(df, celltype, target, field):
        dat = df[(df['cell_type'] == celltype) & (df['target'] == target)]
        ligand = dat['ligand'].values[0]
        receptor = dat['receptor'].values[0]
        gene_type = 'Gene'
        if ligand and receptor:
            gene_type = 'Ligand/Receptor'
            ligand_partners_df = inter_data[(inter_data['cell_type_ligand'] == celltype) & (inter_data['ligand'] == target)]
            receptor_partners_df = inter_data[(inter_data['cell_type_receptor'] == celltype) & (inter_data['receptor'] == target)]
            ligand_partners = list(zip(ligand_partners_df['cell_type_receptor'].tolist(), ligand_partners_df['receptor'].str.lower().tolist()))
            receptor_partners = list(zip(receptor_partners_df['cell_type_ligand'].tolist(), receptor_partners_df['ligand'].str.lower().tolist()))
            partners = ligand_partners + receptor_partners
        elif ligand:
            gene_type = 'Ligand'
            partners_df = inter_data[(inter_data['cell_type_ligand'] == celltype) & (inter_data['ligand'] == target)]
            partners = list(zip(partners_df['cell_type_receptor'].tolist(), partners_df['receptor'].str.lower().tolist()))
        elif receptor:
            gene_type = 'Receptor'
            partners_df = inter_data[(inter_data['cell_type_receptor'] == celltype) & (inter_data['receptor'] == target)]
            partners = list(zip(partners_df['cell_type_ligand'].tolist(), partners_df['ligand'].str.lower().tolist()))
        partners_dict = dict()
        for ct, g in partners:
            if ct not in partners_dict:
                partners_dict[ct] = []
            partners_dict[ct].append(g)
        return dat[field].values[0], partners_dict, gene_type

    # Build next ring for DC1 heatmap
    dc1_colormap = ColorTransformer(-1, 1, 'cividis')
    diffusion_data = []
    for celltype in celltypes:
        id = celltype2id[celltype]
        targets = celltype2targets[celltype]
        for i, t in enumerate(targets):
            dc1, partners, gene_type = get_data(outer_data, celltype, t, 'cell_type_dc1')
            # Normalize to 0-1
            #dc1 = (dc1 - min_dc1) / (max_dc1 - min_dc1)
            diffusion_data.append({
                'block_id': id,
                'start': i,
                'end': (i + 1),
                'value': dc1,
                'value_text': f"<h4>Cell Type: {celltype}<br>Target ({gene_type}): {t}<br>DC1: {dc1:.2f}</h4>",
                'color': dc1_colormap(dc1),
                'celltype': celltype,
                'target': t,
                'partners': partners
            })
    if set_progress is not None:
        set_progress((4+(progress_offset*7), 7*(progress_offset+1)))

    # Next ring for Differential abundance
    da_colormap = ColorTransformer(-0.5, 0.5, 'RdBu_r')
    da_data = []
    for celltype in celltypes:
        id = celltype2id[celltype]
        targets = celltype2targets[celltype]
        for i, t in enumerate(targets):
            da, partners, gene_type = get_data(outer_data, celltype, t, 'DA_score')
            da_data.append({
                'block_id': id,
                'start': i,
                'end': (i + 1),
                'value': da,
                'value_text': f"<h4>Cell Type: {celltype}<br>Target ({gene_type}): {t}<br>DA: {da:.2f}</h4>",
                'color': da_colormap(da),
                'celltype': celltype,
                'target': t,
                'partners': partners
            })
    if set_progress is not None:
        set_progress((5+(progress_offset*7), 7*(progress_offset+1)))

    # Next ring for magnitude of CIN-dependent effect
    max_receptor_numSigI1 = inter_data['numSigI1'].max()
    min_receptor_numSigI1 = inter_data['numSigI1'].min()
    numSigI1_data = []
    for celltype in celltypes:
        id = celltype2id[celltype]
        color = saturate_color(celltype2color[celltype], .85)
        targets = celltype2targets[celltype]
        receptor_info = inter_data[inter_data['cell_type_receptor'] == celltype]
        for i, t in enumerate(targets):
            # Set to 0 if target is not a receptor
            if t not in receptor_info['receptor'].values:
                numSigI1_data.append({
                    'block_id': id,
                    'start': i,
                    'end': (i + 1),
                    'value': 0,
                    'value_text': "",
                    'color': color,
                    'celltype': celltype,
                    'target': t,
                    'partners': dict()
                })
            else:
                numSigI1 = receptor_info[(receptor_info['receptor'] == t) & (receptor_info['cell_type_receptor'] == celltype)]['numSigI1'].values[0]
                partners_df = inter_data[(inter_data['cell_type_receptor'] == celltype) & (inter_data['receptor'] == t)]
                partners = list(zip(partners_df['cell_type_ligand'].tolist(), partners_df['ligand'].str.lower().tolist()))
                # Normalize to 0 minimum
                lig_effect = numSigI1 - min_receptor_numSigI1
                numSigI1_data.append({
                    'block_id': id,
                    'start': i,
                    'end': (i + 1),
                    'value': max(lig_effect, 0),
                    'value_text': f"<h4>Cell Type: {celltype}<br>Target (Receptor): {t}<br>numSigI: {numSigI1:d}</h4>",
                    'color': color,
                    'celltype': celltype,
                    'target': t,
                    'partners': partners
                })
    if set_progress is not None:
        set_progress((6+(progress_offset*7), 7*(progress_offset+1)))

    # Next ring for chords connecting ligands to receptors
    log2fc_colormap = ColorTransformer(-0.2, 0.2, 'bwr', alpha=0.8)
    chord_data = []
    text_data = dict()
    for i, inter_row in inter_data.iterrows():
        lig = inter_row['ligand']
        rec = inter_row['receptor']
        source_celltype = inter_row['cell_type_ligand']
        target_celltype = inter_row['cell_type_receptor']
        source_position = celltype2targets[source_celltype].index(lig)
        target_position = celltype2targets[target_celltype].index(rec)

        # Max thickness of ribbons on either end of the target position
        thickness = smooth_step(inter_row['numSigI1'] / max_receptor_numSigI1, 0.33, 3)

        # Allow for overwriting if highlighting
        highlight_chord = (not should_highlight) or (lig in highlighted_genes or rec in highlighted_genes)
        if should_highlight and highlight_chord and (source_celltype, lig) in text_data:
            del text_data[(source_celltype, lig)]
        if should_highlight and highlight_chord and (target_celltype, rec) in text_data:
            del text_data[(target_celltype, rec)]

        if (source_celltype, lig) not in text_data:
            text_data[(source_celltype, lig)] = {
                'block_id': celltype2id[source_celltype],
                'position': source_position+.5,
                'value': lig if highlight_chord else "",
                "celltype": source_celltype,
                'partners': {target_celltype: [rec.lower()]},
            }
        else:
            if target_celltype not in text_data[(source_celltype, lig)]['partners']:
                text_data[(source_celltype, lig)]['partners'][target_celltype] = [rec.lower()]
            else:
                text_data[(source_celltype, lig)]['partners'][target_celltype].append(rec.lower())

        if (target_celltype, rec) not in text_data:
            text_data[(target_celltype, rec)] = {
                'block_id': celltype2id[target_celltype],
                'position': target_position+.5,
                'value': rec if highlight_chord else "",
                "celltype": target_celltype,
                'partners': {source_celltype: [lig.lower()]},
            }
        else:
            if source_celltype not in text_data[(target_celltype, rec)]['partners']:
                text_data[(target_celltype, rec)]['partners'][source_celltype] = [lig.lower()]
            else:
                text_data[(target_celltype, rec)]['partners'][source_celltype].append(lig.lower())

        chord_data.append({
            'color': maybe_brighten(log2fc_colormap, inter_row['MAST_log2FC_ligand'], lig, rec),
            'logfc': inter_row['MAST_log2FC_ligand'],
            'highlighted': highlight_chord,
            'source_celltype': source_celltype,
            'target_celltype': target_celltype,
            'ligand': lig,
            'receptor': rec,
            'value_text': f"<h4>Source: {source_celltype}<br>Ligand: {lig}<br>Target: {target_celltype}<br>Receptor: {rec}<br>Ligand log2FC: {inter_row['MAST_log2FC_ligand']:.2f}<br>Interactions: {inter_row['numSigI1']}</h4>",
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
    chord_data = sorted(chord_data, key=lambda c: (c['highlighted'], c['logfc']), reverse=False)

    ring_width = 15

    legend_group = make_circos_legend(min_receptor_numSigI1, max_receptor_numSigI1,
                                      cin_only,
                                      log2fc_colormap, da_colormap, dc1_colormap)

    return dbc.Row([
        dbc.Col(dbc.Container(
            dashbio.Circos(
                enableDownloadSVG=False,
                enableZoomPan=True,
                layout=layout,
                selectEvent={
                    "0": "both",
                    "1": "both",
                    "2": "both",
                    "3": "both",
                    "4": "both"
                },
                tracks=[{
                    'type': 'TEXT',
                    'data': list(text_data.values()),
                    'config': {
                        'innerRadius': 17.5*ring_width,
                        'outerRadius': 20*ring_width,
                        'style': {
                            'font-size': 7
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
                        'radius': 17*ring_width,
                        'tooltipContent': {
                            'name': 'value_text'
                        }
                    }
                }, {
                    'type': 'HISTOGRAM',
                    'data': numSigI1_data,
                    'config': {
                        'color': {
                            'name': 'color',
                        },
                        'innerRadius': 19*ring_width,
                        'outerRadius': 22*ring_width,
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
                        'innerRadius': 22*ring_width,
                        'outerRadius': 23*ring_width,
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
                        'innerRadius': 23*ring_width,
                        'outerRadius': 24*ring_width,
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
                    "innerRadius": 24*ring_width,
                    "outerRadius": 25*ring_width
                    #"cornerRadius": 4,
                },
                size=800,
                style={
                    'viewBox': '400 400 800 800',
                    'display': 'block',
                    'width': '100%',
                    'height': '100%'
                })
        ), width=9), dbc.Col(
                dcc.Graph(
                    id='circos-legend',
                    figure=legend_group,
                    config={
                        'displaylogo': False,
                        'showTips': False,
                        'displayModeBar': False,
                        'scrollZoom': False,
                        'responsive': False,
                        'showLink': False,
                        'watermark': False
                    }
                ), className='circosLegends'
        )], className='g-0')


def make_circos_legend(min_numSigI1, max_numSigI1,
                       cin_only: bool,
                       logfc_transformer: ColorTransformer,
                       da_transformer: ColorTransformer,
                       dc1_transformer: ColorTransformer) -> go.Figure:
    line_size_min_legend = go.Scatter(
        name="Few",
        x=[None], y=[None],
        mode='lines',
        showlegend=True,
        legendgroup='numSigI1',
        legendgrouptitle=dict(text='Receptor Interactions'),
        line=dict(
            color='black',
            width=1,
            dash='solid'
        )
    )

    line_size_max_legend = go.Scatter(
        name="Many",
        x=[None], y=[None],
        mode='lines',
        showlegend=True,
        legendgroup='numSigI1',
        legendgrouptitle=dict(text='Receptor<br>Interactions'),
        line=dict(
            color='black',
            width=6,
            dash='solid'
        )
    )

    logfc_legend = go.Scatter(
        x=[None], y=[None],
        text=[],
        mode='markers',
        hoverinfo='text',
        showlegend=False,
        marker=dict(
            opacity=0,
            cmin=logfc_transformer.min,
            cmax=logfc_transformer.max,
            color=[],
            colorbar=dict(
                title="Ligand<br>LogFC",
                thickness=25,
                titleside='right',
                bgcolor='rgba(0,0,0,0)',
                len=0.75,
                lenmode='fraction',
                xanchor='right',
                yanchor='top',
                x=1.25,
                y=0,
                tickmode='array',
                tickvals=[logfc_transformer.min, logfc_transformer.max],
                ticktext=["CIN<sup>low</sup>", "CIN<sup>high</sup>" + ("" if cin_only else "/<br>STING<sup>KO</sup>")],
            ),
            colorscale=logfc_transformer.make_plotly_colorscale()
        )
    )

    da_legend = go.Scatter(
        x=[None], y=[None],
        text=[],
        mode='markers',
        hoverinfo='text',
        showlegend=False,
        marker=dict(
            opacity=0,
            cmin=da_transformer.min,
            cmax=da_transformer.max,
            color=[],
            colorbar=dict(
                title="Differential<br>Abundance",
                thickness=25,
                titleside='right',
                bgcolor='rgba(0,0,0,0)',
                len=0.75,
                lenmode='fraction',
                xanchor='right',
                yanchor='top',
                x=1.25,
                y=1.5,
                tickmode='array',
                tickvals=[da_transformer.min, da_transformer.max],
                ticktext=["CIN<sup>low</sup>", "CIN<sup>high</sup>"]
            ),
            colorscale=da_transformer.make_plotly_colorscale()
        )
    )

    dc1_legend = go.Scatter(
        x=[None], y=[None],
        text=[],
        mode='markers',
        hoverinfo='text',
        showlegend=False,
        marker=dict(
            opacity=0,
            cmin=dc1_transformer.min,
            cmax=dc1_transformer.max,
            color=[],
            colorbar=dict(
                title="Diffusion<br>Component 1",
                thickness=25,
                titleside='right',
                bgcolor='rgba(0,0,0,0)',
                len=0.75,
                lenmode='fraction',
                xanchor='right',
                yanchor='top',
                x=1.25,
                y=3,
                tickmode='array',
                tickvals=[],
                ticktext=[]
            ),
            colorscale=dc1_transformer.make_plotly_colorscale()
        )
    )

    fig = go.Figure(
        data=[line_size_max_legend, line_size_min_legend, dc1_legend, da_legend, logfc_legend],
        layout=go.Layout(
            title=dict(text=''),
            titlefont_size=14,
            showlegend=True,
            hovermode='closest',
            autosize=False,
            width=250,
            height=800,  # Match circos
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor="rgba(0,0,0,0)",
            legend=dict(
                xanchor='right',
                yanchor='bottom',
                x=-1.6,
                y=-.5
            ),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, scaleratio=1, scaleanchor='x')
        )
    )
    return fig
