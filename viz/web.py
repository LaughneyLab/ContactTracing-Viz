from typing import List, Tuple, Optional
from uuid import uuid4

import dash
import dash_bootstrap_components as dbc
import numpy as np
import plotly.graph_objects as go
import pandas as pd
from dash import html, Output, Input, dcc, callback, State
import dash_bio as dashbio
from dash.exceptions import PreventUpdate

from viz.util import ColorTransformer, celltype_to_colors, saturate_color, smooth_step, \
    brighten_color, brighten_and_saturate_color


def make_tooltip(content, tooltip_content):
    # uuid4 is used to generate a unique id for the tooltip
    tooltip_id = str(uuid4())
    content_span = html.Span(content,
                             id=tooltip_id,
                             style={'display': 'inline', 'text-decoration': 'underline dotted', 'cursor': 'help'})
    tooltip = dbc.Tooltip(tooltip_content,
                          target=tooltip_id,
                          placement='auto')

    return content_span, tooltip


def make_fdr_slider(id: str, value) -> html.Div:
    if isinstance(value, str):
        value = float("." + value[3:])
    return make_custom_slider(
        id, 0.01, 0.25, value, 0.01,
        lambda x: "fdr" + f"{x:.2f}".split(".")[1]
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
        className='form-range'
    )
    input = dbc.Input(id=id + '-input',
                      type='number',
                      min=min,
                      max=max,
                      step=step,
                      value=value)
    value_store = dcc.Store(id=id, data=transform(value) if transform else value)
    slider_element = html.Div(dbc.Container([
        value_store,
        dbc.Row([
            dbc.Col(slider, align="end"),
            dbc.Col(input, width=2, align="start", class_name='slider_text_input')
        ])
    ], fluid=True))

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
            return transform(slider_store) if transform else slider_store, slider_store, slider_store

        return transform(value) if transform else value, value, value

    return slider_element


def jumbotron(title, main_content, sub_content, *additional_content, dark=False):
    return html.Div(dbc.Container([
        html.H1(title, className='display-3'),
        html.P(main_content, className='lead'),
        html.Hr(className='my-2'),
        html.P(sub_content),
        *[html.P(content) for content in additional_content]
    ], fluid=True, className='py-3'), className="p-3 rounded-3 " + ('dark text-white' if dark else 'light'))


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


def control_panel(submit_btn_id: str, *element_rows: List[dbc.Card]) -> html.Div:
    rows = [dbc.Row(dbc.CardGroup(cols)) for cols in element_rows]
    # Add all rows to a container that will be within an Accordion except for the
    # last row (which contains the submit button!)

    options_accordion = dbc.Accordion(
        dbc.AccordionItem(dbc.Container(rows[:-1], fluid=True), title="Figure Options", item_id='options'),
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
        options_accordion,
        rows[-1]
    ],
        className='mb-3'
    )


def figure_output(title, footer, element, help_info, outline=True) -> html.Div:
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
                dbc.Col(title),
                dbc.Col(help_button, width=1)
            ])),
            dbc.CardBody([
                html.P([
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
        html.P(subtitle, className='lead'),
        html.Hr(className='my-2'),
        html.Br(),
        html.Div(dbc.Container(content), className='card-text')
    ], fluid=True, className='py-3'), className="p-3 rounded-3 " + ('dark text-white' if dark else 'light'))


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
        dbc.Col(dbc.Button(wrap_icon('fa-circle-dot', 'View Circos'), id='circos-button', href='/circos'), width="auto"),
        dbc.Col(dbc.Button(wrap_icon('fa-arrows-left-right-to-line', 'View Cell Type Interactions'), id='ct-interaction-button', href='/interactions'), width="auto"),
        dbc.Col(dbc.Button(wrap_icon('fa-maximize', 'View Ligand Effects'), id='ligand-effect-button', href='/ligand-effects'), width="auto"),
    ])


def make_circos_figure(set_progress, progress_offset: int,
                       outer_data: pd.DataFrame,
                       inter_data: pd.DataFrame,
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
            'celltype': celltype
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

    # Build next ring for DC1 heatmap
    dc1_colormap = ColorTransformer(-1, 1, 'cividis')
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
                'value_text': f"<h3>{celltype}<br>{t} DC1: {dc1:.2f}</h3>",
                'color': dc1_colormap(dc1),
                'target': t
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
            da = outer_data[(outer_data['cell_type'] == celltype) & (outer_data['target'] == t)]['DA_score'].values[0]
            da_data.append({
                'block_id': id,
                'start': i,
                'end': (i + 1),
                'value': da,
                'value_text': f"<h3>{celltype}<br>{t} DA: {da:.2f}</h3>",
                'color': da_colormap(da),
                'target': t
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
                    'target': t
                })
            else:
                lig_effect = receptor_info[(receptor_info['receptor'] == t) & (receptor_info['cell_type_receptor'] == celltype)]['numSigI1'].values[0]
                # Normalize to 0 minimum
                lig_effect -= min_receptor_numSigI1
                numSigI1_data.append({
                    'block_id': id,
                    'start': i,
                    'end': (i + 1),
                    'value': max(lig_effect, 0),
                    'value_text': f"<h3>{celltype}<br>{t} numSigI: {lig_effect:d}</h3>",
                    'color': color,
                    'target': t
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

        text_data[(source_celltype, lig)] = {
            'block_id': celltype2id[source_celltype],
            'position': source_position+.5,
            'value': lig if highlight_chord else ""
        }
        text_data[(target_celltype, rec)] = {
            'block_id': celltype2id[target_celltype],
            'position': target_position+.5,
            'value': rec if highlight_chord else ""
        }

        chord_data.append({
            'color': maybe_brighten(log2fc_colormap, inter_row['MAST_log2FC_ligand'], lig, rec),
            'logfc': inter_row['MAST_log2FC_ligand'],
            'highlighted': highlight_chord,
            'source': {
                'id': celltype2id[source_celltype],
                'start': source_position - thickness,
                'end': source_position + thickness,
                'info': f"{lig}+ {source_celltype}",
            },
            'target': {
                'id': celltype2id[target_celltype],
                'start': target_position - thickness,
                'end': target_position + thickness,
                'info': f"{rec}+ {target_celltype}" if highlight_chord else None,
                'value_text': f"<br>log2FC: {inter_row['MAST_log2FC_ligand']:.2f}<br>Interactions: {inter_row['numSigI1']}",
            },
        })
    # Sort to place red chords on top
    chord_data = sorted(chord_data, key=lambda c: (c['highlighted'], c['logfc']), reverse=False)

    ring_width = 15

    legend_group = make_circos_legend(min_receptor_numSigI1, max_receptor_numSigI1,
                                      log2fc_colormap, da_colormap, dc1_colormap)

    return dbc.Row([
        dbc.Col(dbc.Container(html.Div(
            dashbio.Circos(
                enableDownloadSVG=True,
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
                            'source': 'source',
                            'sourceID': 'info',
                            'target': 'target',
                            'targetID': 'info',
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
                    "outerRadius": 25*ring_width,
                    'size': 800
                    #"cornerRadius": 4,
                }
            )
        , style={'height': 800, 'width': 800})), width=9), dbc.Col(
                dcc.Graph(
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
                )
        )], className='g-0')


def make_circos_legend(min_numSigI1, max_numSigI1,
                       logfc_transformer: ColorTransformer,
                       da_transformer: ColorTransformer,
                       dc1_transformer: ColorTransformer) -> go.Figure:
    line_size_min_legend = go.Scatter(
        name=str(min_numSigI1),
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
        name=str(max_numSigI1),
        x=[None], y=[None],
        mode='lines',
        showlegend=True,
        legendgroup='numSigI1',
        legendgrouptitle=dict(text='Receptor Interactions'),
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
                title="Ligand LogFC",
                thickness=25,
                titleside='right',
                bgcolor='rgba(0,0,0,0)',
                len=0.75,
                lenmode='fraction',
                xanchor='right',
                yanchor='top',
                x=1.75,
                y=0
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
                title="Differential Abundance",
                thickness=25,
                titleside='right',
                bgcolor='rgba(0,0,0,0)',
                len=0.75,
                lenmode='fraction',
                xanchor='right',
                yanchor='top',
                x=1.75,
                y=1.5
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
                title="Diffusion Component 1",
                thickness=25,
                titleside='right',
                bgcolor='rgba(0,0,0,0)',
                len=0.75,
                lenmode='fraction',
                xanchor='right',
                yanchor='top',
                x=1.75,
                y=3
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
                x=-2,
                y=-.5
            ),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, scaleratio=1, scaleanchor='x')
        )
    )
    return fig
