import base64
import os

import dash
import dash_bootstrap_components as dbc
import humanize
from dash import html, callback, Output, Input, dcc, State
from dash.exceptions import PreventUpdate

from data import config
from viz.data import prune_ct_data
from viz.web import jumbotron, wrap_icon, make_data_redirect_buttons

dash.register_page(__name__,
                   path='/',
                   name='Home',
                   order=0)


layout = [
    jumbotron(
        "ContactTracing Interactive Visualizer",
        "This is a web application that allows for ContactTracing outputs to be easily interpreted.",
        "Please select a pre-built results file or upload your own below:",
        dbc.Row([
            dbc.Col(dbc.DropdownMenu(
                label='Included Experiment Files',
                children=[
                    dbc.DropdownMenuItem(wrap_icon('fa-microscope', 'Mouse CIN', right=True),
                                         id='mouse-cin', n_clicks=0),
                    dbc.DropdownMenuItem(wrap_icon('fa-person', 'Human CIN', right=True),
                                         id='human-cin', n_clicks=0)
                ]
            )),
            dbc.Col(html.Div(dbc.Button(wrap_icon('fa-file-import', 'Upload File (.h5ad)'), id='upload-button', n_clicks=0, disabled=not config.ALLOW_UPLOADS), id='upload-button-container')),
        ]),
        dbc.Row([
            dbc.Container("Nothing Selected", id='main-content-container', fluid=True)
        ])
    ),
    dcc.Store(id='builtin-data-selection', storage_type='memory'),
    dcc.Store(id='custom-data-selection', storage_type='memory'),
]

if not config.ALLOW_UPLOADS:
    layout.append(dbc.Tooltip("File uploads are disabled in the configuration file.", target='upload-button-container'))


def make_data_summary(file_path: str, file_name: str, custom=False):
    # Hardcoded descriptions if not custom
    descriptions = []

    if not custom:
        if file_name == 'Mouse CIN':
            descriptions.append(html.P('Based on scRNA-seq data from a mouse model of CIN.', className="card-text"))
        elif file_name == 'Human CIN':
            descriptions.append(html.P('Human CIN FIXME', className="card-text"))

    buttons = make_data_redirect_buttons()

    return dbc.Container([
            dbc.Card([
                dbc.CardHeader(wrap_icon('fa-database', "Selected File")),
                dbc.CardBody([
                    html.H4(f"Selected ContactTracing Data: {file_name}", className="card-title"),
                    html.P("Basic Data Summary:", className="card-text"),
                    *descriptions,
                    dbc.Table([
                        html.Thead(html.Tr([html.Th("Attribute"), html.Th("Value")])),
                        html.Tbody([
                            html.Tr([html.Td("File Name"), html.Td(file_name)]),
                            html.Tr([html.Td("File Size"), html.Td(humanize.naturalsize(os.path.getsize(file_path)))]),
                            html.Tr([html.Td("User Uploaded?"), html.Td("Yes" if custom else "No")]),
                        ])
                    ], bordered=True, responsive=True),
                    buttons,
                ])
            ], color="information", inverse=False),
        ])


def transition_to_data_summary(file_path: str, file_name: str, custom=False):
    # Redirect to data summary page

    data = {'filename': file_name, 'path': file_path, 'custom': custom}

    # FIXME: Hack to include results tsvs while they are not present in the h5ad
    if file_name == 'Mouse CIN':
        data['tsv'] = "data/mouse_ranked_interactions.tsv"
        data['obs'] = "data/CIN_and_STING_degobs.txt"
    elif file_name == 'Human CIN':
        data['tsv'] = "data/cin_ranked_interactions.tsv"
        data['obs'] = "data/NA.txt"

    assert os.path.exists(file_path), f"File {file_path} does not exist"

    return data, make_data_summary(file_path, file_name, custom)


def transition_to_custom_upload():
    return [
        dbc.Container([dcc.Upload(
            id='upload-data',
            children=html.Div(
                wrap_icon('fa-upload',
                          html.Br(),
                          ' Drag and Drop or ',
                          html.A('Select Files')
            )),
            style={
                'width': '50%',
                'height': '7.5em',
                'lineHeight': '3.75em',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '.6em'
            },
            accept='.h5ad',
            multiple=False
        ),
        html.Div(id='file-upload-message')], fluid=True)
    ]


# Wire callback for the reset button
@callback(
    Output('reset-data-indicator', 'data'),
    State('reset-data-indicator', 'data'),
    Input('unselect-button', 'n_clicks'),
    prevent_initial_call=True
)
def reset_data(reset_data, unselect_button):
    if dash.ctx.triggered_id != 'unselect-button':
        raise PreventUpdate

    return 0 if not reset_data else reset_data + 1


@callback(
    Output('data-session', 'data'),
    Input('builtin-data-selection', 'data'),
    Input('custom-data-selection', 'data'),
    Input('reset-data-indicator', 'data'),
    prevent_initial_call=True
)
def set_selected_data(builtin_data, custom_data, reset_data):
    triggered = dash.ctx.triggered_id
    if triggered == 'builtin-data-selection':
        return builtin_data
    elif triggered == 'custom-data-selection':
        return custom_data
    else:
        return None


@callback(
    Output('builtin-data-selection', 'data'),
    Output(component_id='main-content-container', component_property='children'),
    Input(component_id='mouse-cin', component_property='n_clicks'),
    Input(component_id='human-cin', component_property='n_clicks'),
    Input(component_id='upload-button', component_property='n_clicks'),
    Input(component_id='reset-data-indicator', component_property='data'),
    prevent_initial_call=True
)
def update_data_selected(mouse_cin, human_cin, upload_button, reset_data):
    clicked = dash.ctx.triggered_id
    if clicked == 'mouse-cin':
        return transition_to_data_summary('data/mouse_deg_fixed_pruned.h5ad', 'Mouse CIN')
    elif clicked == 'human-cin':
        return transition_to_data_summary('data/human_degboth_pruned.h5ad', 'Human CIN')
    elif clicked == 'upload-button':
        if config.ALLOW_UPLOADS:
            return None, transition_to_custom_upload()
        else:
            return None, dbc.Alert(wrap_icon('fa-triangle-exclamation', "File uploads are disabled on this server"), color='warning')
    elif clicked == 'reset-data-indicator':
        return None, []
    else:
        raise PreventUpdate


@callback(
    Output('custom-data-selection', 'data'),
    Output(component_id='file-upload-message', component_property='children'),
    Input(component_id='upload-data', component_property='contents'),
    Input(component_id='upload-data', component_property='filename'),
    background=True
)
def update_file_upload(contents, filename):
    if contents is None:
        raise dash.exceptions.PreventUpdate

    if filename is None:
        filename = "uploaded_file.h5ad"

    if not filename.endswith('.h5ad') or '/' in filename or '\\' in filename:
        return None, dbc.Alert(wrap_icon('fa-triangle-exclamation', "File must be of type .h5ad"), color="danger")

    file_path = os.path.join('data', filename)

    if os.path.exists(file_path):
        return None, dbc.Alert(wrap_icon('fa-triangle-exclamation', "File already exists"), color="danger")

    with open(file_path, 'wb') as f:
        f.write(base64.decodebytes(contents))

    pruned_file_path = prune_ct_data(file_path)
    os.remove(file_path)

    return transition_to_data_summary(pruned_file_path, filename, custom=True)

