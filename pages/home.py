import base64
import os

import dash
import dash_bootstrap_components as dbc
from dash import html, callback, Output, Input, dcc
from dash.exceptions import PreventUpdate

from data import config
from viz.web import jumbotron, wrap_icon

dash.register_page(__name__,
                   path='/',
                   name='Home',
                   order=0)


layout = [
    jumbotron(
        "ContactTracing Interactive Visualizer",
        "This is a web application that allows for ContactTracing outputs to easily interpreted.",
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
    dcc.Store(id='custom-data-selection', storage_type='memory')
]

if not config.ALLOW_UPLOADS:
    layout.append(dbc.Tooltip("File uploads are disabled in the configuration file.", target='upload-button-container'))


def transition_to_data_summary(file_path: str, file_name: str, custom=False):
    # TODO: Display basic data summary
    # TODO: Set state to point to this file for analyses
    # Redirect to data summary page

    data = {'filename': file_name, 'path': file_path, 'custom': custom}

    summary_content = file_name

    return data, summary_content


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
                'height': '120px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            accept='.h5ad',
            multiple=False
        ),
        html.Div(id='file-upload-message')], fluid=True),
    ]


@callback(
    Output('data-session', 'data'),
    Input('builtin-data-selection', 'data'),
    Input('custom-data-selection', 'data')
)
def set_selected_data(builtin_data, custom_data):
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
)
def update_data_selected(mouse_cin, human_cin, upload_button):
    clicked = dash.ctx.triggered_id
    if clicked == 'mouse-cin':
        return transition_to_data_summary('data/mouse_degboth.h5ad', 'Mouse CIN')
    elif clicked == 'human-cin':
        return transition_to_data_summary('data/human_degboth.h5ad', 'Human CIN')
    elif clicked == 'upload-button':
        if config.ALLOW_UPLOADS:
            return None, transition_to_custom_upload()
        else:
            return None, dbc.Alert(wrap_icon('fa-triangle-exclamation', "File uploads are disabled on this server"), color='warning')
    else:
        raise PreventUpdate


@callback(
    Output('custom-data-selection', 'data'),
    Output(component_id='file-upload-message', component_property='children'),
    Input(component_id='upload-data', component_property='contents'),
    Input(component_id='upload-data', component_property='filename'),
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

    return transition_to_data_summary(contents, filename, custom=True)

