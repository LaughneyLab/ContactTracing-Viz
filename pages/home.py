import os

import dash
import dash_bootstrap_components as dbc
from dash import html, callback, Output, Input, dcc, State
from dash.exceptions import PreventUpdate

from viz.web import jumbotron, wrap_icon, make_data_redirect_buttons

dash.register_page(__name__,
                   path='/',
                   name='Home',
                   order=0)


def make_data_summary(file_path: str, file_name: str, custom=False):
    # Hardcoded descriptions if not custom
    descriptions = []

    if not custom:
        if file_name == 'Mouse CIN':
            descriptions.append(html.P('Based on scRNA-seq data from a mouse model of CIN.', className="card-text"))
        elif file_name == 'Mouse STING':
            descriptions.append(html.P('Based on scRNA-seq data from a mouse model of STING expression.', className="card-text"))
        elif file_name == 'Mouse CIN & STING Intersection':
            descriptions.append(html.P('Based on interactions that are shared between CIN and STING comparisons.', className="card-text"))
        #elif file_name == 'Human CIN':
        #    descriptions.append(html.P('Human CIN FIXME', className="card-text"))

    buttons = make_data_redirect_buttons()

    return dbc.Container([
            dbc.Card([
                dbc.CardHeader(wrap_icon('fa-database', "Selected File")),
                dbc.CardBody([
                    html.H4(f"Selected ContactTracing Data: {file_name}", className="card-title"),
                    html.P("Basic Data Summary:", className="card-text"),
                    *descriptions,
                  #  dbc.Table([
                  #      html.Thead(html.Tr([html.Th("Attribute"), html.Th("Value")])),
                  #      html.Tbody([
                  #          html.Tr([html.Td("File Name"), html.Td(file_name)]),
                  #          html.Tr([html.Td("File Size"), html.Td(humanize.naturalsize(os.path.getsize(file_path)))]),
                  #          html.Tr([html.Td("User Uploaded?"), html.Td("Yes" if custom else "No")]),
                  #      ])
                  #  ], bordered=True, responsive=True),
                    buttons,
                ])
            ], color="information", inverse=False),
        ])


layout = [
    jumbotron(
        "ContactTracing Interactive Visualizer",
        "This is a web application that allows for ContactTracing outputs to be easily interpreted.",
        "Please select a pre-built results file:", # or upload your own below:",
        dbc.Row([
            dbc.Container(make_data_summary('builtin', 'Mouse CIN'), id='main-content-container', fluid=True)
        ])
    ),
]

