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
    return dbc.Container([
            dbc.Card([
                dbc.CardHeader(
                    dbc.Tabs([
                        dbc.Tab(label='Welcome', tab_id='default', labelClassName="text-muted"),
                        dbc.Tab(label='Dataset Information', tab_id='dataset', labelClassName="text-muted"),
                        dbc.Tab(label='Visualization Information', tab_id='plot', labelClassName="text-muted"),
                        dbc.Tab(label='Miscellaneous', tab_id='misc', labelClassName="text-muted"),
                    ], id='home-tabs', active_tab='default', persistence=False)
                ),
                dbc.CardBody([
                    html.H4(id='home-tab-title', className='card-title'),
                    html.P(id='home-tab-content', className='card-text'),
                ])
            ], color="information", inverse=False),
        ])


DEFAULT_TAB = [
    # TODO: Images of example plots
    html.P("By leveraging the intrinsic variability present in the Tumor Microenvironment (TME), we can identify "
           "transcriptional response profiles that are associated with particularly experimental conditions. This web "
           "tool produces visualizations of the ContactTracing toolkit on a mouse model of chromosomally unstable "
           "cancer as described in our associated publication."),
    html.P(["This interactive web tool for exploring the ContactTracing-generated data is built using the ",
            html.A("Plotly Dash framework", href="https://dash.plotly.com/"),
            ". For additional information, please navigate the tabs above."]),
    html.P([
        html.H5("Included Visualizations"),
        make_data_redirect_buttons()
    ]),
    html.P([html.H5("Citation"), dbc.Card("TBA", body=True)])
]
DATASET_DESCRIPTIONS = [
    html.P("The included dataset is a mouse model of chromosomally unstable cancer (CIN) as described in our associated "
           "publication..."),
]
PLOT_DESCRIPTIONS = [
    html.P("The ContactTracing toolkit produces a variety of plots that are used to interpret the results of the "
              "transcriptional profiling. These plots are described below."),
    html.H5(html.A("Circos Plot", href="/circos")),
    html.P("...."),
    html.H5(html.A("Cell Type Interactions Plot", href="/interactions")),
    html.P("...."),
    html.H5(html.A("Downstream Ligand Effects Plot", href="/ligand-effects")),
    html.P("....")
]
MISC_INFO = [
    html.P("Additional information about ContactTracing, the data presented, and contacts can be found below."),
    html.H5("Contact Information")
]


layout = [
    jumbotron(
        "ContactTracing Interactive Visualizer",
        "This is a web application that allows for ContactTracing outputs to be easily interpreted.",
        dbc.Row([
            dbc.Container(make_data_summary('builtin', 'Mouse CIN'), id='main-content-container', fluid=True)
        ])
    ),
]


@callback(
    Output('home-tab-title', 'children'),
    Output('home-tab-content', 'children'),
    Input('home-tabs', 'active_tab')
)
def update_tab_content(active_tab):
    if active_tab == 'default':
        return 'ContactTracing is a tool for profiling the TME', DEFAULT_TAB
    elif active_tab == 'dataset':
        return 'Included datasets for visualization', DATASET_DESCRIPTIONS
    elif active_tab == 'plot':
        return 'Interpreting ContactTracing plots', PLOT_DESCRIPTIONS
    elif active_tab == 'misc':
        return 'Additional Information', MISC_INFO
    else:
        raise PreventUpdate
