import sys

import dash
from dash import Dash, html, dcc, Output, Input
import dash_bootstrap_components as dbc

from viz.web import wrap_icon

app = Dash(__name__,
           suppress_callback_exceptions=True,
           compress=True,
           meta_tags=[
                {'name': 'description', 'content': 'An interactive dashboard to visualize transcriptional responses detected by ContactTracing.'},
                {'name': 'viewport', 'content': 'width=device-width, initial-scale=1'},
                {'name': 'robots', 'content': 'index,follow'},
           ],
           title='ContactTracing',
           external_stylesheets=[dbc.themes.MINTY, dbc.icons.FONT_AWESOME],
           use_pages=True)

layout = html.Div([
    dbc.Navbar(
        dbc.Container(
            [
                html.A(
                    dbc.Row([
                        dbc.Col(dbc.NavbarBrand(wrap_icon('fa-circle-nodes', html.Strong('ContactTracing')),)),
                    ], align='center', className="g-0"),
                    href="/",
                    style={'textDecoration': 'none'}  # Hide hyperlink underline
                ),
                dbc.Nav(
                    [dbc.NavItem(dbc.NavLink('', id='data-selection', active=True))] +
                    [dbc.NavItem(dbc.NavLink(page['name'], href=page['relative_path'])) for page in dash.page_registry.values()],
                    navbar=True,
                    pills=True
                )
            ]
        ),
        color='primary',
        fixed='top',
        dark=True,
        className='mb-5'
    ),
    dcc.Store(id='data-session', storage_type='session'),  # Temporary session storage
    html.Div(id='page-content', children=dbc.Container(dash.page_container, fluid=True, class_name='m-5 pt-5'))
])

app.layout = dbc.Container(fluid=True,
                           children=layout)


@app.callback(
    Output('data-selection', 'children'),
    Input('data-session', 'data')
)
def update_data_selection(data):
    prefix = 'Selected Data: '
    if data is None:
        return prefix + "Nothing"
    else:
        return prefix + data['filename']


if __name__ == '__main__':
    debug = sys.argv[1] == 'debug'
    app.run(port=8050,
            debug=debug,
            dev_tools_ui=debug,
            dev_tools_props_check=debug,
            dev_tools_serve_dev_bundles=debug,
            dev_tools_hot_reload=debug,
            dev_tools_hot_reload_interval=15)
