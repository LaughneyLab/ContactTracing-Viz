import os
import shutil
import sys
from uuid import uuid4

import dash
from dash import html, dcc, Output, Input, State, Dash
import dash_bootstrap_components as dbc

from viz.config import LONG_CALLBACK_EXPIRY

# Workaround for dill crashes: Ref: https://github.com/uqfoundation/dill/issues/332#issuecomment-908826972
try:
    import dill
    import abc
    @dill.register(abc.ABCMeta)
    def save_abc(pickler, obj):
        dill.StockPickler.save_type(pickler, obj)
except:
    pass


# Workaround for recursion error in dbc by removing the deprecated wrappers
# Ref: https://github.com/facultyai/dash-bootstrap-components/issues/892
try:
    def setstate_wrapper(self, state):
        # ensure deprecated & wrapped fields are set to avoid recursive stack explosion in __getattr__
        self.deprecated = state.get("deprecated", None)
        self.wrapped = state.get("wrapped", None)

    setattr(dbc._V1DeprecationWarningWrapper, '__setstate__', setstate_wrapper)
except:
    pass


# Each server instance gets a unique cache
launch_uuid = uuid4()

callback_manager = None
if 'REDIS_URL' in os.environ:
    try:
        from dash import CeleryManager
        from celery import Celery
        celery_app = Celery(__name__, broker=os.environ['REDIS_URL'] + "/0", backend=os.environ['REDIS_URL'] + "/1")
        callback_manager = CeleryManager(celery_app, cache_by=[lambda: launch_uuid], expire=LONG_CALLBACK_EXPIRY)
    except:
        callback_manager = None
if callback_manager is None:
    celery_app = None
    print("WARNING: Celery not available, falling back to diskcache", file=sys.stderr)
    from dash import DiskcacheManager
    from diskcache import Cache
    # Don't remove unless its the main process
    if __name__ == "__main__" and os.path.exists('./ct_viz_cache'):
        shutil.rmtree('./ct_viz_cache')
    callback_manager = DiskcacheManager(Cache('./ct_viz_cache'), cache_by=[lambda: launch_uuid], expire=LONG_CALLBACK_EXPIRY)


# Custom style made with https://bootstrap.build/ based on FLATLY
STYLESHEET = "bootstrap.min.css"

debug = len(sys.argv) > 1 and sys.argv[1] == 'debug'

app = Dash(__name__,
           suppress_callback_exceptions=True,
           compress=True,
           update_title="Please wait...",
           meta_tags=[
                {'name': 'description', 'content': 'ContactTracing, a new method to interrogate cell-cell interactions and their biological effects from single-cell data.'},
                {'name': 'viewport', 'content': 'width=device-width, initial-scale=1'},
                {'name': 'robots', 'content': 'index,follow'},
                {'property': 'og:title', 'content': 'ContactTracing Visualizer'},
                {'property': 'og:type', 'content': 'website'},
                {'itemprop': 'name', 'content': 'ContactTracing Visualizer'},
                {'itemprop': 'description', 'content': 'ContactTracing, a new method to interrogate cell-cell interactions and their biological effects from single-cell data.'},
                {'property': 'og:description', 'content': 'ContactTracing, a new method to interrogate cell-cell interactions and their biological effects from single-cell data.'},
                {'itemprop': 'image', 'content': '/assets/ct_cover.jpg'},
                {'property': 'og:image', 'content': '/assets/ct_cover.jpg'},
               {'property': 'og:url', 'content': 'http://contacttracing.laughneylab.com/'},
           ],
           title='ContactTracing',
           external_scripts=['https://sharonchoong.github.io/svg-exportJS/svg-export.min.js',
                             "/assets/circos_hooks.js", '/assets/plotly_hooks.js'],
           external_stylesheets=[
               # Fonts
               dict(rel="preconnect",
                    href="https://fonts.googleapis.com"),
               dict(rel="preconnect",
                    href="https://fonts.gstatic.com",
                    crossorigin=""),
               dict(rel="stylesheet",
                    href="https://fonts.googleapis.com/css2?family=Courier+Prime:wght@400;700&family=IBM+Plex+Mono:wght@400;700&family=Noto+Sans:wght@400;600&family=Source+Code+Pro:wght@400;700&display=swap"),
               # Stylesheets
               dbc.themes.FLATLY,
               dbc.icons.BOOTSTRAP,
               dbc.icons.FONT_AWESOME,
               STYLESHEET,
               "custom.css"  # Custom overrides
           ],
           use_pages=True,
           background_callback_manager=callback_manager,
           show_undo_redo=debug
           )
server = app.server


def access_code_page():
    #page_container = dash.page_container
    CODE = 'cintme123'

    access_modal = dbc.Modal([
        dbc.ModalHeader("Enter Access Code", close_button=False),
        dbc.ModalBody([
            dbc.Input(id='access-code', type='password', placeholder='Enter access code', className='mb-3', autofocus=True),
        ]),
        dbc.ModalFooter(dbc.Button("Enter", id="access-modal-close", className="ms-auto", n_clicks=0))
    ], is_open=True, size='xl', backdrop='static', keyboard=False, centered=True)

    container = dbc.Container(fluid=True, children=[
        html.Div(id='access-wrapper', children=access_modal)
    ])

    @dash.callback(
        Output('access-wrapper', 'children'),
        Output('access-code-entered', 'data'),
        Input('access-modal-close', 'n_clicks'),
        Input('access-code', 'n_submit'),
        State('access-code', 'value'),
        State('access-code-entered', 'data'),
    )
    def access_modal_close(n_clicks, n_submit, code, was_entered):
        if (n_clicks is not None and n_clicks < 1) and (n_submit is not None and n_submit < 1) and not was_entered:
            return access_modal, False
        if code == CODE or was_entered:
            return None, True
        return access_modal, False

    return container


layout = [
    dbc.Navbar(
        dbc.Container(
            [
                html.A(
                    dbc.Row([
                        dbc.Col(html.Img(src=app.get_asset_url('ct_logo.png'), height="30px"), width='auto'),
                        dbc.Col(dbc.NavbarBrand(html.Strong('ContactTracing', style={'fontSize': '1.75rem'}))),
                    ], align='center', className="g-0"),
                    href="/",
                    style={'textDecoration': 'none'}  # Hide hyperlink underline
                ),
                dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
                dbc.Collapse(dbc.Nav(
                    [dbc.NavItem(dbc.NavLink(page['name'], href=page['relative_path'])) for page in dash.page_registry.values()] +
                    [dbc.NavItem(dbc.NavLink('GitHub', href='https://github.com/LaughneyLab/ContactTracing_tutorial'))],
                    navbar=True,
                    pills=True
                ), navbar=True, id='navbar-collapse', className='justify-content-end'),
            ]
        ),
        color='light',
        fixed='top',
        dark=False,
        light=True,
        className='mb-5',
    ),
    # dcc.Store(id='access-code-entered', storage_type='local', data='debug' in sys.argv),  # Temporary session storage
    # dcc.Store(id='is_mobile_device', storage_type='memory', data=False),  # Temporary session storage
    html.Br(), html.Br(),
    html.Div(id='page-content', children=dbc.Container([
        dash.page_container,
        # access_code_page() if 'debug' not in sys.argv else None
    ], fluid=True))
]


app.layout = dbc.Container(fluid=True,
                           children=layout)


# app.clientside_callback("""
# function() {
#     return (navigator.maxTouchPoints || 'ontouchstart' in document.documentElement);
# }
# """, Output('is_mobile_device', 'data'))


@app.callback(
    Output("navbar-collapse", "is_open"),
    [Input("navbar-toggler", "n_clicks")],
    [State("navbar-collapse", "is_open")],
)
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


if __name__ == '__main__':
    # debug = False
    app.run(port=8000,
            debug=debug,
            dev_tools_ui=debug,
            dev_tools_props_check=debug,
            dev_tools_serve_dev_bundles=debug,
            dev_tools_hot_reload=debug,
            dev_tools_prune_errors=False,
            use_reloader=False)
