import os
import shutil
import sys

import dash
from dash import Dash, html, dcc, Output, Input
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate

from data.config import LONG_CALLBACK_EXPIRY
from viz.data import compile_data
from viz.web import wrap_icon


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


# Compile our data for optimized queries
compile_data()


callback_manager = None
if 'REDIS_URL' in os.environ:
    try:
        from dash import CeleryManager
        from celery import Celery
        celery_app = Celery(__name__, broker=os.environ['REDIS_URL'] + "/0", backend=os.environ['REDIS_URL'] + "/1")
        callback_manager = CeleryManager(celery_app, expire=LONG_CALLBACK_EXPIRY)
    except:
        callback_manager = None
if callback_manager is None:
    print("WARNING: Celery not available, falling back to diskcache", file=sys.stderr)
    from dash import DiskcacheManager
    from diskcache import Cache
    # Don't remove unless its the main process
    if __name__ == "__main__" and os.path.exists('./ct_viz_cache'):
        shutil.rmtree('./ct_viz_cache')
    callback_manager = DiskcacheManager(Cache('./ct_viz_cache'), expire=LONG_CALLBACK_EXPIRY)

"""
font-family: 'Courier Prime', monospace;
font-family: 'IBM Plex Mono', monospace;
font-family: 'Noto Sans', sans-serif;
font-family: 'Source Code Pro', monospace;
"""


# Custom style made with https://bootstrap.build/ based on FLATLY
STYLESHEET = "bootstrap.min.css"

app = Dash(__name__,
           suppress_callback_exceptions=True,
           compress=True,
           meta_tags=[
                {'name': 'description', 'content': 'An interactive dashboard to visualize transcriptional responses detected by ContactTracing.'},
                {'name': 'viewport', 'content': 'width=device-width, initial-scale=1'},
                {'name': 'robots', 'content': 'index,follow'},
           ],
           title='ContactTracing',
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
           ],
           use_pages=True,
           background_callback_manager=callback_manager
           )
server = app.server


layout = html.Div([
    dbc.Navbar(
        dbc.Container(
            [
                html.A(
                    dbc.Row([
                        dbc.Col(html.Img(src=app.get_asset_url('ct_logo.png'), height="30px"), width='auto'),
                        dbc.Col(dbc.NavbarBrand(html.Strong('ContactTracing'))),
                    ], align='center', className="g-0"),
                    href="/",
                    style={'textDecoration': 'none'}  # Hide hyperlink underline
                ),
                dbc.Nav(
                    [dbc.NavItem(dbc.NavLink(page['name'], href=page['relative_path'])) for page in dash.page_registry.values()],
                    navbar=True,
                    pills=True
                )
            ]
        ),
        color='light',
        fixed='top',
        dark=False,
        light=True,
        className='mb-5',
    ),
    dcc.Store(id='data-session', storage_type='session'),  # Temporary session storage
    html.Div(id='page-content', children=dbc.Container(dash.page_container, fluid=True, class_name='mt-5 pt-5'))
])


app.layout = dbc.Container(fluid=True,
                           children=layout)


if __name__ == '__main__':
    debug = sys.argv[1] == 'debug'
   # debug = False
    app.run(port=8000,
            debug=debug,
            dev_tools_ui=debug,
            dev_tools_props_check=debug,
            dev_tools_serve_dev_bundles=debug,
            dev_tools_hot_reload=debug,
            dev_tools_prune_errors=False,
            use_reloader=False)
