import dash
import dash_bootstrap_components as dbc
from dash_extensions.enrich import html, callback, Output, Input
from dash.exceptions import PreventUpdate

from viz.docs import home_welcome_info, home_dataset_descriptions, home_plot_descriptions, home_misc_info
from viz.web import jumbotron

dash.register_page(__name__,
                   path='/',
                   name='Home',
                   order=0)


layout = [
    jumbotron(
        "ContactTracing Interactive Visualizer",
        "This is a web application that allows for ContactTracing outputs to be easily interpreted.",
        dbc.Row([
            dbc.Container(dbc.Container([
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
                        html.Hr(),
                        html.P(id='home-tab-content', className='card-text'),
                    ])
                ], color="information", inverse=False),
            ]), id='main-content-container', fluid=True)
        ])
    ),
]


@callback(
    Output('contact-info', "children"),
    Input('contact-info-button', "n_clicks")
)
def show_contact_info(n_clicks):
    if n_clicks > 0:
        return dbc.Card([
            html.P(['Dr. Ashley Laughney: ',
                   html.A('Ashley.laughney@gmail.com', href='mailto:Ashley.laughney@gmail.com')]),
            html.P(["Austin Varela: ", html.A('aav4003@med.cornell.edu', href='mailto:aav4003@med.cornell.edu')])
        ], body=True)


@callback(
    Output("circos-plot-help-offcanvas", "is_open"),
    Output("cell-type-plot-help-offcanvas", "is_open"),
    Output("ligand-effect-plot-help-offcanvas", "is_open"),
    Input("full-circos-help-button", "n_clicks"),
    Input("full-cell-type-help-button", "n_clicks"),
    Input("full-ligand-effect-help-button", "n_clicks"),
)
def show_offcanvas_content(circos_n_clicks, celltype_n_clicks, ligand_effect_n_clicks):
    triggered = dash.callback_context.triggered_id
    if triggered == "full-circos-help-button":
        return circos_n_clicks > 0, False, False
    elif triggered == "full-cell-type-help-button":
        return False, celltype_n_clicks > 0, False
    elif triggered == "full-ligand-effect-help-button":
        return False, False, ligand_effect_n_clicks > 0
    else:
        raise PreventUpdate


@callback(
    Output('home-tab-title', 'children'),
    Output('home-tab-content', 'children'),
    Input('home-tabs', 'active_tab')
)
def update_tab_content(active_tab):
    if active_tab == 'default':
        return [html.I("ContactTracing"), ", studying the impact of chromosomal instability (CIN) & STING on the tumor ecosystem."], \
            home_welcome_info()
    elif active_tab == 'dataset':
        return [html.I("ContactTracing"), " is a unique approach for profiling cellular responses to the tumor microenvironment."], \
            home_dataset_descriptions()
    elif active_tab == 'plot':
        return [html.I("ContactTracing"), " produces information-rich data on condition-specific genetic responses."], \
            home_plot_descriptions()
    elif active_tab == 'misc':
        return ['Additional information on ', html.I("ContactTracing"), ' and this website.'], \
            home_misc_info()
    else:
        raise PreventUpdate
