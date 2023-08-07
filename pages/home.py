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
        html.Div(),  # "ContactTracing Interactive Visualizer",
        html.Div(),   # "This is a web application that allows for ContactTracing outputs to be easily interpreted.",
        dbc.Row([
            dbc.Accordion([
                dbc.AccordionItem([
                    html.H4([html.I("ContactTracing"), ", studying the impact of chromosomal instability (CIN) & STING on the tumor ecosystem."], className='card-title'),
                    html.Hr(),
                    html.P(),
                    html.Div(home_welcome_info(), className='card-text'),
                ], title='Welcome'),
                dbc.AccordionItem([
                    html.H4([html.I("ContactTracing"), " is a unique approach for profiling cellular responses to the tumor microenvironment."], className='card-title'),
                    html.Hr(),
                    html.P(),
                    html.Div(home_dataset_descriptions(), className='card-text'),
                ], title='Dataset Information'),
                dbc.AccordionItem([
                    html.H4([html.I("ContactTracing"), " produces information-rich data on condition-specific genetic responses."], className='card-title'),
                    html.Hr(),
                    html.P(),
                    html.Div(home_plot_descriptions(), className='card-text'),
                ], title='Visualization Information'),
                dbc.AccordionItem([
                    html.H4(['Additional information on ', html.I("ContactTracing"), ' and this website.'], className='card-title'),
                    html.Hr(),
                    html.P(),
                    html.Div(home_misc_info(), className='card-text'),
                ], title='Miscellaneous'),
            ], active_item=['item-0', 'item-1', 'item-2', 'item-3'], always_open=True)
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
            html.Div(["Austin Varela: ", html.A('aav4003@med.cornell.edu', href='mailto:aav4003@med.cornell.edu')]),
            html.P(),
            html.Div(["Dr. Melissa Hubisz: ", html.A('mt269@cornell.edu', href='mailto:mt269@cornell.edu')])
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
