import dash_bootstrap_components as dbc
from dash import html, Output, Input


def jumbotron(title, main_content, sub_content, *additional_content, dark=False):
    return html.Div(dbc.Container([
        html.H1(title, className='display-3'),
        html.P(main_content, className='lead'),
        html.Hr(className='my-2'),
        html.P(sub_content),
        *[html.P(content) for content in additional_content]
    ], fluid=True, className='py-3'), className="p-3 rounded-3 " + ('bg-dark text-white' if dark else 'bg-light'))


def interactive_panel(title, *content, dark=False):
    return dbc.Card([
        dbc.CardHeader(title),
        dbc.CardBody([
            html.Div(content)
        ])
    ], className='mb-3 ' + ('bg-dark text-white' if dark else 'bg-light'))


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
        dbc.Col(dbc.Button(wrap_icon('fa-arrows-left-right-to-line', 'View Cell Type Interactions'), id='ct-interaction-button', href='/interactions')),
        dbc.Col(dbc.Button(wrap_icon('fa-maximize', 'View Ligand Effects'), id='ligand-effect-button', href='/ligand-effects')),
        dbc.Col(dbc.Button(wrap_icon('fa-circle-dot', 'View Circos'), id='circos-button', href='/circos')),
        dbc.Col(dbc.Button(wrap_icon('fa-link-slash', 'Unselect Data'), id='unselect-button', n_clicks=0))  # Callback in pages/home
    ])
