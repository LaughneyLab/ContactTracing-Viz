import dash

from viz.web import interactive_panel, wrap_icon

dash.register_page(__name__,
                   path='/circos',
                   name='Circos',
                   order=2)

layout = [
    interactive_panel(wrap_icon('fa-circle-dot', 'Cell Type to Cell Type Interactions'), "WIP")
]

