import dash

from viz.web import interactive_panel, wrap_icon

dash.register_page(__name__,
                   path='/interactions',
                   name='Interactions',
                   order=1)

layout = [
    interactive_panel(wrap_icon('fa-arrows-left-right-to-line', 'Cell Type to Cell Type Interactions'), "WIP")
]
