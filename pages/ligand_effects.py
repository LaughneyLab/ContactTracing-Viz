import dash

from viz.web import interactive_panel, wrap_icon

dash.register_page(__name__,
                   path='/ligand-effects',
                   name='Ligand Effects',
                   order=1)

layout = [
    interactive_panel(wrap_icon('fa-maximize', 'Cell Type to Cell Type Interactions'), "WIP")
]
