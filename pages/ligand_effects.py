import dash

dash.register_page(__name__,
                   path='/ligand-effects',
                   name='Ligand Effects',
                   order=1)

layout = []
