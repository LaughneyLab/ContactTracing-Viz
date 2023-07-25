import anndata as ad
import pandas as pd

__name__ = "jupyter_support"

from dash import jupyter_dash

# Configure the viewing mode
# external = display a link to click
# tab = open in a new browser tab
# jupyterlab = display in tab if run in jupyterlab
jupyter_dash.default_mode = "external"

# Support notebooks behind a proxy
try:
    jupyter_dash.infer_jupyter_proxy_config()
except: pass

from viz import data


def run_app_in_jupyter(
        interactions: pd.DataFrame,
        adata: ad.AnnData,
        deg_adata: ad.AnnData,
        target_stats: pd.DataFrame,
        condition_name: str,
        skip_compile: bool,
        port: int
):
    """
    Run the visualization within a Jupyter notebook.
    :param interactions: The interactions file run with ContactTracing.
    :param adata: The contactTracing results output data.
    :param port: The port to host the app in.
    """
    # Remove cin-tme specific references
    data.EXCLUDED_CELL_TYPES = []

    if not skip_compile:
        print("Compiling the provided dataset...This may take awhile.", flush=True)
        data.compile_custom_dataset(adata, deg_adata, target_stats, interactions, condition_name)

    from app import app

    return app.run(
        port=port,
        jupyter_mode='external',
        jupyter_width='100%',
        jupyter_height=1000,
    )
