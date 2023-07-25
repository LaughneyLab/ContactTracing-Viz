# ContactTracing-Viz
[![DOI](https://zenodo.org/badge/569471421.svg)](https://zenodo.org/badge/latestdoi/569471421)

A repository for visualizing ContactTracing Results in a web browser.


## Local Deployment for exploring custom results

A limited version of the site can be used to explore custom-run ContactTracing results following the tutorial located 
here: https://github.com/LaughneyLab/ContactTracing_tutorial

This requires the additional installation of jupyter, diskcache, dill, and multiprocess

Then you may follow the steps in the `interactive_notebook.ipynb` notebook to generate the results and visualize results.


## Deployment for the official CIN-TME Results Website

NOTE: For local deployment, this uses the diskcache library for job queueing.

When deployed, it is recommended to use celery + redis instead. To enable this:
1. Install redis and set the `REDIS_URL` environment variable to the redis url
2. Install the celery python package.
