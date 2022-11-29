# ContactTracing-Viz

A repository for visualizing ContactTracing Results in a web browser.

NOTE: For local deployment, this uses the diskcache library for job queueing.
When deployed, it is recommended to use celery + redis instead. To enable this:
1. Install redis and set the `REDIS_URL` environment variable to the redis url
2. Install the celery python package.
