FROM continuumio/miniconda3 AS build

# Use the faster libmamba solver for conda
RUN conda install --yes -n base conda-libmamba-solver
RUN conda config --system --set solver libmamba

RUN conda create -n condaenv python=3.9 pip

ADD requirements.txt .

# Install requirements with conda (with fallback to pip)
# https://gist.github.com/kobybibas/a0d126a8bc007999ae0bcf8b9980fafa
RUN while read requirement; do conda install -n condaenv -c conda-forge --yes $requirement || (conda run -n condaenv pip install --no-cache-dir $requirement); done < requirements.txt

# Prepare to use redis
RUN conda run -n condaenv pip install --no-cache-dir celery[redis]

# Webserver
RUN conda install -n condaenv -c conda-forge --yes gunicorn

# Lets us move the conda environment into a smaller image
RUN conda install -c conda-forge --yes conda-pack
RUN conda-pack -n condaenv -o /tmp/condaenv.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/condaenv.tar && \
    rm /tmp/condaenv.tar
RUN /venv/bin/conda-unpack


FROM redis:bullseye AS runtime

# Copy environment
COPY --from=build /venv /venv

ENV REDIS_URL="redis://localhost:6379"
EXPOSE 8000

SHELL ["/bin/bash", "-c"]

ADD . .

ENTRYPOINT redis-server --daemonize yes && \
    source /venv/bin/activate && \
    celery -A app:celery_app worker --loglevel=info --concurrency=2 --detach  && \
    gunicorn app:server --workers 4 --bind '0.0.0.0:8000'
