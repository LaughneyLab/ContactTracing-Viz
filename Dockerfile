FROM continuumio/miniconda3 AS build

RUN conda create -n condaenv python=3.9 pip

ADD requirements.txt .

# Install requirements with conda (with fallback to pip)
# https://gist.github.com/kobybibas/a0d126a8bc007999ae0bcf8b9980fafa
RUN while read requirement; do conda install -n condaenv -c conda-forge --yes $requirement || (conda activate condaenv && pip install --no-cache-dir $requirement); done < requirements.txt

# Prepare to use redis
RUN conda activate condaenv && pip install --no-cache-dir celery[redis]

# Webserver
RUN conda install -n condaenv -c conda-forge --yes gunicorn

# Lets us move the conda environment into a smaller image
RUN conda install -c conda-forge conda-pack
RUN conda-pack -n condaenv -o /tmp/condaenv.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/condaenv.tar && \
    rm /tmp/condaenv.tar
RUN /venv/bin/conda-unpack


FROM redis:7-bullseye AS runtime

# Copy environment
COPY --from=build /venv /venv

ENV REDIS_URL="127.0.0.1:6379"
EXPOSE 8000

ADD . .

SHELL ["/bin/sh", "-c"]

ENTRYPOINT redis-server && \
    source /venv/bin/activate && \
    gunicorn app:server --workers 3 --daemon --bind '127.0.0.1:8000' && \
    celery -A app:celery_app worker --loglevel=info --concurrency=2