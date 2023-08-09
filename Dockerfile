FROM mambaorg/micromamba:alpine AS build

# Use the faster libmamba solver for conda
#RUN conda install --yes -n base conda-libmamba-solver
#RUN conda config --system --set solver libmamba

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml .
RUN micromamba install -y -n base -f environment.yml

#RUN conda create -n condaenv python=3.9 pip

#ADD requirements.txt .

# Install requirements with conda (with fallback to pip)
# https://gist.github.com/kobybibas/a0d126a8bc007999ae0bcf8b9980fafa
#RUN while read requirement; do conda install -n condaenv -c conda-forge --yes $requirement || (conda run -n condaenv pip install --no-cache-dir $requirement); done < requirements.txt
#RUN conda run -n condaenv pip install --no-cache-dir -r requirements.txt

# Prepare to use redis
#RUN conda run -n condaenv pip install --no-cache-dir celery[redis]

# Webserver
#RUN conda install -n condaenv -c conda-forge --yes gunicorn

# Lets us move the conda environment into a smaller image
RUN micromamba install -y -n base -c conda-forge conda-pack && \
    micromamba clean --all --yes

WORKDIR /

USER root

RUN micromamba run -n base conda-pack --prefix /opt/conda -o /tmp/condaenv.tar.gz && \
    mkdir venv && cd venv && tar xf /tmp/condaenv.tar.gz && \
    rm /tmp/condaenv.tar.gz

WORKDIR /venv

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN bin/conda-unpack

# Compile data and default figures
#WORKDIR /tmp
#COPY . .
# Pre-compile the data if it's not available
#RUN micromamba run -n base \
#    python viz/data.py

# Pre-compile default figures
#RUN PYTHONPATH=./ micromamba run -n base \
#    python pages/circos.py
#RUN PYTHONPATH=./ micromamba run -n base \
#    python pages/interactions.py
#RUN PYTHONPATH=./ micromamba run -n base \
#    python pages/ligand_effects.py
#RUN mkdir -p /figures && mv data/compiled/*.pkl /figures/

#RUN rm -rf /tmp/*


FROM redis:bullseye AS runtime

# Copy environment
COPY --from=build /venv /venv
#COPY --from=build /figures /figures

WORKDIR /app

ENV REDIS_URL="redis://localhost:6379"
EXPOSE 8000

SHELL ["/bin/bash", "-c"]

ADD . .
#RUN rm -rf data/compiled

# File system data
#ADD my_efs .

#RUN mv my_efs/data/compiled data/ || true

VOLUME /app/data

ENTRYPOINT redis-server --daemonize yes --maxmemory 4g --latency-tracking no && \
    source /venv/bin/activate && \
    celery -A app:celery_app worker --loglevel=info --detach && \
    sleep 5 && \
    PYTHONPATH=./ fastwsgi app:server --port 8000
#    PYTHONPATH=./ gunicorn app:server --workers 1 --bind '0.0.0.0:8000'
