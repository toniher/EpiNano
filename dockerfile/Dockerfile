# https://micromamba-docker.readthedocs.io/en/latest/quick_start.html
FROM mambaorg/micromamba:2.5
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
  micromamba clean --all --yes
USER root
RUN apt update && \
    apt upgrade -y && \
    apt install -y unzip procps && \
    rm -rf /var/lib/apt/lists/*


RUN mkdir -p /opt/epinano
WORKDIR /opt/epinano
COPY epinano* /opt/epinano/
COPY Epinano* /opt/epinano/
COPY LICENSE /opt/epinano/
COPY misc/ /opt/epinano/misc/
COPY models/ /opt/epinano/models/
#COPY Reference_sequences/ /opt/epinano/Reference_sequences/
#COPY test_data /opt/epinano/test_data/

# Workaround to make it work with Apptainer/Singularity
RUN unlink /opt/conda/bin/python
RUN echo "micromamba run -n base /opt/conda/bin/python3.11 \"\$@\"" > /opt/conda/bin/python
RUN mv /opt/conda/bin/Rscript /opt/conda/bin/Rscript.orig
RUN echo "micromamba run -n base /opt/conda/bin/Rscript.orig \"\$@\"" > /opt/conda/bin/Rscript

RUN chmod a+rx *py *sh *R
RUN chmod a+rx /opt/conda/bin/python /opt/conda/bin/Rscript

USER mambauser
ENV PATH="$PATH:/opt/epinano:/opt/conda/bin"
