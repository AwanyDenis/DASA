FROM continuumio/miniconda3
WORKDIR /work
VOLUME /work

VOLUME /data


# RUN mkdir ref
# RUN mkdir InputFastqDir

COPY Snakefile config.yaml environment.yaml ./
COPY envs/ ./envs/

RUN conda env create -f environment.yaml 

ENV PATH /opt/conda/envs/pipeline_env/bin:$PATH
RUN /bin/bash -c "source activate pipeline_env"
