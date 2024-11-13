FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-python:1.1.5

USER root

WORKDIR /build

# Install R for Purity Reviewer rpy2 dependency
RUN apt-get update && \
    apt-get install -y --no-install-recommends r-base && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

# Install Purity Reviewer
RUN git clone https://github.com/getzlab/PurityReviewer.git && \
    pip install ./PurityReviewer && \
    rm -rf /build/PurityReviewer

USER $USER