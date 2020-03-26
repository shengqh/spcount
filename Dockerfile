FROM shengqh/bioinfo:base

ARG BOWTIE_VERSION="1.2.3"
RUN cd /opt; \
    wget https://github.com/BenLangmead/bowtie/releases/download/v${BOWTIE_VERSION}/bowtie-${BOWTIE_VERSION}-linux-x86_64.zip; \
    unzip bowtie-${BOWTIE_VERSION}-linux-x86_64.zip; \
    rm bowtie-${BOWTIE_VERSION}-linux-x86_64.zip
ENV PATH=$PATH:/opt/bowtie-${BOWTIE_VERSION}-linux-x86_64
RUN bowtie --version

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir --upgrade \
    pysam \
    biopython \
    pytabix

RUN apt-get install -y tabix
    
ENV CPDSEQER_VERSION="0.0.1"
RUN pip3 install git+git@github.com:shengqh/spcount.git
RUN spcount -h
