FROM sd2e/apps:python3-miniconda

RUN conda config --add channels r \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda install cutadapt --yes \
    && conda install bowtie2 --yes \
    && conda install samtools --yes \
    && conda install bbmap --yes

# create env?

# ENV PATH "/opt/bin/:$PATH"
# ADD config.yml /config.yml
ADD src /opt/src
