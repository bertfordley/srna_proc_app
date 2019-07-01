FROM sd2e/apps:python2-miniconda

RUN conda config --add channels r \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda install cutadapt --yes \
    && conda install bowtie2 --yes \
    && conda install samtools --yes \
    && conda install bbmap --yes

ADD /src/srna_fastq_processor.py /src/srna_fastq_processor.py
ADD /tests /tests
ADD /indices /indices

