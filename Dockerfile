FROM sd2e/apps:python3-miniconda

RUN conda config --add channels r \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda install cutadapt --yes \
    && conda install bowtie2 --yes \
    && conda install samtools --yes \
    && conda install bbmap --yes

RUN mkdir /opt/scripts/

#ADD /src/srna_fastq_processor.py /opt/scripts/srna_fastq_processor.py
#RUN chmod 777 /opt/scripts/srna_fastq_processor.py

# ENV PATH "/opt/bin/:$PATH"
# ADD config.yml /config.yml


