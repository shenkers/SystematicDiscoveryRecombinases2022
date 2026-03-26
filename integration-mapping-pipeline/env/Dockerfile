FROM continuumio/miniconda3

RUN apt-get -y update && apt-get -y install gcc make libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev g++

# DEBUG
# RUN apt-get -y install vim less htop
# DEBUG

RUN conda install -n base -c conda-forge conda-libmamba-solver && \
    conda config --set solver libmamba

RUN conda install -y -c conda-forge -c bioconda -c r \
    python=3.8 snakemake samtools bowtie2 htslib emboss cd-hit bwa fastp fastqc r r-essentials

RUN conda install -y -c conda-forge -c bioconda \
    click pandas biopython pysam editdistance scipy networkx tqdm umi_tools \
    pybedtools bioconductor-genomicranges bioconductor-helloranges \
    bioconductor-biostrings r-bedr r-ggseqlogo

RUN pip install git+https://github.com/hsulab-arc/MGEfinder-integration-mapping-pipeline.git

WORKDIR /integration-mapping-pipeline/snakemake
