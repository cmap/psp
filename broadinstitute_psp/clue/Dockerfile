FROM continuumio/anaconda

RUN mkdir -p /usr/src/psp && \
mkdir -p /proteomics/qcnorm_dir && \
mkdir -p /proteomics/sim_dir

WORKDIR /usr/src/psp/proteomics-signature-pipeline/broadinstitute_psp

COPY . /usr/src/psp

RUN conda create -y -n psp -c bioconda pandas=0.20.3 scipy=0.19.0 h5py=2.7.0 cmapPy=2.0.1
RUN ["chmod", "+x", "/usr/src/psp/proteomics-signature-pipeline/setup_psp.sh"]
RUN ["/bin/bash", "/usr/src/psp/proteomics-signature-pipeline/setup_psp.sh"]
RUN ["chmod", "+x", "/usr/src/psp/proteomics-signature-pipeline/broadinstitute_psp/clue/prot_query.sh"]

ENTRYPOINT ["/bin/bash","clue/prot_query.sh"]
