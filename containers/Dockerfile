# Same environment as containers/subfam.def. Build from the repository root:
#   docker build -f containers/Dockerfile -t subfam .
FROM mambaorg/micromamba:2.0.5
USER root
COPY environment.yml /opt/subfam/environment.yml
RUN micromamba create -y -n subfam -f /opt/subfam/environment.yml && micromamba clean -a -y
COPY SubFam.sh kmer_order.c /opt/subfam/
RUN chmod 755 /opt/subfam/SubFam.sh && ln -s /opt/subfam/SubFam.sh /usr/local/bin/SubFam.sh
ENV PATH=/opt/conda/envs/subfam/bin:$PATH LC_ALL=C
WORKDIR /data
ENTRYPOINT ["/opt/subfam/SubFam.sh"]
