# pull base image
FROM cancerbits/dockr:3.6.3-c

# docker build -t cancerbits/dockr:meth_imprint -f meth_imprint.Dockerfile .

LABEL maintainer Florian Halbritter "florian.halbritter@ccri.at"
LABEL version meth_imprint

# e.g.:
# RUN installGithub.r chris-mcginnis-ucsf/MULTI-seq@a969bc453376f47a6a578a3d2e103f318f13e388
# RUN R -e "BiocManager::install(c('doParallel'))"
# RUN install2.r --error qlcMatrix


RUN R -e "BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene')"
RUN R -e "BiocManager::install('dmrseq')"
RUN R -e "BiocManager::install('ComplexHeatmap')"
RUN R -e "BiocManager::install('circlize')"
RUN R -e "BiocManager::install('ggstance')"

