FROM ghcr.io/prefix-dev/pixi:latest

RUN pixi config set --global run-post-link-scripts insecure

RUN pixi global install -c conda-forge -c bioconda r-base

RUN pixi global install -c conda-forge -c bioconda htslib

RUN pixi global install -c conda-forge -c bioconda bedtools

RUN pixi global install --environment r-base -c conda-forge -c bioconda \
    r-argparser \
    r-data.table \
    r-tidyverse \
    r-magrittr \
    r-pbmcapply \
    r-bedr \
    r-cowplot \
    r-mclust \
    r-mvtnorm \
    r-tibble \
    r-matrixstats \
    r-ggrepel \
    bioconductor-limma \
    bioconductor-edger \
    bioconductor-pcatools \
    bioconductor-fastseg \
    bioconductor-biostrings \
    bioconductor-plyranges \
    bioconductor-genomicranges

RUN find ${HOME}/.pixi/envs/r-base/bin -name '*bioconductor-*-post-link.sh' | \
   xargs -I % bash -c "PREFIX=${HOME}/.pixi/envs/r-base PATH=${HOME}/.pixi/envs/r-base/bin:${PATH} %"

CMD ["bash"]
