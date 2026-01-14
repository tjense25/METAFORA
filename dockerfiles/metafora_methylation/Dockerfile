FROM ghcr.io/prefix-dev/pixi:latest

RUN pixi config set --global run-post-link-scripts insecure

RUN pixi global install -c conda-forge -c bioconda htslib

RUN pixi global install -c conda-forge -c bioconda bedtools

RUN pixi global install -c bioconda -c conda-forge -c epi2melabs modbam2bed

CMD ["bash"]
