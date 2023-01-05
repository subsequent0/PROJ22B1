# hash:sha256:db82f0110e9d241394214d077e4420a4bfe0c6f6b264889117e4474428a9ad73
FROM registry.codeocean.com/codeocean/r-studio:2022.07.0-548-r4.2.1-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        clusttree \
        data.table \
        ggsci \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("celldex")' \
    && Rscript -e 'remotes::install_version("ggplot2")' \
    && Rscript -e 'remotes::install_version("monocle")'

RUN Rscript -e 'remotes::install_github( \
        "Cellchat")' \
    && Rscript -e 'remotes::install_github( \
        "Rcolorbrewer")' \
    && Rscript -e 'remotes::install_github( \
        "Seurat")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "SingleR", \
        "cowplot", \
        "tidyverse" \
    ))' # Original versions: latest latest latest
