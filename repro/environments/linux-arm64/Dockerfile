FROM julia:1.12.6-bookworm AS julia

FROM mambaorg/micromamba:2.3.3

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash ca-certificates coreutils libxml-libxml-perl patch rsync \
    && rm -rf /var/lib/apt/lists/*

COPY --from=julia /usr/local/julia /usr/local/julia

USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge \
        python=3.11.16 \
        cmake=3.31.8 \
        make=4.4.1 \
        git=2.55.0 \
        gcc_linux-aarch64=14.3.0 \
        gfortran_linux-aarch64=14.3.0 \
        esmf=8.9.1=mpi_mpich_hc00c11e_100 \
        mpich=4.3.2 \
        libblas=3.11.0 \
        liblapack=3.11.0 \
        libnetcdf=4.9.3 \
        netcdf-fortran=4.6.2 \
    && micromamba clean --all --yes

USER root
RUN ln -sf /usr/bin/perl /opt/conda/bin/perl \
    && ln -sf /opt/conda/bin/aarch64-conda-linux-gnu-ar /opt/conda/bin/ar \
    && ln -sf /opt/conda/bin/aarch64-conda-linux-gnu-ranlib /opt/conda/bin/ranlib \
    && ln -sf /opt/conda/bin/aarch64-conda-linux-gnu-gcc /opt/conda/bin/gcc \
    && ln -sf /opt/conda/bin/aarch64-conda-linux-gnu-g++ /opt/conda/bin/g++ \
    && ln -sf /opt/conda/bin/aarch64-conda-linux-gnu-gfortran /opt/conda/bin/gfortran
USER $MAMBA_USER

# Keep Debian Perl ahead of the conda dependency so CIME sees XML::LibXML;
# compiler, MPI, Python, and numerical-library tools still resolve from /opt/conda.
ENV PATH=/usr/local/julia/bin:/usr/local/bin:/usr/bin:/bin:/opt/conda/bin \
    NETCDF=/opt/conda \
    NETCDF_PATH=/opt/conda \
    PNETCDF=/opt/conda \
    PNETCDF_PATH=/opt/conda \
    ESMFMKFILE=/opt/conda/lib/esmf.mk \
    JULIA_DEPOT_PATH=/opt/gmd/julia-depot

WORKDIR /work
