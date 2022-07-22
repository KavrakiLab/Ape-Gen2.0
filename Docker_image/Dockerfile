FROM ubuntu:bionic
ENV DEBIAN_FRONTEND=noninteractive

# The following line is important for singularity
#RUN mkdir /scratch /work /home1 /gpfs /corral-repl /corral-tacc /data

RUN apt-get update && \
    apt-get -y install \
        cmake \
        cython \
        doxygen \
        freeglut3-dev \
        git \
        gfortran \
        libboost-date-time1.62-dev \
        libboost-filesystem1.62-dev \
        libboost-iostreams1.62-dev \
        libboost-program-options1.62-dev \
        libboost-serialization1.62-dev \
        libboost-system1.62-dev \
        libboost-thread1.62-dev \
        libboost-timer1.62-dev \
        libeigen3-dev \
        libfftw3-dev \
        libfreetype6-dev \
        libglew-dev \
        libglm-dev \
        libpng-dev \
        libpython-dev \
        libxml2-dev \
        python-pip \
        python-pmw \
        swig \
        wget \
        vim


RUN ["/bin/bash", "-c", "wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"]
RUN chmod 0755 Miniconda3-latest-Linux-x86_64.sh
RUN ./Miniconda3-latest-Linux-x86_64.sh -b -p /conda
ENV PATH "/conda/bin:$PATH"

# environment.yml has pdbfixer, numpy, mdtraj and openmm
COPY ./Docker_image/environment.yml /environment.yml
RUN conda env create -f /environment.yml
RUN echo "source activate apegen" >> ~/.bashrc
ENV PATH "/conda/envs/apegen/bin:$PATH"


# This is manual, but won't work
#RUN conda create --name apegen python=3.6 \
#    && echo "conda activate apegen" >> ~/.bashrc \
#    && conda init bash \
#    && conda activate apegen \
#    && conda update conda \
#    && conda install -c schrodinger pymol \ 
#    && conda install -c conda-forge/label/cf202003 requests \
#    && conda install -c conda-forge pdbfixer \
#    && conda install -c conda-forge numpy \
#    && conda install -c conda-forge mdtraj \
#    && conda install -c conda-forge openmm \
#    && conda install -c conda-forge biopython \
#    && conda install -c salilab modeller \
#    && conda install -c conda-forge beautifulsoup4
#ENV PATH "/conda/envs/apegen/bin:$PATH"

# These are okay outside environment.yml because they install to PATH
# Perhaps better if inside yml?
RUN conda install -c bioconda smina
RUN conda install -c bioconda autodock-vina

# This for some reason creates conflicts in the conda .yml file, not sure as to why
RUN pip install mpire
RUN pip install pdb-tools

# RCD
RUN conda install -c intel mkl
ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/conda/lib/"
RUN wget -qO- http://chaconlab.org`wget -qO- http://chaconlab.org/modeling/rcd/rcd-download|grep txz|grep uk-button| cut -f4 -d\"` | tar Jxf - && \
    cp RCD_v1.40_Linux_20190228/bin/rcd /usr/local/bin

# Probably you should be cloning the directory if you want to build the image yourself?
# Makes sense for other people to download I guess
# RUN git clone https://github.com/KavrakiLab/APE-Gen.git
ADD . $HOME/data/

# replace files that cause warnings and errors
# COPY /replacement_files/replace_files.sh /replace_files.sh
# RUN echo "sh /replace_files.sh" >> ~/.bashrc

WORKDIR /data

RUN wget -qO- http://chaconlab.org`wget -qO- http://chaconlab.org/modeling/korp/down-korp | grep Korp6Dv1.txz |grep uk-button | cut -f4 -d\"` | tar Jxf - && \
    cp ./Korp6Dv1/bin/korpe ./RCD_required_files/korpe && \
    cp ./Korp6Dv1/korp6Dv1.bin ./RCD_required_files/korp6Dv1.bin && \
    rm -r ./Korp6Dv1

ENTRYPOINT ["bash"]

# # Build this Dockerfile into an image by getting to the Dockerfile dir and executing:
# # docker build . -t kavrakilab/apegen2.0
