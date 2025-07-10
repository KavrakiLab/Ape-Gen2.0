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

RUN ["/bin/bash", "-c", "wget https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Miniforge-pypy3-Linux-x86_64.sh"]
RUN chmod 0755 Miniforge-pypy3-Linux-x86_64.sh
RUN bash Miniforge-pypy3-Linux-x86_64.sh -b -p /mamba
ENV PATH "/mamba/bin:$PATH"
RUN rm Miniforge-pypy3-Linux-x86_64.sh


# environment.yml has pdbfixer, numpy, mdtraj and openmm
COPY ./Docker_image/environment.yml /environment.yml
RUN mamba env create -f /environment.yml
RUN echo "source activate apegen" >> ~/.bashrc
ENV PATH "/mamba/envs/apegen/bin:$PATH"

# These are okay outside environment.yml because they install to PATH
# Perhaps better if inside yml?
RUN mamba install -c bioconda smina
RUN mamba install -c bioconda autodock-vina

# This for some reason creates issues with smina, I am not sure why
RUN pip install mpire # 
# RUN pip install pdb-tools

# Probably you should be cloning the directory if you want to build the image yourself?
# Makes sense for other people to download I guess
# RUN git clone https://github.com/KavrakiLab/APE-Gen.git
ADD . $HOME/data/

# replace files that cause warnings and errors
# COPY /replacement_files/replace_files.sh /replace_files.sh
# RUN echo "sh /replace_files.sh" >> ~/.bashrc

# RCD
# RUN mamba install -c intel mkl # This is here, libraries are missing if I leave it in the .yml
RUN mamba install -y -c conda-forge mkl

ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/mamba/lib/"
RUN wget -qO- http://chaconlab.org`wget -qO- http://chaconlab.org/modeling/rcd/rcd-download|grep txz|grep uk-button| cut -f4 -d\"` | tar Jxf - && \
    cp RCD_v1.40_Linux_20190228/bin/rcd /usr/local/bin && \
    cp RCD_v1.40_Linux_20190228/dunbrack.bin /data/RCD_required_files/dunbrack.bin

# KORPE
RUN wget -qO- http://chaconlab.org`wget -qO- http://chaconlab.org/modeling/korp/down-korp|grep Korp6Dv1.txz|grep uk-button| cut -f4 -d\"` | tar Jxf - && \
    cp Korp6Dv1/korp6Dv1.bin /data/RCD_required_files/korp6Dv1.bin && \
    cp Korp6Dv1/bin/korpe /data/RCD_required_files/korpe

# Change dir
WORKDIR /data

RUN cp ./helper_scripts/prepare_ligand4.py /mamba/envs/apegen/MGLToolsPckgs/AutoDockTools/Utilities24/ # Functional
RUN grep -rl 'simtk.openmm' /mamba/ | xargs sed -i 's/simtk.openmm/openmm/g' # Avoiding warnings when importing

ENTRYPOINT ["bash"]

# # Build this Dockerfile into an image by getting to the Dockerfile dir and executing:
# # docker build . -t kavrakilab/apegen2.0