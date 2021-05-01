FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
WORKDIR /opt

# Install dependencies
RUN apt-get -qq update \
&&  apt-get -qq --no-install-recommends install \
       ca-certificates git cmake ninja-build g++ zlib1g-dev \
       libmkl-dev libmkl-interface-dev libmkl-threading-dev libmkl-computational-dev \
       python3 python-is-python3 python3-pip \
&&  python -m pip install numpy meshio

# Fetch & patch FEBio3
RUN git clone --depth 1 --branch v3.3.2 https://github.com/febiosoftware/FEBio.git febio_src
RUN sed -i 's|sprintf(ops.szcnf, "%sfebio.xml", szpath);|sprintf(ops.szcnf, "%s", "/etc/febio/febio.xml");|g' febio_src/FEBio3/FEBioApp.cpp
RUN sed -i '/char szpath/d' febio_src/FEBio3/FEBioApp.cpp
RUN sed -i '/febio::get_app_path(szpath/d' febio_src/FEBio3/FEBioApp.cpp

# Build & install FEBio3
RUN cmake -S febio_src \
          -B febio_bld \
          -GNinja \
          -DCMAKE_BUILD_TYPE=Release \
          -DUSE_MKL=ON \
          -DMKL_LIB_DIR=/usr/lib/x86_64-linux-gnu \
          -DOMP_LIB=/usr/lib/x86_64-linux-gnu/libiomp5.so \
          -DMKL_INC=/usr/include/mkl \
&& cmake --build febio_bld \
&& mkdir -p /etc/febio \
&& mv febio_bld/bin/febio3 /usr/bin/febio3 \
&& mv febio_bld/bin/febio.xml /etc/febio/febio.xml

# Remove built files
RUN rm -rf febio_src febio_bld \
&& apt-get -qq remove git cmake ninja-build g++ libmkl-dev libmkl-interface-dev libmkl-threading-dev libmkl-computational-dev \
&& apt-get clean \
&& rm -rf /var/cache/apt/lists
