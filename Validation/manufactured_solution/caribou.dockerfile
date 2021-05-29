FROM jnbrunet/caribou-ubuntu-1804-builder:v20.12 as build
ENV DEBIAN_FRONTEND noninteractive
ENV SOFA_ROOT /opt/sofa
ENV CARIBOU_ROOT /opt/sofa/plugins/SofaCaribou
WORKDIR /opt

# Install SOFA
RUN curl --output sofa.zip -L  "https://github.com/sofa-framework/sofa/releases/download/v20.12.02/SOFA_v20.12.02_Linux.zip" \
&&  unzip sofa.zip -d temp \
&&  mv temp/`ls temp` $SOFA_ROOT \
&&  rm -rf temp

# Install pybind11
RUN git clone --depth 1 -b v2.4  https://github.com/pybind/pybind11.git /tmp/pybind11 \
&&  cmake -S/tmp/pybind11 -B/tmp/pybind11/build -DPYBIND11_TEST=OFF -DCMAKE_BUILD_TYPE=Release \
&& cmake --install /tmp/pybind11/build \
&& rm -rf /tmp/pybind11

# Compile Caribou
RUN git clone --depth 1 https://github.com/jnbrunet/caribou.git /tmp/caribou \
&&  cmake \
        -DCARIBOU_BUILD_TESTS=ON \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=$CARIBOU_ROOT \
        -B/tmp/caribou/build \
        -S/tmp/caribou \
&&  cmake --build /tmp/caribou/build -- -j8  && cmake --install /tmp/caribou/build \
&& rm -rf /tmp/caribou

# Install bindings
RUN ln -s $CARIBOU_ROOT/lib/python3/site-packages/* /usr/local/lib/python3.7/dist-packages/

FROM ubuntu:20.04

ENV SOFA_ROOT /opt/sofa
ENV CARIBOU_ROOT /opt/sofa/plugins/SofaCaribou
ENV LD_LIBRARY_PATH /opt/sofa/lib

# Install dependencies
RUN apt-get update && \
    apt-get -qq --no-install-recommends install software-properties-common && \
    add-apt-repository -y ppa:deadsnakes && \
    apt-get -qq update && \
    apt-get -qq --no-install-recommends install \
        libpython3.7 python3.7 python3-pip libgl1-mesa-glx libglib2.0-0 libglx0 libopengl0 libharfbuzz0b libgomp1 curl ca-certificates unzip libnss3 libfontconfig libxss1 libasound2 libopus0 libxslt1.1  && \
    apt-get clean && \
    rm /usr/bin/python3 && ln -s /usr/bin/python3.7 /usr/bin/python3 && ln -s /usr/bin/python3.7 /usr/bin/python && \
    python -m pip install numpy meshio sympy

# Install SOFA + SofaCaribou
COPY --from=build /opt/sofa /opt/sofa

# Install bindings
RUN ln -s $CARIBOU_ROOT/lib/python3/site-packages/* /usr/local/lib/python3.7/dist-packages/ && \
    ln -s $SOFA_ROOT/plugins/SofaPython3/lib/python3/site-packages/* /usr/local/lib/python3.7/dist-packages/ && \
    ln -s $SOFA_ROOT/bin/runSofa /usr/local/bin/runSofa
