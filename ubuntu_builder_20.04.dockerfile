FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get -qq update && apt-get -qq upgrade

RUN apt-get -qq --no-install-recommends install libpython3.8 python3.8 python3-pip libpython3.8-dev \
&&  apt-get -qq --no-install-recommends install qtbase5-dev libqt5charts5-dev libqt5opengl5-dev qtwebengine5-dev libopengl0 libeigen3-dev libglew-dev zlib1g-dev \
&&  apt-get -qq --no-install-recommends install libboost1.71-dev libboost-system1.71-dev libboost-filesystem1.71-dev libboost-program-options1.71-dev libboost-thread1.71-dev \
&&  apt-get -qq --no-install-recommends install libmkl-dev libmkl-interface-dev libmkl-threading-dev libmkl-computational-dev \
&&  apt-get -qq --no-install-recommends install g++ ccache ninja-build curl unzip git \
&&  python3.8 -m pip install numpy \
&&  apt-get clean

# Install CMake
ADD https://github.com/Kitware/CMake/releases/download/v3.20.1/cmake-3.20.1-linux-x86_64.sh /tmp
RUN chmod a+x /tmp/cmake-3.20.1-linux-x86_64.sh \
&&  /tmp/cmake-3.20.1-linux-x86_64.sh --skip-license --prefix=/usr/local \
&&  rm /tmp/cmake-3.20.1-linux-x86_64.sh

# Install pybind11
RUN git clone --depth 1 -b v2.4  https://github.com/pybind/pybind11.git /tmp/pybind11 \
&&  cmake -GNinja -S/tmp/pybind11 -B/tmp/pybind11/build -DPYBIND11_TEST=OFF -DCMAKE_BUILD_TYPE=Release \
&& cmake --install /tmp/pybind11/build \
&& rm -rf /tmp/pybind11

# Install VTK
RUN git clone --depth 1 --branch v9.0.1 https://gitlab.kitware.com/vtk/vtk.git /tmp/vtk \
&&  cd /tmp/vtk \
&&  git submodule update --init \
&&  mkdir build && cd build \
&&  cmake \
      -GNinja \
      -DCMAKE_C_COMPILER_LAUNCHER=ccache \
      -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
      -DCMAKE_BUILD_TYPE=Release \
      -DVTK_GROUP_ENABLE_Imaging=NO \
      -DVTK_GROUP_ENABLE_MPI=NO \
      -DVTK_GROUP_ENABLE_Qt=NO \
      -DVTK_GROUP_ENABLE_Rendering=NO \
      -DVTK_GROUP_ENABLE_StandAlone=DEFAULT \
      -DVTK_GROUP_ENABLE_Views=NO \
      -DVTK_GROUP_ENABLE_Web=NO \
      -DVTK_RELOCATABLE_INSTALL=ON \
      -DVTK_MODULE_ENABLE_VTK_IOLegacy=YES \
      -DBUILD_SHARED_LIBS=OFF \
      .. \
&&   cmake --build . \
&&   cmake --install . \
&& cd / \
&& rm -rf /tmp/vtk

WORKDIR /opt
ENV SOFA_ROOT=/opt/sofa
ENV PYTHONPATH=/opt/sofa/plugins/SofaPython3/lib/python3/site-packages