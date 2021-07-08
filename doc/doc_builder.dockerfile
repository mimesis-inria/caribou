FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get -qq update \
&&  apt-get -qq --no-install-recommends install doxygen python3 python3-pip texlive-base texlive-font-utils \
&&  apt-get clean

RUN pip3 install -r https://raw.githubusercontent.com/jnbrunet/caribou/master/doc/sphinx/source/requirements.txt