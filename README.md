# Caribou

[![Documentation Status](https://readthedocs.org/projects/caribou/badge/?version=latest)](https://caribou.readthedocs.io/en/latest/?badge=latest)
[![MacOS](https://github.com/jnbrunet/caribou/actions/workflows/macos.yml/badge.svg)](https://github.com/jnbrunet/caribou/actions/workflows/macos.yml)
[![Linux](https://github.com/jnbrunet/caribou/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/jnbrunet/caribou/actions/workflows/ubuntu.yml)

The caribou project is aimed at multiphysics computation. 
It brings a plugin that complements [SOFA multiphysics framework](https://www.sofa-framework.org/). 
It also provides generic c++ utilities, and SOFA components such as solvers and forcefields.

The project is composed of two modules:
1. **Caribou library** brings multiple geometric, linear analysis and topological tools that are designed to 
be as independent as possible from external projects.
2. **Sofa caribou library** is built on top of the **caribou library**, but brings new components to
the SOFA project as a plugin.

### Quick installation

You can get the latest binaries compatible with your SOFA version by following these links:

| Sofa version | Linux                                                                                                      |MacOS | Windows |
| ------------ | -----------------------------------------------------------------------------------------------------------|------|---------|
| v20.06       | [SofaCaribou_latest.tar.gz](https://caribou.jnbrunet.com/artifacts/linux/v20.06/SofaCaribou_latest.tar.gz) | N/A  | N/A     |
| v20.12       | [SofaCaribou_latest.tar.gz](https://caribou.jnbrunet.com/artifacts/linux/v20.12/SofaCaribou_latest.tar.gz) | N/A  | N/A     |
| v21.06       | [SofaCaribou_latest.tar.gz](https://caribou.jnbrunet.com/artifacts/linux/v21.06/SofaCaribou_latest.tar.gz) | N/A  | N/A     |
| master       | [SofaCaribou_latest.tar.gz](https://caribou.jnbrunet.com/artifacts/linux/master/SofaCaribou_latest.tar.gz) | N/A  | N/A     |

To install them, simply expand the archive content directly into the `plugins` directory of your 
SOFA installation. For example, the following commands will install SOFA v21.06 and install Caribou plugin and its python bindings afterwards:

```console
$ # Set the environment variable that will contain the path to SOFA's installation
$ export SOFA_ROOT=$PWD/SOFA_V21.06

$ # Install SOFA V21.06 into SOFA_ROOT 
$ # (SKIP THIS IF YOU ALREADY HAVE INSTALLED SOFA)
$ curl --progress-bar --output sofa.zip -L "https://github.com/sofa-framework/sofa/releases/download/v21.06.00/SOFA_v21.06.00_Linux.zip"  && unzip -qq sofa.zip -d temp && mv temp/`ls temp` $SOFA_ROOT && rm -rf temp && rm -rf sofa.zip

$ # Install Caribou into the SOFA_ROOT/plugins directory
$ curl --progress-bar --output SofaCaribou.tar.gz -L "https://caribou.jnbrunet.com/artifacts/linux/v21.06/SofaCaribou_latest.tar.gz"
$ tar zxf SofaCaribou.tar.gz --directory $SOFA_ROOT/plugins

$ # Tell python where it should find SP3 and Caribou
$ export PYTHONPATH="$SOFA_ROOT/plugins/SofaPython3/lib/python3/site-packages:$SOFA_ROOT/plugins/SofaCaribou/lib/python3/site-packages"

$ # Test the installation
$ python3.7 -c "import SofaRuntime;print(SofaRuntime);import SofaCaribou;print(SofaCaribou)"
```

### Documentation

You can read the documentation for building and using Caribou at the following link:
[https://caribou.readthedocs.io](https://caribou.readthedocs.io)