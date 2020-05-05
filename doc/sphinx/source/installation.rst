 .. _installation:

Installation
============

Prerequisites - Caribou
-----------------------

At the moment, the only way to install Caribou is by compiling it on your computer. You can look at the following table
for the dependencies that you will have to install before starting the compilation.

+------------------------------------------------------------------------------------------+--------------+------------------------------------+--------------------------------+---------------------------------------------------------------------+
| Package                                                                                  | Type         | Ubuntu                             | Mac OSX                        | Description                                                         |
+------------------------------------------------------------------------------------------+--------------+------------------------------------+--------------------------------+---------------------------------------------------------------------+
| `Eigen <http://eigen.tuxfamily.org/dox/>`__                                              | **Required** | ``sudo apt install libeigen3-dev`` | ``brew install eigen``         | Used everywhere inside Caribou has the main linear algebra library. |
+------------------------------------------------------------------------------------------+--------------+------------------------------------+--------------------------------+---------------------------------------------------------------------+
| Python 3 libs                                                                            | Optional     | ``sudo apt install python3-dev``   | ``brew install python3``       | Used for the python bindings of Caribou objects.                    |
+------------------------------------------------------------------------------------------+--------------+------------------------------------+--------------------------------+---------------------------------------------------------------------+
| `pybind11 <https://pybind11.readthedocs.io/en/stable/>`__                                | Optional     | ``sudo apt install pybind11-dev``  | ``brew install pybind11``      | Used for the python bindings of Caribou objects.                    |
+------------------------------------------------------------------------------------------+--------------+------------------------------------+--------------------------------+---------------------------------------------------------------------+
| OpenMP                                                                                   | Optional     | ``sudo apt install libomp-dev``    | ``brew install libomp``        | Used to parallelize the computation of some components.             |
+------------------------------------------------------------------------------------------+--------------+------------------------------------+--------------------------------+---------------------------------------------------------------------+
| Intel MKL                                                                                | Optional     | ``sudo apt install libmkl-dev``    | ``conda install -c intel mkl`` | Used by Eigen for optimization of some math operations.             |
+------------------------------------------------------------------------------------------+--------------+------------------------------------+--------------------------------+---------------------------------------------------------------------+
| `SOFA Framework <https://www.sofa-framework.org/community/doc/>`__                       | Optional     | See below for more information                                                                                                            |
+------------------------------------------------------------------------------------------+--------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| `SofaPython3 <https://github.com/sofa-framework/plugin.SofaPython3#pluginsofapython3>`__ | Optional     | See below for more information                                                                                                            |
+------------------------------------------------------------------------------------------+--------------+-------------------------------------------------------------------------------------------------------------------------------------------+

Prerequisites - Sofa
--------------------

If you want to compile the Caribou's SOFA plugin (alias **SofaCaribou**), you will need to either:

* [**Method 1**] Install the `SOFA binaries and headers <https://www.sofa-framework.org/download/>`__
  somewhere on your computer and note its installation directory. Or,
* [**Method 2**] Compile it following the `SOFA build documentation <https://www.sofa-framework.org/community/doc/getting-started/build/linux/>`__.
  Once it is built, execute the installation by going into the build directory of SOFA (for example,
  */home/user/sofa/build/master/*), and using the command ``cmake --install .``

Once done, export the installation path of SOFA inside the ``SOFA_INSTALL`` environment variable. For example,

.. code-block:: bash

    $ export SOFA_INSTALL="/home/user/sofa/build/master/install"


.. note::

   To make sure your ``SOFA_INSTALL`` is well defined, you can verify that the following file path exists:

   .. code-block:: bash

        $SOFA_INSTALL/lib/cmake/SofaFramework/SofaFrameworkTargets.cmake

Prerequisites - SofaPython3
---------------------------
If you want the **SofaCaribou** python bindings, you will also need to compile the `SofaPython3 plugin <https://github.com/sofa-framework/plugin.SofaPython3>`__.
This plugin allows you to compile it using two different methods.

* `"in-tree" build type <https://github.com/sofa-framework/plugin.SofaPython3#in-tree-build>`__, which is, building the plugin
  along with the compilation of SOFA using the ``CMAKE_EXTERNAL_DIRECTORIES`` CMake variable of SOFA. This means that when
  you compile SOFA, you also compile the SofaPython3 plugin at the same time, using the same build directory. The plugin
  binaries will be found at the same place as the SOFA ones.
* `"out-of-tree" build type <https://github.com/sofa-framework/plugin.SofaPython3#out-of-tree-build>`__, which is,
  building the plugin in its own build directory, outside of SOFA.

If you used `"in-tree" build type <https://github.com/sofa-framework/plugin.SofaPython3#in-tree-build>`__, nothing more has to be done.
The installation files of the SofaPython3 plugin have been installed alongside the SOFA ones.

If instead you used the `"out-of-tree" build type <https://github.com/sofa-framework/plugin.SofaPython3#out-of-tree-build>`__,
you will need to install the built files by using the command ``cmake --install .``
in the build directory of the plugin (similarly to what you have done with SOFA in the last section).

Once done, export the installation path of SofaPython3 inside the ``SP3_INSTALL`` environment variable. For example, for
an out-of-tree build in the */home/user/plugin.SofaPython3* directory:

.. code-block:: bash

    $ export SP3_INSTALL="/home/user/plugin.SofaPython3/build/master/install"

For an "in-tree" build type, the ``SP3_INSTALL`` environment variable will be exactly the same as the ``SOFA_INSTALL``
environment variable defined earlier.


.. note::

   To make sure your ``SP3_INSTALL`` is well defined, you can verify that the following file path exists:

   .. code-block:: bash

        $SP3_INSTALL/lib/cmake/SofaPython3/PluginTargets.cmake


Compiling Caribou
-----------------
All right, at this point you should have everything needed to compile Caribou. If you are also building SofaCaribou and
its python bindings, you also have defined the ``SOFA_INSTALL`` and ``SP3_INSTALL`` environment variables.