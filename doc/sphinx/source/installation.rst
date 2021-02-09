 .. _installation:

Installation
============

Prerequisites - Caribou
-----------------------

At the moment, the only way to install Caribou is by compiling it on your computer. You can look at the following table
for the dependencies that you will have to install before starting the compilation.

+------------------------------------------------------------------------------------------+--------------+-----------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| Package                                                                                  | Type         | Ubuntu                                  | Mac OSX                                                                                                                           | Description                                                         |
+==========================================================================================+==============+=========================================+===================================================================================================================================+=====================================================================+
| `Eigen <http://eigen.tuxfamily.org/dox/>`__                                              | **Required** | ``sudo apt install libeigen3-dev``      | ``brew install eigen``                                                                                                            | Used everywhere inside Caribou has the main linear algebra library. |
+------------------------------------------------------------------------------------------+--------------+-----------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| Python 3 libs                                                                            | Optional     | ``sudo apt install python3-dev``        | ``brew install python3``                                                                                                          | Used for the python bindings of Caribou objects.                    |
+------------------------------------------------------------------------------------------+--------------+-----------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| `pybind11 <https://pybind11.readthedocs.io/en/stable/>`__                                | Optional     | ``sudo apt install pybind11-dev``       | ``brew install pybind11``                                                                                                         | Used for the python bindings of Caribou objects.                    |
+------------------------------------------------------------------------------------------+--------------+-----------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| OpenMP                                                                                   | Optional     | ``sudo apt install libomp-dev``         | ``brew install libomp``                                                                                                           | Used to parallelize the computation of some components.             |
+------------------------------------------------------------------------------------------+--------------+-----------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| Intel MKL                                                                                | Optional     | ``sudo apt install libmkl-full-dev``    | `instructions here <https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library/choose-download/macos.html>`__ | Used by Eigen for optimization of some math operations.             |
+------------------------------------------------------------------------------------------+--------------+-----------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| Google test suite                                                                        | Optional     | ``sudo apt install libgtest-dev``       | `instructions here <https://stackoverflow.com/questions/15852631/how-to-install-gtest-on-mac-os-x-with-homebrew>`__               | Used for the unit tests of Caribou.                                 |
+------------------------------------------------------------------------------------------+--------------+-----------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------+
| `SOFA Framework <https://www.sofa-framework.org/community/doc/>`__                       | Optional     | See below for more information                                                                                                                                                                                                                    |
+------------------------------------------------------------------------------------------+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `SofaPython3 <https://github.com/sofa-framework/plugin.SofaPython3#pluginsofapython3>`__ | Optional     | See below for more information                                                                                                                                                                                                                    |
+------------------------------------------------------------------------------------------+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

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


Compiling
---------
All right, at this point you should have everything needed to compile Caribou. If you are also building SofaCaribou and
its python bindings, you also have defined the ``SOFA_INSTALL`` and ``SP3_INSTALL`` environment variables.

Start by cloning the Caribou source code and create a build directory inside of it.

.. code-block:: bash

    $ git clone https://github.com/jnbrunet/caribou.git
    $ cd caribou
    $ mkdir build
    $ cd build

Next, cmake will be use to configure the build option. It is used with the following format: ``cmake -DVAR=VALUE ..``
where **VAR** is the name of a configuration variable and **VALUE** is the value assigned to the variable. Caribou provides
the following configuration variables:

+-----------------------------+--------+----------+-------------------------------------------------------------------------------------------+
| Var                         | Value  | Default  | Description                                                                               |
+=============================+========+==========+===========================================================================================+
| CARIBOU_USE_FLOAT           | ON/OFF | OFF      | Specify if the floating point type should be float (OFF) or double(ON).                   |
+-----------------------------+--------+----------+-------------------------------------------------------------------------------------------+
| CARIBOU_BUILD_TESTS         | ON/OFF | OFF      | Whether or not the test suite of Caribou should be build.                                 |
+-----------------------------+--------+----------+-------------------------------------------------------------------------------------------+
| CARIBOU_WITH_SOFA           | ON/OFF | ON       | Compile the Caribou's SOFA plugin (SofaCaribou).                                          |
+-----------------------------+--------+----------+-------------------------------------------------------------------------------------------+
| CARIBOU_OPTIMIZE_FOR_NATIVE | ON/OFF | ON       | Tell the compiler to optimize Caribou following the architecture of your computer.        |
+-----------------------------+--------+----------+-------------------------------------------------------------------------------------------+
| CARIBOU_WITH_PYTHON_3       | ON/OFF | ON       | Compile Caribou's python bindings.                                                        |
+-----------------------------+--------+----------+-------------------------------------------------------------------------------------------+
| CARIBOU_WITH_MKL            | ON/OFF | ON       | Compile Caribou with Intel® Math Kernel Library (MKL) support.                            |
+-----------------------------+--------+----------+-------------------------------------------------------------------------------------------+
| CARIBOU_WITH_OPENMP         | ON/OFF | ON       | Compile Caribou with OpenMP support.                                                      |
+-----------------------------+--------+----------+-------------------------------------------------------------------------------------------+
| CMAKE_INSTALL_PREFIX        | Path   | install/ | Specify where the built files (following the `make install` command) should be installed. |
+-----------------------------+--------+----------+-------------------------------------------------------------------------------------------+

If you are compiling the Caribou's SOFA plugin, you will also need to tell cmake where it should find it. This can be
done by setting the cmake variable ``CMAKE_PREFIX_PATH`` to ``$SOFA_INSTALL/lib/install/cmake``. The same
thing needs to be done with SofaPython3 if you are also compiling Caribou's python bindings. In this case, you can set
``CMAKE_PREFIX_PATH`` to ``$SOFA_INSTALL/lib/install/cmake;$SP3_INSTALL/lib/install/cmake`` (note the semicolon ``;``
between the two paths).

For example, if you want to compile Caribou with MKL support and python bindings:

.. code-block:: bash

    $ cmake -DCARIBOU_WITH_MKL=ON -DCARIBOU_WITH_PYTHON_3=ON ..

If you want to compile Caribou with SOFA and python bindings:

.. code-block:: bash

    $ cmake -DCARIBOU_WITH_PYTHON_3=ON -DCMAKE_PREFIX_PATH="$SOFA_INSTALL/lib/cmake;$SP3_INSTALL/lib/cmake" ..

You can now start the compilation.

.. code-block:: bash

    $ cmake --build . -j4
    $ cmake --install .

The last command (``cmake --install .``) installed all the built files inside the directory ``install`` (or the directory
specified by the cmake variable ``CMAKE_INSTALL_PREFIX`` if you changed it). Export this path to the environment variable
``CARIBOU_INSTALL``:

.. code-block:: bash

    $ export CARIBOU_INSTALL="${PWD}/install"

.. note::

   To make sure your ``CARIBOU_INSTALL`` is well defined, you can verify that the following file path exists:

   .. code-block:: bash

        $CARIBOU_INSTALL/lib/cmake/Caribou/CaribouTargets.cmake


Installing python bindings
--------------------------

If you compiled the Caribou's python bindings, and you want them to be found automatically by your python scripts,
you can create a symbolic link to the binding directories inside Python's site-package path:

For linux, this can be done with the following command:

.. code-block:: bash

    $ ln -sFfv $(find $CARIBOU_INSTALL/lib/python3.7/site-packages -maxdepth 1 -mindepth 1) $(python3 -m site --user-site)

And for Mac OSX:

.. code-block:: bash

    $ ln -sFfv $(find $CARIBOU_INSTALL/lib/python3.7/site-packages -d 1) $(python3 -m site --user-site)

You can test that the bindings have been correctly installed by starting a python shell and import Caribou:

.. code-block:: python

    import Caribou

    # Do the following only if you compiled the Caribou's SOFA plugin
    import SofaRuntime
    import SofaCaribou
