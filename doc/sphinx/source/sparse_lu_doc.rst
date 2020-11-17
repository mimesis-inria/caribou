.. _sparse_lu_doc:
.. role:: important
.. role:: warning

<LUSolver />
==================

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::solver::LUSolver`

Implementation of a sparse LU linear solver.


.. list-table::
    :widths: 1 1 1 100
    :header-rows: 1
    :stub-columns: 0

    * - Attribute
      - Format
      - Default
      - Description
    * - backend
      - option
      - Eigen
      - Solver backend to use.
            * **Eigen**
                | Eigen LU solver (SparseLU).
                | **[default]**

            * **Pardiso**
                Pardiso LU solver. :warning:`Requires Intel® Math Kernel Library (MKL) installed.`
    * - symmetric
      - bool
      - False
      - Allows to explicitly state that the system matrix will be symmetric. This will in turn enable various optimizations.
        This option is only used by the Eigen backend.

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <StaticODESolver newton_iterations="10" correction_tolerance_threshold="1e-8" residual_tolerance_threshold="1e-8" printLog="1" />
                <LUSolver backend="Pardiso" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject('StaticODESolver', newton_iterations=10, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
            node.addObject('LUSolver', backend="Pardiso")


Available python bindings
*************************

None at the moment.
