.. _sparse_llt_doc:
.. role:: important
.. role:: warning

<SparseLLTSolver />
===================

Implementation of a sparse Cholesky (:math:`LL^T`) linear solver.


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
                | Eigen LLT solver (SimplicialLLT).
                | **[default]**

            * **Pardiso**
                Pardiso LLT solver. :warning:`Requires Intel® Math Kernel Library (MKL) installed.`

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <StaticODESolver newton_iterations="10" correction_tolerance_threshold="1e-8" residual_tolerance_threshold="1e-8" printLog="1" />
                <SparseLLTSolver backend="Pardiso" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject('StaticODESolver', newton_iterations=10, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
            node.addObject('SparseLLTSolver', backend="Pardiso")


Available python bindings
*************************

None at the moment.
