.. _sparse_llt_doc:
.. role:: important
.. role:: warning

<LLTSolver />
===================

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::solver::LLTSolver`

Implementation of a sparse Cholesky (:math:`LL^T`) direct linear solver.

This class provides a :math:`LL^T` Cholesky factorizations of sparse matrices that are selfadjoint and positive definite.
In order to reduce the fill-in, a symmetric permutation P is applied prior to the factorization such that the
factorized matrix is :math:`P A P^{-1}`.

The component uses the Eigen SimplicialLLT class as the solver backend.


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
                Pardiso LLT solver.

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <StaticODESolver newton_iterations="10" correction_tolerance_threshold="1e-8" residual_tolerance_threshold="1e-8" printLog="1" />
                <LLTSolver backend="Pardiso" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject('StaticODESolver', newton_iterations=10, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
            node.addObject('LLTSolver', backend="Pardiso")


Available python bindings
*************************

None at the moment.
