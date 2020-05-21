.. _cg_solver_doc:
.. role:: important
.. role:: warning

<ConjugateGradientSolver />
===========================

Implementation of a Conjugate Gradient linear solver for selfadjoint (symmetric) positive-definite matrices.

.. list-table::
    :widths: 1 1 1 100
    :header-rows: 1
    :stub-columns: 0

    * - Attribute
      - Format
      - Default
      - Description
    * - printLog
      - bool
      - false
      - Output informative messages at the initialization and during the simulation.
    * - verbose
      - bool
      - false
      - Output convergence status at each iterations of the CG.
    * - maximum_number_of_iterations
      - int
      - 25
      - Maximum number of iterations before diverging.
    * - residual_tolerance_threshold
      - float
      - 1e-5
      - Convergence criterion: The CG iterations will stop when the residual norm of the residual
        :math:`\frac{|r_{k}|}{|r_0|} = \frac{|r_{k} - a_k A p_k|}{|b|}` at iteration k is lower than
        this threshold (here :math:`b` is the right-hand side vector).
    * - preconditioning_method
      - option
      - None
      - Preconditioning method used.

            * **None**: No preconditioning, hence the complete matrix won't be built. **(default)**
            * **Identity**: A naive preconditioner which approximates any matrix as the identity matrix. This is
              equivalent as using **None**, except the system matrix is built.
            * **Diagonal**: Preconditioning using an approximation of A.x = b by ignoring all off-diagonal entries of A.
              Also called Jacobi preconditioner, work very well on diagonal dominant matrices.
              This preconditioner is suitable for both selfadjoint and general problems. The diagonal entries are pre-inverted and stored into a dense vector.
              See `here <https://eigen.tuxfamily.org/dox/classEigen_1_1DiagonalPreconditioner.html>`__ for more details.
            * **IncompleteCholesky**: Preconditioning based on an modified incomplete Cholesky with dual threshold.
              See `here <https://eigen.tuxfamily.org/dox/classEigen_1_1IncompleteCholesky.html>`__ for more details.
            * **IncompleteLU**: Preconditioning based on the incomplete LU factorization.
              See `here <https://eigen.tuxfamily.org/dox/classEigen_1_1IncompleteLUT.html>`__ for more details.

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <StaticODESolver newton_iterations="10" correction_tolerance_threshold="1e-8" residual_tolerance_threshold="1e-8" printLog="1" />
                <ConjugateGradientSolver maximum_number_of_iterations="2500" residual_tolerance_threshold="1e-12" preconditioning_method="Diagonal" printLog="0" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject('StaticODESolver', newton_iterations=10, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
            node.addObject('ConjugateGradientSolver', maximum_number_of_iterations=2500, residual_tolerance_threshold=1e-12, preconditioning_method="Diagonal", printLog=False)


Available python bindings
*************************

.. py:class:: ConjugateGradientSolver

    .. py:function:: K()

        :return: Reference to the system matrix
        :rtype: :class:`scipy.sparse.csc_matrix`
        :note: No copy involved.

        Get the system matrix A = (mM + bB + kK) as a compressed sparse column major matrix.

    .. py:function:: assemble(m, b, k)

        :param m: Factor for the mass (M) matrix.
        :type m: float
        :param b: Factor for the damping (b) matrix.
        :type b: float
        :param k: Factor for the stiffness (K) matrix.
        :type k: float
        :return: Reference to the system matrix
        :rtype: :class:`scipy.sparse.csc_matrix`
        :note: No copy involved.

        Get the system matrix A = (mM + bB + kK) as a compressed sparse column major matrix.