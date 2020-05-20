 .. _static_ode_doc:
 .. role:: important

<StaticODESolver />
===================

Implementation of a Newton-Raphson static ODE solver.


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
    * - newton_iterations
      - int
      - 1
      - Number of newton iterations between each load increments (normally, one load increment per simulation time-step.
    * - correction_tolerance_threshold
      - float
      - 1e-5
      - Convergence criterion: The newton iterations will stop when the norm of correction \|du| reach this threshold.
    * - residual_tolerance_threshold
      - float
      - 1e-5
      - Convergence criterion: The newton iterations will stop when the relative norm of the residual
        :math:`\frac{|R_k|}{|R_0|} = \frac{|f_k - Ku_k|}{|f_0 - Ku_0|}` at iteration k is lower than this threshold.
        Use a negative value to disable this
    * - shoud_diverge_when_residual_is_growing
      - bool
      - false
      - Divergence criterion: The newton iterations will stop when the residual is greater than the one from the
        previous iteration.
    * - warm_start
      - bool
      - false
      - For iterative linear solvers, use the previous solution has a warm start. Note that for the first newton step,
        the current position is used as the warm start.

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

.. py:class:: StaticODESolver

    :var iteration_times: List of times (in nanoseconds) that each Newton-Raphson iteration took to compute in the last call to Solve().
    :vartype iteration_times: list [int]

    :var squared_residuals: The list of squared residual norms (:math:`|r|^2`) of every newton iterations of the last solve call.
    :vartype squared_residuals: list [:class:`numpy.double`]

    :var squared_initial_residual: The initial squared residual (:math:`|r_0|^2`) of the last solve call.
    :vartype squared_initial_residual: :class:`numpy.double`

    :var iterative_linear_solver_squared_residuals: The list of squared residual norms (:math:`|r|^2`) of every iterative linear solver iterations, for each newton iterations of the last solve call.
    :vartype iterative_linear_solver_squared_residuals: list [ list [:class:`numpy.double`] ]

    :var iterative_linear_solver_squared_rhs_norms: List of squared right-hand side norms (:math:`|b|^2`) of every newton iterations before the call to the solve method of the iterative linear solver.
    :vartype iterative_linear_solver_squared_rhs_norms: list [:class:`numpy.double`]


