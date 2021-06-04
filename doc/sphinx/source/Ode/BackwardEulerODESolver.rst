 .. _backward_euler_ode_doc:
 .. role:: important

<BackwardEulerODESolver />
==========================

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::ode::BackwardEulerODESolver`

Implementation of an implicit backward euler solver compatible with non-linear materials.

We are trying to solve to following

.. math::
    \boldsymbol{M} \ddot{\boldsymbol{x}} + \boldsymbol{C} \dot{\boldsymbol{x}} + \boldsymbol{R}(\boldsymbol{x}) = \boldsymbol{P}

where :math:`\boldsymbol{M}` is the mass matrix, :math:`\boldsymbol{C}` is the damping matrix, :math:`\boldsymbol{R}`
is the (possibly non-linear) internal elastic force residual and :math:`\boldsymbol{P}` is the external
force vector (for example, gravitation force or surface traction).

We first transform this second-order differential equation to a first one by introducing two independant variables:

.. math::
    \boldsymbol{v} &= \dot{\boldsymbol{x}} \\
    \boldsymbol{a} &= \ddot{\boldsymbol{x}}

Using the `backward Euler scheme <https://en.wikipedia.org/wiki/Backward_Euler_method>`_, we pose the following approximations:

.. math::
     \boldsymbol{x}_{n+1} &= \boldsymbol{x}_{n} + h \boldsymbol{v}_{n+1} ~~~~ (1) \\
     \boldsymbol{v}_{n+1} &= \boldsymbol{v}_{n} + h \boldsymbol{a}_{n+1} ~~~~ (2)

where :math:`h` is the delta time between the steps :math:`n` and :math:`n+1`.

Substituting :math:`(2)` inside :math:`(1)` gives

.. math::
    \boldsymbol{x}_{n+1} &= \boldsymbol{x}_{n} + h \left[ \boldsymbol{v}_{n} + h \boldsymbol{a}_{n+1} \right] \\
                   &= \boldsymbol{x}_{n} + h \boldsymbol{v}_{n} + h^2 \boldsymbol{a}_{n+1}

And the problem becomes:

.. math::
    \boldsymbol{M} \boldsymbol{a}_{n+1} + \boldsymbol{C} \left[ \boldsymbol{v}_{n} + h \boldsymbol{a}_{n+1} \right] + \boldsymbol{R}(\boldsymbol{x}_{n} + h \boldsymbol{v}_{n} + h^2 \boldsymbol{a}_{n+1}) = \boldsymbol{P}_n

where :math:`\boldsymbol{a}_{n+1}` is the vector of unknown accelerations.

Finally, assuming  :math:`\boldsymbol{R}` is non-linear in :math:`\boldsymbol{x}_{n+1}`, we iteratively solve for :math:`\boldsymbol{a}_{n+1}`
using the `Newton-Raphson method <https://en.wikipedia.org/wiki/Newton's_method#Nonlinear_systems_of_equations>`_ and
using the previous approximations to back propagate it inside the velocity and position vectors.

Let :math:`i` be the Newton iteration number for a given time step :math:`n`, we pose

.. math::
     \boldsymbol{F}(\boldsymbol{a}_{n+1}^i) &= \boldsymbol{M} \boldsymbol{a}_{n+1}^i + \boldsymbol{C} \left[ \boldsymbol{v}_{n} + h \boldsymbol{a}_{n+1}^i \right] + \boldsymbol{R}(\boldsymbol{x}_{n} + h \boldsymbol{v}_{n} + h^2 \boldsymbol{a}_{n+1}^i) - \boldsymbol{P}_n \\
     \boldsymbol{J} = \frac{\partial \boldsymbol{F}}{\partial \boldsymbol{a}_{n+1}} \bigg\rvert_{\boldsymbol{a}_{n+1}^i} &= \boldsymbol{M} + h \boldsymbol{C} + h^2 \boldsymbol{K}(\boldsymbol{a}_{n+1}^i)

where :math:`\boldsymbol{x}_{n}` and :math:`\boldsymbol{x}_{n}` are the position and velocity vectors at the beginning of the time
step, respectively, and remains constant throughout the Newton iterations.

We then solve for :math:`\boldsymbol{a}_{n+1}^{i+1}` with

.. math::
     \boldsymbol{J} \left [ \Delta \boldsymbol{a}_{n+1}^{i+1} \right ] &= - \boldsymbol{F}(\boldsymbol{a}_{n+1}^i) \\
     \boldsymbol{a}_{n+1}^{i+1} &= \boldsymbol{a}_{n+1}^{i} + \Delta \boldsymbol{a}_{n+1}^{i+1}

And propagate back the new acceleration using :math:`(1)` and :math:`(2)`.

In addition, this component implicitly adds a Rayleigh's damping matrix :math:`\boldsymbol{C}_r = r_m \boldsymbol{M} + r_k \boldsymbol{K}(\boldsymbol{x}_{n+1})`.
We therefore have

.. math::
     \boldsymbol{F}(\boldsymbol{a}_{n+1}^i) &= \boldsymbol{M} \boldsymbol{a}_{n+1}^i + \boldsymbol{C} \left[ \boldsymbol{v}_{n} + h \boldsymbol{a}_{n+1}^i \right] + \boldsymbol{R}(\boldsymbol{x}_{n} + h \boldsymbol{v}_{n} + h^2 \boldsymbol{a}_{n+1}^i) - \boldsymbol{P}_n \\
                                &= \boldsymbol{M} \boldsymbol{a}_{n+1}^i + (\boldsymbol{C}_r+\boldsymbol{C}) \left[ \boldsymbol{v}_{n} + h \boldsymbol{a}_{n+1}^i \right] + \boldsymbol{R}(\boldsymbol{x}_{n} + h \boldsymbol{v}_{n} + h^2 \boldsymbol{a}_{n+1}^i) - \boldsymbol{P}_n \\
                                &= \boldsymbol{M} \boldsymbol{a}_{n+1}^i + (r_m\boldsymbol{M}+r_k\boldsymbol{K}) \left[ \boldsymbol{v}_{n} + h \boldsymbol{a}_{n+1}^i \right] + \boldsymbol{C} \left[ \boldsymbol{v}_{n} + h \boldsymbol{a}_{n+1}^i \right] + \boldsymbol{R}(\boldsymbol{x}_{n} + h \boldsymbol{v}_{n} + h^2 \boldsymbol{a}_{n+1}^i) - \boldsymbol{P}_n \\
                                &= \left[ (1 + hr_m)\boldsymbol{M} + h\boldsymbol{C} + hr_k\boldsymbol{K} \right] \boldsymbol{a}_{n+1}^i
                                 + \left[ r_m\boldsymbol{M} + \boldsymbol{C} + r_k \boldsymbol{K} \right] \boldsymbol{v}_n
                                 + \left[ \boldsymbol{R}(\boldsymbol{x}_{n} + h \boldsymbol{v}_{n} + h^2 \boldsymbol{a}_{n+1}^i) - \boldsymbol{P}_n \right] \\
     \boldsymbol{J} = \frac{\partial \boldsymbol{F}}{\partial \boldsymbol{a}_{n+1}} \bigg\rvert_{\boldsymbol{a}_{n+1}^i} &= (1 + hr_m)\boldsymbol{M} + h \boldsymbol{C} + h(h+r_k) \boldsymbol{K}(\boldsymbol{a}_{n+1}^i)

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
    * - rayleigh_stiffness
      - double
      - 0.0
      - The stiffness factor :math:`r_k` used in the Rayleigh's damping matrix :math:`\boldsymbol{D} = r_m \boldsymbol{M} + r_k \boldsymbol{K}`.
    * - rayleigh_mass
      - double
      - 0.0
      - The mass factor :math:`r_m` used in the Rayleigh's damping matrix :math:`\boldsymbol{D} = r_m \boldsymbol{M} + r_k \boldsymbol{K}`.
    * - newton_iterations
      - int
      - 1
      - Number of newton iterations between each load increments (normally, one load increment per simulation time-step).
    * - correction_tolerance_threshold
      - float
      - 1e-5
      - Relative convergence criterion: The newton iterations will stop when the norm of correction
        :math:`\frac{|\delta \boldsymbol{u}^{k}|}{\sum_{i=0}^k|\delta \boldsymbol{u}^{i}|}` reaches this threshold.
    * - residual_tolerance_threshold
      - float
      - 1e-5
      - Relative convergence criterion: The newton iterations will stop when the relative norm of the residual
        :math:`\frac{|\boldsymbol{R}_k|}{|\boldsymbol{R}_0|}` at iteration k is lower than this threshold.
        Use a negative value to disable this criterion.
    * - absolute_residual_tolerance_threshold
      - float
      - 1e-15
      - Absolute convergence criterion: The newton iterations will stop when the absolute norm of the residual
        :math:`|\boldsymbol{R}_k|` at iteration k is lower than this threshold. This criterion is also used to
        detect the absence of external forces and skip useless Newton iterations.
        Use a negative value to disable this criterion.
    * - pattern_analysis_strategy
      - option
      - BEGINNING_OF_THE_TIME_STEP
      - Define when the pattern of the system matrix should be analyzed to extract a permutation matrix. If the sparsity and
        location of the coefficients of the system matrix doesn't change much during the simulation, then this analysis can
        be avoided altogether, or computed only one time at the beginning of the simulation. Else, it can be done at the
        beginning of the time step, or even at each reformation of the system matrix if necessary. The default is to
        analyze the pattern at each time step.

        **Options:**
            * NEVER
            * BEGINNING_OF_THE_SIMULATION
            * BEGINNING_OF_THE_TIME_STEP **(default)**
            * ALWAYS
    * - linear_solver
      - LinearSolver
      - None
      - Linear solver used for the resolution of the system. Will be automatically found in the current context node if
        none is supplied.
    * - converged
      - bool
      - N/A
      - Whether or not the last call to solve converged.

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <BackwardEulerODESolver rayleigh_stiffness="0.1" rayleigh_mass="0.1" newton_iterations="10" correction_tolerance_threshold="1e-8" residual_tolerance_threshold="1e-8" printLog="1" />
                <LLTSolver backend="Pardiso" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject('BackwardEulerODESolver', rayleigh_stiffness=0.1, rayleigh_mass=0.1, newton_iterations=10, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
            node.addObject('LLTSolver', backend='Pardiso')


Available python bindings
*************************

None at the moment.
