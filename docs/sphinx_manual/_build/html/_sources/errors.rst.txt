Common Errors and How to Fix Them
==================================

**Undefined function** ``ColumnGrid``

The ``src/`` directory is not on the MATLAB path. Ensure that
``addpath(fullfile(script_dir, '..', '..', 'src'))`` appears at the top of
the script and that ``script_dir`` resolves to the correct location.

**Negative concentrations appear during the run**

The time step is too large for the current sinking speed or grid spacing.
Evaluate the CFL condition for the largest bin. With ``n_sections = 30``,
``dt = 0.25`` day, and ``dz = 50`` m, the CFL evaluates to 0.48 which is
safely stable. Reducing ``dt`` or coarsening the size grid will improve
stability.

**Transfer efficiency is unrealistically high (above 10%)**

Check ``n_sections``. With ``n = 20``, the largest model bin does not reach
the critical size for disaggregation in deep water, so particles pile up and
the deep flux is greatly overestimated. Always use ``n_sections = 30`` with
``dt = 0.25`` day for production runs.

**The spinup loop runs all 80 cycles without converging**

First, check that the flux boundary condition array ``phi_bc_daily`` contains
non-zero values. If the UVP data file is not found or returns an empty array,
the model column drains to zero. Second, check that ``enable_disagg = true``
and that ``disagg_dmax_A`` is set.

**Error:** ``operator_split`` **requires** ``disagg_dmax_cm`` **or** ``disagg_epsilon``

Set either ``cfg.disagg_dmax_cm = 1.0`` for a fixed :math:`D_\mathrm{max}`
(in cm) or ``cfg.disagg_dmax_A = 9.39e-6 * 5`` for the Parker
:math:`\varepsilon`-scaling. The error is issued by
``SimulationConfig.validate()``.
