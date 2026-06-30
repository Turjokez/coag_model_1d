Configuration Parameters
========================

All model options are controlled through the ``SimulationConfig`` class::

    cfg = SimulationConfig();
    cfg.n_sections     = 30;
    cfg.alpha          = 0.10;
    cfg.sinking_law    = 'kriest_8';
    cfg.enable_disagg  = true;
    sim = ColumnSimulation(cfg, col_grid, prof);

.. list-table:: Commonly used parameters in ``SimulationConfig``
   :widths: 28 18 40 14
   :header-rows: 1

   * - Parameter
     - Default
     - Description
     - Units
   * - **Size grid**
     -
     -
     -
   * - ``n_sections``
     - 20†
     - Number of size bins. Use 30 for production.
     - --
   * - ``d0``
     - 20×10\ :sup:`-4`
     - Diameter of the smallest particle.
     - cm
   * - ``fr_dim``
     - 2.33
     - Particle fractal dimension :math:`D_f`.
     - --
   * - ``r_to_rg``
     - 1.6
     - Collision radius to radius of gyration.
     - --
   * - **Physics switches**
     -
     -
     -
   * - ``enable_coag``
     - true
     - Activates coagulation.
     - bool
   * - ``enable_disagg``
     - false
     - Activates disaggregation.
     - bool
   * - ``enable_zoo``
     - false
     - Activates zooplankton grazing.
     - bool
   * - ``enable_microbe``
     - false
     - Activates microbial remineralization.
     - bool
   * - ``enable_mining``
     - false
     - Activates micro-zooplankton mining.
     - bool
   * - **Sinking**
     -
     -
     -
   * - ``sinking_law``
     - ``'current'``
     - Sinking speed law; use ``'kriest_8'``.
     - --
   * - ``ds_kernel_mode``
     - ``'sinking_law'``
     - DS kernel mode; always use ``'sinking_law'``.
     - --
   * - **Coagulation**
     -
     -
     -
   * - ``alpha``
     - 1.0
     - Stickiness :math:`\alpha`. Typical range 0.01–0.5.
     - --
   * - **Disaggregation**
     -
     -
     -
   * - ``disagg_mode``
     - ``'legacy'``
     - Use ``'operator_split'`` in production.
     - --
   * - ``disagg_dmax_A``
     - 9.39×10\ :sup:`-6`
     - Parker constant for :math:`D_\mathrm{max}`. Multiply by 5 for EXPORTS.
     - m
   * - **Zooplankton**
     -
     -
     -
   * - ``zoo_c``
     - 10\ :sup:`-4`
     - Filter feeder clearance rate.
     - m³ ind⁻¹ d⁻¹
   * - ``zoo_s``
     - 10\ :sup:`-4`
     - Flux feeder capture cross-section.
     - m² ind⁻¹
   * - ``zoo_p``
     - 0.3
     - Egestion fraction.
     - --
   * - ``zoo_ic``
     - 7
     - Fecal pellet target bin (bin 8, ~115 μm).
     - --
   * - **Microbial remineralization**
     -
     -
     -
   * - ``microbe_r0``
     - 0.03
     - Base remineralization rate :math:`r_0`.
     - day⁻¹
   * - ``microbe_use_temp``
     - false
     - Apply Q10 temperature scaling.
     - bool
   * - ``microbe_q10``
     - 2.0
     - Q10 factor.
     - --

† Defaults marked † differ from the class default in production runs.
