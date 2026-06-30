Model Physics and Biology
==========================

This section describes the governing equations for each process in the model,
the physical or biological motivation for including each process, and the
meaning of every parameter that appears.

Size Sectional Method
---------------------

The particle size distribution spans roughly four orders of magnitude in
diameter. The model divides the size range into discrete size bins (sections),
with bins spaced logarithmically so that each successive bin spans twice the
mass of the previous one:

.. math::
   m_{i+1} = 2\,m_i.

The representative diameter of bin :math:`i` is

.. math::
   d_i = d_0 \cdot 2^{(i-1)/3},

where :math:`d_0 = 20\,\mu\mathrm{m}` is the unit particle diameter. With 30 sections,
bin diameters range from approximately 20 μm (bin 1) to 10 mm (bin 30).

Marine aggregates are fractal objects. The fractal dimension :math:`D_f`
relates the number of primary particles :math:`N` in an aggregate to its
outer radius :math:`r_g`:

.. math::
   N \propto \left(\frac{r_g}{r_0}\right)^{D_f}.

The default value :math:`D_f = 2.33` is consistent with Jackson & Burd (1998).

Gravitational Settling
----------------------

For a solid sphere in a viscous fluid at low Reynolds number, the Stokes
terminal velocity is

.. math::
   w_\mathrm{Stokes}(d) = \frac{(\rho_p - \rho_f)\,g\,d^2}{18\,\mu}.

Marine aggregates are porous and fractal, so their sinking speed follows the
empirical Kriest (2002) power law:

.. math::
   w(d) = 66\,d_\mathrm{cm}^{0.62} \quad [\mathrm{m\,day^{-1}}],

where :math:`d_\mathrm{cm}` is diameter in centimeters. This is the
``kriest_8`` law used by default. The exponent 0.62 is much smaller than
the Stokes exponent of 2, reflecting that large porous aggregates sink more
slowly than solid spheres.

Fecal pellets use Stokes law with excess density
:math:`\Delta\rho = 0.15\,\mathrm{g\,cm^{-3}}`:

.. math::
   w_\mathrm{fp}(d) = \frac{\Delta\rho\,g\,d^2}{18\,\mu}.

For a 115 μm pellet (bin 8), this gives :math:`w_\mathrm{fp} \approx 69`
m day\ :sup:`-1`, roughly 16 times faster than marine snow at the same size.

Coagulation
-----------

The Smoluchowski equation gives the rate of change of biovolume in bin
:math:`k` due to coagulation:

.. math::
   \left.\frac{d\phi_k}{dt}\right|_\mathrm{coag}
   = \frac{1}{2}\sum_{i+j \to k} \alpha\,\beta_{ij}\,n_i\,n_j\,v_k
   - \phi_k \sum_{j=1}^{N} \alpha\,\beta_{kj}\,n_j,

where :math:`n_i = \phi_i / v_i` is number concentration, :math:`v_i` is
bin volume, and :math:`\alpha` is the stickiness coefficient. Three mechanisms
drive collisions:

**Brownian motion** (dominant below ~1 μm):

.. math::
   \beta_{ij}^\mathrm{Br}
   = \frac{2\,k_B\,T}{3\,\mu}
     \left(\frac{1}{r_i} + \frac{1}{r_j}\right)(r_i + r_j).

**Fluid shear** (dominant 10–100 μm):

.. math::
   \beta_{ij}^\mathrm{sh}
   = \frac{4}{3}\,\dot{\gamma}\,(r_i + r_j)^3,
   \qquad \dot{\gamma} = \left(\frac{\varepsilon}{\nu}\right)^{1/2}.

**Differential settling** (dominant above ~100 μm):

.. math::
   \beta_{ij}^\mathrm{DS}
   = \pi\,(r_i^\mathrm{eff} + r_j^\mathrm{eff})^2\,|w_i - w_j|.

The total kernel is :math:`\beta_{ij} = \beta_{ij}^\mathrm{Br} + \beta_{ij}^\mathrm{sh} + \beta_{ij}^\mathrm{DS}`.

Disaggregation
--------------

Turbulence limits the maximum stable aggregate size (Parker et al. 1972):

.. math::
   D_\mathrm{max}(\varepsilon) = A\,\varepsilon^{-1/2},

where :math:`A = 9.39 \times 10^{-6}` m (multiplied by 5 for EXPORTS-NA).
Particles exceeding :math:`D_\mathrm{max}` are fragmented in an operator
split, conserving total biovolume:

.. math::
   \Delta\phi_{i-1} = +\tfrac{2}{3}\,\phi_i^\mathrm{excess},\quad
   \Delta\phi_{i-2} = +\tfrac{1}{3}\,\phi_i^\mathrm{excess},\quad
   \Delta\phi_i     = -\phi_i^\mathrm{excess}.

One-Dimensional Vertical Transport
-----------------------------------

Each size bin is transported downward by its sinking speed:

.. math::
   \frac{\partial \phi_i}{\partial t}
   + \frac{\partial (w_i\,\phi_i)}{\partial z}
   = S_i(z, t).

Discretized with the first-order upwind scheme:

.. math::
   \phi_i^{k,\,n+1}
   = \phi_i^{k,\,n}
   - \frac{w_i\,\Delta t}{\Delta z}\bigl(\phi_i^{k,\,n} - \phi_i^{k-1,\,n}\bigr)
   + \Delta t\,S_i^{k,\,n}.

The CFL stability condition :math:`w_\mathrm{max}\,\Delta t / \Delta z < 1`
evaluates to :math:`96 \times 0.25 / 50 = 0.48` with the default settings,
safely below the stability limit.

Zooplankton Grazing
-------------------

Following Stemmann et al. (2004), two feeding strategies are modelled.
Filter feeders remove particles at a rate proportional to concentration:

.. math::
   \left.\frac{d\phi_i}{dt}\right|_\mathrm{ff} = -c\,Z_c(z)\,\phi_i.

Flux feeders intercept sinking particles; their rate scales with sinking speed:

.. math::
   \left.\frac{d\phi_i}{dt}\right|_\mathrm{sf} = -s\,w_i\,Z_f(z)\,\phi_i.

A fraction :math:`p` of ingested biovolume is egested as fecal pellets into
bin 8 (~115 μm).

Microbial Remineralization
--------------------------

Applied as an operator-split exact exponential decay after each time step:

.. math::
   \phi_i \leftarrow \phi_i \cdot e^{-r\,\Delta t}.

With optional Q10 temperature scaling:

.. math::
   r = r_0 \cdot Q_{10}^{(T - T_\mathrm{ref})/10},

where :math:`Q_{10} = 2.0` and :math:`T_\mathrm{ref} = 20^\circ\mathrm{C}`.

Micro-Zooplankton Mining
------------------------

Small copepods (*Oncaea*) bite pieces from large aggregates (Stemmann et al.
2004, Eq. 25):

.. math::
   \left.\frac{d\phi_i}{dt}\right|_\mathrm{mine}
   = -s_m\,w_i\,Z_m\,\phi_i\;\mathbf{1}(d_i \geq d_\mathrm{min}),

where :math:`Z_m = 250` ind m\ :sup:`-3` and mining is restricted to bins
≥ 12 (d ≈ 254 μm).
