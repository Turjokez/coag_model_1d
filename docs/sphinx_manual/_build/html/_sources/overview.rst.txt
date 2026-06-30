Overview of the Model
=====================

This manual describes a one-dimensional numerical model for the transport,
aggregation, and fragmentation of marine particles in the ocean water column.
The model is designed to study the biological carbon pump: the process by
which particles formed at the ocean surface sink to depth, exporting carbon
from the upper ocean.

The state variable is the biovolume concentration :math:`\phi(z, d, t)`
[m\ :sup:`3` m\ :sup:`-3`], representing the volume of particles per unit
volume of seawater at depth :math:`z`, particle diameter :math:`d`, and
time :math:`t`. Biovolume is used rather than mass because it is the quantity
measured directly by the Underwater Vision Profiler (UVP) and because it is
exactly conserved when two particles coagulate.

The model domain is a single vertical column of ocean, 1000 m deep,
discretized into 20 layers of 50 m thickness each. The particle size
distribution is resolved into 30 logarithmically spaced size bins spanning
approximately 2 μm to 10 mm. Time integration uses a first-order upwind
scheme with a time step of :math:`\Delta t = 0.25` day. The model is
implemented in MATLAB using a set of object-oriented classes described in
:doc:`structure`.

The physical and biological processes represented are gravitational settling,
coagulation (aggregation), turbulent disaggregation (breakup), zooplankton
grazing, fecal pellet production and transport, microbial remineralization,
and micro-zooplankton mining of large aggregates. Each process is described
in detail in :doc:`physics`.
