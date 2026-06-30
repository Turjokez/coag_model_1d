Running a Data-Driven Comparison with EXPORTS-NA
=================================================

The production workflow compares the model directly to observations from
the EXPORTS-NA campaign (North Atlantic, May 2021). The procedure consists
of loading daily turbulence profiles from ``keps_for_dave.mat``,
constructing an observation-based flux boundary condition from the UVP data,
spinning up the model to quasi-steady state, and accumulating model output
on the same days and depths as the UVP casts.

All UVP data must be filtered to 100–2000 μm before any comparison or use
as boundary forcing:

.. code-block:: matlab

    mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;

Without this filter the raw UVP data contains 85–93% zooplankton by biovolume
at all depths, which would bias the boundary condition and the model-data
comparison.

The flux boundary condition at 100 m is constructed each day from the UVP
spectrum at that depth. The daily turbulence profile is updated at the start
of each model day from ``keps_day``, so the shear coagulation and
disaggregation rates vary both with depth and time through the cruise.
The model is compared to UVP observations at 125, 325, and 475 m by
averaging the model state over the subset of days on which UVP casts were
collected at each depth.
