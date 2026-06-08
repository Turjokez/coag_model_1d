# EXPORTS North Atlantic Data

This folder keeps the downloaded EXPORTS North Atlantic data used for the
1-D particle model.

Main use in the model:
- CTD and float data give temperature and salinity for `DepthProfile`.
- UVP, Coulter, LISST, ViewSizer, and IFCB give particle size distributions.
- Sediment traps and 234Th give particle flux and transfer efficiency.
- Microbial and particle respiration data constrain microbial loss `r`.
- TEP and sinking velocity data help constrain stickiness and sinking speed.
- Protistan grazing data can help constrain the mining term.
- NPP and POC data set or check surface production and surface POC.

First model-data step:
1. Parse CTD into clean `T(z,t)` and `S(z,t)`.
2. Parse UVP top 5 m into daily surface particle forcing.
3. Map observed particle size bins into model section bins.
4. Compare model and observed size spectra with depth.

Important remaining gaps:
- VMP turbulence data for `epsilon(z)` is still missing.
- Zooplankton depth profiles `Zc(z)` and `Zf(z)` still need PI contact.
- EXPORTSNA fecal pellet production is still not available here.
