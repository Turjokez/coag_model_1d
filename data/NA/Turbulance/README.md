# VMP / Turbulence Data

Source: `keps_for_dave.mat`.

Main use:
- This gives real turbulent dissipation `epsilon(z,t)`.
- The model uses `epsilon(z)` to set the disaggregation limit `D_max(z)`.

Data in this folder:
- `S.z`: depth, negative ocean convention.
- `S.eps`: turbulent dissipation, units `m^2/s^3`.
- `S.T`, `S.S`, `S.rho`: temperature, salinity, density.
- `S.kappa_T`: thermal diffusivity.

Model use:
- Load this with `scripts/data/load_keps.m`.
- Convert `eps` from `m^2/s^3` to `cm^2/s^3` by multiplying by `1e4`.
- Use the time-mean profile for the first model-data run.
- Data only reaches about 300 m, so the loader holds the deepest value below that.

First processing step:
- Run `scripts/data/test_load_keps.m`.
- Check that `epsilon(z)` is high near the surface and low at depth.
