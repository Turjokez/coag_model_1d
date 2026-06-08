# CTD Profiles

Source: SeaBASS EXPORTSNA, UCSB/CRSEO.

Main author / PI: Norman Nelson group.

Data in this folder:
- Binned CTD profiles from DY131 and JC214.
- Variables include depth, temperature, salinity, oxygen, chlorophyll, and optical variables.

Model use:
- Use temperature and salinity to build `DepthProfile`.
- Temperature changes viscosity, so it affects sinking speed and differential settling.
- Oxygen is useful as background context for microbial respiration.
- Surface POC files in this archive can help set or check surface particle carbon.

First processing step:
- Make one clean table with `date, time, lat, lon, depth, temperature, salinity, oxygen`.
- Then interpolate `T(z)` and `S(z)` to the model depth grid.
