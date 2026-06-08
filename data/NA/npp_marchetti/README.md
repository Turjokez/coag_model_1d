# Net Primary Production and Chlorophyll

Source: SeaBASS EXPORTSNA, UNC Marchetti archive.

Main author / PI: Adrian Marchetti.

Data in this folder:
- Chlorophyll file.
- 13C/15N POC/PON production file.

Model use:
- Use NPP to set or check the surface production term.
- Use chlorophyll to check bloom timing and surface biomass.
- POC/PON can help connect model biovolume to carbon units.

First processing step:
- Parse NPP and chlorophyll by date and depth.
- Estimate `surface_pp_mu` or daily surface source strength.
