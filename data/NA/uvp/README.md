# Processed UVP Particle Size Distributions

Source: SeaBASS EXPORTSNA, UAF McDonnell archive.

Main author / PI: Andrew McDonnell / Lee Karp-Boss.

Data in this folder:
- Processed UVP5 Level-2 particle size distributions.
- Files include `PSD_DNSD_*` and `PSD_DVSD_*` columns.
- Size range starts near 57 um and extends to large aggregates.

Model use:
- This is the main observed particle size distribution for model comparison.
- Use top 5 m UVP data as daily surface boundary forcing.
- Use deeper UVP profiles to compare model particle size spectra with depth.
- `PSD_DVSD_*` is most useful for model biovolume `phi`.
- `PSD_DNSD_*` is useful for particle number checks.

First processing step:
- Parse the differential UVP file.
- Map UVP sizes to model section bins.
- Save a clean `uvp_model_bins.csv` file.
