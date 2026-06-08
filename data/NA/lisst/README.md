# LISST / ViewSizer Particle Data

Source: SeaBASS EXPORTSNA, USM Zhang archive.

Main author / PI: Xiaodong Zhang.

Data in this folder:
- LISST-VSF particle files.
- ViewSizer particle file.
- These are particle and scattering measurements.

Model use:
- Use as an extra check on small to intermediate particle size spectra.
- This can help fill the gap between Coulter/IFCB and UVP.
- It is not the first comparison target, but it is useful for checking the shape of the PSD.

First processing step:
- Read the `PSD_DNSD_*` fields.
- Map sizes to model bins.
- Compare with Coulter and UVP in overlapping size ranges.

Note:
- The curation table listed LISST-DEEP under UAF/McDonnell, but that archive did not expose a separate LISST-DEEP file. These USM LISST/ViewSizer files are the available particle-size support data.
