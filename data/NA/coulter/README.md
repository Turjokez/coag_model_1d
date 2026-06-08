# Coulter Counter

Source: SeaBASS EXPORTSNA, USM Zhang archive.

Main author / PI: Xiaodong Zhang.

Data in this folder:
- Coulter Counter particle size distribution.
- Size range is about 2-60 um.

Model use:
- Use this to constrain the lowest model particle bins.
- This helps because UVP starts at larger sizes.
- It is useful for checking whether bin 1 and nearby small bins are realistic.

First processing step:
- Parse `PSD_DNSD_*` fields.
- Convert number spectrum to volume spectrum if needed.
- Map the sizes into model bins.
