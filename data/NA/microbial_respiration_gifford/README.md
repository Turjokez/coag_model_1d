# Microbial Respiration

Source: SeaBASS EXPORTSNA, UNC Gifford archive.

Main author / PI: Scott Gifford.

Data in this folder:
- Bottle and flow-through oxygen respiration measurements.
- These give community and bacterial respiration rates.

Model use:
- Use to constrain microbial loss rate `microbe_r0`.
- This is important for remineralization of sinking particle material.
- It can help tune how much carbon is lost before particles reach depth.

First processing step:
- Convert respiration rates to a daily loss rate where possible.
- Compare with model microbial loss `r * particle_biovolume`.
