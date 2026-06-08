# Imaging FlowCytobot

Source: SeaBASS EXPORTSNA, WHOI Sosik archive.

Main author / PI: Heidi Sosik / Lee Karp-Boss.

Data in this folder:
- Imaging FlowCytobot particle and plankton files.
- Many small files (`2980` files).
- These cover smaller particles than UVP.

Model use:
- Use this to constrain the small size bins that UVP misses.
- This is most useful for bins below about 100 um.
- It can help set the surface condition for small particles and phytoplankton-like particles.

First processing step:
- Do this after UVP parsing.
- Build a small-particle spectrum and merge it with UVP by size.
- Keep the first model run simple if this takes too long.
