# Passow TEP and Sinking Velocity

Source: SeaBASS EXPORTSNA, MU Newfoundland Passow archive.

Main author / PI: Uta Passow.

Data in this folder:
- TEP and biogeochemical files.
- Marine Snow Catcher sinking velocity file.
- Marine snow volume file.

Model use:
- TEP can help constrain particle stickiness `alpha`.
- Sinking velocity can validate the `kriest_8` sinking law.
- Marine snow volume helps check aggregate size and volume assumptions.
- Biogeochemical particle type data can help check fecal/aggregate density assumptions.

First processing step:
- Parse sinking velocity versus particle size or particle type.
- Compare observed sinking speeds with model `SettlingVelocityService`.
- Parse TEP profile and use it as a prior for `alpha`.
