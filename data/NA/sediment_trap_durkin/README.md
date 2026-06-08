# Durkin Sediment Trap Flux

Source: SeaBASS EXPORTSNA, MBARI Durkin archive.

Main author / PI: Colleen Durkin.

Data in this folder:
- Sediment trap particle flux files.
- Includes classified gel trap particle fluxes.
- May include fecal pellet and aggregate categories.

Model use:
- This is one of the best datasets for checking model particle flux.
- Use it to compare aggregate flux and fecal pellet flux separately.
- It can validate the separate fecal pellet pool `Y_fp`.

First processing step:
- Parse flux by depth and particle class.
- Compare model `Y` flux and `Y_fp` flux at matching trap depths.
