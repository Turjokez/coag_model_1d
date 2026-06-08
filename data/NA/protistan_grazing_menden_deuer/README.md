# Protistan Grazing

Source: SeaBASS EXPORTSNA, URI Menden-Deuer archive.

Main author / PI: Susanne Menden-Deuer.

Data in this folder:
- Protistan biomass and grazing files.
- Includes grazing measurements and flow cytometry grazing data.

Model use:
- Use as a first constraint for the mining term.
- The model mining parameter `mining_Zm` needs a depth profile.
- These data help decide whether mining should be weak, strong, surface-focused, or depth-focused.

First processing step:
- Parse grazing rate and protist biomass by depth.
- Build a rough `Zm(z)` profile for the mining term.
