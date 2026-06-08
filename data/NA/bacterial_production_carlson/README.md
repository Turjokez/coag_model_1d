# Bacterial Production

Source: SeaBASS EXPORTSNA, UCSB Carlson archive.

Main author / PI: Craig Carlson.

Data in this folder:
- Bacterial abundance.
- Bacterial production / remineralization context.
- DOC files.

Model use:
- Use as supporting data for microbial loss.
- It helps estimate bacterial growth efficiency and check if `microbe_r0` is reasonable.
- This is not the direct model forcing, but it supports the microbial term.

First processing step:
- Read production and DOC fields.
- Use with respiration data to estimate bacterial growth efficiency.
