# Report -- June 9, 2026
## First Data-Driven 1-D Column Run (EXPORTS North Atlantic)

---

We ran the 1-D model driven by real EXPORTS-NA data for the first time. Surface phi comes from UVP observations each day. The model predicts what happens below. We compare model output to the observed UVP profile at depth.

---

## Contents

1. Setup
2. What we fixed before running
3. Results
4. What the figures show
5. Next steps

---

## 1. Setup

**Data:** EXPORTS North Atlantic cruise, May 4--29, 2021. 22 days with UVP casts.

**Surface forcing:** Each day, the top model layer is set to the UVP particle size distribution from the top 5 m, filtered to sizes < 2000 µm. Sizes above 2000 µm are excluded because raw UVP at depth is 85--93% zooplankton-sized objects, not aggregates.

**Model config:**
- n_sections = 30, dt = 0.25 day (4 steps per day)
- sinking: kriest\_8
- all physics on: coagulation, disaggregation (depth-varying D\_max from ε(z)), zooplankton (Stemmann 2004 depth profiles), fecal pellets (cross-coagulation), microbial remineralization, mining
- alpha = 0.50, microbe\_r0 = 0.01 day⁻¹ (best-fit from 2D grid search)

**Real physics profile:** ε(z), T(z), S(z) from keps\_for\_dave.mat (VMP data, same cruise). Depth range: surface to 500 m (where VMP data ends).

Script: `scripts/data/run_data_column_daily.m`

---

## 2. What We Fixed Before Running

The previous comparison used raw UVP phi as the target. This was wrong. Raw UVP includes large objects (> 2000 µm) that are mostly zooplankton, not aggregates. The model only tracks aggregates and fecal pellets. Comparing model phi to raw UVP phi gave a false "36× mismatch" at depth.

After filtering UVP to < 2000 µm, the mismatch disappeared:

| Depth | Model | UVP < 2000 µm |
|-------|-------|----------------|
| surface | 2.062e-05 cm³ cm⁻³ | 2.060e-05 cm³ cm⁻³ |
| 175 m | 3.240e-06 cm³ cm⁻³ | 3.311e-06 cm³ cm⁻³ |

The same filter was applied to the 2D optimizer. Best fit: alpha = 0.50, r0 = 0.01 day⁻¹, loss = 0.88. The loss surface is now meaningful.

---

## 3. Results

The model tracks the UVP surface forcing well. The cruise-mean depth profile matches UVP at the two checked depths. The time-depth picture shows the model responding to the daily surface forcing signal.

**Figures produced:**

| Figure | File | What it shows |
|--------|------|----------------|
| Cruise-mean depth profile | `data_daily_depth_profile.png` | Model vs UVP < 2000 µm, 0--500 m |
| Surface time series | `data_daily_surface_time.png` | Model surface phi vs UVP forcing each day |
| Time-depth heatmap | `data_daily_timedepth.png` | Model phi(depth, day) over the full cruise |
| Selected-day profiles | `data_daily_profiles_selected.png` | Model on day 1 / mid / end vs UVP mean |

---

## 4. What the Figures Show

**Cruise-mean depth profile** (`data_daily_depth_profile.png`): Model and UVP agree well in the upper 200 m. Both decrease roughly with depth.

**Surface time series** (`data_daily_surface_time.png`): The surface layer follows the UVP forcing. The Dirichlet BC is working correctly -- surface is reset to UVP data before and after each substep.

**Time-depth heatmap** (`data_daily_timedepth.png`): Shows the column evolving over 26 days. Phi is highest in the top 100 m and decreases with depth. Day-to-day variability at the surface (driven by the UVP forcing) propagates downward over several days via sinking.

**Selected-day profiles** (`data_daily_profiles_selected.png`): Early, mid, and late cruise model profiles overlaid on the UVP cruise mean. Day-to-day variability is visible in the upper 100 m. Below 200 m the profiles are more stable.

---

## 5. Next Steps

The model is now running with real data. The immediate gaps are:

1. **Date-resolved UVP comparison.** The current comparison uses UVP cruise mean at depth. To do "May 18 model vs May 18 UVP" we need to parse the UVP .sb file by date at each depth. The date information is in the file -- this is the next script to write (`parse_uvp_daily.m`).

2. **Finer grid search or conjugate gradient** near alpha = 0.50, r0 = 0.01. The 6×6 grid has step size 0.1 in alpha; a finer pass around the best point would tell us how sharp the minimum is.

3. **Contact Amy and Debbie** for net tow zooplankton data to constrain Zc, Zf depth profiles. Current profiles are from Stemmann 2004 Atlantic data; EXPORTS-specific profiles would improve the fit.

4. **EXPORTS data curation table** -- build the full model variable → dataset → SeaBASS link → contact table that Adrian asked for.
