# Forward Plan -- June 11, 2026
## Based on Adrian meeting + Jackson (1997) + Zhang et al. (2022)

---

## What the papers say

**Jackson et al. (1997) — Monterey Bay, 1 µm to 1 cm:**
Particle size spectrum is well-fitted by a power law across six instruments from 1 µm to 10 mm. Number spectrum slope b ≈ 3.0 (aperture instrument), image-based instruments give b ≈ 2.96–3.00. Most particle volume sits in the 0.1–3 mm range. Result is consistent with coagulation-disaggregation theory.

**Zhang et al. (2022) — Ocean Station Papa, 0.02 µm to 2000 µm:**
Seven instruments combined. Junge slope (number spectrum) ≈ −3.8 for Coulter Counter, LISST, IFCB, ViewSizer. UVP gives slope −2.6 — significantly shallower than all other instruments. This means UVP alone underestimates small particles. Slope shows little variation in the upper 75 m. Particles 1–100 µm account for 70–90% of solid volume.

**Key implication for us:** UVP starts at 100 µm. If we use raw UVP as the surface boundary condition, we are missing 70–90% of the solid volume sitting in the 1–100 µm range. We need to extend the spectrum downward using a power law, with slope ~ −3.8 (number) as a starting point.

---

## Plan — 5 steps in order

---

### Step 1 — Find the mass overproduction bug

**What:** Run the mass budget diagnostic with real UVP surface forcing. Print how much mass enters at surface, how much is produced by coagulation/biology/zooplankton, how much leaves at the bottom, and what the net is vs UVP.

**Why Adrian said this:** Without the overlap bin fix, model gives 30–40× more mass than UVP. This should NOT depend on bin mapping. The real cause is likely that bins below 100 µm (no UVP data → currently set to zero or something arbitrary) are generating or destroying mass incorrectly when the model runs coagulation.

**Script:** Modify `run_compare_spectrum.m` or a new `run_mass_budget_check.m` to print a budget table per depth layer: surface flux in, coagulation production, disagg loss, zoo loss, microbe loss, bottom flux out.

**Success check:** Identify which term is 30–40× too big compared to UVP.

---

### Step 2 — Fix surface boundary condition: power-law extension to 1 µm

**What:** At each surface cast, take the UVP volume spectrum φ(d) from 100 µm upward, fit a power law in log–log space, and extrapolate it downward to 1 µm (model bin 1). Use this extended spectrum as the surface BC for all bins.

**Why:** UVP misses the 1–100 µm range. Jackson and Zhang both show the spectrum continues as a power law below 100 µm with slope ~ −3.8 (number). If we set those bins to zero, the model invents particles through coagulation and over-predicts mass.

**How to fit:** In log–log space, fit a line to log(φ) vs log(d) using the available UVP bins (say 100–500 µm range for cleaner fit, away from the large-particle noise). Extrapolate that line to model bin centers below 100 µm.

**Code location:** `get_daily_surface_phi.m` — add a power-law fill step after reading UVP bins.

**Expected behavior:** Model starts with more realistic small-bin mass → coagulation produces aggregates more like UVP → no more 30–40× overproduction.

---

### Step 3 — ε floor at 10⁻⁸ m²/s³ below 100 m

**What:** In `load_keps_daily.m` (or wherever ε(z) is passed to the model), clamp:

```matlab
eps_floor = 1e-8;   % m^2/s^3
keps_day.eps = max(keps_day.eps, eps_floor);
```

**Why:** Below 100 m, ε from keps_for_dave drops very low. With the Parker formula, D_max = Da × ε^(−1/4) becomes unrealistically large (> 10 mm or more) in the deep. This lets particles grow without any fragmentation check. Adrian said use whatever ε is at 100 m and hold it constant below — the data shows it is ~ 10⁻⁸ there.

**One line change** in `load_keps_daily.m` after interpolation.

---

### Step 4 — Start model at 100 m, predict 100–500 m

**What:** Write a new run script `run_100m_start.m`. Instead of surface UVP, use UVP at 100 m as the top boundary condition. Run the model from 100 m down to 500 m. Compare model output to UVP at 150, 200, 300, 400, 500 m.

**Why:** Below 100 m there is no primary production, minimal zooplankton biology, and the main driver is what sinks from above. It is a much cleaner test of sinking + coagulation + disaggregation. If we can match 100–500 m, then we know the deep physics is working and we can later tackle the complicated surface layer separately.

**Adrian's prediction:** Model will likely under-predict small particles at depth. Hypothesis: zooplankton swimming disaggregates large sinking particles into small ones — this is not in the current disagg formula (which uses ε only).

**Figure to make:** Same spectrum comparison as Figure 2 in report_june11, but top BC at 100 m and depth range 100–500 m.

---

### Step 5 — Invest in Alldredge (logistic) disagg, not Parker

**What:** Stop trying to tune D_a in the Parker formula. Instead, work on making the logistic (Alldredge) mode give physically correct behavior. Specifically:
- The redistribution bug is already fixed (p=0, uniform). ✓
- Next: check that r_max = C₀ × ε^(−B) with C₀ = 2×10⁻³ cm, B = 0.45 gives D_max in the right range vs UVP.
- Sweep C₀ the same way we swept D_a (×1, ×3, ×5) and compare to UVP at 75 m and 200 m.

**Why:** Parker (1972) was derived for wastewater treatment — inorganic flocs in high-shear lab conditions. It is not right for marine snow. Alldredge formula was measured on actual marine snow aggregates. Even if not perfect, it is the right physical analogy. Adrian said: "invest time in trying to get the Alldredge model to work the way we think it should."

**Script:** Copy `run_dmax_sensitivity.m` → `run_alldredge_sensitivity.m`, replace D_a sweep with C₀ sweep.

---

## Order of priority

| Step | What | Estimated effort |
|------|------|-----------------|
| 1 | Mass budget diagnostic | 1–2 hours |
| 2 | Power-law surface BC extension | 2–3 hours |
| 3 | ε floor (one line) | 15 min |
| 4 | 100 m start run + figure | 1–2 hours |
| 5 | Alldredge C₀ sweep | 2–3 hours |

Step 3 is trivial — do it immediately alongside Step 1. Steps 1 and 2 together fix the surface forcing problem, which is the root cause of the mass bug. Step 4 is the cleaner test Adrian wants. Step 5 is the longer-term disagg direction.

---

## What we are NOT doing yet

- Grid search on α × r₀ — not valid until surface forcing is fixed (Step 2).
- Zooplankton-only disaggregation term — wait for Adrian's input after Step 4 results.
- 51 µm two-class split — Adrian's separate analysis; not urgent for the column model.

---

## Key numbers from the papers

| Quantity | Value | Source |
|----------|-------|--------|
| Number spectrum slope | −3.8 (most instruments) | Zhang et al. 2022 |
| UVP slope alone | −2.6 (shallower — misses small particles) | Zhang et al. 2022 |
| Volume distribution peak | 0.1–3 mm (most mass) | Jackson et al. 1997 |
| Small particle volume fraction | 70–90% of solid volume in 1–100 µm | Zhang et al. 2022 |
| Fractal dimension D (small particles) | 2.8–3.5 (mean 3.3, nearly solid) | Zhang et al. 2022 |
| Fractal dimension D (large particles) | < 3 (porous, must be accounted for) | Zhang et al. 2022 |
