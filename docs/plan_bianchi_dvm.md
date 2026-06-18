# Plan — Bianchi-Style DVM Implementation
# Status: PENDING Adrian review. Do not implement until BC test results reviewed.

---

## Why This Is Different from What We Already Tried

The Archibald (2019) approach we implemented is a **static rerouting fraction**: a fixed fraction of surface grazing fecal flux is redirected to a target depth band each timestep. It has no memory, no gut pool, no time lag.

The reason it failed: 115 µm fecal pellets sink at ~69 m/day. They transit the 300-500 m band in ~2 days. Rerouting a small flux into fast-sinking particles cannot build standing stock.

Bianchi et al. (2013, *Global Biogeochemical Cycles*) is different in three ways:

1. **Explicit migrant gut pool** — fecal is not released instantly; it builds up in the gut during surface feeding and clears over hours at depth. This produces a spatially spread release.
2. **Migrant mortality at depth** — a fraction of migrant biomass dies at depth, contributing directly to POC (marine snow, not just fecal).
3. **DOC release** — dissolved organic carbon from migrant respiration and sloppy feeding at depth is not tracked here, but the POC part is.

The key difference: mortality at depth produces **marine snow aggregates (Y), not just fecal pellets (Y_fp)**. These sink slower and aggregate. This is the mechanism that could populate the 200-1200 µm range.

---

## Equations

Let $Z_m$ [ind m$^{-3}$] be the migrant zooplankton concentration. In the current model, zoo profiles are prescribed via Stemmann Fig 1. For DVM, we need to distinguish:

- $Z_{\rm res}(z)$: resident (non-migrating) zooplankton — stays at prescribed depth
- $Z_{\rm mig}$: migrating fraction, which moves between surface and deep daily

**Feeding zone (surface, nighttime, $z \leq z_{\rm feed}$):**

$$\frac{dG}{dt} = I(z,t) - \frac{G}{\tau_{\rm gut}} \tag{1}$$

where $G$ [volume individual$^{-1}$] is the gut content and $I$ is the ingestion rate from the normal grazing formulation. At the start of the descent, $G_0$ is the accumulated gut content.

**Descent and gut clearance (daytime, $z_{\rm feed} < z \leq z_{\rm deep}$):**

During descent the migrant is at depth. Gut clears with time constant $\tau_{\rm gut}$:

$$G(t) = G_0\, e^{-t/\tau_{\rm gut}} \tag{2}$$

The fecal production rate at depth is:

$$\phi_{\rm fp}^{\rm dvm}(z,t) = Z_{\rm mig}(z)\, \frac{G_0}{\tau_{\rm gut}}\, e^{-t/\tau_{\rm gut}} \cdot p_{\rm eg} \tag{3}$$

where $p_{\rm eg}$ is the egestion fraction (same as `zoo_p` currently).

In a 1-D daily-timestep model, equations (2-3) simplify. On each day, migrant grazers in the feeding zone accumulate $G_0$ over one night. On the following day, $G_0$ is released in the deep band over timescale $\tau_{\rm gut}$.

**Migrant mortality at depth:**

$$\phi_{\rm mort}(z) = Z_{\rm mig}(z)\, m_{\rm dvm} \tag{4}$$

where $m_{\rm dvm}$ [day$^{-1}$] is the daily mortality rate of migrants at depth. This POC goes into the marine snow array $Y$ (not fecal), at a bin corresponding to zooplankton body size.

**Practical 1-D discrete version:**

For each timestep in the feeding zone ($z_k \leq z_{\rm feed}$):

$$G_0 \mathrel{+}= I_k \cdot \Delta t \quad \text{(accumulated overnight)} \tag{5}$$

On the descent day, distribute into deep layers with exponential weighting:

$$\Delta Y_{\rm fp}^{(k)} = p_{\rm eg} \cdot G_0 \cdot Z_{\rm mig} \cdot w_k, \quad w_k \propto e^{-z_k / (w_{\rm mig} \cdot \tau_{\rm gut})} \tag{6}$$

where $w_{\rm mig}$ [m day$^{-1}$] is the descent speed (~100 m day$^{-1}$ for large copepods).

---

## State Variables Needed

| Variable | Location | Description |
|----------|----------|-------------|
| `G_gut` | scalar or per-layer | accumulated gut content from previous night |
| `Z_mig(z)` | new config param or derived from Stemmann profile × dvm_p | migrant concentration profile |
| `w_mig` | config param | migrant descent speed [m/day] |
| `tau_gut` | config param | gut clearance time constant [day] |
| `m_dvm` | config param | migrant mortality rate at depth [day⁻¹] |
| `dvm_mort_bin` | config param | which Y bin gets mortality POC |

---

## Code Changes Required

### SimulationConfig.m — new params

```matlab
% Bianchi-style DVM (gut pool + mortality)
enable_dvm_bianchi = false;
dvm_tau_gut        = 0.25;   % gut clearance time [day] (~6 hours)
dvm_w_mig          = 100;    % migrant descent speed [m/day]
dvm_m_mort         = 0.01;   % migrant daily mortality rate [day^-1]
dvm_mort_bin       = 15;     % Y bin for mortality POC (~400 um)
dvm_mig_frac       = 0.5;    % fraction of Stemmann zoo that migrates
```

### ColumnRHS.m — stepY() changes

Three additions to the grazing block:

1. **Gut accumulation** (in feeding zone layers): add ingested BV to `G_gut` pool instead of (or in addition to) normal fecal.

2. **Gut release** (in deep layers): at each deep layer, add $\Delta Y_{\rm fp}$ from equation (6) based on exponential weight and current `G_gut`.

3. **Mortality source** (in deep layers): add `dvm_m_mort * Z_mig(k) * dt` to `Y(k, dvm_mort_bin)` each timestep.

A new persistent or passed variable `G_gut` is needed between timesteps — this is the key structural change that Archibald-style does not have.

### New helper: dvm_mig_profile.m

Derives $Z_{\rm mig}(z)$ from Stemmann profile × `dvm_mig_frac`. Returns layer-indexed vector.

---

## Parameters to Discuss with Adrian

| Parameter | Value range | Source |
|-----------|-------------|--------|
| $\tau_{\rm gut}$ | 4-12 hours (0.17-0.5 day) | Dagg & Walser 1987; Paffenhöfer & Knowles 1980 |
| $w_{\rm mig}$ | 50-200 m day$^{-1}$ | Wiebe et al. 1992; Childress et al. 1980 |
| $m_{\rm dvm}$ | 0.005-0.02 day$^{-1}$ | Bianchi et al. 2013 Table 1 |
| `dvm_mig_frac` | 0.2-0.8 | Archibald et al. 2019; Steinberg et al. 2000 |
| `dvm_mort_bin` | bin 12-18 (~250-700 µm) | copepod body size at EXPORTS site |

---

## Why This Might Actually Work (Unlike Archibald)

The Archibald-style test failed because fast-sinking small fecal at 115 µm transits the deep layer in ~2 days → negligible time-mean standing stock.

Bianchi-style has two differences:

1. **Mortality POC goes into Y (marine snow), not Y_fp.** Marine snow aggregates are larger and slower sinking. They also coagulate with other particles, which can shift the size distribution toward larger bins — exactly the 200-1200 µm range that is missing.

2. **Gut clearance spreads release over $\tau_{\rm gut}$ = 6-12 hours** at depths that the migrant passes through. With descent speed 100 m/day, gut clears over ~40-100 m of descent. This distributes fecal over a broader depth band than our single-bin injection.

Whether these differences are large enough to explain a 66% deficit is unknown. That is an empirical question requiring the simulation.

---

## Decision Needed from Adrian

1. Is $m_{\rm dvm}$ (mortality at depth) significant at the EXPORTS site?
2. What is the appropriate $Z_{\rm mig}$ — fraction of Stemmann profile, or a separate published value?
3. Is the gut clearance time approach appropriate, or does Adrian prefer a simpler daily-mean active flux approach (e.g., from Steinberg et al. 2000 EXPORTS data)?

**Do not implement until BC test result and Adrian's input on these three questions.**
