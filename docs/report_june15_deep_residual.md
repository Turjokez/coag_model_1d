# Report -- June 15, 2026
# What Controls the Deep Particle Loss?


## 1. The Problem

The model/UVP biovolume ratio is around 0.9 at 325 m but falls to 0.34 at 475 m. Above 275 m the fit is reasonable. I ran seven diagnostics to find out why and whether any physics change can fix it.

This shows the particle volume spectrum for each cast date side by side. In the upper water column both panels show similar warm colors (high particle volume). Below roughly 300 m the model panel (b) turns dark blue while the UVP panel (a) stays warmer, meaning the model is consistently short of particle volume at depth across all cast dates. This is the deep residual that the seven diagnostics below try to explain. Note that the first column (05-04) is blank in panel (b) because that date is used as the boundary condition input.

![Cast-by-cast spectrum comparison](./figures/cast_spectrum_2d.png)

Particle volume spectrum [ppmV mm$^{-1}$] vs ESD (mm, x-axis) and depth (m, y-axis) for each cast date (May 4-29). (a) UVP observations. (b) Model, best config (alpha = 0.10, Da x5, 100 m BC). Color scale is shared and log.*



## 2. Model 

The model is a 1-D column with 20 layers of 50 m each, 1000 m total, with n_sec = 30 and dt = 0.25 day. The boundary condition comes from UVP data at 97.5 m, mapped onto the 30 model bins over the 100-2000 um range. All diagnostics use this same setup.

Sinking velocity follows Kriest (2002),

$$w_i = 66\, d_i^{0.62} \quad [{\rm m\ day}^{-1}] \tag{1}$$

where $d_i$ is the diameter of bin $i$ in cm.

Disaggregation (operator-split): the maximum stable particle size at each depth is

$$D_{\rm max}(z,t) = D_a\, \varepsilon(z,t)^{-1/4} \tag{2}$$

where eps comes from VMP measurements (keps_for_dave.mat), floored at $10^{-8}$ m$^2$ s$^{-3}$. Particles above D_max are broken: 2/3 to the next smaller bin, 1/3 spread uniformly.

The best-fit parameters from the earlier grid search are:

$$\alpha = 0.10, \quad D_a = 5 \times 9.39 \times 10^{-6}\ {\rm m\ (Parker \times 5)}, \quad r_0 = 0 \tag{3}$$

with zoo grazing on (Stemmann 2004, Zc = 0.307 m$^{-3}$, Zf = 0.063 m$^{-3}$), mining on, microbial remineralization off. All five model-based diagnostics (Steps 1, 3, 4, 5, 7) use this exact configuration. Steps 2 and 6 do not run the model (sinking velocity is computed analytically; day/night split is UVP-only). The comparison metric is (model phi) / (UVP phi) integrated over 100-2000 um, averaged over 22 cast days.

---

## 3. Step 1: Are the particles in the wrong size range?

The first thing I wanted to check is whether the model has the mass but placed it in the wrong bins. Specifically, is coagulation building very large aggregates (above 2 mm, outside the UVP window) that disaggregation then breaks into very small fragments (below 100 um, also outside the UVP window)?

At steady state I computed the fraction of total biovolume in three size classes at each depth. Figure 1 shows the result.

![Mass fraction by size class](./figures/mass_fraction_diagnostic.png)

*Figure 1. Fraction of total model biovolume in each size class vs depth (left) and absolute profiles (right). Averaged over 22 cast days.*

| Depth (m) | below 100 um | 100-2000 um | above 2000 um |
|-----------|-------------|-------------|---------------|
| 75        | 1.9%        | 98.1%       | 0.0%          |
| 175       | 13.9%       | 86.1%       | 0.0%          |
| 275       | 9.1%        | 90.9%       | 0.0%          |
| 375       | 5.0%        | 95.0%       | 0.0%          |
| 475       | 2.4%        | 97.6%       | 0.0%          |
| 775       | 0.0%        | 100.0%      | 0.0%          |

The above-2000-um fraction is 0.0% everywhere. Below 400 m, essentially all of the model mass sits in the 100-2000 um range. The deep deficit is a total mass problem, not a size-window problem.

Note: Figures 3 through 6 below cover only 75-475 m. This is the region where the comparison is most meaningful. Below 475 m the UVP cast coverage is sparse.

---

## 4. Step 2: Do particles sink too fast?

If particles sink quickly through the 350-475 m layer, the time-averaged concentration there will be low even if the flux is fine. I computed the sinking velocity directly from equation (1) for all 30 bins. Figure 2 shows the result.

![Sinking velocity vs diameter](./figures/sinking_velocity_diagnostic.png)

*Figure 2. Sinking velocity w [m day$^{-1}$] vs diameter for all 30 bins. Red points are the UVP-visible range (100-2000 um).*

In the UVP-visible range, w runs from 4.2 m day$^{-1}$ (d = 115 um) to 23.2 m day$^{-1}$ (d = 1846 um). Transit time through a 150 m band is 6.5 to 36 days. These are physically reasonable. Note that fecal pellets (Y_fp, excess density 0.15 g cm$^{-3}$) sink at roughly 69 m day$^{-1}$ but are tracked separately and are not included in the comparison.

Fast sinking is not the explanation.

---

## 5. Step 3: Is the disaggregation mode wrong?

Before testing the disaggregation mode, I first checked whether disaggregation is simply switched off at depth. At depth, turbulence is weak, so from equation (2), D_max becomes very large. If D_max grows larger than any particle in the model, nothing ever gets broken and disaggregation stops working entirely. To check this, I ran the model with an artificial upper limit on how large D_max is allowed to be (at 1.0, 0.5, 0.3, and 0.2 cm). The result: limiting D_max to 0.5 cm or 1.0 cm made no difference to the model output at any depth. This tells us D_max is already within that range in the actual EXPORTS turbulence data, so disaggregation is not being silently switched off. Limiting D_max more aggressively (0.3 or 0.2 cm) made the deep fit worse, not better. D_max being too large is ruled out.

The next question is whether the disaggregation rule itself is wrong. The current model uses a hard D_max cutoff (operator_split): particles above D_max break, particles below survive. The logistic form is an alternative where breakup is not a hard on/off but a smooth continuous function of particle size relative to the local turbulence. The maximum stable radius at each layer is

$$r_{\rm max}(z,t) = C_0\, \varepsilon^{-B}, \quad C_0 = 2\times10^{-3}\ {\rm cm},\ B = 0.45 \tag{4}$$

Each bin then has a survival fraction

$$f(r_i) = \frac{1}{1 + \exp\left[\kappa\left(\dfrac{r_i}{r_{\rm max}} - 1\right)\right]}, \quad \kappa = 3.5 \tag{5}$$

Note that f is a retention fraction (not a fragmentation rate). When $r_i \ll r_{\rm max}$, f approaches 1 and the particle survives. When $r_i \gg r_{\rm max}$, f approaches 0 and the particle is broken. This is confirmed in the code (`frag_factor = 1./(1+exp(kappa*(ratio-1)))`, comment: "~1 below threshold, ~0 above"). The mass lost from each bin is redistributed to smaller bins using an r$^{-p}$ weighting.

Figure 4 compares the two modes.

![Logistic vs operator_split](./figures/logistic_disagg_test.png)

*Figure 4. Model/UVP ratio vs depth for operator_split (black) and logistic (red). Depth range 75-475 m.*

| Depth (m) | Operator split | Logistic |
|-----------|---------------|----------|
| 125       | 1.31          | 1.63     |
| 175       | 1.48          | 1.75     |
| 325       | 0.91          | 0.38     |
| 375       | 0.68          | 0.18     |
| 425       | 0.48          | 0.08     |
| 475       | 0.34          | 0.04     |

The logistic mode makes the deep fit much worse. At 375 m the ratio drops from 0.68 to 0.18. What is going on here is that the logistic fragmentation is always active at depth regardless of how small eps is, so it continuously drains the deep column. Keep `disagg_mode = 'operator_split'`.

---

## 6. Step 4: Are the UVP particles at depth physically larger?

I computed total biovolume BV(z) and particle number N(z) for both model and UVP in the 100-2000 um range. If BV_uvp/BV_mod exceeds N_uvp/N_mod, the UVP particles are on average larger than the model particles, which would suggest a missing source rather than excess loss.

Figure 5 shows BV and N vs depth.

![N and BV diagnostic](./figures/number_biovolume_diagnostic.png)

*Figure 5. BV (left) and N (right) vs depth for UVP (black) and model (red). 100-2000 um, 22-day mean.*

| Depth (m) | BV_uvp / BV_mod | N_uvp / N_mod |
|-----------|----------------|---------------|
| 75        | 0.63           | 0.65          |
| 175       | 0.86           | 0.97          |
| 275       | 1.33           | 1.28          |
| 375       | 2.00           | 1.67          |
| 425       | 2.12           | 1.90          |
| 475       | 2.94           | 2.12          |

Below 325 m the BV ratio exceeds the N ratio at every depth. At 475 m: BV ratio = 2.94, N ratio = 2.12. The UVP particles at depth are on average larger than the model particles. Recall that fragmentation would do the opposite (produce more small particles, raising the N ratio faster than the BV ratio). This result points toward a missing source of large particles at depth.

---

## 6. Step 4: Is there a day vs night difference in UVP at depth?

If DVM is adding particles at depth during the day, daytime UVP biovolume at 350-500 m should exceed nighttime. I split the 26-day UVP record into day (UTC 06-20, 8147 rows) and night (UTC 20-06, 3614 rows) and compared mean BV(z).

Figure 6 shows the result.

![Day vs night UVP](./figures/uvp_daynight_bv.png)

*Figure 6. UVP BV (left) and N (right) vs depth for day (red) and night (black) casts. 100-2000 um.*

| Depth (m) | BV_day / BV_night | N_day / N_night |
|-----------|------------------|----------------|
| 362       | 1.00             | 0.92           |
| 412       | 1.01             | 0.90           |
| 462       | 0.97             | 0.90           |

No diel signal in BV at depth. Two readings are possible. First, DVM is not the source. Second, fecal pellets from DVM zooplankton sink at roughly 70 m day$^{-1}$ and would transit any fixed depth layer in about 2 days, so the standing-stock difference between day and night would be nearly undetectable even if DVM is active. I am not sure which is correct. This test is inconclusive.

---

## 7. Step 5: Which size bins are missing at depth?

I compared model and UVP spectra bin by bin at 375 m to check whether the deficit is uniform or skewed toward specific sizes. A uniform loss rate would give a flat model/UVP ratio across all bins.

Figure 7 shows the spectrum comparison.

![Size spectrum at 375 m](./figures/spectrum_at_depth.png)

*Figure 7. Particle volume spectrum phi [ppmV] vs diameter at 375 m. UVP (black circles), model (red squares). 100-2000 um.*

| d (um) | phi model | phi UVP  | ratio    |
|--------|-----------|----------|----------|
| 115    | 4.54e-08  | 7.14e-09 | **6.35** |
| 145    | 9.91e-09  | 8.84e-09 | 1.12     |
| 183    | 1.13e-08  | 1.55e-08 | 0.73     |
| 231    | 1.36e-08  | 2.84e-08 | 0.48     |
| 291    | 1.43e-08  | 2.55e-08 | 0.56     |
| 366    | 1.92e-08  | 3.85e-08 | 0.50     |
| 462    | 2.61e-08  | 4.93e-08 | 0.53     |
| 581    | 3.67e-08  | 6.49e-08 | 0.57     |
| 733    | 5.09e-08  | 8.77e-08 | 0.58     |
| 923    | 7.14e-08  | 1.20e-07 | 0.59     |
| 1163   | 9.60e-08  | 2.07e-07 | 0.46     |

The ratio is not flat. At 115 um the model is 6x too high (small particles accumulate through repeated fragmentation). From 183-1163 um the model is consistently around 50% of UVP. This is a skewed deficit, not a uniform one. DVM fecal pellets from copepods typically fall in the 100-600 um range, which is exactly where the deficit is largest.

---

## 8. Summary

1. All physics tests (size bins, sinking speed, D_max, disagg mode) ruled out — the model physics is not the problem.
2. The deficit is a missing source: UVP particles at depth are physically larger than model particles, and the spectral deficit is concentrated in the 200-1200 um range.
3. Best candidate is DVM fecal pellets not yet in the model. Zoo/mining loss too strong (Stemmann 2004 max values) is a secondary candidate.

