# Report — May 06, 2026
## 1-D Depth-Dependent Coagulation Model: Build, Verification, and Fragmentation Diagnosis

---

## What this covers


- Clearing up the DS kernel question — whether it was actually using the right sinking law or not
- Re-running the full April step chain on the May 06 code to confirm nothing changed
- Running a process isolation matrix to find what was causing a huge particle number explosion with fragmentation on
- Building a new set of OOP classes for clean 1-D column runs with depth-varying kernel scaling


---

## 1. Sorting out the DS kernel

Earlier figures showed all four sinking laws looking almost the same in the DS kernel panels. I needed to check whether the code was actually picking up the law choice, or silently using the same formula no matter what.

The DS kernel depends on the *speed difference* between two particles, not their absolute speeds:

$$\beta_{DS}(d_1, d_2) \;=\; \frac{\pi}{4}\,(d_1 + d_2)^2\,\bigl|w(d_1) - w(d_2)\bigr|$$

This means two laws that look similar in absolute speed can still give very different collision rates — especially when one particle is much larger than the other. The right test is to fix the small particle at $d_1 = 1\,\mu\text{m}$ and scan the larger partner size $d_2$. At large $d_2$, the small particle's speed becomes negligible and the kernel is driven almost entirely by $w(d_2)$, so the law differences become visible.

![DS settling speed laws](figures/kernel_ds_settling_speed_laws.png)

The figure above shows how the four laws compare in absolute sinking speed. The `kriest_8` and `kriest_9` laws share the same exponent ($b = 0.62$) but `kriest_9` is twice as fast. The `current` law is the slowest at small sizes but has a steeper slope and overtakes the others at large sizes. `siegel_2025` sits in between.

![DS linecuts fixed sizes](figures/kernel_ds_linecuts_fixed_sizes.png)

When I fix $d_1 = 1\,\mu\text{m}$ and scan $d_2$, the laws clearly separate at large sizes. At $d_2 = 10{,}000\,\mu\text{m}$ the kernel spans roughly two orders of magnitude between `kriest_8` and `current`. The earlier figures that looked flat were plotted with a normalization that forced all four laws to meet at the same reference diameter — that forced meeting is what made them look similar, not the physics. Once I removed that constraint, the correct separation came back.

The code path is also confirmed: `ds_kernel_mode = 'sinking_law'` routes the DS kernel through the named law. The old `legacy` mode was not updating DS values when the law changed — that is why the earlier panels were misleading.

![DS size-size map](figures/kernel_ds_size_size_repo_laws.png)

The full size-size beta map shows the expected pattern — near-zero along the diagonal where the two particles have similar speeds, and large values off-diagonal where the speeds differ. The law choice clearly changes the shape of this map, which affects which size pairs dominate the collision process.

To see where DS and shear dominate relative to each other, I mapped $\log_{10}(\beta_{DS}/\beta_\text{shear})$ at two dissipation levels. The shear kernel is:

$$\beta_\text{shear}(d_1, d_2) \;=\; \frac{\pi}{6}\,\left(\frac{\varepsilon}{\nu}\right)^{1/2}\,\left(\frac{d_1 + d_2}{2}\right)^3$$

where $\varepsilon$ is turbulent dissipation and $\nu$ is kinematic viscosity.

![DS vs shear ratio low eps](figures/apr08_diffsed_vs_shear_ratio_eps_1em08.png)

At $\varepsilon = 10^{-8}\,\text{W\,kg}^{-1}$ (calm interior), DS dominates almost everywhere — about 97% of the size-size space for `kriest_8`.

![DS vs shear ratio high eps](figures/apr08_diffsed_vs_shear_ratio_eps_1em04.png)

At $\varepsilon = 10^{-4}\,\text{W\,kg}^{-1}$ (energetic surface), shear expands and DS drops to about 63%. At the default conditions used in the 1-D tests (`kriest_8`, calm interior), DS is the dominant collision mechanism.

---

## 2. Re-running the April 

After sorting out the DS kernel, I re-ran the full April 09 step chain on the May 06 code. The point was not to find new results but to have a confirmed record showing the model still behaves the same way after the code updates.

### Step 1 — Sinking speed laws

I compared the four laws in speed and travel time to 1000 m. The laws are:

$$w_\text{kriest8} = 66\,d^{0.62}, \quad w_\text{kriest9} = 132\,d^{0.62} \quad (d \text{ in cm},\; w \text{ in m/day})$$

$$w_\text{siegel} = 20.2\,D_\text{mm}^{0.67} \quad (D_\text{mm} \text{ in mm})$$

and `current`, which uses the image-to-volume relation in the repo.

![rerun step1 speed laws](figures/may06_rerun_step1_sinking_speed_laws.png)

![rerun step1 travel time](figures/may06_rerun_step1_travel_time_laws.png)

All four laws separate cleanly. At $d = 1\,\text{mm}$, travel times to 1000 m are about 63 days (`kriest_8`), 32 days (`kriest_9`), 50 days (`siegel`), and 9 days (`current`). These match the April numbers. Mean timing error is about 2.6%, which comes from discretization, not the laws themselves.

### Step 2 — Transport scheme

I tested upwind vs Lax-Wendroff with advection only (no diffusion or process terms).

![rerun step2 transport error](figures/may06_rerun_step2_transport_travel_error.png)

The figure shows arrival-time error vs particle size for both schemes. The result:

- **Upwind**: `neg_count = 0`, mean error ~ 0.09%
- **Lax-Wendroff**: 550,638 negative values, larger timing error

Upwind is the safe choice for any run where concentrations must stay physically non-negative. That settled the scheme selection.

### Step 3 — Diffusion

I added flux-form vertical diffusion. The full form I used is:

$$\frac{\partial}{\partial z}\!\left(K_z\,\frac{\partial N}{\partial z}\right) = K_z\,\frac{\partial^2 N}{\partial z^2} + \frac{\partial K_z}{\partial z}\,\frac{\partial N}{\partial z}$$

Both terms are needed — the drift term on the right becomes important where $K_z$ changes steeply with depth. Cell-face values are arithmetic means: $K_{z,\text{face}} = \frac{1}{2}(K_z(k) + K_z(k+1))$, and zero-flux conditions are applied at both the top and bottom.

![rerun step3 diffusion signal width](figures/may06_rerun_step3_diffusion_signal_width.png)

![rerun step3 diffusion conservation](figures/may06_rerun_step3_diffusion_conservation.png)

The pulse at 1000 m is broader with diffusion on — for the 100 µm bin, the width grows from 17.8 to 23.2 days. Conservation error stays at $10^{-13}\%$ (machine precision). The signal spreads but no volume is created or lost.

### Step 4 — Coagulation

The combined collision kernel is built from three mechanisms:

$$\beta_{ij} \;=\; \alpha\,\bigl(\beta^\text{Br}_{ij} + \beta^\text{sh}_{ij} + \beta^\text{DS}_{ij}\bigr)$$

where $\alpha$ is the stickiness coefficient. The Brownian kernel is:

$$\beta^\text{Br}_{ij} \;=\; \frac{2\,k_B\,T}{3\,\mu}\,\frac{(d_i + d_j)^2}{d_i\,d_j}$$

The shear kernel is given in Section 1, and the DS kernel is also from Section 1. In the code, `BetaAssembler.combineAndScale` builds the combined matrix, applying the correct physical prefactors for each mechanism and converting units to m³/day.

The sectional gain-loss update for bin $k$ is:

$$\frac{dN_k}{dt}\bigg|_\text{coag} \;=\; \underbrace{\sum_{\substack{i \leq j \\ v_i + v_j \in [v_k,\, v_{k+1})}} \beta_{ij}\,N_i N_j}_{\text{gain from mergers}} \;-\; \underbrace{N_k \sum_j \beta_{kj}\,N_j}_{\text{loss from collisions}}$$

In the code this is handled by the `b25` and `b1` matrices: `b25` carries the loss term, `b1` carries the gain from adjacent bin mergers.

![rerun step4 coag conservation](figures/may06_rerun_step4_coagulation_conservation.png)

![rerun step4 coag psd](figures/may06_rerun_step4_coagulation_column_psd.png)

Conservation: $1.67 \times 10^{-5}\%$, `neg_count = 0`. The PSD figure shows coagulation shifting material into larger bins (red above black at large sizes) — moderate effect because sinking is also removing particles throughout the run.

### Step 5 — Fragmentation

The legacy fragmentation rate per bin $k$ grows exponentially with bin index:

$$\frac{dv_k}{dt}\bigg|_\text{frag} \;=\; -r_k\,(v_k - c_4\,v_{k+1}), \qquad r_k = \frac{c_3\,c_4^k}{86400}$$

where $v_k$ is volume in bin $k$ and $r_k$ is in s$^{-1}$. The factor $(v_k - c_4 v_{k+1})$ can change sign, which means the term is not strictly a loss — it can become a source when $v_k < c_4 v_{k+1}$. That sign reversal is one reason the legacy term is not mass-conserving.

![rerun step5 frag conservation](figures/may06_rerun_step5_fragmentation_conservation.png)

![rerun step5 frag small-size volume](figures/may06_rerun_step5_fragmentation_small_size_volume.png)

In the isolated step test: conservation error $3.83 \times 10^{-5}\%$, `neg_count = 0`. The time to reach 80% small-size volume fraction decreases from 83.7 to 68.2 days when fragmentation is on — which is the right direction.

### Step 6 — Depth structure with variable $K_z(z)$

I switched from a constant water column to four depth-dependent profiles: $K_z(z)$, $T(z)$, $S(z)$, $\rho(z)$. At this step, only $K_z(z)$ feeds into the transport. The others are stored for later use.

Profile ranges:

- $K_z(z)$: $1.36 \times 10^{-6}$ to $1.50 \times 10^{-3}$ m²/s
- Temperature: 4.02 to 18.00 °C
- Salinity: 34.20 to 34.97 psu
- Density: 1024.50 to 1026.66 kg/m³

![rerun step6 depth profiles](figures/may06_rerun_step6_depth_profiles.png)

![rerun step6 depth conservation](figures/may06_rerun_step6_depth_conservation.png)

Conservation: $3.56 \times 10^{-5}\%$. The PSD shifts slightly when the depth structure is added, but not erratically — it reflects the varying diffusivity.

### Step 7 — Depth-dependent sinking speed

The last step in the chain uses the depth profiles to correct the sinking speed at each layer. Kinematic viscosity at depth:

$$\nu(z) \;=\; \frac{\mu(T(z))}{\rho(z)}$$

Sinking speed correction:

$$w(z) \;=\; w_\text{ref} \cdot \frac{\nu_\text{ref}}{\nu(z)}$$

This comes from Stokes' law — sinking speed is inversely proportional to viscosity. Colder, denser water at depth is more viscous, so particles sink more slowly as they go deeper. In this run: speed scale goes from 1.000 at the surface to 0.682 at 2000 m, a 32% slowdown.

![rerun step7 sinking scale](figures/may06_rerun_step7_sinking_scale_profile.png)

![rerun step7 sinking conservation](figures/may06_rerun_step7_sinking_conservation.png)

![rerun step7 sinking psd](figures/may06_rerun_step7_sinking_final_column_psd.png)

Conservation: $4.54 \times 10^{-5}\%$. The PSD shifts upward because slower deep sinking keeps more particles in the column. And with fragmentation at `c3 = 0.02`, total number change jumps to **+869.8%** — that is the red flag.

### Budget across all steps

![budget split](figures/apr08_1d_panel_budget_split.png)

Each panel shows three curves: volume in the column (black), volume exported out the bottom (red), and the tracked total — in + out (blue). The blue line stays near 100% in all panels, meaning no volume is being created or destroyed as each process is added. In the step-7 panel the black curve stays higher for longer and the red rises more slowly, which matches the slower deep sinking.

---

## 3. New base 1-D depth run

After the step chain confirmed things were still working, I ran a fresh base case as the starting point for all future work. Settings:

- Sinking law: `kriest_8`
- Column depth: 1000 m
- Size sections: 5 (the full particle size range is split into 5 logarithmically spaced bins, each holding the total concentration for particles of that size class)
- Both $K_z(z)$ and $w(z)$ depth-dependent
- Coagulation and fragmentation both off

The point of this run is just to have a clean transport baseline.

![may06 1d base pulse profiles](figures/may06_1d_base_depth_pulse_profiles.png)

![may06 1d base depth size snapshots](figures/may06_1d_base_depth_size_snapshots.png)

The pulse center moves down steadily over time. Larger particles reach greater depth faster, as expected. Trust checks:

- `max CFL = 0.054`
- `neg_count = 0`
- Biovolume change: $-0.000\%$ (pulse has not reached the bottom in 60 days)

![may06 1d base conservation](figures/may06_1d_base_depth_conservation.png)

This is the reference all process-addition runs are compared against.

---

## 4. Isolating the fragmentation problem

I ran a four-case matrix to find exactly what causes the +869% number jump. The same depth-dependent setup ran for 60 days, with coagulation and fragmentation toggled independently:

1. Transport + sinking only
2. + Coagulation
3. + Fragmentation
4. + Both

![may06 matrix total number](figures/may06_step7_matrix_total_number.png)

![may06 matrix conservation](figures/may06_step7_matrix_conservation.png)

![may06 matrix small-size volume](figures/may06_step7_matrix_small_size_volume.png)

Results:

| Case | Total number change | neg_count |
|------|---------------------|-----------|
| Sink only | −28.4% | 0 |
| Sink + coag | −28.4% | 0 |
| Sink + frag | +869.8% | 0 |
| Sink + coag + frag | +869.8% | 0 |

Coagulation does nothing to the number inventory on this timescale. Fragmentation alone causes the explosion — and it is the same whether coagulation is on or off. All four cases have `neg_count = 0` and conservation errors below $10^{-4}\%$. This is not a numerical problem. It is entirely the fragmentation parameter.

---

## 5. Why fragmentation explodes: the rate imbalance

To understand the +869%, I computed the fragmentation rate coefficient per bin and compared it to the sinking removal rate — no simulation needed, just grid math (`run_may06_frag_rate_analysis.m`).

The fragmentation rate at bin $k$ is $r_\text{frag}(k) = c_3 \cdot c_4^k$. With `c3 = 0.02` and `c4 = 1.45`, this grows exponentially. The sinking removal rate is $r_\text{sink}(k) = w_k / H$.

![frag rate analysis](figures/may06_frag_rate_analysis.png)

The figure shows the sinking rate (solid black) and four fragmentation rate curves for different `c3` values. At `c3 = 0.02`, fragmentation is between 40 and 2900 times faster than sinking depending on the bin. At the middle bins where most particles actually sit, the ratio is around 100 to 200. Fragmentation is creating small particles far faster than sinking can remove them — that is why the inventory explodes.

To bring the two rates into balance, `c3` needs to be in the range $10^{-4}$ to $2 \times 10^{-4}$. Even at `c3 = 0.0002`, the fragmentation curve still runs near or above the sinking curve at the large-size end. There is no clean calibration with `c4 = 1.45` that makes the legacy term physically consistent across the full size range.

---

## 6. c3 parameter sweep

To find a safe working range numerically, I ran the four-case matrix at eight `c3` values from $10^{-4}$ to $0.02$, each for 60 days in the 0-D slab model.

![c3 sweep matrix](figures/may06_c3_sweep_matrix.png)

The figure shows total number inventory change as a function of `c3`. Key readings:

- Below `c3 = 0.002`: `sink + frag` stays within about 2% of `sink only` — fragmentation is not overwhelming sinking
- Above `c3 = 0.002`: the combined coag + frag case becomes non-monotonic — it goes positive at `c3 = 0.005` and `0.01`, then collapses at `c3 = 0.02`. This is not a physically meaningful regime.
- `c3 ≤ 0.002` is the safe range; `c3 ≤ 0.001` is safer if both coagulation and fragmentation are running together

The negative concentration counts (930–1400) are at $10^{-30}$ level — floating-point noise in empty bins, not real instability.

---

## 7. Operator-split fragmentation

In parallel with the `c3` sweep, I tested the operator-split fragmentation mode. This is based on Alldredge et al. (1990), where the critical aggregate size depends on local turbulence:

$$d_\text{max} \;=\; C \cdot \varepsilon^{-\gamma}$$

Aggregates larger than $d_\text{max}$ are broken apart by turbulent shear; smaller ones are stable. The constants I used are $C = 3.0$ mm and $\gamma = 0.15$, which sits between Alldredge's measurements ($\gamma = 0.11$ for Chaetoceros, $\gamma = 0.29$ for Nitzschia). This mode only *redistributes* volume into smaller bins — it never creates or destroys mass — so biovolume should be conserved exactly.

I ran the four-case matrix with `disagg_mode = 'operator_split'` and swept $d_\text{max}$ from 0.5 cm to 5.0 cm.

![operator split matrix](figures/may06_operator_split_matrix.png)

Results across all four $d_\text{max}$ values:

- `neg_count = 0` everywhere (compared to 931–1388 with the legacy term)
- Biovolume: `sink + frag` matches `sink only` to within 0.002% — truly mass-conserving

All four $d_\text{max}$ values give the same result because the initial spectrum has no particles large enough to cross even the smallest threshold ($d_\text{max} = 0.5$ cm = 5 mm). The fragmentation operator checks each bin against $d_\text{max}$ and only acts on bins where $d > d_\text{max}$. With no large aggregates present, it has nothing to do. To actually see the $d_\text{max}$ sensitivity, I need to first run coagulation long enough to build large aggregates, then add fragmentation. That is the next experiment.

The operator-split mode is the right choice going forward — it conserves mass, produces no negatives, and the fragmentation threshold is directly tied to a physical quantity ($\varepsilon$).

---

## 8. New OOP classes for 1-D column runs (Phase 1)

The column work up to this point used a mix of scripts and the existing `CoagulationSimulation` class. I wrote five new OOP classes in `1d-model-testing/src/` to make column runs cleaner and easier to extend.

**`ColumnGrid`** — holds the depth grid. Takes total depth $H$ and number of cells $n_z$, gives back cell midpoints, face positions, and spacing $\Delta z$.

**`DepthProfile`** — stores the four depth-dependent ocean fields: $K_z(z)$, $T(z)$, $S(z)$, $\rho(z)$, and kinematic viscosity $\nu(z) = \mu(T(z))/\rho(z)$. Two constructors:

- `DepthProfile.flat(n_z)` — constant water column
- `DepthProfile.typical(n_z)` — warm mixed layer + thermocline + cold deep water

It also provides three scaling functions (`brownianScale`, `shearScale`, `dsScale`) that return per-layer correction factors for the kernel scaling in `ColumnRHS`.

**`ColumnTransport`** — one explicit Euler transport step. The upwind advection flux at the downward face of cell $k$:

$$F_k \;=\; w(k)\,\max(Y_k,\,0)$$

The $\max(Y_k, 0)$ prevents negative concentrations from producing unphysical upward fluxes. The advective tendency:

$$\left.\frac{dY_k}{dt}\right|_\text{adv} \;=\; \frac{F_{k-1} - F_k}{\Delta z}$$

Top boundary: $F_0 = 0$ (no inflow from above). Bottom: open exit. Diffusion in flux form:

$$\left.\frac{dY_k}{dt}\right|_\text{diff} \;=\; \frac{1}{\Delta z}\left(K_{k+1/2}\,\frac{Y_{k+1} - Y_k}{\Delta z} - K_{k-1/2}\,\frac{Y_k - Y_{k-1}}{\Delta z}\right)$$

where $K_{k+1/2} = \frac{1}{2}(K_z(k) + K_z(k+1))$, converted from m²/s to m²/day. Full update:

$$Y_\text{new} \;=\; \max\!\left(Y + \Delta t\,(T_\text{adv} + T_\text{diff}),\; 0\right)$$

**`ColumnRHS`** — one full time step: transport first, then depth-scaled process rates at each layer. The design challenge was how to apply depth-dependent kernel scaling without rebuilding the full beta matrices at every layer and every step. The solution: store the three component matrices pre-scaled in the constructor:

```matlab
f_brown = alpha * conBr  * day_to_sec * scale_brown;
f_shear = alpha * gamma  * day_to_sec * scale_shear;
f_ds    = alpha * setcon * day_to_sec * scale_ds;

b25_brown = f_brown .* raw_brown.b25;
b1_brown  = f_brown .* raw_brown.b1;
```

These factors match exactly what `combineAndScale` uses in the 0-D code, so when all depth scales equal 1, the combined matrix recovers the flat-column result exactly. At depth layer $k$, the combined scaled matrix is assembled as:

$$\mathbf{B}_{25}^{(k)} \;=\; s_\text{Br}(k)\,\mathbf{B}_{25}^\text{Br} \;+\; s_\text{sh}(k)\,\mathbf{B}_{25}^\text{sh} \;+\; s_\text{DS}(k)\,\mathbf{B}_{25}^\text{DS}$$

This costs only three multiplies and two adds per layer instead of recomputing from scratch. Process rates use sub-stepping (`proc_substeps = 10` by default) for explicit Euler stability, with non-negativity clipping at each sub-step.

**`ColumnSimulation`** — wraps everything: owns the grid, profile, and RHS, initializes the state, and runs the main time loop. Scripts only interact with this class and the config.

Two existing files also needed changes:

- `SimulationConfig`: added `proc_substeps = 10` property
- `CoagulationRHS`: added `evaluateScaled(t, v, b25_scaled, b1_scaled)` — same logic as `evaluate()` but takes caller-supplied beta matrices, which is how `ColumnRHS` passes per-layer matrices without a separate `CoagulationRHS` object per depth

**Pulse test (Phase 1 check).** A concentration pulse is injected at the top and tracked as it sinks. Results:

- `neg_count = 0`
- Biovolume change: $-0.000\%$ (pulse still in column at 60 days)
- Max CFL: 0.054

![phase 1 pulse test](figures/may06_column_pulse_test.png)

Phase 1 is complete.

---

## 9. Depth-dependent kernel scaling (Phase 2)

With the classes working, the next thing to verify was whether the depth-scaling was correct. The key question: when all depth scale factors are 1, does the column model recover exactly the same result as the 0-D model?

I ran `ColumnSimulation` with `DepthProfile.flat()` (constant $T$, $S$, $\rho$, $K_z$) and compared it to a 0-D `CoagulationSimulation` run with the same config and no transport. The relative L2 difference came out at $3.40 \times 10^{-4}$ — well within the tolerance expected from different time-stepping methods. The scaling is correct.

For the ocean-profile test, I ran flat vs `DepthProfile.typical()` side by side. In the ocean profile, the midwater layer sits in the thermocline where the water is warmer and less viscous. This should increase the Brownian collision rate there relative to a cold flat column. The midwater biovolume in the ocean profile came out **41.89% higher** than in the flat-profile run — in the right direction and physically reasonable. The depth-dependent kernel scaling is active and doing what it should.

![phase 2 depth scaling](figures/may06_phase2_depth_scaling.png)

---

## Where things stand

The base 1-D column is in good shape. Transport is clean, the step chain re-verifies correctly on the May 06 code, and the depth-dependent sinking and coagulation kernels are both active and cross-validated. The fragmentation problem is understood: the legacy term at `c3 = 0.02` is unstable because the exponential rate coefficient grows far faster than the sinking removal rate. Safe range is `c3 ≤ 0.002`, and the operator-split mode is preferred going forward because it conserves mass, has no negatives, and uses a physically grounded breakup criterion.

