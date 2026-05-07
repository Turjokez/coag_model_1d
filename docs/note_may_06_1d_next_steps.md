# Note - May 06, 2026 - 1-D next steps

The base 1-D depth-dependent model is now working as a stable starting point. In this base, depth dependence is active through `Kz(z)` in transport and through `w(z)` in sinking speed. The pulse test moves down with time in a physically reasonable way, there are no negative concentrations, and tracked-volume conservation error stays very small.

To keep progress clear, the latest matrix run used four simple cases on the same setup: `transport + sink`, `transport + sink + coag`, `transport + sink + frag`, and `transport + sink + coag + frag`. This was done to isolate which process changes the system most. The main result is that fragmentation controls the large jump in total number, while transport and coagulation cases stay close to each other in this test.

So the next work is one-by-one. First, keep the same stable depth-dependent base and tune fragmentation parameters (`c3`, `c4`) to reduce unrealistic number growth. Second, re-run the same matrix after each tuning change and keep the same checks (`neg_count`, conservation error, total number trend). Third, after fragmentation is in a reasonable range, add the next biology terms on top of this same base model.
