# opti_rev7.py — conservative-physics fixes

Working checklist from the `solar-uav-design` cross-review (Adam Yang) + the
"How the Solar UAV Optimizer Works" writeup. Ordered by impact. The whole
review points at one conclusion: **night-bus watts and night duration are the
only binding quantity** — mass has ~3 kg slack and span is pinned at the wall,
so fixes that raise night power honestly matter most.

Legend: [x] done · [ ] todo · 🔴 critical · 🟠 high · 🟡 medium/pending

| # | Area | Problem | Fix | Pri | Status |
|---|------|---------|-----|-----|--------|
| 1 | Motor I₀ | Night torque is tiny → no-load current I₀ is a large fraction of night bus; omitting it understates watts | Add `motor_i0` to `i_cruise` | 🔴 | [x] |
| 2 | Prop operating point | Fixed CoP 0.80 + guessed advance ratio; risk of an invented (unphysical) cruise point | FLAG both; tie to real prop map (T,P vs J) later | 🔴 | [ ] flagged |
| 3 | ESC efficiency | 0.95 optimistic below ~50% duty; night cruise ~27% duty | Named FLAGGED `esc_efficiency`; conservative value | 🟠 | [x] |
| 4 | Cell temp derate | Cells ~75 °C at desert noon; −0.32%/°C ⇒ ~16% midday loss uncredited | Per-timestep efficiency `η·[1+γ(T_cell−25)]` | 🟠 | [x] |
| 5 | Battery charge-rate cap | Model banks unlimited afternoon surplus; real pack limited by charge C-rate → refill looks too easy | Per-step charge-power cap (FLAGGED `charge_c_rate`) | 🟠 | [x] |
| 6 | Battery usable Wh derate | Datasheet Wh may not hold under real load/temp (Si-anode) | 8% usable derate (`battery_usable_derate=0.92`); rate-derate N/A at 0.06C | 🟠 | [x] |
| 7 | Battery spec | 6 A is *their* pack, not ours; need the real cell | Amprius SA03 pack locked (483 Wh, 0.25C) | 🟠 | [x] |
| 8 | Objective | Mass has slack, span at wall — minimizing span optimizes a non-binding quantity | Consider max night-energy-margin objective | 🟠 | [ ] design |
| 9 | Battery-true night | Pack must cover every hour solar<load (~14 h), not civil night | Verify low-sun flux honesty; consider pvlib Ineichen | 🟡 | [ ] |
| 10 | Site consistency | We use lat 37.4°; their site El Mirage is 34.65° | Confirm actual site; reconcile lat/alt/temps | 🟡 | [ ] |

## NOT worth fixing (per their doc)
- Chasing airfoil c_d (tenths of a watt vs ~39 W prop/motor/ESC bite)
- Adding more solar cells (noon already spills; adds mass/drag, hurts night)
- Washout (taper λ≈0.6 already puts Oswald e at the 0.98 cap)
- T-tail downwash / flight-cooling / endplate credits (conservative-physics violations)

## Already ahead of their model — keep
Integrated lateral-directional stability, aileron roll-rate sizing,
dihedral/weathercock constraints, smooth gradient (IPOPT) solve.
