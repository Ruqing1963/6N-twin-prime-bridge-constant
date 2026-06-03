# 6N Twin-Prime Bridge Constant (Part VIII)

The arithmetic structure of the Part VI bridge constant C(d), and the
identification of its one global constant as the inverse twin density.

**Background.** Part VI defined the bridge constant C(d) = r(d|ω)/P(N+d twin|ω),
an ω-independent number relating the conditional gap preference to the
right-centre survival, but left its dependence on the separation d unexplained
(C(42)=41.2, C(210)=27.7 — no pattern given).

**The arithmetic structure (S₁₀, 23,988,173 twin centres, d=1..30).**

```
    C(d) = C0 / S(d),     C0 ≈ 62.75
```

where **S(d)** is the Hardy–Littlewood admissibility of the separation — the
standard singular-series factor for a twin pair at centre-step d:

```
    S(d) = prod_{q>3} [ (#{r mod q : r, r+d both q-safe}/q) / ((q-2)/q)^2 ],
    q-safe: r not in dead(q) = {±6⁻¹ mod q}.
```

The product **C(d)·S(d) is constant to CV 0.41%** across all 30 separations,
while C(d) alone ranges over a factor of four (35 to 158). S(d) is cutoff-stable.

**The constant C0 identified: C0 = 1/ρ (inverse twin density).** At the ω-merged
level the definitions give P_merged(d) = ρ·S(d) — the standard two-point
correlation form (density × admissibility) — holding to CV 0.39%. So

```
    C0 = 1/ρ,   ρ = twin-centre line density.
```

Verified across shells: 1/ρ = 50.3 vs fitted C0 = 50.1 on S₉; 62.5 vs 62.8 on
S₁₀ (each to 0.4%). The shell-growth of C0 is exactly the thinning of ρ. **No
free constant remains** in r(d|ω): ρ is itself the Hardy–Littlewood twin density.

**Consequence — the fully resolved gap preference:**

```
    r(d|ω) = S(d) · ρ · P(N+d twin|ω)
```

separation-dependence = classical admissibility S(d); factor-count dependence =
the closed-form survival P (Part VII, K·∏ f_q); the two meet through the twin
density ρ. The conditional 6N theory shares its d-dependence with the classical
singular series, and the bridge constant is the density factor of the two-point
correlation.

> **Scope.** Experimental / computational number theory; results for d=1..30 on
> shells S₉, S₁₀. No claim about the infinitude of twin primes or any k-tuple
> conjecture.

Part I: doi:10.5281/zenodo.20470367 · VI: doi:10.5281/zenodo.20517990 ·
VII: doi:10.5281/zenodo.20518470

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   └── Cd_arith_S10.csv     d, 6d, C_meas, S_adm, C_times_S, C0_over_S  (S10)
├── code/
│   ├── bridge_constant.py   measures C(d), computes S(d), verifies C(d)S(d)=C0
│   │                        and C0=1/rho; emits Cd_arith_S{K}.csv
│   └── make_arith_fig.py    builds the 3-panel figure from ../data
├── figures/                fig_paper8_arith.{pdf,png}
└── paper/                  Chen_6N_Paper8.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Measure C(d), S(d), verify C(d)S(d)=C0 and C0=1/rho. Default S10 (~12 min).
python code/bridge_constant.py            # S10
MAXK=9 python code/bridge_constant.py     # S9 (faster; shows C0 shrinks with rho)

# 2. Figure (reads ../data/Cd_arith_S10.csv).
cd code && python make_arith_fig.py
```

### Conventions (same as Parts I–VII)

- Twin centre: N with 6N−1, 6N+1 both prime. Centre-step d; physical gap 6d.
- dead(q) = {±6⁻¹ mod q}; a centre is q-safe iff its residue ∉ dead(q).
- C(d) = r(d|ω)/P(N+d twin|ω), ω-averaged (ω-independent per Part VI).
- S(d): Hardy–Littlewood two-centre admissibility (singular-series factor), q≤5000.
- ρ = (#twin centres)/(#N in shell) = twin-centre line density.
- Engine: complete segmented-sieve factorisation + deterministic interval-sieve
  primality; S₁₀ twin count 23,988,173 matches Part I; "N+d twin" by binary search.

## License

MIT — see `LICENSE`.
