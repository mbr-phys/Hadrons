# HQE operator insertions and a combined heavy-meson 3pt contraction

## Summary

This adds the code for Tiers 1 and 2 of `lifetimes_round2.md` §8:

- **Tier 1:** three small modules that insert HQE operators on a propagator (covariant derivatives,
  clover field strength, gluonic insertions), plus one shared header that defines the lattice conventions.
- **Tier 2:** `MContraction::HeavyMeson3pt`, one module per (source wall, sink wall, flow time). It computes:
  - the two-quark HQE three-point functions ($\bar QQ$, $\bar Q\gamma^0Q$, $\mu_\pi^2$, $\mu_G^2$, $\rho_D^3$, $\rho_{LS}^3$ and their remainders);
  - four-quark three-point functions with arbitrary $\Gamma_1\otimes\Gamma_2$ and derivative legs. This covers the dim-6 bag
    parameters, the dim-7 $\Delta F=2$ operators ($R_2$, $R_3$, $\tilde R$, index-coupled) and the $\Delta F=0$ dim-7 operators.

No existing module is modified. **Nothing has been compiled or run yet.**

## Motivation

The existing modules can only build these correlators as chains: `RHQInsertionI` (repeatedly), then
`Gamma3pt` or `WeakNonEye3pt`. That route:
- holds a full propagator per intermediate step;
- recomputes the same derivatives for every operator;
- needs thousands of modules per configuration;
- cannot express $\Gamma_1\neq\Gamma_2$ or index-coupled four-quark operators.

`HeavyMeson3pt` differentiates each propagator leg once per direction. It streams over directions and shares
the result across all two- and four-quark operators. The standalone Tier-1 modules remain as the reference
implementation and for ad-hoc use.

## New files

| File | Content |
|---|---|
| `Hadrons/Modules/MContraction/HQEOperators.hpp` | Header-only helpers in `Hadrons::HQE` (conventions below) |
| `Hadrons/Modules/MUtilities/CovariantDerivative.{hpp,cpp}` | `MUtilities::CovariantDerivative{Propagator,Fermion}` |
| `Hadrons/Modules/MGauge/FieldStrength.{hpp,cpp}` | `MGauge::FieldStrength` |
| `Hadrons/Modules/MUtilities/GaugeInsertion.{hpp,cpp}` | `MUtilities::GaugeInsertion{Propagator,Fermion}` |
| `Hadrons/Modules/MContraction/HeavyMeson3pt.{hpp,cpp}` | `MContraction::HeavyMeson3pt` |
| `tests/Test_heavy_meson_3pt.cpp` (+ `tests/Makefile.am`) | Validation workflow |

`Modules.hpp`/`modules.inc` are git-ignored. Re-run `bootstrap.sh` (or `Hadrons/make_module_list.sh`) before building.

## Conventions (`HQEOperators.hpp`)

All are Euclidean, in lattice units, with Grid's gamma matrices.

| Object | Definition |
|---|---|
| symmetric $D_\mu f$ | $\tfrac12[U_\mu(x)f(x+\hat\mu)-U_\mu^\dagger(x-\hat\mu)f(x-\hat\mu)]$ |
| forward / backward | $U_\mu(x)f(x+\hat\mu)-f(x)$ / $f(x)-U^\dagger_\mu(x-\hat\mu)f(x-\hat\mu)$ |
| Laplacian over `dirs` | $\sum_\mu[U_\mu f(x+\hat\mu)+U^\dagger_\mu f(x-\hat\mu)-2f(x)]$ (nearest neighbour) |
| temporal BC | for $\mu=t$, the link on the last time slice is multiplied by `bc` (default $-1$) |
| $F_{\mu\nu}$ | Grid clover `WilsonLoops::FieldStrength`, made traceless; anti-Hermitian, $\approx i a^2 gG_{\mu\nu}$ up to the sign convention for $A_\mu$ |
| adjoint $D_\mu M$ | $\tfrac12[U_\mu M(x+\hat\mu)U_\mu^\dagger - U^\dagger_\mu(x-\hat\mu)M(x-\hat\mu)U_\mu(x-\hat\mu)]$ |
| $\sigma_{\mu\nu}$ | $\tfrac i2[\gamma_\mu,\gamma_\nu]=i\gamma_\mu\gamma_\nu$ ($\mu\neq\nu$) |
| $\sigma\!\cdot\!B$, $\sigma\!\cdot\!E$ | $\sum_{i<j<3}\sigma_{ij}F_{ij}$, $\sum_{j<3}\sigma_{j3}F_{j3}$ |

Minkowski phases ($D_0=iD_4$, $G_{0i}=iG_{4i}$, the factor $i$ in $F$) are **not** applied anywhere. They belong to the
analysis and should be fixed by the sign checks in `lifetimes_round2.md` §7.

## Tier 1 modules

- **`CovariantDerivative{Propagator,Fermion}`**: `q`, `gauge`, `type` = `symmetric|forward|backward|laplacian`, `dirs`, `bc`.
  - Derivative types take exactly one direction.
  - `laplacian` with `dirs = "0 1 2"` gives the spatial $\mathbf D^2$.
- **`FieldStrength`** (outputs a `GaugeLinkField`): `gauge`, `mode` = `clover|divergence`, `mu`, `nu`, `dirs`.
  - `divergence` gives $\sum_{\rho\in{\rm dirs}}[D_\rho,F_{\rho\nu}]$. For example `nu = 3`, `dirs = "0 1 2"` is the
    Hermitian Darwin insertion.
- **`GaugeInsertion{Propagator,Fermion}`**: `q`, `gauge`, `matrix`, `mode` = `matrix|sigmaB|sigmaE`, `gamma`.
  - Output is $\Gamma\,M\,q$, $\Gamma\,\sigma\!\cdot\!B\,q$ or $\Gamma\,\sigma\!\cdot\!E\,q$.

All three can be chained with each other and with `RHQInsertion*`, and fed to `Gamma3pt`.

## `MContraction::HeavyMeson3pt`

### Parameters
```
gauge      flowed gauge field (links and F at the operator flow time)
bc         temporal fermion BC for time derivatives (+1/-1, default -1)
flowTime   label written to the result
hqe        { qSrc, qSnk, spectator, sinkNoise, tSnk, gammaSnk, gammaSrc, gammas, ops }   (skipped if qSnk empty)
fourQuark  [ { name, qBarLeft, qLeft, qBarRight, qRight, gammaIn, gammaOut, gammaPairs,
               derivatives: [ { legs, dirs, couple } ] } ]
output     result file stem
```

### Two-quark part
$$
C[X](t)=\sum_{\mathbf x}\mathrm{tr}\big[W(x)\,X\,q_3(x)\big],\quad
W=\Gamma_A K\Gamma_B\,\gamma_5q_2^\dagger\gamma_5,\quad
K=\sum_{\mathbf z,\,z_t=t_{\rm snk}}\eta^\dagger(z)S_l(z)
$$
- $q_2$ = `qSrc` (flowed heavy quark from the source wall), $q_3$ = `qSnk` (flowed heavy quark from the sink wall),
  $S_l$ = `spectator` (*unflowed* light quark from the source wall), $\eta$ = `sinkNoise` (the sink-wall source field,
  smeared if appropriate). With an empty `sinkNoise` the sink is a plain wall sum.
- With $\Gamma_A=\Gamma_B=\gamma_5$ and $X=\Gamma$ this is exactly `Gamma3pt`. The conjugated noise sink is built in, so
  Tier-0 item C1 is not needed for this module.
- `gammas`: plain $\Gamma$ insertions (`all` = 16), stored with `op = "gamma"`.
- `ops` (or `all`):

| op | $X$ | notes reference |
|---|---|---|
| `D2` | $\sum_i D^-_iD^+_i$ | $\mu_\pi^2$ |
| `D4`, `D4L` | $D_4$ on $Q$; $D_4$ on $\bar Q$ (left-acting) | mass check, Re parts |
| `D4D4`, `D4D4D4` | $D^-_4D^+_4$, $D_4D^-_4D^+_4$ | full-$D$ appendix |
| `D2D4` | $\sum_iD^-_iD^+_iD_4$ | $\mathcal R_D$ |
| `sigmaB`, `sigmaE` | $\sigma\!\cdot\!B$, $\sigma\!\cdot\!E$ | $\mu_G^2$, $\bar Q\mathcal EQ$ |
| `sigmaBD4`, `sigmaED4` | $\sigma\!\cdot\!B\,D_4$, $\sigma\!\cdot\!E\,D_4$ | $\mathcal R_{LS}$, full-$D$ $\rho_{LS}$ |
| `DiE` | $\sum_iD_i(F_{i3}\,\cdot)$ | $\rho_D^3$ invariant part |
| `divE` | $\sum_i[D_i,F_{i3}]$ | $\rho_D^3$, Hermitian form |
| `DiSigmaE` | $\sum_{ij}D_i(\sigma_{ij}F_{j3}\,\cdot)$ | $\rho_{LS}^3$ invariant part |

Spatial derivatives use exact summation by parts on the time slice, e.g. $\sum_{\mathbf x}\mathrm{tr}[W D^-_iD^+_iq_3]=-\sum_{\mathbf x}\mathrm{tr}[(W\overleftarrow D^+_i)(D^+_iq_3)]$,
so each leg is differentiated at most once per direction and only nearest-neighbour stencils appear.
This assumes periodic spatial boundaries (no twist).

### Four-quark part
Here $L=q_L\,\Gamma_{\rm in}^\dagger\gamma_5q_{\bar L}^\dagger\gamma_5$ and $R=q_{\bar R}\,\Gamma_{\rm out}\gamma_5q_R^\dagger\gamma_5$.
The module outputs `trace 1` $=\mathrm{tr}[L\Gamma_1R\Gamma_2]$ and `trace 2` $=\mathrm{tr}[L\Gamma_1]\,\mathrm{tr}[R\Gamma_2]$.
With $\Gamma_1=\Gamma_2$ this is exactly `WeakNonEye3pt`.
- `gammaPairs`: `diag` (16 $G\otimes G$, default), `all` (256), or `"(G1 G2)(G1 G2)..."`.
- Each `derivatives` entry gives a sum over $\rho\in$ `dirs` (default `0 1 2 3`, Euclidean) of the contraction with a
  symmetric $D_\rho$ on every leg in `legs`. On the adjoint legs (`qBarLeft`, `qRight`) the derivative acts on the antiquark to the left.
  - `couple = G1|G2` additionally replaces $\Gamma_{1,2}\to\gamma_\rho\Gamma_{1,2}$ (index-coupled operators).
  - An empty `legs` gives the underived contraction. No entries also means underived only.
- Results are labelled `name`, `deriv` (e.g. `qBarLeft qLeft@0123`, `qLeft@3`, `qLeft@0123:G2`, `none`), `gamma1`,
  `gamma2`, `trace`. Assembling $R_i$, $\tilde R_i$ (Fierz, Wick assignments, colour structures) is left to the analysis.
  `lifetimes_round2.md` §4.4 gives the leg/topology table for $R_{2,3}$.

### Implementation notes
- **Streaming.** The loop over $\rho$ builds $D_\rho$ of each *distinct* leg propagator once, across all four-quark entries,
  uses it for every specification that needs $\rho$, then reuses the buffer. The HQE part does the same over spatial $i$.
- **Accumulation.** Correlators are accumulated with `sliceSum` per term, so no per-operator lattice accumulators are kept
  for the four-quark part.
- **Memory.** Temporaries (≈6 propagators in the HQE part, one per distinct differentiated leg plus `L`, `R` in the
  four-quark part, and $F_{\mu\nu}$) are allocated locally in `execute` and are not tracked by the environment's memory profiler.

## Testing (`tests/Test_heavy_meson_3pt.cpp`)

The test uses point sources at both ends, so every check is an exact identity on a single configuration. It needs $N_t>4$.
It runs three chains (random gauge, a random gauge transformation of it, unit gauge) and writes everything to
`heavy_meson_3pt.*`. Comparisons to make offline:

1. `hm3pt` `gamma` vs `ref_gamma3pt` (`Gamma3pt` through `MSink::SourceSink`).
2. `D4` vs `ref_D4` (`CovariantDerivative`); vs `0.5*ref_rhq_D4` (`RHQInsertionI`) except at operator times $0$, $N_t-1$.
3. `D2` vs `ref_D2`. This tests the summation by parts against the directly applied Laplacian.
4. `sigmaB`, `sigmaE`, `divE` vs `ref_sigmaB`, `ref_sigmaE`, `ref_divE` (`GaugeInsertion` and `FieldStrength`).
5. `DeltaF2`, `deriv = none` vs `ref_weakNonEye` (both traces, 16 $G$).
6. `hm3pt_g` == `hm3pt` for **all** entries. This is gauge covariance, including every derivative, field-strength and four-quark derivative entry.
7. `hm3pt_unit`: all field-strength entries vanish.

Not covered by an exact identity: the reality and hermiticity properties, $D_4+D_4^L\approx\partial_t C$, and the physical
sign checks ($\mu_\pi^2>0$, $C[iD_0]/C[\gamma^0]\to m_Q$). See `lifetimes_round2.md` §7.

## Known limitations / follow-ups

- Not compiled. Expect small template or expression-template fixes on the first build.
- Spatial summation by parts assumes periodic spatial BCs. Twisted boundary conditions are not supported.
- `bc` only affects time-direction derivative links; the field strength always uses the periodic gauge field.
- The clover $F_{\mu\nu}$ is unimproved. At finite flow time its $O(a^2)$ effects are $O(a^2/\tau)$.
- Not included: the gluon-EOM Darwin route (needs quark loops), eye diagrams, baryons, and a general `MSink::SourceSink`
  adjoint option (C1, only needed for the Tier-0 `Gamma3pt` route with Z2 noise).
- The `HeavyMesonLifetimes` application has not been changed yet.

## Rebasing onto the current `origin/develop`

Branch state when this was written (checked with read-only git commands; local `origin/develop` = `1d64626a`):

- `test/AdjointFlow` tracks `mkb/test/AdjointFlow` (your fork). `origin` is upstream `aportelli/Hadrons`.
- 20 branch commits since the merge base `cc1a68a2`, including one merge commit (`d149698c`, a merge of the fork's copy of this branch). Upstream has 29 commits since the same base.
- **Uncommitted:** `tests/Makefile.am` (modified) and the new files from this PR (untracked). `notes/`, `lifetimes_round2.md` and `newCode.md` are also untracked.
- **Overlap with upstream:**
  - Only `tests/Makefile.am` is changed on both sides. Upstream adds `Test_gamma3pt`; keep both entries.
  - Upstream also rewrote `MContraction::Gamma3pt`, which `tests/Test_heavy_meson_3pt.cpp` uses (see step 5).
  - No upstream path clashes with the new untracked files.

### 1. Safety net
```bash
git branch backup/AdjointFlow-pre-rebase          # current HEAD, restore with: git reset --hard backup/AdjointFlow-pre-rebase
```

### 2. Commit this work (rebase needs a clean tree)
```bash
git add Hadrons/Modules/MContraction/HQEOperators.hpp \
        Hadrons/Modules/MUtilities/CovariantDerivative.{hpp,cpp} \
        Hadrons/Modules/MGauge/FieldStrength.{hpp,cpp} \
        Hadrons/Modules/MUtilities/GaugeInsertion.{hpp,cpp}
git commit -m "HQE insertion modules: CovariantDerivative, FieldStrength, GaugeInsertion"
git add Hadrons/Modules/MContraction/HeavyMeson3pt.{hpp,cpp} tests/Test_heavy_meson_3pt.cpp tests/Makefile.am newCode.md
git commit -m "HeavyMeson3pt: combined HQE two-quark and four-quark 3pt contraction + test"
```
Leave `notes/` and `lifetimes_round2.md` untracked, or commit them as you prefer. Untracked files that upstream does not touch are left alone by the rebase.
If you would rather not commit yet, use `git stash push -u` before step 4 and `git stash pop` after it.

### 3. Update the upstream ref
```bash
git fetch origin
git log --oneline HEAD..origin/develop             # what you are rebasing onto
```

### 4. Rebase
```bash
git rebase origin/develop
```
- **Merge commit `d149698c`.** The default rebase linearises history and drops it, replaying the non-merge commits of both
  sides in order. If both sides of that merge contained the same change, git skips the duplicate as empty. A conflict between
  two of your own commits is resolved like any other. The two "fix stochastic meson" commits (`67b1d829`, `bcf85fea`) are
  consecutive, distinct changes, not duplicates; keep both. To squash or reorder first, run `git rebase -i origin/develop`
  in your own terminal (interactive rebase is not available from Claude Code).
- **Conflict in `tests/Makefile.am`.** Keep both sides: upstream's `Test_gamma3pt` entries and this branch's
  `Test_propagator_flow`, `Test_fermion_flow`, `Test_adjoint_flow`, `Test_stochastic_meson` and `Test_heavy_meson_3pt` entries,
  in both `EXTRA_PROGRAMS` and the `_SOURCES`/`_LDADD` blocks.
- After each resolution: `git add <file> && git rebase --continue`. To give up: `git rebase --abort`.

### 5. Required follow-up: new `Gamma3pt` interface
Upstream changed `Gamma3pt::Par::gamma` from `std::string` to `std::vector<std::string>` of `"gSnk gV gSrc"` triads.
It also added `momProjector`, and the result is now `corr[gamma][momentum][time]`. In `tests/Test_heavy_meson_3pt.cpp`, inside the
`gamma3pt` lambda, replace
```cpp
g3Par.gamma = gamma;
```
with
```cpp
g3Par.gamma        = {"Gamma5 " + gamma + " Gamma5"};
g3Par.momProjector = "";
```
Upstream's contraction is $\mathrm{tr}[g_{\rm Snk}\,q_1\,(g_{\rm Src}^\dagger\gamma_5)\,q_2^\dagger\,(\gamma_5g_V)\,q_3]$, which for
$\gamma_5$ vertices is identical to before. `HeavyMeson3pt` therefore still matches it with `gammaSnk = gammaSrc = Gamma5`.
For non-$\gamma_5$ interpolators, `HeavyMeson3pt`'s `gammaSrc` corresponds to upstream's $g_{\rm Src}^\dagger$. Consider
switching to upstream's convention (apply `adj` to `gammaSrc`) so the two modules agree for all vertices. The offline comparison
of check 1 must read the new `Gamma3pt` result layout (momentum index 0).

### 6. Rebuild and check
```bash
./bootstrap.sh                                      # regenerates Hadrons/Modules.hpp and modules.inc (git-ignored)
# reconfigure/build as usual, then build the tests:
make -C <build>/tests Test_heavy_meson_3pt Test_gamma3pt Test_fermion_flow Test_adjoint_flow
```
Rerun the branch tests, then the new checks listed under **Testing**.

### 7. Publish
The rebase rewrites the branch's history, so updating your fork needs a force push:
```bash
git push --force-with-lease mkb test/AdjointFlow
```
Anyone else with a checkout of `mkb/test/AdjointFlow` must reset to the new branch rather than pull/merge.
