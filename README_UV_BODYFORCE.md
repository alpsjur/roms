# UV_BODYFORCE

A ROMS option that lets you apply a prescribed acceleration (optionally varying with depth) to the
`u` and `v` momentum equations, instead of relying on winds, pressure
gradients, or Coriolis to drive the flow.

It exists mainly to build **idealized 1D/single-column test cases**, used for validating turbulence closures and, in particular, the
**`STRUCTURE_MIXING`** drag/mixing parametrization (see
[`README_STRUCTURE_MIXING.md`](README_STRUCTURE_MIXING.md)). The
[`mixtest_1d`](https://github.com/alpsjur/mixtest_1d) project is built entirely around this option.

---

## 1. What it does

`UV_BODYFORCE` adds two things to the momentum equations at every interior
level, in both the `u` and `v` directions:

1. **A prescribed body force** `bfrc_u`/`bfrc_v` (m/s²) — a 3D field
   (rho-vertical-levels, U/V horizontal points) read from the grid NetCDF
   file. It can be uniform or vary with depth, so you can prescribe
   whatever forcing profile (i.e. shear) you want.
2. **A linear (Rayleigh) damping term** `-bfrc_cd * u` (and same for `v`) —
   applied at every level.

Put together, the momentum tendency added at each level is:

$$\frac{du}{dt}\Big|_{\text{bodyforce}} = \text{ramp}(t)\cdot\text{bfrc\_u} - \text{bfrc\_cd}\cdot u$$

The damping term turns a constant forcing into a **target
velocity**: if you set `bfrc_u = bfrc_cd * u_target`, then in steady state
(ignoring other dynamics) $u \to u_{\text{target}}$:

$$\frac{du}{dt} = \text{bfrc\_cd}\,(u_{\text{target}} - u) \;\;\Rightarrow\;\; u(t) = u_{\text{target}}\,(1-e^{-\text{bfrc\_cd}\,t})$$

If instead you set `bfrc_cd = 0`, the damping term vanishes and you get a
pure constant acceleration: $u(t) = \text{bfrc\_u}\cdot t$.  This is handy for
testing the response of another drag term (e.g. `STRUCTURE_MIXING`'s
quadratic drag) in isolation, since then the only thing balancing the
forcing is whatever other drag you've switched on.

A linear ramp from 0 to full strength between `bfrc_tstr` and `bfrc_tend`
(days) is applied on top, so the forcing doesn't switch on with a shock at
`t=0`.

---

## 2. Setting it up

### a) Turn it on

In your application header:
```c
#define UV_BODYFORCE
```
You'll typically pair it with `#undef UV_ADV` and `#undef UV_COR` in
idealized single-column setups, so the body force is the only thing
driving the flow (as done in `mixtest_1d`'s `mixtest_1d.h`).

### b) Three new scalar parameters in your `.in` file

```
    BFRC_TSTR == 0.0d0     ! days, ramp start
    BFRC_TEND == 0.0d0     ! days, ramp end
     BFRC_CD == 0.005d0    ! 1/s, linear damping rate
```
- If `BFRC_TEND <= BFRC_TSTR`, the force is applied at full strength from
  the start (no ramp).
- `bfrc_cd` is a
  linear rate (1/s). Because it's
  applied explicitly, it must satisfy `bfrc_cd * DT << 1` for numerical
  stability (e.g. with `DT = 40 s`, keep `bfrc_cd` well below ~0.02 s⁻¹).
  Set it to `0.0` if you just want a pure constant acceleration with no
  damping.

### c) A new 3D field in your grid file: `bfrc_u`/`bfrc_v`

| | |
|---|---|
| Variable names | `bfrc_u`, `bfrc_v` |
| Dimensions | `(s_rho, eta_u, xi_u)` and `(s_rho, eta_v, xi_v)` |
| Units | m/s² |
| Meaning | prescribed acceleration at U/V-points, per vertical level |

Like `STRUCTURE_MIXING`'s `str_a`, both fields are optional. If
missing from the grid file, ROMS sets them to zero and prints a warning.

The easiest way to build them is from a target steady-state velocity
profile $u_{\text{target}}(z)$ (uniform or depth-varying) using the
relation above:

```python
bfrc_u = bfrc_cd * u_target   # (s_rho, eta_u, xi_u), m/s2
bfrc_v = bfrc_cd * v_target
```

`mixtest_1d`'s `tools/make_grd.py` does exactly this (function
`build_bfrc`), with two modes:
- **`uniform`/`profile`**: derive `bfrc` from a target velocity (constant,
  or read from a `depth  u  [v]` text file and interpolated onto the
  vertical grid), combined with `bfrc_cd` as above.
- **`force`**: bypass the velocity-target relation entirely and write a
  prescribed constant acceleration directly (paired with `BFRC_CD = 0`),
  useful for isolating the response of another drag term.

---


## 3. What happens under the hood

- Everything lives in `rhs3d.F`, just before the main `K_LOOP` (i.e.
  applied once per time step, at every interior level, in both
  directions).
- The prescribed force and linear damping are added to `ru`/`rv` at
  U/V-points, scaled by layer thickness `Hz` and cell area (`1/(pm·pn)`),
  the same way other body forces (bottom stress, vegetation drag) are
  applied. It folds automatically into the barotropic forcing; no
  changes are needed elsewhere (e.g. `step2d.F`).
- The ramp factor is computed once per time step from `tdays`,
  `bfrc_tstr`, `bfrc_tend`.
- `bfrc_u`/`bfrc_v` are read in `get_grid.F` (optional, defaults to zero
  with a warning if absent. Both the classic `nf90` and `PIO` I/O paths
  are supported).
- `DIAGNOSTICS_UV`, when active, accumulates this term into the existing
  `M3vvis` diagnostic slot (there's no dedicated diagnostic index for it).

### Coexisting with `STRUCTURE_MIXING`

Inside a structure zone (`str_a > 0`), the linear `UV_BODYFORCE` damping
and the quadratic `STRUCTURE_MIXING` drag both act on the flow at once.
That means the exact target `u_target = bfrc_u/bfrc_cd` is only reached
where `str_a = 0`; inside a structure zone the extra quadratic drag pulls
the steady state below that target. 

---

## 4. Testing

The `mixtest_1d` test suite validates this option against closed-form
solutions:

| Test | Setup | Analytical check |
|---|---|---|
| `test_UV_BODYFORCE` | `mode: force`, `BFRC_CD = 0`, uniform, `str_a = 0` | $u(t) = F\cdot t$, rtol 1e-5 |
| `test_UV_BODYFORCE_PROFILE` | `mode: profile`, depth-varying target, `str_a = 0` | $u_k(t) = U_0(z_k)(1-e^{-\text{bfrc\_cd}\,t})$ per level, rtol 2e-2 (mixing smooths the prescribed shear a bit, hence the looser tolerance) |
| `test_STRUCTURE_DRAG` | `mode: force`, `BFRC_CD = 0`, `str_a > 0` | $u(t) = \sqrt{F/\alpha}\tanh(t\sqrt{F\alpha})$, i.e. constant force balanced by quadratic structure drag |


