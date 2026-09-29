# STRUCTURE_MIXING

A ROMS option that adds drag from submerged structures (think turbine
foundations, bridge piers, floating platforms) and feeds the energy that
drag removes from the mean flow back into the GLS turbulence closure as
extra mixing. It's based on **Rennau, Schimmels & Burchard (2012)**.

---

## 1. The physics

Structures are represented as vertical cylinders. Two things
happen wherever cylinders are present in a grid cell:

**Step 1 — momentum drag.** The structures slow the flow down:

$$G_d^u = -\tfrac{1}{2} C_D\, a\, u\, \sqrt{u^2+v^2}, \qquad
  G_d^v = -\tfrac{1}{2} C_D\, a\, v\, \sqrt{u^2+v^2}$$

- $C_D$ — a drag coefficient (dimensionless)
- $a$ — frontal area density of the cylinders, in m⁻¹: $a = Xd/A_{\text{cell}}$,
  where $X$ is the number of cylinders in a grid cell, $d$ their diameter,
  and $A_{\text{cell}}$ the horizontal area of the grid cell.

**Step 2 — turbulence injection.** The kinetic energy removed from the mean
flow, $P_d = \tfrac{1}{2} C_D\, a\,(u^2+v^2)^{3/2}$, is injected into the GLS turbulence closure, boosting
both TKE ($k$) and the generic length-scale variable ($\psi$):

$$\frac{\partial k}{\partial t} = \ldots + P + P_d + G - \varepsilon$$
$$\frac{\partial \psi}{\partial t} = \ldots + \frac{\psi}{k}\left(c_{\psi1}P + c_{\psi4}P_d + c_{\psi3}G - c_{\psi2}\varepsilon\right)$$

The only new closure coefficient is $c_{\psi4}$ (`gls_c4`), which plays the
same role for structure drag that $c_{\psi1}$ plays for shear production.

`a` is a fully general 3D field, so the same code handles
bottom-mounted monopiles, floating turbine near the surface, or
structures at any depth in between — wherever `a > 0`, drag and mixing are
applied, no assumptions about orientation.

---

## 2. Requirements

- `STRUCTURE_MIXING` **requires** `GLS_MIXING` . It is not implemented for for example `MY25_MIXING`. This is enforced
  at compile-check time. ROMS will refuse to run otherwise.
- Forward (nonlinear) model only. There's no tangent-linear/adjoint support, so it can't currently be used with 4D-Var.

---

## 3. Setting it up

### a) Turn it on

In your application header (e.g. `mycase.h`):
```c
#define STRUCTURE_MIXING
#define GLS_MIXING
```

### b) Two new scalar parameters in your `.in` file

Add these next to the other `GLS_*` coefficients:
```
    STR_CD == 0.63d0        ! cylinder drag coefficient C_D
    GLS_C4 == 1.4d0         ! structure-drag GLS coefficient c_psi4
```
`GLS_C4 = 1.4` is what Rennau et al. (2012) used as realistic mixing for k-ε. Since GLS
coefficients are closure-dependent, so if you switch closures it should be calibrated.

### c) A new 3D field in your grid file: `str_a`

This is where the structures live. You need to add a variable
called `str_a` to your ROMS grid NetCDF file:

| | |
|---|---|
| Dimensions | `(s_rho, eta_rho, xi_rho)` |
| Units | m⁻¹ |
| Meaning | frontal area density, `a = X·d / A_cell` |
| No structures value | `0.0` |
| Land cells | should also be `0.0` (masked out) |

If `str_a` isn't found in the grid file, ROMS just
sets it to zero everywhere and print a warning. The run then behaves exactly
like `STRUCTURE_MIXING` was off.

There's no built-in tool in ROMS to build this field. You'll need
to compute it yourself and add it to the grid file. The vertical levels you populate with nonzero values determine where the structure "lives" (all levels for a monopile, top levels for a floating structure, etc.).

---

## 4. Worked example

Say you have a single monopile foundation, diameter $d = 8$ m, sitting in
one grid cell of horizontal area $A_{\text{cell}} = 1/(pm \cdot pn) = 250000$ m²
(a 500 m × 500 m cell), and it only occupies the top 4 sigma layers out of
20. Then, for that column:

$$a = \frac{X d}{A_{\text{cell}}} = \frac{1 \times 8}{250000} = 3.2\times10^{-5}\ \text{m}^{-1}$$

You'd set `str_a(i,j,k) = 3.2e-3` for `k = 17,..., 20` (top four rho-levels) and
`str_a(i,j,k) = 0.0` for `k = 1,...,16` at that `(i,j)`, and `0.0` everywhere
else in the domain — including at masked/land cells.

A minimal Python snippet to add this to an existing grid file with
`xarray`, given a boolean 3D mask `structure_mask(s_rho, eta_rho, xi_rho)`
of where the cylinder occupies each cell, cylinder diameter `d`, and count
`X` per cell:

```python
import xarray as xr

ds = xr.open_dataset("my_grid.nc", mode="a")

A_cell = 1.0 / (ds["pm"] * ds["pn"])          # cell area [m^2]
a = xr.where(structure_mask, X * d / A_cell, 0.0)
a = a.where(ds["mask_rho"] > 0, 0.0)           # zero out land cells

ds["str_a"] = a.transpose("s_rho", "eta_rho", "xi_rho")
ds["str_a"].attrs = {"long_name": "structure frontal area density",
                      "units": "meter-1"}
ds.to_netcdf("my_grid_with_structures.nc")
```

Then in your `.in` file:
```
    STR_CD == 0.63d0
    GLS_C4 == 1.4d0
```

The rest (drag + mixing) happens automatically
wherever `str_a > 0`.

---

## 5. What happens under the hood

- **Momentum drag** is added in `rhs3d.F`, inside the main vertical loop,
  at every interior rho-level. It's folded automatically into the
  barotropic forcing (`rufrc`/`rvfrc`) — no changes needed to the 2D mode.
- **TKE/GLS production** (`Pd_struct`) is added in `gls_corstep.F`, in the
  same block that already computes shear and buoyancy production, at
  W-points.
- For efficiency, `str_a` is combined once with the grid metrics into a
  precomputed field `GRID(ng)%str_a_omn = str_a / (pm·pn)`, computed a
  single time in `metrics.F` at initialization (the same way ROMS
  precomputes `fomn = f/(pm·pn)`). This keeps the per-timestep drag
  calculation in `rhs3d.F` from re-doing that division every step.
- `str_cd` and `gls_c4` are written to your history/restart/average output
  files (so runs stay self-documenting/reproducible), and both the classic
  NetCDF (`nf90`) and parallel I/O (`PIO`) code paths are supported.

A quick note on energy bookkeeping: the drag term (`rhs3d.F`) and the TKE
injection (`gls_corstep.F`) are computed at slightly different points
(rho- vs. W-points) with slightly different velocity averaging, so the
energy removed from the mean flow and the energy injected into turbulence
will be *close* but not bit-identical. This is the same kind of
approximation ROMS already makes between shear production and viscous
dissipation. So don't expect exact energy
conservation in a budget analysis.

---


## 6. Good to know

- **Nested grids:** `str_cd` and `gls_c4` are one value per grid; `str_a`
  must be supplied separately in each nested grid's NetCDF file.
- **Diagnostics:** there's no dedicated `DIAGNOSTICS_UV` output for the
  structure drag term yet (e.g. no `M3fstr`/history-file diagnostic).
- **Numerical stability:** the drag term is treated explicitly/semi-
  implicitly (current code uses the flow speed and velocities both at the
  same, most recent time level, `nrhs`, inside `rhs3d.F`). This is stable
  for realistic drag coefficients and structure densities. If you use an
  unusually large `str_a` or `str_cd` and see instability, this is a good
  first place to look.
- **Rennau et al. (2012)** studied bottom-mounted structures; using this
  option for floating structures (e.g. floating wind turbine foundations) assumes the
  same drag + mixing physics apply, which is a reasonable
  but unvalidated extrapolation.

---

## Reference

Rennau, H., Schimmels, S., & Burchard, H. (2012). *On the effect of structure-
induced resistance and mixing on inflows into the Baltic Sea: A numerical
model study.* Ocean Dynamics.
