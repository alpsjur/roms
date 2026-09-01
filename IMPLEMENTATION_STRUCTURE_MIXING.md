# Implementation Plan: Parametrization of Structural Mixing in ROMS

**Goal:** Implement drag from submerged structures (e.g. turbine foundations) and the
corresponding turbulent mixing enhancement into ROMS, following Rennau, Schimmels &
Burchard (2012).

---

## Table of Contents

1. [Overview and Physics](#1-overview-and-physics)
2. [Naming Conventions](#2-naming-conventions)
3. [Key Design Decisions](#3-key-design-decisions)
4. [ROMS Code Architecture — Relevant Parts](#4-roms-code-architecture--relevant-parts)
5. [New CPP Define](#5-new-cpp-define)
6. [New Scalar Parameters](#6-new-scalar-parameters)
7. [New 3D Field: Frontal Area Density `str_a`](#7-new-3d-field-frontal-area-density-str_a)
8. [Complete List of Files to Modify](#8-complete-list-of-files-to-modify)
9. [Detailed Change-by-File Guide](#9-detailed-change-by-file-guide)
   - [9.1 `ROMS/Include/cppdefs.h`](#91-romsincludecppdefinesh)
   - [9.2 `ROMS/Modules/mod_scalars.F`](#92-romsmodulesmod_scalarsf)
   - [9.3 `ROMS/Modules/mod_mixing.F`](#93-romsmodulesmod_mixingf)
   - [9.4 `ROMS/Utility/checkdefs.F`](#94-romsutilitycheckdefsf)
   - [9.5 `ROMS/Utility/read_phypar.F`](#95-romsutilityread_phparf)
   - [9.6 `ROMS/Utility/get_grid.F`](#96-romsutilityget_gridf)
   - [9.7 `ROMS/Utility/def_info.F`](#97-romsutilitydef_infof)
   - [9.8 `ROMS/Utility/wrt_info.F`](#98-romsutilitywrt_infof)
   - [9.9 `ROMS/Nonlinear/rhs3d.F`](#99-romsnonlinearrhs3df)
   - [9.9a `ROMS/Modules/mod_grid.F` and `ROMS/Utility/metrics.F` (efficiency optimization)](#99a-romsmodulesmod_gridf-and-romsutilitymetricsf-efficiency-optimization)
   - [9.10 `ROMS/Nonlinear/gls_corstep.F`](#910-romsnonlineargls_corstepf)
   - [9.11 `ROMS/Nonlinear/step3d_uv.F` (optional)](#911-romsnonlinearstep3d_uvf-optional)
10. [User Input File Changes](#10-user-input-file-changes)
11. [NetCDF Grid File Changes](#11-netcdf-grid-file-changes)
12. [Important Numerical Considerations](#12-important-numerical-considerations)
13. [Testing Strategy](#13-testing-strategy)
14. [Things to Keep in Mind](#14-things-to-keep-in-mind)

---

## 1. Overview and Physics

The parametrization adds two coupled effects of submerged structures.

### Step 1 — Momentum drag

Modelled as drag on an array of cylinders:

$$G_d^u = -\tfrac{1}{2} C_D\, a\, u\, \sqrt{u^2+v^2}, \qquad
  G_d^v = -\tfrac{1}{2} C_D\, a\, v\, \sqrt{u^2+v^2}$$

where:
- $C_D$ — drag coefficient (dimensionless);
- $a$ \[m$^{-1}$\] — frontal area density of the cylinders: $a = X d / A_\text{cell}$  
  ($X$ = number of cylinders, $d$ = cylinder diameter, $A_\text{cell}$ = horizontal cell area).

Energy extracted from the mean flow:
$$P_d = \tfrac{1}{2} C_D\, a\, (u^2+v^2)^{3/2} \quad [\text{m}^2\,\text{s}^{-3}]$$

### Step 2 — Injection into the GLS turbulence closure

The extracted energy $P_d$ is injected into both TKE ($k$) and the generic length-scale variable ($\psi$):

$$\frac{\partial k}{\partial t} + \ldots = \mathcal{D}_k + P + P_d + G - \varepsilon$$

$$\frac{\partial \psi}{\partial t} + \ldots = \mathcal{D}_\psi + \frac{\psi}{k}\!\left(c_{\psi1}\,P + c_{\psi4}\,P_d + c_{\psi3}\,G - c_{\psi2}\,\varepsilon\right)$$

$P$ is shear production, $G$ buoyancy production. The only new GLS coefficient is $c_{\psi4}$
(role analogous to $c_{\psi1}$ for shear production but applied to the structure drag source).

---

## 2. Naming Conventions

Consistent naming is critical. All names below follow the patterns already used in ROMS:

| Quantity | Symbol | Fortran name | `.in` keyword | NetCDF name |
|---|---|---|---|---|
| Cylinder drag coefficient | $C_D$ | `str_cd(ng)` | `STR_CD` | `str_cd` (scalar) |
| GLS structure-drag coeff. | $c_{\psi4}$ | `gls_c4(ng)` | `GLS_C4` | `gls_c4` (scalar) |
| Frontal area density field | $a(i,j,k)$ | `MIXING(ng)%str_a` | — | `str_a` |
| CPP define | — | `STRUCTURE_MIXING` | — | — |

**Rationale:**
- `gls_c4` follows exactly the pattern of `gls_c1`, `gls_c2`, `gls_c3m`, `gls_c3p` already in `mod_scalars.F`. These correspond to $c_{\psi1}$–$c_{\psi3}$ in Umlauf & Burchard (2003);
  `gls_c4` is the natural next entry.
- `str_cd` follows the lowercase short-name style of `rdrg`, `rdrg2`.
- `str_a` as a member of `MIXING(ng)` follows the short names of other fields in that struct
  (`Akv`, `bvf`, `tke`, `gls`).

---

## 3. Key Design Decisions

| Decision                      | Choice                                                                           | Rationale                                                                                      |
| ----------------------------- | -------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| CPP guard                     | `STRUCTURE_MIXING`                                                                 | Must co-exist with `GLS_MIXING`; enforced in `checkdefs.F`                                     |
| `str_cd`, `gls_c4`            | Scalars per grid in `mod_scalars.F`                                              | Same pattern as all GLS coefficients                                                           |
| `str_a(i,j,k)`                | 3D pointer in `mod_mixing.F`                                                     | Static, time-invariant; zero outside structure zone                                            |
| Input of `str_a`              | Read from ROMS grid NetCDF                                                       | Analogous to `rdrag2`; constructed externally                                                  |
| Vertical extent of structures | `str_a=0` outside the structure depth range                                      | Works for any orientation: bottom-up (monopile), surface-down (floating turbine), or mid-water |
| Momentum drag location        | `rhs3d.F`                                                                        | All interior 3D body forces live here                                                          |
| TKE/GLS injection             | `gls_corstep.F`                                                                  | Corrector step where all production terms are added                                            |
| Time-stepping                 | Explicit (same as shear production)                                              | Consistent with existing GLS code                                                              |
| Stability                     | Semi-implicit on drag rate (Option C); fully implicit in `step3d_uv.F` if needed | See §12                                                                                        |
| Tangent/Adjoint               | Out of scope                                                                     | Forward model only                                                                             |

---

## 4. ROMS Code Architecture — Relevant Parts

- **`ru(i,j,k,nrhs)` / `rv(i,j,k,nrhs)`** — 3D momentum tendency arrays. All interior body
  forces accumulate here in `rhs3d.F`. At the end of `rhs3d.F` these are summed into the
  barotropic forcing `rufrc`/`rvfrc` automatically — no separate 2D modification needed.
- **`tke(i,j,k,3)` / `gls(i,j,k,3)`** — TKE and GLS variable at W-points (k=0..N,
  interfaces between sigma layers). Produced and dissipated in `gls_corstep.F`.
- **`gls_c1(ng)`..`gls_c3p(ng)`** — Existing GLS closure coefficients in `mod_scalars.F`.
  New `gls_c4(ng)` follows this pattern exactly.
- **`MIXING(ng)%Akv`**, **`tke`**, **`gls`** etc. — Fields in the `T_MIXING` derived type
  in `mod_mixing.F`. New `str_a` is added here.
- **`GRID(ng)%rdrag2`** — Existing 2D quadratic drag coefficient read from the grid NetCDF.
  Template for reading the 3D `str_a` field in `get_grid.F`.
- **`gls_corstep.F` line ~752–779** — The production/dissipation block where `Kprod` and
  `Pprod` (shear + buoyancy) are computed and added to `tke`/`gls`. The structure drag
  injection `Pd_struct` is added immediately after this block.
- **`rhs3d.F` K_LOOP** — Starting around line 498, this loop computes all baroclinic body
  forces at each layer k=1..N. Structure drag is added here.

---

## 5. New CPP Define

Add to `ROMS/Include/cppdefs.h` in the vertical mixing options section (near `GLS_MIXING`):

```
** STRUCTURE_MIXING        if structure-induced drag and mixing (Rennau et al. 2012)    **
```

In the user's application header file (e.g. `mycase.h`):
```c
#define STRUCTURE_MIXING
```

`STRUCTURE_MIXING` requires `GLS_MIXING`. This is enforced in `checkdefs.F` (§9.4).

---

## 6. New Scalar Parameters

Two new scalars per grid, following the `gls_c1`/`gls_c2` pattern:

| Parameter | Fortran name | `.in` keyword | Typical value | Description |
|---|---|---|---|---|
| Cylinder drag coefficient | `str_cd` | `STR_CD` | 1.0 | $C_D$ in the drag formula |
| GLS structure-drag coeff. | `gls_c4` | `GLS_C4` | 1.0 | $c_{\psi4}$ in the $\psi$ equation |

`gls_c4 = 1.0` is recommended for k-ε by Rennau et al. (2012). It controls how strongly
the structure drag modifies the turbulent length scale.

---

## 7. New 3D Field: Frontal Area Density `str_a`

- **Fortran:** `MIXING(ng)%str_a(LBi:UBi, LBj:UBj, N(ng))`
- **Grid location:** rho-points horizontally, rho-levels vertically (k=1..N)
- **Units:** m$^{-1}$; zero where no structures exist
- **Time dependence:** static (this implementation)
- **Source:** read from the ROMS grid NetCDF file (variable name `str_a`)
- **Construction:** by a separate external Python script (see `a_struct_script_notes.md`)

`str_a` is a **general 3D field** — it can be non-zero at any combination of sigma
levels. The ROMS code applies drag and injects TKE wherever `str_a > 0`, with no
assumption about whether structures are bottom-mounted, surface-attached, or at an
intermediate depth. This makes the implementation equally valid for:
- **Bottom-fixed structures** (monopile/jacket wind turbines, bridge piers): non-zero
  `str_a` in the bottom sigma levels.
- **Floating structures** (semi-submersible or spar wind turbines): non-zero `str_a` in
  the upper sigma levels, representing the submerged hull/column.
- **Arbitrary depth range:** any set of sigma levels can be activated.

The field is stored in `mod_mixing.F`. When computing GLS production, `str_a` is averaged
vertically to W-points inside `gls_corstep.F`. When computing momentum drag, `str_a` is
averaged horizontally to u/v-points inside `rhs3d.F`. The averaging is done in the ROMS
code, not in the external script.

> **Update (efficiency optimization):** the momentum drag term in `rhs3d.F` requires
> `str_a/(pm*pn)` (see §9.9 and §12.8 for why this normalization is needed). Dividing by
> `pm*pn` inside the `i,j,k` loop, every time step, is wasteful because `str_a`, `pm`,
> and `pn` are all time-invariant. Instead, the compound field
> `GRID(ng)%str_a_omn(i,j,k) = str_a(i,j,k) * omn(i,j)` (where `omn = 1/(pm*pn)`,
> already available in `mod_grid.F`) is precomputed **once**, at initialization, in
> `ROMS/Utility/metrics.F` — the same routine that computes the analogous compound field
> `fomn = f*omn`. `rhs3d.F` then only needs to average `str_a_omn` to u/v-points; it no
> longer divides by `pm*pn` or accesses `MIXING(ng)%str_a` at all. See §9.9a.
>
> Note: `gls_corstep.F` does **not** use `str_a/(pm*pn)` — it uses `MIXING(ng)%str_a`
> directly (averaged vertically to W-points) for the TKE/GLS production term. This
> optimization only applies to the momentum drag term in `rhs3d.F`.

---

## 8. Complete List of Files to Modify

### Core (required)

| # | File | Change |
|---|---|---|
| 1 | `ROMS/Include/cppdefs.h` | Document `STRUCTURE_MIXING` option |
| 2 | `ROMS/Modules/mod_scalars.F` | Declare, allocate, deallocate `str_cd`, `gls_c4` |
| 3 | `ROMS/Modules/mod_mixing.F` | Declare, allocate, initialise, deallocate `str_a` |
| 4 | `ROMS/Utility/checkdefs.F` | Enforce `GLS_MIXING` co-requirement; print option |
| 5 | `ROMS/Utility/read_phypar.F` | Parse `STR_CD`, `GLS_C4` from `.in` file and print |
| 6 | `ROMS/Utility/get_grid.F` | Read 3D `str_a` from grid NetCDF (nf90 **and** PIO) |
| 6a | `ROMS/Modules/mod_grid.F` | Declare, allocate, initialise, deallocate precomputed `str_a_omn = str_a/(pm*pn)` |
| 6b | `ROMS/Utility/metrics.F` | Compute `GRID(ng)%str_a_omn = MIXING(ng)%str_a * omn` once, at initialization |
| 7 | `ROMS/Utility/def_info.F` | Define `str_cd`, `gls_c4` as variables in output files (nf90 **and** PIO) |
| 8 | `ROMS/Utility/wrt_info.F` | Write `str_cd`, `gls_c4` values to output files (nf90 **and** PIO) |
| 9 | `ROMS/Nonlinear/rhs3d.F` | Add drag body force to `ru`/`rv` at all interior levels, using precomputed `GRID(ng)%str_a_omn` |
| 10 | `ROMS/Nonlinear/gls_corstep.F` | Add `Pd_struct` to tke and gls equations |

### Optional (recommended for numerical robustness or analysis)

| # | File | Change |
|---|---|---|
| 11 | `ROMS/Nonlinear/step3d_uv.F` | Add implicit drag to tridiagonal diagonal (if explicit form causes instability) |
| 12 | `ROMS/Modules/mod_scalars.F` | Add diagnostic indices `M3fstr`, `M2fstr` |
| 13 | `ROMS/Nonlinear/rhs3d.F` | Accumulate `DiaRU`/`DiaRV` for structure drag (`DIAGNOSTICS_UV`) |
| 14 | `ROMS/Utility/def_his.F` | Define `str_drag_u`, `str_drag_v` in history file |
| 15 | `ROMS/Utility/wrt_his.F` | Write structure drag diagnostic to history file |

---

## 9. Detailed Change-by-File Guide

### 9.1 `ROMS/Include/cppdefs.h`

Add one documentation line in the vertical mixing section (near line 237, close to `GLS_MIXING`):

```
** STRUCTURE_MIXING        if structure-induced drag and mixing (Rennau et al. 2012)    **
```

No functional code lives here — this file documents available CPP options only.

---

### 9.2 `ROMS/Modules/mod_scalars.F`

Follow the pattern of the `gls_c1`/`gls_c2` block (around lines 1726–1751 for declarations,
~3396 for allocation, ~4177 for deallocation).

#### Declarations (inside the `GLS_MIXING` block)
```fortran
# ifdef STRUCTURE_MIXING
!    str_cd        Structure cylinder drag coefficient (non-dimensional).
!    gls_c4        Structure-drag GLS coefficient (non-dimensional). Analogous
!                  to gls_c1 for shear production (Rennau et al. 2012, c_psi4).

        real(r8), allocatable :: str_cd(:)
        real(r8), allocatable :: gls_c4(:)
# endif
```

#### Allocation
```fortran
# ifdef STRUCTURE_MIXING
      IF (.not.allocated(str_cd)) THEN
        allocate ( str_cd(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(gls_c4)) THEN
        allocate ( gls_c4(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
# endif
```

#### Deallocation
```fortran
# ifdef STRUCTURE_MIXING
      IF (allocated(str_cd))   deallocate ( str_cd )
      IF (allocated(gls_c4))   deallocate ( gls_c4 )
# endif
```

---

### 9.3 `ROMS/Modules/mod_mixing.F`

Follow the pattern of `Akv`, `bvf` (around lines 237–248 for declaration, allocation, init,
deallocation).

#### Declaration (inside `T_MIXING`, under `GLS_MIXING` block)
```fortran
# ifdef STRUCTURE_MIXING
!  str_a      Frontal area density of structures [m-1] at rho-points.

          real(r8), pointer :: str_a(:,:,:)
# endif
```

#### Allocation
```fortran
# ifdef STRUCTURE_MIXING
      allocate ( MIXING(ng) % str_a(LBi:UBi,LBj:UBj,N(ng)) )
      Dmem(ng)=Dmem(ng)+REAL((UBi-LBi+1)*(UBj-LBj+1)*N(ng),r8)
# endif
```

#### Initialisation
```fortran
# ifdef STRUCTURE_MIXING
      MIXING(ng) % str_a = 0.0_r8
# endif
```

#### Deallocation
```fortran
# ifdef STRUCTURE_MIXING
      IF (.not.destroy(ng, MIXING(ng)%str_a, MyFile,                     &
     &                 __LINE__, 'MIXING(ng)%str_a')) RETURN
# endif
```

---

### 9.4 `ROMS/Utility/checkdefs.F`

Add two blocks. Follow the style of the existing `GLS_MIXING` check (around line 1272).

```fortran
! Require GLS_MIXING when STRUCTURE_MIXING is active.
#if defined STRUCTURE_MIXING && !defined GLS_MIXING
      IF (Master) WRITE (stdout,10)                                       &
     &  'STRUCTURE_MIXING requires GLS_MIXING to be defined.'
      exit_flag=5
#endif

! Print option if active.
#ifdef STRUCTURE_MIXING
      IF (Master) WRITE (stdout,20) 'STRUCTURE_MIXING',                    &
     &  'Structure-induced drag and GLS mixing (Rennau et al. 2012)'
      is=LEN_TRIM(Coptions)+2
      Coptions(is:is+16)=' STRUCTURE_MIXING,'
#endif
```

---

### 9.5 `ROMS/Utility/read_phypar.F`

Follow the pattern of `GLS_C1`/`GLS_C2` (around line 982).

#### CASE parsing
```fortran
# ifdef STRUCTURE_MIXING
            CASE ('STR_CD')
              Npts=load_r(Nval, Rval, Ngrids, str_cd)
            CASE ('GLS_C4')
              Npts=load_r(Nval, Rval, Ngrids, gls_c4)
# endif
```

#### Output summary (in the parameter print section near the other GLS parameters)
```fortran
# ifdef STRUCTURE_MIXING
          WRITE (out,200) str_cd(ng), 'STR_CD',                          &
     &      'Structure cylinder drag coefficient (nondimensional).'
          WRITE (out,200) gls_c4(ng), 'GLS_C4',                          &
     &      'Structure-drag GLS coefficient c_psi4 (nondimensional).'
# endif
```

---

### 9.6 `ROMS/Utility/get_grid.F`

Template: the `rdrag2` reading block (around line 2199 in `get_grid_nf90`, ~4911 in
`get_grid_pio`). The key difference is that `str_a` is 3D: dimensions
`(xi_rho, eta_rho, s_rho)` in the NetCDF file and stored in `MIXING(ng)%str_a`.

Add in **both** `get_grid_nf90` and `get_grid_pio`:

```fortran
# ifdef STRUCTURE_MIXING
!
!  Read in structure frontal area density (str_a) [m-1].
!  This field is 3D (xi_rho x eta_rho x s_rho). It is zero where
!  no structures exist. The variable is optional: if absent from the
!  grid file, str_a remains zero (no structural drag).
!
      CALL netcdf_inq_var (ng, iNLM, GRDname(ng),                        &
     &                     MyFile, __LINE__,                              &
     &                     ncid=ncGRDid(ng),                             &
     &                     VarName='str_a',                              &
     &                     VarId=varid)
      IF (exit_flag.eq.NoError) THEN
        CALL nf_fread3d (ng, iNLM, GRDname(ng), ncGRDid(ng),             &
     &                   'str_a', varid,                                 &
     &                   0, nf90_double, 1, Vsize,                       &
     &                   scale, Fmin, Fmax,                              &
     &                   LBi, UBi, LBj, UBj, 1, N(ng),                  &
     &                   MIXING(ng) % str_a,                             &
     &                   MIXING(ng) % str_a)
        IF (exit_flag.ne.NoError) RETURN
        IF (Master) THEN
          WRITE (stdout,30) 'structure frontal area density: str_a',     &
     &                      Fmin, Fmax
        END IF
      ELSE
!       Variable absent from grid file: str_a remains zero.
        exit_flag=NoError
        IF (Master) WRITE (stdout,*)                                      &
     &    ' WARNING: str_a not found in grid file; str_a set to zero.'
      END IF
# endif
```

> **Note:** Verify the exact calling signature of `nf_fread3d` by cross-checking with
> other 3D variable reads in `get_grid_nf90`. The PIO variant uses `pio_netcdf_*` calls
> instead — follow the same approach used for other 3D fields in `get_grid_pio`.

---

### 9.7 `ROMS/Utility/def_info.F`

Defines the scalar parameters as NetCDF variables in history/restart/average output files.
Follow the `gls_c1` block (around line 1861). Apply to **both** the nf90 and PIO subroutines.

```fortran
# ifdef STRUCTURE_MIXING
      Vinfo( 1)='str_cd'
      Vinfo( 2)='structure cylinder drag coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                    &
     &               1, (/0/), Aval, Vinfo, ncname,                      &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='gls_c4'
      Vinfo( 2)='structure-drag GLS coefficient (c_psi4)'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                    &
     &               1, (/0/), Aval, Vinfo, ncname,                      &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
```

---

### 9.8 `ROMS/Utility/wrt_info.F`

Writes the scalar parameter values into output files. Follow the `gls_c1` write block
(around line 620). Apply to **both** `wrt_info_nf90` and `wrt_info_pio`.

```fortran
# ifdef STRUCTURE_MIXING
      CALL netcdf_put_fvar (ng, model, ncname, 'str_cd',                 &
     &                      str_cd(ng), (/0/), (/0/),                    &
     &                      ncid=ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c4',                 &
     &                      gls_c4(ng), (/0/), (/0/),                    &
     &                      ncid=ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
```

The PIO variant uses `pio_netcdf_put_fvar` — follow the same pattern as `gls_c1` in
`wrt_info_pio`.

> **Why this matters:** Without these two files, `str_cd` and `gls_c4` are silently absent
> from all output files. Runs cannot be reproduced from their outputs.

---

### 9.9 `ROMS/Nonlinear/rhs3d.F`

Adds the structural drag body force to `ru`/`rv` at each interior rho-level.

#### Where to place the code
Inside `K_LOOP : DO k=1,N(ng)` (starting around line 498), after the Coriolis and
pressure-gradient terms, near the end of the loop body. Use a clearly labelled block.

#### Required module access
`GRID(ng)%str_a_omn` (the precomputed `str_a/(pm*pn)` field — see §9.9a) lives in
`mod_grid`, which `rhs3d.F` already `USE`s unconditionally. **Pass `str_a_omn` as an
explicit subroutine argument** (consistent with ROMS style), analogous to how `fomn` is
passed. Add `GRID(ng)%str_a_omn` to the `CALL rhs3d_tile(...)` in the outer `rhs3d`
wrapper, and add the corresponding dummy argument declaration in `rhs3d_tile`. Since
`str_a` itself is no longer referenced in this file, `rhs3d.F` does not need
`USE mod_mixing` for `STRUCTURE_MIXING` (it is still needed when `WEC` is defined, for
`rustr3d`/`rvstr3d`).

#### Local variable declarations (add to `rhs3d_tile`)
```fortran
# ifdef STRUCTURE_MIXING
      real(r8) :: str_a_omn_u, str_a_omn_v, v_at_u, u_at_v, spd
# endif
```

#### Code block inside K_LOOP
```fortran
# ifdef STRUCTURE_MIXING
!
!-----------------------------------------------------------------------
!  Structure-induced drag body force (Rennau, Schimmels & Burchard 2012)
!  G_d = -0.5 * str_cd * str_a * u * |V|
!  str_a_omn = str_a/(pm*pn) is precomputed in metrics.F; only averaging
!  to u/v-points is needed here (no division in the time-step loop).
!-----------------------------------------------------------------------
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              str_a_omn_u=0.5_r8*(str_a_omn(i-1,j,k)+str_a_omn(i,j,k))
              IF (str_a_omn_u.gt.0.0_r8) THEN
                v_at_u=0.25_r8*(v(i-1,j  ,k,nstp)+v(i,j  ,k,nstp)+      &
     &                          v(i-1,j+1,k,nstp)+v(i,j+1,k,nstp))
                spd=SQRT(u(i,j,k,nstp)**2+v_at_u**2)
                cff=0.5_r8*str_cd(ng)*str_a_omn_u*spd*                   &
     &              0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
                ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff*u(i,j,k,nstp)
              END IF
            END DO
          END DO
!
          DO j=JstrV,Jend
            DO i=Istr,Iend
              str_a_omn_v=0.5_r8*(str_a_omn(i,j-1,k)+str_a_omn(i,j,k))
              IF (str_a_omn_v.gt.0.0_r8) THEN
                u_at_v=0.25_r8*(u(i  ,j-1,k,nstp)+u(i+1,j-1,k,nstp)+   &
     &                          u(i  ,j  ,k,nstp)+u(i+1,j  ,k,nstp))
                spd=SQRT(u_at_v**2+v(i,j,k,nstp)**2)
                cff=0.5_r8*str_cd(ng)*str_a_omn_v*spd*                   &
     &              0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
                rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff*v(i,j,k,nstp)
              END IF
            END DO
          END DO
# endif
```

#### Units check
`ru(i,j,k,nrhs)` accumulates volume-integrated momentum tendencies. The drag acceleration
$G_d^u$ [m/s²] is applied over layer thickness `0.5*(Hz(i-1,j,k)+Hz(i,j,k))` [m], giving a
contribution with the same units as how surface/bottom body forces are applied. ✓

#### Barotropic propagation — no extra action needed
At the bottom of `rhs3d.F`, `ru(i,j,k,nrhs)` is summed over all k into `rufrc(i,j)`,
which is the depth-integrated baroclinic forcing passed to the barotropic 2D mode.
The structural drag is automatically included. No changes to `step2d.F` are needed.

---

### 9.9a `ROMS/Modules/mod_grid.F` and `ROMS/Utility/metrics.F` (efficiency optimization)

`mod_grid.F` declares a new 3D pointer, guarded by `STRUCTURE_MIXING`:
```fortran
# ifdef STRUCTURE_MIXING
          real(r8), pointer :: str_a_omn(:,:,:)
# endif
```
allocated as `GRID(ng)%str_a_omn(LBi:UBi,LBj:UBj,N(ng))`, initialised to zero, and
destroyed alongside `GRID(ng)%Hz`, following the same pattern as the existing
`fomn = f*omn` compound field.

`metrics.F` computes it once, immediately after `omn(i,j) = 1/(pm(i,j)*pn(i,j))` is
computed (the same place `fomn` is computed):
```fortran
# ifdef STRUCTURE_MIXING
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            str_a_omn(i,j,k)=str_a(i,j,k)*omn(i,j)
          END DO
        END DO
      END DO
# endif
```
This requires passing `MIXING(ng)%str_a` (input) and `GRID(ng)%str_a_omn` (output)
through the `metrics`/`metrics_tile` argument lists, and extending the existing
`# if defined DIFF_3DCOEF || defined VISC_3DCOEF` guard on `USE mod_mixing` to also
include `STRUCTURE_MIXING`. `str_a_omn` is exchanged the same way as `fomn`/`omn`
(`exchange_r3d_tile`/`mp_exchange3d` for periodic/distributed boundaries), since
`metrics.F` is a one-time initialization routine (called from `set_grid.F`), this
computation and exchange happens only once per run, not every time step.

`get_grid.F` is unchanged: it still reads `str_a` into `MIXING(ng)%str_a`, which remains
the sole source field; `str_a_omn` is purely a derived, precomputed convenience field.

---

### 9.10 `ROMS/Nonlinear/gls_corstep.F`

Adds the structure drag production $P_d$ to both the TKE and GLS equations.

#### Where to place the code
Inside the production/dissipation loop
(`DO j=Jstr,Jend` → `DO i=Istr,Iend` → `DO k=1,N(ng)-1`), immediately after the
shear+buoyancy production block (around lines 772–779). Note: k here indexes W-points
(interfaces); at interface k, the adjacent rho-layers are k (below) and k+1 (above).

#### Passing `str_a` into the subroutine
Add `MIXING(ng)%str_a` as an explicit argument to both the outer `gls_corstep` wrapper
call and the `gls_corstep_tile` subroutine signature, following the same pattern as `Akv`,
`Lscale` etc. `u` and `v` are already passed.

Dummy argument declaration in `gls_corstep_tile`:
```fortran
# ifdef STRUCTURE_MIXING
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: str_a(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: str_a(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
```

#### `str_cd` and `gls_c4` access
Both are in `mod_scalars`, which `gls_corstep_tile` already `USE`s. No additional `USE`
statements needed.

#### Ghost cells
The code accesses `u(i+1,j,k,nstp)` and `v(i,j+1,k,nstp)`. This already occurs in the
existing shear computation (around line 355), so halo exchange is already done before
`gls_corstep` is called. No additional `mp_exchange` calls are needed.

#### Local variable declarations
```fortran
# ifdef STRUCTURE_MIXING
      real(r8) :: str_a_wpt, spd2_k, spd2_k1, spd2_wpt, Pd_struct
# endif
```

#### Code block (immediately after the existing shear/buoyancy production lines ~772–779)
```fortran
# ifdef STRUCTURE_MIXING
!
!  Structure-drag TKE and GLS production (Rennau et al. 2012).
!  Pd = 0.5 * str_cd * str_a * (u^2+v^2)^(3/2)   [m2/s3]
!  str_a is averaged from rho-levels k and k+1 to W-point k.
!  Speed is averaged from rho-levels k and k+1 to W-point k.
!
            str_a_wpt=0.5_r8*(str_a(i,j,k)+str_a(i,j,k+1))
            IF (str_a_wpt.gt.0.0_r8) THEN
!             Speed squared at rho-level k (below interface)
              spd2_k =0.25_r8*(u(i  ,j,k  ,nstp)+u(i+1,j,k  ,nstp))**2 &
     &                +0.25_r8*(v(i,j  ,k  ,nstp)+v(i,j+1,k  ,nstp))**2
!             Speed squared at rho-level k+1 (above interface)
              spd2_k1=0.25_r8*(u(i  ,j,k+1,nstp)+u(i+1,j,k+1,nstp))**2 &
     &                +0.25_r8*(v(i,j  ,k+1,nstp)+v(i,j+1,k+1,nstp))**2
!             Average to W-point
              spd2_wpt=0.5_r8*(spd2_k+spd2_k1)
              Pd_struct=0.5_r8*str_cd(ng)*str_a_wpt*                     &
     &                  spd2_wpt*SQRT(MAX(spd2_wpt,0.0_r8))
!
!             Add to TKE equation (k equation)
              tke(i,j,k,nnew)=tke(i,j,k,nnew)+dt(ng)*cff*Pd_struct
!
!             Add to GLS equation (psi equation), coefficient gls_c4
              gls(i,j,k,nnew)=gls(i,j,k,nnew)+                          &
     &                        dt(ng)*cff*gls_c4(ng)*Pd_struct*           &
     &                        gls(i,j,k,nstp)/                           &
     &                        MAX(tke(i,j,k,nstp),gls_Kmin(ng))
            END IF
# endif
```

> `cff = 0.5*(Hz(i,j,k)+Hz(i,j,k+1))` — already computed immediately above this block
> for the shear/buoyancy production.

---

### 9.11 `ROMS/Nonlinear/step3d_uv.F` (optional)

This change is only needed if the explicit drag in §9.9 causes numerical instability
(see §12 for when this occurs).

The vertical viscosity tridiagonal system in `step3d_uv.F` solves for `u(nnew)`.
The diagonal coefficient at layer k is:
```fortran
BC(i,k) = Hzk(i,k) - FC(i,k) - FC(i,k-1)
```

To add **fully implicit** structural drag, augment the diagonal before the forward
substitution. This replaces — not supplements — the explicit form in `rhs3d.F`. Choose
only one approach.

```fortran
# ifdef STRUCTURE_MIXING
!  Implicit structural drag: increase diagonal by drag rate * Hz * dt.
            str_a_u=0.5_r8*(MIXING(ng)%str_a(i-1,j,k)+                  &
     &                      MIXING(ng)%str_a(i  ,j,k))
            IF (str_a_u.gt.0.0_r8) THEN
              v_at_u=0.25_r8*(v(i-1,j,k,nstp)+v(i,j,k,nstp)+            &
     &                        v(i-1,j+1,k,nstp)+v(i,j+1,k,nstp))
              spd=SQRT(u(i,j,k,nstp)**2+v_at_u**2)
              BC(i,k)=BC(i,k)+0.5_r8*str_cd(ng)*str_a_u*spd*Hzk(i,k)
            END IF
# endif
```

This requires `USE mod_mixing` at the top of `step3d_uv.F`. It is unconditionally stable
but slightly more complex to implement.

---

## 10. User Input File Changes

Add the following to the ROMS `.in` file, in the GLS parameter block next to `GLS_C1` etc.:

```
! STR_CD       Structure cylinder drag coefficient (nondimensional).
! GLS_C4       Structure-drag GLS coefficient (c_psi4, nondimensional).
!              GLS_C4 = 1.0 is recommended for k-epsilon (Rennau et al. 2012).
!              The coefficient appears in the psi equation:
!              d(psi)/dt = ... + (psi/k) * (gls_c1*P + gls_c4*Pd + ...) + ...

    STR_CD == 1.0d0
    GLS_C4 == 1.0d0
```

---

## 11. NetCDF Grid File Changes

The field `str_a` must be added to the ROMS grid NetCDF file (`GRDNAME`).

### Variable specification

```
Variable name:  str_a
Dimensions:     (s_rho, eta_rho, xi_rho)   ! confirm against get_grid.F implementation
Units:          m-1
long_name:      "structure frontal area density"
valid_min:      0.0
_FillValue:     0.0
```

- `s_rho` index 0 (Python) = k=1 (Fortran) = **bottom** rho-level.
- `s_rho` index N-1 (Python) = k=N (Fortran) = **top** rho-level.
- Cells without structures: `str_a = 0`.

> **Dimension order:** Confirm the exact expected dimension order by looking at how
> `nf_fread3d` is called for other 3D fields in `get_grid.F`. Document it explicitly once
> the ROMS reading code is written, and keep the script consistent.

The construction of `str_a` is handled by a separate Python script — see
`a_struct_script_notes.md`.

---

## 12. Important Numerical Considerations

### Stability of the momentum drag

The structural drag is $G_d^u = -\lambda u$ where $\lambda = \frac{1}{2} C_D\, a\, |V|$
[s$^{-1}$]. Three treatment options exist:

| Option | Where | Stability | Complexity |
|---|---|---|---|
| A — Explicit | `rhs3d.F`, using `u(nstp)` | Requires $\lambda\Delta t < 1$ | Simplest |
| B — Implicit | `step3d_uv.F`, added to diagonal `BC` | Unconditionally stable | Moderate |
| C — Semi-implicit (recommended for start) | `rhs3d.F`, speed $|V|$ from `nstp`, `u` stepped forward | Stable for typical drag | Simple |

The plan as written implements **Option C**. Instability is only likely for very dense
structures (large `str_a`) or coarse time steps. If it occurs, upgrade to **Option B** in
`step3d_uv.F`. If Option B is used, remove the drag from `rhs3d.F` — do not add it in
both places.

### Approximate energy consistency

The momentum drag (`rhs3d.F`) and GLS injection (`gls_corstep.F`) both compute $P_d$, but
at **different grid points** using **slightly different velocity interpolations**:
- `rhs3d.F`: drag at rho-levels k=1..N, speed averaged to u/v-points.
- `gls_corstep.F`: $P_d$ at W-points k=1..N-1, speed averaged over adjacent rho-layers.

The energy removed from the mean flow and the energy injected into TKE will be close but
not exactly equal. This is the same type of approximation that exists between shear
production of TKE and the dissipation from vertical viscosity in the standard code.
It does not need to be "fixed" but should be acknowledged in energy budget analyses.

### Vertical indexing

- In `gls_corstep.F`, the production loop runs k=1..N-1 over **W-points** (interfaces).
- At W-point k, rho-layer k is below and rho-layer k+1 is above.
- `str_a(i,j,k)` is at rho-level k, so `str_a_wpt = 0.5*(str_a(i,j,k)+str_a(i,j,k+1))`.

### Minimum value guard

`Pd_struct` is always non-negative. The `MAX(spd2_wpt, 0.0)` inside the square root
prevents floating-point issues from any near-zero negative values of `spd2_wpt`.

### `N2S2_HORAVG` — no changes needed

When `N2S2_HORAVG` is defined, ROMS smooths `shear2` and `buoy2` horizontally.
`Pd_struct` does not need analogous smoothing — the velocity field is already smooth at
grid scale and `str_a` is a user-specified static field.

---

## 13. Testing Strategy

1. **Zero-drag test:** Set `str_a = 0` everywhere. Results must be bit-identical to a run
   without `STRUCTURE_MIXING`. Verifies no spurious changes.

2. **Single-column test:** Non-zero `str_a` in one column only. Verify:
   - Current is reduced in layers with non-zero `str_a`.
   - TKE and mixing coefficients are enhanced in those layers.
   - All other columns are unaffected.

3. **Blockage-ratio sanity check:** Compute `str_a * Hz` at the structure location.
   This is the dimensionless frontal blockage (fraction of column blocked). For the
   drag parametrization to be valid, this should typically be well below 0.5.

4. **Energy budget test:** Confirm that the domain-integrated drag power extracted by the
   momentum equation (sum of $P_d \cdot Hz \cdot dx \cdot dy$) approximately equals the
   domain-integrated TKE injection (within the approximation noted in §12).

5. **Idealized channel test:** Run to steady state with a structure zone at mid-water.
   Verify velocity reduction in the structure zone and enhanced mixing above/below.

6. **Qualitative comparison with Rennau et al. (2012):** Their Figure 4 shows reduced
   inflow in an idealized Baltic inlet. Reproduce the qualitative behaviour.

---

## 14. Things to Keep in Mind

1. **All GLS coefficient values are closure-dependent.**
   `gls_c4 = 1.0` is for k-ε. For k-ω, a different value may be appropriate. Exposing
   `gls_c4` as a free parameter is exactly what is needed to explore this.

2. **`MY25_MIXING` is out of scope.**
   The Mellor-Yamada 2.5 closure has its own corrector step in
   `ROMS/Nonlinear/my25_corstep.F`. An analogous injection would be needed there if the
   user ever switches closures. The `checkdefs.F` guard prevents accidental use of
   `STRUCTURE_MIXING` with `MY25_MIXING`.

3. **No Tangent/Adjoint/Representer coverage.**
   This plan covers the forward nonlinear model only. Tangent-linear
   (`Tangent/tl_rhs3d.F`) and adjoint files would need updating if 4D-Var is ever used.

4. **Land masking.**
   Zero velocity at land cells means the drag is automatically zero there even if `str_a`
   is non-zero. For cleanliness, the external script should still enforce `str_a = 0` at
   land cells using `mask_rho`.

5. **Multiple nested grids.**
   `str_cd`, `gls_c4` are arrays dimensioned `(Ngrids)` — one value per grid. `str_a`
   is per-grid via `MIXING(ng)%str_a`. If using nesting, each grid file must contain its
   own `str_a` variable.

6. **Both PIO and nf90 I/O paths must be updated.**
   ROMS supports two parallel I/O backends. Changes to `get_grid.F`, `def_info.F`, and
   `wrt_info.F` must be applied in **both** the `_nf90` and `_pio` subroutines.

7. **The diagnostic output system.**
   If `DIAGNOSTICS_UV` is defined, adding a `M3fstr` diagnostic index (analogous to
   `M3fveg` in `mod_scalars.F`) allows the structure drag force to be output in history
   files. This requires additional changes to `def_his.F` and `wrt_his.F` but is optional
   for the initial implementation.

8. **Physical meaning of the area A in `a = X*d/A`.**
   `A` is the horizontal plan-view area of the ROMS grid cell: `A = 1/(pm*pn)`. See the
   external script notes for how to compute `str_a` correctly.

9. **Code style.**
   Follow ROMS formatting: 72-character lines with `&` continuation, `! ` comments,
   CPP guards with one space (`# ifdef`, `# endif`), lowercase variable names.

10. **Version control.**
    Keep the implementation on a dedicated feature branch. Commit logically by file
    (one commit per modified file) to simplify review and debugging.

11. **Floating offshore wind turbines — surface proximity.**
    When structures are near the surface (e.g. semi-submersible floaters), the structure
    drag production `Pd_struct` is injected into TKE/GLS at W-points near k=N. This is
    physically correct and the code handles it without any changes. However, note that the
    GLS surface boundary conditions set `tke(i,j,N,nnew)` and `gls(i,j,N,nnew)` as
    Dirichlet conditions based on wind stress. The production loop runs only over interior
    W-points k=1..N-1, so the surface BC itself is not overridden. Structure production
    in the top rho-layer (k=N) contributes via W-point k=N-1, which is correct.
    The net effect is that structure-induced mixing adds to (not replaces) wind-driven
    near-surface turbulence.

12. **Rennau et al. (2012) is for bottom-mounted structures.**
    The reference paper deals with bridge piers in the Fehmarn Belt. For floating wind
    turbine applications the physics is the same (cylinder drag + GLS injection), but the
    depth range is different. The coefficient `gls_c4` may need recalibration for floating
    structures; using `gls_c4 = 1.0` as a starting point is reasonable.
