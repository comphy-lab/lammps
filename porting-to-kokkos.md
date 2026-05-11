# Porting Pair Styles to KOKKOS

This document provides a structured plan and reference for porting some
pair styles to the KOKKOS package.  It covers the required code changes,
documentation updates, and groups the work into Copilot-session-sized
batches ordered from simplest to most complex.

---

## Background and Motivation

The KOKKOS package provides GPU-accelerated (CUDA/HIP) and many-core CPU
(OpenMP/threads) variants of LAMMPS styles using Kokkos abstractions.
Porting a pair style to KOKKOS means:

1. Creating `src/KOKKOS/pair_<name>_kokkos.h` and
   `src/KOKKOS/pair_<name>_kokkos.cpp`.
2. Registering the new files in `src/KOKKOS/Install.sh`.
3. Making small changes to the **base class** header (adding `virtual` to
   `allocate()` and adding `if (copymode) return;` to the destructor of the
   base class if not already present).
4. Updating the RST documentation in `doc/src/pair_<name>.rst`.
5. Updating `doc/src/Commands_pair.rst`.

---

## Key Rules for KOKKOS Porting

### 1.  `if (copymode) return;` in the base-class destructor

When a Kokkos pair style is copied for use on-device (via `copymode = 1`),
the copy's destructor must not free memory that belongs to the original object.
Every **base class** that a `_kokkos` style will inherit from must have:

```cpp
PairFoo::~PairFoo()
{
  if (copymode) return;   // MUST be the first line
  // ... rest of destructor ...
}
```

If the base-class destructor already has this line, no change is needed.
Check with:

```bash
grep -n "copymode" src/<PKG>/pair_<name>.cpp
```

### 2.  `virtual void allocate()` in the base-class header

The KOKKOS style overrides `allocate()` to replace plain arrays with dual
views.  The override only works correctly in C++ when the base-class
declaration is `virtual`.  Change:

```cpp
void allocate();          // old
virtual void allocate();  // required
```

Check with:

```bash
grep -n "allocate" src/<PKG>/pair_<name>.h
```

### 3.  Per-atom view members: `typename AT::` not `DAT::`

Every member that stores a per-atom view obtained from
`atomKK->k_<field>.view<DeviceType>()` **must** be declared with the
template-parameterized alias `typename AT::t_<type>`, **never** with the
hardcoded alias `DAT::t_<type>`.

```cpp
// WRONG — kk/host instantiation fails with Kokkos ViewMapping assertion:
DAT::t_tagint_1d_randomread molecule;

// CORRECT — adapts to DeviceType (LMPDeviceType or LMPHostType):
typename AT::t_tagint_1d_randomread molecule;
```

`AT` is `ArrayTypes<DeviceType>`.  `DAT` is `ArrayTypes<LMPDeviceType>`.
In GPU builds, the pair style is instantiated for *both* `LMPDeviceType` and
`LMPHostType`.  Views with `DAT::` types live in GPU memory; assigning a
host-side `atomKK` view to them triggers a static assertion in Kokkos.

**Safe to use `DAT::`:** dual views that the pair style allocates itself and
never assigns from an `atomKK` view:
`k_eatom`, `k_vatom`, `k_cutsq`, `k_cut_ljsq`, `k_cut_coulsq`.

**Must use `typename AT::`:** all per-atom views set via
`atomKK->k_<field>.view<DeviceType>()`:
`x`, `c_x`, `f`, `type`, `q`, `molecule`, `radius`, etc.

Quick audit command:
```bash
grep "DAT::t_" src/KOKKOS/pair_<name>_kokkos.h
```
Every hit should be a dual-view (type contains `tdual` or `ttransform`), not
a plain per-atom view (type contains `t_`).

### 4.  `src/KOKKOS/Install.sh` entries

Add two lines per new file:

```sh
action pair_<name>_kokkos.cpp pair_<name>.cpp   # if base is in another package
action pair_<name>_kokkos.h   pair_<name>.h     # idem
```

If the base class is in the **core** (`src/`), omit the second argument:

```sh
action pair_<name>_kokkos.cpp
action pair_<name>_kokkos.h
```

### 5.  Documentation changes

**In `doc/src/pair_<name>.rst`:**

* Add `.. index:: pair_style <name>/kk` near the top.
* Add `*<name>/kk*` to the `Accelerator Variants:` line.

**In `doc/src/Commands_pair.rst`:**

* Add `k` to the accelerator letter string in the entry for this style
  (letters: `g`=GPU, `i`=INTEL, `k`=KOKKOS, `o`=OPENMP, `t`=OPT).

---

## Standard KOKKOS Template (pair_kokkos.h framework)

Most simple pairwise styles (no three-body, no many-body embedding) fit the
`pair_kokkos.h` template.  The KOKKOS class must implement:

| Function | COUL_FLAG=0 | COUL_FLAG=1 |
|---|---|---|
| `compute_fpair` | required | required (VdW part) |
| `compute_evdwl` | required | required |
| `compute_fcoul` | return 0 | required |
| `compute_ecoul` | return 0 | required |

The template automatically handles half/half-thread/full neighbor list
dispatch, reduction of forces, and energy/virial accumulation.

See `src/KOKKOS/pair_born_kokkos.{h,cpp}` and
`src/KOKKOS/pair_lj_cut_kokkos.{h,cpp}` as canonical examples.

### Typical header skeleton (COUL_FLAG=0)

```cpp
#ifdef PAIR_CLASS
PairStyle(<name>/kk,     Pair<Name>Kokkos<LMPDeviceType>);
PairStyle(<name>/kk/device,Pair<Name>Kokkos<LMPDeviceType>);
PairStyle(<name>/kk/host,  Pair<Name>Kokkos<LMPHostType>);
#else
#ifndef LMP_PAIR_<NAME>_KOKKOS_H
#define LMP_PAIR_<NAME>_KOKKOS_H
#include "pair_kokkos.h"
#include "pair_<name>.h"
#include "neigh_list_kokkos.h"
namespace LAMMPS_NS {
template<class DeviceType>
class Pair<Name>Kokkos : public Pair<Name> {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  Pair<Name>Kokkos(class LAMMPS *);
  ~Pair<Name>Kokkos() override;
  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_<name> {
    KOKKOS_INLINE_FUNCTION params_<name>() { /* zero all */ }
    KOKKOS_INLINE_FUNCTION params_<name>(int) { /* zero all */ }
    KK_FLOAT cutsq, /* ... pair-specific params ... */;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(...) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(...) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(...) const { return 0; }

  Kokkos::DualView<params_<name>**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_<name>**,...>::t_dev_const_um params;
  params_<name> m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  Kokkos::DualView<KK_FLOAT**,DeviceType> k_cutsq;
  typename Kokkos::DualView<KK_FLOAT**,...>::t_dev d_cutsq;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;
  Kokkos::DualView<KK_FLOAT*,DeviceType> k_eatom, k_vatom;
  int neighflag, newton_pair;
  int nlocal, nall, eflag, vflag;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,FULL,true>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,FULL,false>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALF,true>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALF,false>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALFTHREAD,true>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<Pair<Name>Kokkos<DeviceType>,
    FULL,CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<Pair<Name>Kokkos<DeviceType>,
    HALF,CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<Pair<Name>Kokkos<DeviceType>,
    HALFTHREAD,CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<Pair<Name>Kokkos<DeviceType>,
    CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<Pair<Name>Kokkos<DeviceType>>(
    Pair<Name>Kokkos<DeviceType>*);
};
}
#endif
#endif
```

For `COUL_FLAG=1`, add `compute_fcoul` and `compute_ecoul`, add
`special_coul`, `qqrd2e`, `q` array, and change `CoulLongTable<0>` to
`CoulLongTable<1>` in the friend declarations.

### Typical `.cpp` skeleton

```cpp
template<class DeviceType>
Pair<Name>Kokkos<DeviceType>::Pair<Name>Kokkos(LAMMPS *lmp) : Pair<Name>(lmp)
{
  respa_enable = 0;
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

template<class DeviceType>
Pair<Name>Kokkos<DeviceType>::~Pair<Name>Kokkos()
{
  if (copymode) return;
  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->destroy_kokkos(k_cutsq, cutsq);
  }
}

template<class DeviceType>
void Pair<Name>Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in; vflag = vflag_in;
  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  ev_init(eflag,vflag,0);
  if (eflag_atom) { /* reallocate k_eatom */ }
  if (vflag_atom) { /* reallocate k_vatom */ }
  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);
  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal; nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0..3] = ...;
  NeighListKokkos<DeviceType>* k_list = ...;
  ...
  EV_FLOAT ev = pair_compute<Pair<Name>Kokkos<DeviceType>,
    CoulLongTable<COUL_FLAG>>(this, k_list);
  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) { virial[0..5] += ev.v[0..5]; }
  if (eflag_atom) { k_eatom.template modify<DeviceType>(); k_eatom.sync_host(); }
  if (vflag_atom) { k_vatom.template modify<DeviceType>(); k_vatom.sync_host(); }
  if (vflag_fdotr) pair_virial_fdotr_compute(this);
  copymode = 0;
}

template<class DeviceType>
void Pair<Name>Kokkos<DeviceType>::allocate()
{
  Pair<Name>::allocate();
  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_<name>**,...>("Pair<Name>::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

template<class DeviceType>
void Pair<Name>Kokkos<DeviceType>::init_style()
{
  Pair<Name>::init_style();
  // Request full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
}

template<class DeviceType>
double Pair<Name>Kokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = Pair<Name>::init_one(i,j);
  k_params.view_host()(i,j).<field> = static_cast<KK_FLOAT>(<field>[i][j]);
  /* ... fill all params fields ... */
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i < MAX_TYPES_STACKPARAMS+1 && j < MAX_TYPES_STACKPARAMS+1)
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
  m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();
  return cutone;
}

// Explicit instantiation at the bottom of the .cpp file
namespace LAMMPS_NS {
template class Pair<Name>Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class Pair<Name>Kokkos<LMPHostType>;
#endif
}
```

---

## Candidate Pair Styles to Port

52 pair styles are possible candidates for porting to KOKKOS.  They are
grouped below by porting complexity, from easiest to hardest.

---

## Group 1 — Simple Pairwise, No Coulomb (EXTRA-PAIR; ~13 styles)

**Status: COMPLETED**

**Complexity:** Low.  Use the `pair_kokkos.h` template with `COUL_FLAG=0`.
Implement `compute_fpair` and `compute_evdwl`; `compute_ecoul` returns 0.
Pattern to follow: `pair_born_kokkos` or `pair_beck_kokkos`.

**Session size recommendation:** 4–6 styles per session.

| Pair style | Package | Base header | Status |
|---|---|---|---|
| `mie/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_mie_cut.h` | **done** |
| `nm/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_nm_cut.h` | **done** |
| `gauss/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_gauss_cut.h` | **done** |
| `harmonic/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_harmonic_cut.h` | **done** |
| `cosine/squared` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_cosine_squared.h` | **done** |
| `lj/smooth` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_lj_smooth.h` | **done** |
| `lj/smooth/linear` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_lj_smooth_linear.h` | **done** |
| `morse/smooth/linear` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_morse_smooth_linear.h` | **done** |
| `ufm` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_ufm.h` | **done** |
| `wf/cut` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_wf_cut.h` | **done** |
| `born/gauss` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_born_gauss.h` | **done** |
| `lj/mdf` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_lj_mdf.h` | **done** |
| `lennard/mdf` | EXTRA-PAIR | `src/EXTRA-PAIR/pair_lennard_mdf.h` | **done** |

**Required base-class changes:**
- Add `virtual void allocate();` in each header (if not already `virtual`).
- Add `if (copymode) return;` as first line of each destructor.

**Documentation:** Add `.. index:: pair_style <name>/kk` and add `*<name>/kk*`
to `Accelerator Variants`.  Add `k` to each entry in `Commands_pair.rst`.

---

## Group 2 — Simple Pairwise, No Coulomb (continued; ~5 styles)

**Status: COMPLETED**

**Complexity:** Low.  Same as Group 1 but with more parameters per pair-type.

| Pair style | Package | Notes | Status |
|---|---|---|---|
| `momb` | EXTRA-PAIR | London + Buckingham, extra exponent params | **done** |
| `buck/mdf` | EXTRA-PAIR | Buckingham with switching, see `pair_buck_kokkos` | **done** |
| `pedone` | EXTRA-PAIR | Morse-like with r^-6 correction, no Coulomb | **done** |
| `lj/pirani` | EXTRA-PAIR | Generalized LJ with Pirani potential, per-type m,n | **done** |
| `ylz` | ASPHERE | Orientation-dependent potential, requires custom kernel | **done** |

**Note:** `ylz` is in the ASPHERE package but computes only short-range VdW
with fixed charge; it should use `COUL_FLAG=0`.

---

## Group 3 — Pairwise with Short-Range Coulomb (Wolf/DSF/cut; ~8 styles)

**Status: COMPLETED**

**Complexity:** Moderate.  Use the `pair_kokkos.h` template with `COUL_FLAG=1`.
Add `compute_fcoul` and `compute_ecoul`.  The Coulomb part uses no long-range
tables (Wolf damping, DSF cutoff, or plain cut).
Pattern to follow: `pair_lj_cut_coul_dsf_kokkos`.

| Pair style | Package | Coulomb type | Status |
|---|---|---|---|
| `born/coul/wolf` | EXTRA-PAIR | Wolf damped | **done** |
| `lj/cut/coul/wolf` | EXTRA-PAIR | Wolf damped | **done** |
| `nm/cut/coul/cut` | EXTRA-PAIR | cut | **done** |
| `coul/diel` | EXTRA-PAIR | dielectric screening | **done** |
| `coul/shield` | INTERLAYER | exponential screening | **done** |
| `buck6d/coul/gauss/dsf` | MOFFF | Gaussian-damped DSF | **done** |
| `coul/cut/global` | EXTRA-PAIR | tiny wrapper of `coul/cut` | **done** |
| `lj/cut/sphere` | EXTRA-PAIR | radius from per-atom property; no Coulomb (COUL_FLAG=0) | **done** |

**Notes:**
- `coul/diel`: stores scalar dielectric parameters (`a_eps`, `b_eps`, `eps_s`) as `KK_FLOAT`
  members and copies them before each kernel launch.
- `coul/shield`: includes per-atom `molecule` view to skip same-layer pairs; the taper
  function (Tap, dTap) is inlined using Horner's method with hardcoded coefficients.
  **Important:** the `molecule` member must be declared `typename AT::t_tagint_1d_randomread`
  (not `DAT::t_tagint_1d_randomread`), otherwise the `kk/host` instantiation fails to
  compile in GPU builds.  See *Lessons Learned* for the general rule.
- `coul/cut/global`: inherits `PairCoulCutKokkos<DeviceType>` directly; only overrides
  `coeff()` (enforces 2 args) and `extract()`.
- `lj/cut/sphere`: COUL_FLAG=0; the per-atom `radius` view (declared `typename AT::
  t_kkfloat_1d_randomread`) is synced; cutsq stores the maximum possible cutoff per
  type pair; the actual per-atom cutoff and sigma are computed inside
  `compute_fpair`/`compute_evdwl`.

---

## Group 4 — Pairwise with Long-Range Coulomb (Ewald/PPPM; ~10 styles)

**Complexity:** Moderate.  `COUL_FLAG=1`, use `init_tables()` for Coulomb.
Pattern to follow: `pair_born_kokkos` extended with `pair_lj_cut_coul_long_kokkos`.

| Pair style | Package | Notes |
|---|---|---|
| `nm/cut/coul/long` | EXTRA-PAIR | NM with long-range Ewald |
| `lj/switch3/coulgauss/long` | YAFF | YAFF + Gaussian charge Ewald |
| `mm3/switch3/coulgauss/long` | YAFF | YAFF MM3 variant |
| `buck6d/coul/gauss/long` | MOFFF | MOFFF Buckingham + Ewald |

---

## Group 5 — TIP4P Water Models (~4 styles)

**Complexity:** High.  The TIP4P geometry correction (virtual oxygen site M)
requires special handling in the neighbor loop.  Study the existing
`pair_lj_cut_tip4p_long_kokkos` carefully before attempting.

| Pair style | Package | Notes |
|---|---|---|
| `tip4p/cut` | MOLECULE | TIP4P cut |
| `tip4p/long` | KSPACE | TIP4P long Coulomb |
| `lj/cut/tip4p/cut` | MOLECULE | LJ + TIP4P cut |
| `lj/cut/tip4p/long` | KSPACE | LJ + TIP4P long |

---

## Group 6 — Granular, Colloidal, Lubrication (~2 styles)

**Complexity:** High.  These involve contact mechanics (overlap detection,
friction history) or hydrodynamic lubrication requiring per-atom state.

| Pair style | Package | Notes |
|---|---|---|
| `gran/hooke` | GRANULAR | frictionless Hertz; no history arrays needed |
| `gran/hertz/history` | GRANULAR | requires shear history Kokkos view |

---

## Group 7 — MANYBODY Three-body (SW-like, Tersoff tables; ~3 styles)

**Complexity:** High.  These use a short neighbor list in addition to the
regular neighbor list.  Follow `pair_sw_kokkos` and `pair_tersoff_kokkos`
as patterns; the short-neigh build kernel and the three-body loop kernel
must be implemented manually (the `pair_kokkos.h` template does not support
three-body).

| Pair style | Package | Notes |
|---|---|---|
| `sw/mod` | MANYBODY | SW variant; inherits from PairSW; follow `pair_sw_kokkos` |
| `nb3b/harmonic` | MANYBODY | three-body harmonic, independent from SW hierarchy |
| `tersoff/mod/c` | MANYBODY | Tersoff mod + carbon correction |

---

## Group 8 — Moderate Many-Body / Geometry-Dependent (~6 styles)

**Complexity:** High.  These require non-trivial per-atom or angle-dependent
loops that do not map onto the `pair_kokkos.h` template.

| Pair style | Package | Notes |
|---|---|---|
| `atm` | MANYBODY | three-body Axilrod-Teller-Muto; triple loop |
| `local/density` | MANYBODY | local density embedding; two-pass loop |
| `gayberne` | ASPHERE | ellipsoid-ellipsoid; per-atom quaternions on device |
| `resquared` | ASPHERE | re-squared ellipsoid; similar to gayberne |
| `line/lj` | ASPHERE | line-segment LJ; per-atom shape/orientation |
| `tri/lj` | ASPHERE | triangle LJ; per-atom shape/orientation |

---

## Group 9 — H-Bond, Dipole, Peridynamics (~6 styles)

**Complexity:** High.  Angle-dependent (H-bond) or orientation-dependent
(dipole) or continuum mechanics (peri) interactions.

| Pair style | Package | Notes |
|---|---|---|
| `hbond/dreiding/lj` | MOLECULE | angle-dependent H-bond; triple loop |
| `hbond/dreiding/morse` | MOLECULE | Morse variant of above |
| `hbond/dreiding/lj/angleoffset` | EXTRA-MOLECULE | angle-offset extension |
| `hbond/dreiding/morse/angleoffset` | EXTRA-MOLECULE | angle-offset extension |
| `lj/sf/dipole/sf` | DIPOLE | shifted-force dipole + LJ |
| `lj/expand/sphere` | EXTRA-PAIR | per-atom sphere radii required |

---

## Documentation Checklist per Style

For each pair style ported, update the following:

1. **`src/KOKKOS/pair_<name>_kokkos.h`** — new file
2. **`src/KOKKOS/pair_<name>_kokkos.cpp`** — new file
3. **`src/KOKKOS/Install.sh`** — add `action` lines
4. **Base class `.h` in `src/<PKG>/pair_<name>.h`**:
   - Make `allocate()` virtual
5. **Base class `.cpp` in `src/<PKG>/pair_<name>.cpp`**:
   - Add `if (copymode) return;` as first line of destructor
6. **`doc/src/pair_<name>.rst`**:
   - Add `.. index:: pair_style <name>/kk` near the top and so that they are ordered in groups by the base pair style and then for the same base style alphabetically according to the accelerator suffix. So `/kk` follows `/gpu` but comes before `/omp`.
   - Add `*<name>/kk*` to `Accelerator Variants:`
   - Order the list of accelerator variants alphabetically
7. **`doc/src/Commands_pair.rst`**:
   - Add letter `k` to the entry for this style
   - When adding to existing letter use (ko) and *not* (o,k), that is order alphabetically and do not use a comma.

---

## Build and Test Instructions

### Build for testing

```bash
mkdir -p build
cmake -S cmake -B build \
  -C cmake/presets/gcc.cmake \
  -C cmake/presets/most.cmake \
  -D PKG_KOKKOS=on \
  -D Kokkos_ENABLE_OPENMP=on \
  -D DOWNLOAD_POTENTIALS=off \
  -D ENABLE_TESTING=on \
  -G Ninja
cmake --build build -j 4
```

### Style checks before committing

```bash
cd src && make check-whitespace && make check-permissions
```

### Running unit tests

```bash
cd build && ctest -V -R pair
```

---

## Summary of Common Mistakes to Avoid

1. **Forgetting `if (copymode) return;`** — causes double-free crash on GPU.
2. **Forgetting `virtual void allocate()`** — silently calls the wrong allocate,
   missing KK dual-view setup, leading to incorrect results or crashes.
3. **Wrong `CoulLongTable` template parameter** — use `<0>` when
   `COUL_FLAG=0`, `<1>` when `COUL_FLAG=1`.
4. **Missing friend declarations** — `pair_compute*` helpers are templated
   free functions; the class must declare them as friends.
5. **Not using `KK_FLOAT`** — all floating-point fields in `params_*` structs
   must use `KK_FLOAT` (not `double`).
6. **Not adding `// NOLINTNEXTLINE` before each `KOKKOS_INLINE_FUNCTION`** —
   required to suppress linter warnings in LAMMPS.
7. **Wrong Install.sh action format** — when base class is in another package
   (e.g., EXTRA-PAIR), pass both filenames; when base is in `src/`, omit
   second filename.
8. **Not rebuilding with `make purge` after switching build systems.**
9. **Using `DAT::t_<type>` for per-atom views assigned from `atomKK`** —
   every per-atom view assigned via `atomKK->k_<field>.view<DeviceType>()`
   must be declared `typename AT::t_<type>` (not `DAT::t_<type>`).  `AT` is
   `ArrayTypes<DeviceType>`, so the view type adapts to the template parameter.
   `DAT` is hardcoded to `ArrayTypes<LMPDeviceType>`, which compiles fine for
   the `kk/device` instantiation but triggers a Kokkos `ViewMapping` static
   assertion when the class is instantiated for `LMPHostType` (the `kk/host`
   variant required in GPU builds).  Standard views that hold their own data
   (e.g., `DAT::ttransform_kkacc_1d k_eatom`, `DAT::ttransform_kkfloat_2d
   k_cutsq`) are correctly declared with `DAT::` because they are allocated
   locally and are not assigned from an `atomKK->k_*.view<DeviceType>()` call.

---

## Lessons Learned (updated after Group 2)

### Global scalar parameters in pair styles

Some pair styles (e.g., `momb`) have global scalar parameters set by
`settings()` rather than per-pair parameters set by `pair_coeff`.  These
scalars are stored as plain `double` members of the base class (e.g.,
`sscale`, `dscale`).  For the KOKKOS version, copy them to `KK_FLOAT`
members (e.g., `m_sscale`, `m_dscale`) in `compute()` before the kernel
launch.  They can then be captured by the `KOKKOS_INLINE_FUNCTION` methods
through `this` (which is captured by the Kokkos functor).

### MDF (Mei-Davenport-Fernando) tapering

The `buck/mdf` style has a smooth taper applied when the inter-atomic
distance is in the range `[cut_inner, cut_outer]`.  In the standard
`pair_kokkos.h` template, the `compute_fpair` function receives `rsq` and
must compute `r = sqrt(rsq)` itself.  The tapering logic (computing `d`,
`tt`, `dt`) can be inlined directly in `compute_fpair` and `compute_evdwl`
using a conditional `if (rsq > cut_inner_sq)`.  Store `cut_inner`,
`cut_inner_sq`, and the outer `cutsq` in the per-pair `params_*` struct so
the kernel can branch correctly.

### Orientation-dependent pair styles (ylz, gayberne, etc.)

The `ylz` style computes orientation-dependent forces and torques based on
per-atom quaternion data.  It does not fit the `pair_kokkos.h` template
because:
1. It needs to read quaternion data from `AtomVecEllipsoidKokkos`.
2. It accumulates torques as well as forces.
3. The inner loop is not a simple `compute_fpair` scalar.

For such styles a **custom kernel** is required, following the pattern of
`pair_lj_cut_dipole_cut_kokkos`:

- Use `TagPairXXXKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>` tags.
- Override `compute()` to launch the appropriate template instantiation with
  `Kokkos::parallel_for` or `Kokkos::parallel_reduce`.
- Accumulate torques with atomic operations via a view with
  `Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>`.
- Access the ellipsoid bonus data (quaternions) via
  `avecKK->k_bonus.view<DeviceType>()` where `avecKK` is a pointer to
  `AtomVecEllipsoidKokkos` cast from `atom->style_match("ellipsoid")`.
- Use `MathExtraKokkos::quat_to_mat()` to convert quaternion to rotation
  matrix on-device.
- Set `datamask_read` to include `TORQUE_MASK | ELLIPSOID_MASK`.
- Set `datamask_modify` to include `TORQUE_MASK`.
- Store the neighbor-list arrays explicitly (`d_neighbors`, `d_numneigh`,
  `d_ilist`) rather than relying on the `pair_compute` helper.

### lj/pirani: Kokkos::pow in device code

The ILJ/Pirani potential uses `pow(rx, n(rx))` where `n(rx)` is itself a
function of the inter-atomic distance.  `Kokkos::pow(base, exp)` is the
device-portable equivalent of `std::pow`.  Using `std::pow` inside a
`KOKKOS_INLINE_FUNCTION` causes a compiler error on GPU targets.  Always
prefer `Kokkos::` math functions in device kernels:
`Kokkos::sqrt`, `Kokkos::exp`, `Kokkos::pow`, `Kokkos::log`,
`Kokkos::sin`, `Kokkos::cos`, etc.

### Force sign convention

Be careful with force sign conventions.  In the base CPU code, forces are
often accumulated as `f[i][0] += delx * fpair` where `fpair = F/r` and
`delx = xj - xi`.  The KOKKOS kernel must use the same convention.  For the
`ylz` style the force on atom `i` due to atom `j` is `-dU/dr * rhat` where
`rhat` points from `i` to `j` (i.e., `r12 = xj - xi`); the force on `j`
is the Newton-3 partner `-force_on_i`.

---

## Lessons Learned (updated after Group 3, session May 2026)

### Per-atom views: always use `typename AT::` not `DAT::`

The most dangerous portability mistake for COUL_FLAG=1 styles (and any style
accessing extra per-atom arrays) is declaring a per-atom view with the
hardcoded `DAT::` type alias instead of the template-parameterized `AT::`:

```cpp
// WRONG — compiles for kk/device but fails at kk/host instantiation:
DAT::t_tagint_1d_randomread molecule;

// CORRECT — adapts to DeviceType, works for both kk/device and kk/host:
typename AT::t_tagint_1d_randomread molecule;
```

The rule is simple: any member that is assigned from
`atomKK->k_<field>.view<DeviceType>()` must be declared with `typename AT::`.
Views that manage their own data (dual views allocated in `allocate()`, e.g.,
`k_eatom`, `k_vatom`, `k_cutsq`) are correctly declared with `DAT::` because
they do not hold per-atom data obtained from `atomKK`.

This bug manifests as a Kokkos `ViewMapping` static assertion at compile time
when targeting GPU architectures (because GPU builds instantiate both
`Pair<X>Kokkos<LMPDeviceType>` and `Pair<X>Kokkos<LMPHostType>`).  On CPU-only
builds both `LMPDeviceType` and `LMPHostType` share the same memory space, so
the bug is silent.

The `pair_coul_shield_kokkos` style contains a `molecule` array (to skip
same-layer pairs) that was originally declared as `DAT::t_tagint_1d_randomread`
and was fixed to `typename AT::t_tagint_1d_randomread` in this session.

**Checklist:** when writing any new KK pair style, audit every per-atom
array member with `grep "DAT::" pair_<name>_kokkos.h`.  Only these are safe
as `DAT::`:

- `DAT::ttransform_kkacc_1d k_eatom`
- `DAT::ttransform_kkacc_1d_6 k_vatom`
- `DAT::ttransform_kkfloat_2d k_cutsq`
- `DAT::tdual_kkfloat_2d k_cut_ljsq` / `k_cut_coulsq`

Everything else that is set to `atomKK->k_<field>.view<DeviceType>()` must
use `typename AT::`.

### `coul/shield`: per-molecule filtering via the `molecule` view

`pair_coul_shield` skips interactions between atoms in the same layer (same
molecule id).  In the KOKKOS kernel this is implemented with a per-atom
`molecule` view:

```cpp
// in compute_fcoul / compute_ecoul:
if (molecule(i) == molecule(j)) return static_cast<KK_FLOAT>(0.0);
```

The view is synced in `compute()` before the kernel launch:

```cpp
molecule = atomKK->k_molecule.view<DeviceType>();
```

Because the `molecule` view's type must match `DeviceType`'s memory space,
it must be declared `typename AT::t_tagint_1d_randomread`, not
`DAT::t_tagint_1d_randomread`.  This is the concrete example that motivated
the general rule described above.

### `coul/shield`: inline taper function

The `coul/shield` potential uses a polynomial taper function to smooth the
force to zero at the cutoff.  In the KOKKOS kernel, all math must be
device-portable.  The taper and its derivative are computed via Horner's
method with hardcoded coefficients (exactly the same coefficients as the CPU
base class), called directly inside `compute_fcoul` and `compute_ecoul`.
No LAMMPS utility function is needed; the 8 polynomial coefficients are
embedded as literal constants in the kernel.

### `coul/cut/global`: thin wrapper inheriting from `PairCoulCutKokkos`

`pair_coul_cut_global_kokkos` does not need its own full implementation.
It is sufficient to:
1. Create a subclass `PairCoulCutGlobalKokkos<DeviceType>` that inherits
   from `PairCoulCutKokkos<DeviceType>`.
2. Override `coeff()` to enforce exactly two arguments (as in the base CPU
   class).
3. Override `extract()` to return the `epsilon` pointer.
4. Redirect the `PairStyle` macro entries to the new class.

No new `compute_fpair`/`compute_fcoul` is needed; the parent handles all
kernel logic.

### `lj/cut/sphere`: per-atom radius with a global `cutsq` ceiling

`pair_lj_cut_sphere_kokkos` reads per-atom radii and computes per-pair
cutoffs on-the-fly inside `compute_fpair`/`compute_evdwl`.  `cutsq` stores
the **maximum possible** squared cutoff per atom-type pair (set in
`init_one()`), so the neighbor list is built conservatively.  Inside the
kernel, the actual cutoff for atoms `i`–`j` is
`(sigma_ij + radius[i] + radius[j])^2`, and pairs beyond that distance
return zero force regardless of the neighbor list entry.

The `radius` per-atom view is declared `typename AT::t_kkfloat_1d_randomread
radius` (using `typename AT::`, not `DAT::`) and is synced in `compute()`:

```cpp
radius = atomKK->k_radius.view<DeviceType>();
```

The same `typename AT::` rule applies here as to `molecule`.
