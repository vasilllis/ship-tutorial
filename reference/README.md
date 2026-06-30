# LLP momentum-smearing study

Two equivalent implementations — one C++ ROOT macro, one PyROOT script — of a
small [`RDataFrame`](https://root.cern/doc/master/classROOT_1_1RDataFrame.html)
analysis that studies how detector momentum resolution propagates into the
reconstructed observables of a long-lived particle (LLP) decaying to two
daughters. The example sample is a Heavy Neutral Lepton (HNL) decay.

| File | Language | Notes |
|------|----------|-------|
| `llp_simple_analysis.C`  | C++ (ROOT macro) | Run through the ROOT interpreter. |
| `llp_simple_analysis.py` | Python (PyROOT)  | Declares the same C++ helpers via `gInterpreter.Declare`, then drives the analysis from Python. Identical numerical results and (near-)identical performance. |

Both produce the same output: `llp_smearing.png` plus a short console summary.

---

## What it computes

For each event (after requiring exactly two daughters) the daughter transverse
momenta `px`, `py` are smeared, and five observables are reconstructed from the
smeared kinematics:

1. **Invariant mass** `m_inv` of the LLP, from the summed (smeared) daughter
   4-momenta.
2. **Transverse momentum** `pT` of the reconstructed LLP.
3. **Transverse vertex residual** `Δr = r_reco − r_true`, with `r = √(x²+y²)`.
4. **Longitudinal vertex residual** `Δz = z_reco − z_true`.
5. **Vertex impact parameter** `d³ᴰ_LLP`: perpendicular distance from the
   primary vertex (origin) to the reconstructed LLP line.

Each observable is computed at four smearing levels — **σ = 0, 1%, 5%, 10%** —
and the four versions are overlaid (non-stacked) in a `THStack`, coloured along
the `kCool` palette.

---

## How the reconstruction works

**Vertex (DOCA method).** The decay vertex is reconstructed as the mid-point of
the segment of closest approach between the two daughter tracks. Each track is
modelled as a straight line through a fixed **perigee anchor** (`A1`, `A2` — the
foot of the perpendicular from the beamline origin onto each *true* daughter
line) with direction taken from the *smeared* momentum. `docaVertex()` solves
the two-line closest-approach problem in closed form and returns
`{doca, mid_x, mid_y, mid_z}`.

Anchoring at the perigee (rather than at the common true vertex) is what makes
the DOCA respond to the smearing: at **σ = 0** both lines pass through the true
vertex, so the DOCA is zero and the reconstructed vertex equals the truth —
hence all residuals collapse to zero. This is a built-in **closure check**: if
`Δr`, `Δz`, and the vertex IP are not ~0 at σ=0, something is inconsistent.

**Energies.** Daughter energies are recomputed from the smeared 3-momenta and
the *true* masses (`E = √(p² + m²)`); the mass is never smeared.

**Smearing model.** A multiplicative Gaussian factor `(1 + σ·N(0,1))` is applied
independently to `px` and `py`; `pz` is left untouched. The RNG is seeded from
the event entry number (plus a per-level/per-axis "salt"), so results are
**deterministic and reproducible even with multithreading enabled**.

> ⚠️ This smearing is a deliberate pedagogical simplification, **not** a
> realistic detector model. A real detector smears the momentum magnitude with a
> pT-dependent resolution, smears track angles separately, and has finite
> position resolution. Here only transverse momentum is perturbed, and the track
> positions (perigee anchors) are taken as perfectly known.

---

## Requirements

- **ROOT 6.x** with `RDataFrame` (any reasonably recent 6.x release).
- For the Python version: a **PyROOT**-enabled ROOT build (`import ROOT` works).
- No external Python packages beyond ROOT.

Multithreading is enabled via `ROOT::EnableImplicitMT()` and is safe here
because the smearing is entry-seeded.

---

## Input data

Both scripts read the tree **`Events`** from a hard-coded file:

```
./HNL_1.000e+00_1.000e+00_0.000e+00_1.000e+00_0.000e+00_data.root
```

Branches used: `d_px`, `d_py`, `d_pz`, `d_m` (per-event `vector<float>` of the
two daughters) and the true decay vertex `vtx_x`, `vtx_y`, `vtx_z`.

To run on a different file, edit the `RDataFrame(...)` line near the top of
`llp_simple_analysis()` in either file.

---

## Running

**C++ (interpreted):**

```bash
root -l -b -q llp_simple_analysis.C
```

**Python:**

```bash
python llp_simple_analysis.py
# or:  root -l -b -q llp_simple_analysis.py   (via PyROOT)
```

Each run writes `llp_smearing.png` to the current directory and prints a summary
table.

---

## Output

- **`llp_smearing.png`** — a 3×2 canvas with five panels (invariant mass, pT,
  transverse residual `Δr`, longitudinal residual `Δz`, vertex IP). The residual
  panels use a linear y-axis and symmetric ranges centred on zero; the vertex-IP
  panel uses log-y. The sixth pad is intentionally left empty.
- **Console summary** — per smearing level: `⟨m_inv⟩`, and the transverse and
  longitudinal vertex resolutions `σ_r`, `σ_z` (the `StdDev` of `Δr`, `Δz`).

Expected behaviour: `⟨m_inv⟩` should barely move with σ (the smearing is
unbiased on average), while `σ_r` and `σ_z` grow with σ — roughly linearly at
small σ — and are ~0 at σ=0.

---

## Configuration knobs

All near the top of `llp_simple_analysis()`:

| Variable | Meaning |
|----------|---------|
| `sig`    | The list of smearing levels to sweep (e.g. `{0.0, 0.01, 0.05, 0.10}`). |
| `lab`    | Legend labels, one per entry in `sig`. |
| `NB`     | Number of histogram bins (shared by all observables). |
| filename | Input ROOT file (the `RDataFrame` constructor argument). |
| `PV`     | Primary-vertex assumption (origin). Change if your sample has a beamspot offset. |

Histogram axis ranges are determined adaptively from the widest (10%) smearing
level so the four overlays line up.

---

## Implementation notes

- The analysis runs in **two passes**: pass 1 (`Min`/`Max` on the 10% level)
  fixes the common binning; pass 2 books the per-level histograms. Both are lazy
  `RDataFrame` actions triggered on first dereference.
- The C++ and Python `docaVertex`, `perigee`, `smear*`, etc. are byte-for-byte
  the same logic — the Python file embeds them in a single
  `gInterpreter.Declare(...)` string, which is the idiomatic way to keep the
  per-event inner loop compiled while driving everything else from Python.
- `Δr` is the difference of cylindrical *radii*, so a vertex displaced purely
  tangentially (same radius, different position) contributes little to `Δr`. If
  you want the full in-plane miss instead, change the `dr_res` definition to
  `√(Δx² + Δy²)`, or split into signed `Δx`, `Δy`.
