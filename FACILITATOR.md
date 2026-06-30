# Facilitator key -- NOT for students

The student path is `tutorial.ipynb`; `solutions.ipynb` is the completed version.

## TODO 1 -- transverse momentum
```python
reco_pt = reco.Define("pt", "sqrt(px_tot*px_tot + py_tot*py_tot)")
h_pt = reco_pt.Histo1D(("h_pt","p_{T}; p_{T} [GeV]; events",60,0,5), "pt")
```

## TODO 2 -- mass resolution
```python
def mass_resolution(sigma):
    d = (df.Filter("d_px.size() == 2")
           .Define("px_s", f"smear(d_px, {sigma})")
           .Define("py_s", f"smear(d_py, {sigma})")
           .Define("E_s",  "sqrt(px_s*px_s + py_s*py_s + d_pz*d_pz + d_m*d_m)")
           .Define("m_s",  "sqrt(pow(Sum(E_s),2)-pow(Sum(px_s),2)-pow(Sum(py_s),2)-pow(Sum(d_pz),2))"))
    return d.StdDev("m_s").GetValue()
```

## Expected numbers
| sigma | mean mass | resolution (StdDev) |
|------:|:---------:|:-------------------:|
| 0     | 1.000 GeV | ~0.000 (closure)    |
| 1%    | 1.000 GeV | ~0.007 GeV          |
| 5%    | 1.000 GeV | ~0.034 GeV          |
| 10%   | ~1.002 GeV| ~0.067 GeV          |

Mean stays at 1 GeV (unbiased); resolution grows roughly linearly (~0.67 x sigma in GeV).

## Notes / common stumbles
- Run cells **in order** -- later cells reuse `df`, `reco`, and the declared `smear`.
- The notebook reads from EOS (`/eos/user/m/matclim/SHiPsim/SS_2026_data`). If the
  filename there differs, edit the `named = ...` line in the data-path cell. A local
  `data/` copy is the fallback if EOS isn't reachable.
- Smearing uses `gRandom` (single-threaded, seeded 123) so results are reproducible;
  the reference macro uses an entry-seeded scheme that's also reproducible *with* MT.
- At sigma=0 the smear helper still draws (and multiplies by 0) so the peak is the
  truth -- a built-in closure check.
