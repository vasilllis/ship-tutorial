#!/usr/bin/env python3
"""Build tutorial.ipynb (TODOs) and solutions.ipynb (filled) for the LLP study."""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

def md(s):   return new_markdown_cell(s.strip("\n"))
def code(s): return new_code_cell(s.strip("\n"))

intro = md(r"""
# LLP hands-on: reconstructing a Heavy Neutral Lepton

This sample is a simulated **Heavy Neutral Lepton (HNL)** — a hypothetical
long-lived particle. Each HNL flies tens of metres downstream and decays into
**two daughters** at a displaced vertex (just the kind of signature SHiP looks for).

You'll use **ROOT's RDataFrame** to:
1. reconstruct the HNL **invariant mass** from its daughters,
2. study how detector **momentum resolution** smears that measurement.

Run the cells top to bottom (Shift+Enter). The **TODO** cells are yours to fill in.
""")

setup_md = md(r"""
## 0. Setup

Runs on **SWAN** with any recent LCG stack (ROOT + Python included).
`%jsroot on` makes ROOT plots interactive inside the notebook.
""")
setup_code = code(r"""
%jsroot on
import ROOT
print("ROOT", ROOT.gROOT.GetVersion())
""")

data_md = md(r"""
## 1. The data (on EOS)

The samples live on **EOS** — CERN's shared storage, mounted at `/eos/...` on SWAN
and lxplus. We don't copy data into the git repo; we read it straight from EOS.

The tree is called **`Events`**, one row per HNL decay. The branches we'll use:

| branch | meaning |
|--------|---------|
| `LLP_m` | the *true* HNL mass (what we should reconstruct) |
| `vtx_x, vtx_y, vtx_z` | the decay vertex position [m] |
| `d_px, d_py, d_pz, d_m` | the two daughters' momenta & masses (a vector per event) |

First, locate the file on EOS:
""")
path_code = code(r"""
import os, glob

EOS_DIR = "/eos/user/m/matclim/SHiPsim/SS_2026_data"
named   = f"{EOS_DIR}/HNL_1.000e+00_1.000e+00_0.000e+00_1.000e+00_0.000e+00_data.root"

if os.path.exists(named):
    DATA = named                                   # the HNL sample on EOS
elif glob.glob(f"{EOS_DIR}/*.root"):
    DATA = sorted(glob.glob(f"{EOS_DIR}/*.root"))[0]   # any sample in that folder
else:
    DATA = "data/hnl_data.root"                    # fallback: local copy in the repo

print("reading:", DATA)
""")
data_code = code(r"""
df = ROOT.RDataFrame("Events", DATA)

print("events:    ", df.Count().GetValue())
print("true mass: ", round(df.Mean("LLP_m").GetValue(), 4), "GeV")
""")

vtx_md = md(r"""
The HNL is *long-lived*, so it decays far from the collision point. Let's see where
along the beam (`vtx_z`):
""")
vtx_code = code(r"""
c0 = ROOT.TCanvas("c0", "", 700, 500)
h_z = df.Histo1D(("h_z", "decay vertex; z [m]; decays", 100, 0, 100), "vtx_z")
h_z.Draw()
c0.Draw()
""")
vtx_note = md(r"""
The decays sit tens of metres downstream — that displaced vertex is the key
experimental handle for finding long-lived particles.
""")

mass_md = md(r"""
## 2. Reconstruct the invariant mass

The HNL's mass is the invariant mass of its two daughters: add their 4-momenta
and take $m = \sqrt{E^2 - p_x^2 - p_y^2 - p_z^2}$.

In RDataFrame you build new columns with **`Define`** (string expressions, evaluated
in C++). `Sum(...)` adds up the two-element daughter vector.
""")
mass_code = code(r"""
reco = (df
    .Filter("d_px.size() == 2", "exactly two daughters")
    .Define("d_E",    "sqrt(d_px*d_px + d_py*d_py + d_pz*d_pz + d_m*d_m)")  # per-daughter energy
    .Define("px_tot", "Sum(d_px)")
    .Define("py_tot", "Sum(d_py)")
    .Define("pz_tot", "Sum(d_pz)")
    .Define("E_tot",  "Sum(d_E)")
    .Define("m_inv",  "sqrt(E_tot*E_tot - px_tot*px_tot - py_tot*py_tot - pz_tot*pz_tot)")
)

c1 = ROOT.TCanvas("c1", "", 700, 500)
h_m = reco.Histo1D(("h_m", "reconstructed mass; m_{inv} [GeV]; events", 60, 0.9, 1.1), "m_inv")
h_m.Draw()
c1.Draw()
""")
mass_note = md(r"""
A sharp peak at **1 GeV** — you just reconstructed the HNL mass from its decay products.
""")

pt_md = md(r"""
## 3. TODO — reconstruct the transverse momentum

The HNL's transverse momentum is $p_T = \sqrt{p_{x,\mathrm{tot}}^2 + p_{y,\mathrm{tot}}^2}$.
Add one more `Define` for `pt` and histogram it (try a range 0–5 GeV).
""")
pt_todo = code(r"""
# TODO: define "pt" from px_tot and py_tot, then histogram it
# reco_pt = reco.Define("pt", "sqrt(px_tot*px_tot + py_tot*py_tot)")
# h_pt = reco_pt.Histo1D(("h_pt", "p_{T}; p_{T} [GeV]; events", 60, 0, 5), "pt")
# c = ROOT.TCanvas(); h_pt.Draw(); c.Draw()
""")
pt_done = code(r"""
reco_pt = reco.Define("pt", "sqrt(px_tot*px_tot + py_tot*py_tot)")
h_pt = reco_pt.Histo1D(("h_pt", "p_{T}; p_{T} [GeV]; events", 60, 0, 5), "pt")
c_pt = ROOT.TCanvas("c_pt", "", 700, 500)
h_pt.Draw()
c_pt.Draw()
""")

smear_md = md(r"""
## 4. Detector resolution: smear the momenta

A real detector never measures momentum perfectly. We model that by multiplying
each daughter's $p_x, p_y$ by a random factor $(1 + \sigma\cdot\mathcal{N}(0,1))$.

The little C++ helper below does the smearing. **You don't need to read it** — we just
declare it to ROOT once so RDataFrame can use it by name.
""")
smear_helper = code(r"""
ROOT.gInterpreter.Declare(r'''
#include <ROOT/RVec.hxx>
#include "TRandom.h"
// toy momentum smearing: scale each component by (1 + sigma * Gaussian)
ROOT::RVec<float> smear(const ROOT::RVec<float>& v, double sigma){
    ROOT::RVec<float> out(v.size());
    for (std::size_t i = 0; i < v.size(); ++i)
        out[i] = v[i] * (1.0 + sigma * gRandom->Gaus());
    return out;
}
''')
ROOT.gRandom.SetSeed(123)   # reproducible
""")
smear_fn_md = md(r"""
Now a helper that builds the smeared invariant-mass histogram for a given $\sigma$:
""")
smear_fn = code(r"""
def smeared_mass_hist(sigma, name):
    d = (df
        .Filter("d_px.size() == 2")
        .Define("px_s", f"smear(d_px, {sigma})")
        .Define("py_s", f"smear(d_py, {sigma})")
        .Define("E_s",  "sqrt(px_s*px_s + py_s*py_s + d_pz*d_pz + d_m*d_m)")
        .Define("m_s",  "sqrt(pow(Sum(E_s),2) - pow(Sum(px_s),2) - pow(Sum(py_s),2) - pow(Sum(d_pz),2))")
    )
    return d.Histo1D((name, "mass vs smearing; m_{inv} [GeV]; events", 60, 0.6, 1.4), "m_s")
""")

overlay_md = md(r"""
## 5. Overlay several resolutions

Draw the reconstructed mass for $\sigma = 0,\,1\%,\,5\%,\,10\%$ on the same axes.
Watch the peak broaden as the detector gets worse.
""")
overlay_code = code(r"""
sigmas = [0.00, 0.01, 0.05, 0.10]
colors = [ROOT.kCyan+1, ROOT.kBlue, ROOT.kViolet, ROOT.kMagenta]

c2 = ROOT.TCanvas("c2", "", 800, 600)
leg = ROOT.TLegend(0.62, 0.66, 0.88, 0.88)
leg.SetBorderSize(0)

hists = []   # keep references alive
for i, (s, col) in enumerate(zip(sigmas, colors)):
    h = smeared_mass_hist(s, f"h_{int(s*100)}").GetValue()
    h.SetLineColor(col); h.SetLineWidth(2); h.SetDirectory(0)
    hists.append(h)
    h.Draw("hist" if i == 0 else "hist same")
    leg.AddEntry(h, f"#sigma = {int(s*100)}%", "l")

leg.Draw()
c2.Draw()
""")

res_md = md(r"""
## 6. TODO — quantify the resolution

The peak gets wider with $\sigma$. Put a number on it: for each $\sigma$, print the
**mass resolution** = the standard deviation of the smeared mass.

*Hint:* build the same dataframe as in `smeared_mass_hist`, but return `d.StdDev("m_s")`
instead of a histogram.
""")
res_todo = code(r"""
# TODO: for each sigma in [0.0, 0.01, 0.05, 0.10], print StdDev of the smeared mass
# for s in [0.0, 0.01, 0.05, 0.10]:
#     d = (df.Filter("d_px.size() == 2")
#            .Define("px_s", f"smear(d_px, {s})")
#            ... )
#     print(f"sigma={s:>4}:  mass resolution = {d.StdDev('m_s').GetValue():.4f} GeV")
""")
res_done = code(r"""
def mass_resolution(sigma):
    d = (df
        .Filter("d_px.size() == 2")
        .Define("px_s", f"smear(d_px, {sigma})")
        .Define("py_s", f"smear(d_py, {sigma})")
        .Define("E_s",  "sqrt(px_s*px_s + py_s*py_s + d_pz*d_pz + d_m*d_m)")
        .Define("m_s",  "sqrt(pow(Sum(E_s),2) - pow(Sum(px_s),2) - pow(Sum(py_s),2) - pow(Sum(d_pz),2))")
    )
    return d.StdDev("m_s").GetValue()

for s in [0.00, 0.01, 0.05, 0.10]:
    print(f"sigma={s:>4}:  mass resolution = {mass_resolution(s):.4f} GeV")
""")

interp = md(r"""
## 7. What did you find?

- The **mean** mass stays at ~1 GeV — the smearing is unbiased on average.
- The **width** grows with $\sigma$ — that's your detector resolution propagating into
  the measurement. A wider peak makes a signal harder to separate from background.

This is exactly the trade-off a real experiment optimises.
""")

more_md = md(r"""
## 8. (optional) The full study

The complete analysis also reconstructs the **decay vertex** from the two daughter
tracks (the DOCA method) and measures how the *vertex* resolution degrades with
$\sigma$. It's in `reference/llp_simple_analysis.py` (and the C++ `.C`), with a
detailed write-up in `reference/README.md`. Have a look if you're curious.
""")

git_md = md(r"""
## 9. Save your work with git

You changed this notebook — snapshot it so it's safe:
""")
git_code = code(r"""
!git add tutorial.ipynb
!git commit -m "Complete pT and resolution TODOs"
""")

# ---- tutorial ----
tut = new_notebook()
tut.cells = [
    intro, setup_md, setup_code,
    data_md, path_code, data_code, vtx_md, vtx_code, vtx_note,
    mass_md, mass_code, mass_note,
    pt_md, pt_todo,
    smear_md, smear_helper, smear_fn_md, smear_fn,
    overlay_md, overlay_code,
    res_md, res_todo,
    interp, more_md, git_md, git_code,
]

# ---- solutions ----
sol = new_notebook()
sol.cells = [
    md("# LLP hands-on — SOLUTIONS\n\nCompleted version of `tutorial.ipynb`."),
    setup_md, setup_code,
    data_md, path_code, data_code, vtx_md, vtx_code, vtx_note,
    mass_md, mass_code, mass_note,
    md("## 3. pT (solution)"), pt_done,
    smear_md, smear_helper, smear_fn_md, smear_fn,
    overlay_md, overlay_code,
    md("## 6. Resolution (solution)"), res_done,
    interp, more_md,
]

for nb in (tut, sol):
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    nb.metadata["language_info"] = {"name": "python"}

import nbformat
with open("tutorial.ipynb", "w") as fh:  nbformat.write(tut, fh)
with open("solutions.ipynb", "w") as fh: nbformat.write(sol, fh)
for n in ("tutorial.ipynb", "solutions.ipynb"):
    nbformat.validate(nbformat.read(n, as_version=4)); print("valid:", n)
