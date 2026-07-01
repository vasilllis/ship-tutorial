"""SHiP hands-on: reconstruct an LLP step by step — code export of tutorial.ipynb
Run on SWAN (a notebook shows the plots inline with %jsroot on)."""

# # SHiP hands-on: reconstruct an LLP, step by step
#
# We're going to build up a **real analysis** — the momentum-smearing study in
# `reference/llp_simple_analysis.py` — one small piece at a time. A Heavy Neutral
# Lepton (HNL) decays to two daughters at a displaced vertex; we will
#
# 1. reconstruct its **invariant mass**,
# 2. reconstruct its **decay vertex** from the two daughter tracks, and
# 3. watch how detector **momentum resolution** smears both.
#
# **Nothing to solve.** Run each cell (Shift+Enter), read the short note, look at the
# plot. Where you see **▸ Try**, change a number and re-run.

# ## 0. Setup

# (notebook only) %jsroot on
import ROOT
print("ROOT", ROOT.gROOT.GetVersion())

# ## 1. Open the data
#
# The sample lives on **EOS**. Each row of the tree `Events` is one HNL decay, with the
# two daughters stored as short vectors (`d_px`, `d_py`, `d_pz`, `d_m`) and the true
# decay vertex (`vtx_x`, `vtx_y`, `vtx_z`).

import os, glob
EOS_DIR = "/eos/experiment/ship/user/matclim"
named   = f"{EOS_DIR}/HNL_1.000e+00_1.000e+00_0.000e+00_1.000e+00_0.000e+00_data.root"
DATA = named if os.path.exists(named) else (
       sorted(glob.glob(f"{EOS_DIR}/**/*.root", recursive=True)) or ["data/hnl_data.root"])[0]
print("reading:", DATA)

df = ROOT.RDataFrame("Events", DATA)
print("events:", df.Count().GetValue())

# ## 2. Step one of the analysis: require two daughters
#
# The very first line of the reference analysis keeps only events with exactly two
# daughters. `Filter` does that; the label in quotes shows up in a cut-flow report.

base = df.Filter("d_px.size() == 2", "exactly two daughters")

print("kept:", base.Count().GetValue(), "events")
base.Report().Print()      # the cut-flow

# ## 3. The displaced vertex
#
# The HNL is long-lived, so it decays tens of metres downstream. Just histogram the
# truth vertex position `vtx_z`:

c1 = ROOT.TCanvas("c1", "", 700, 500)
h_z = base.Histo1D(("h_z", "decay vertex; z [m]; decays", 100, 0, 100), "vtx_z")
h_z.Draw(); c1.Draw()

# ## 4. Reconstruct the invariant mass
#
# Add the two daughter 4-momenta and take the mass:
# $m=\sqrt{E^2-p_x^2-p_y^2-p_z^2}$. Each `Define` builds one new column; `Sum(...)`
# adds the two entries of a daughter vector. (This is `llpFourMomentum` + `recoEnergy`
# from the reference, written as short expressions.)

reco = (base
    .Define("d_E",   "sqrt(d_px*d_px + d_py*d_py + d_pz*d_pz + d_m*d_m)")  # per-daughter energy
    .Define("m_inv", "sqrt(pow(Sum(d_E),2) - pow(Sum(d_px),2) - pow(Sum(d_py),2) - pow(Sum(d_pz),2))")
    .Define("pt",    "sqrt(pow(Sum(d_px),2) + pow(Sum(d_py),2))")
)

c2 = ROOT.TCanvas("c2", "", 700, 500)
h_m = reco.Histo1D(("h_m", "invariant mass; m_{inv} [GeV]; events", 60, 0.9, 1.1), "m_inv")
h_m.Draw(); c2.Draw()

# A sharp peak at **1 GeV** — the HNL mass. **▸ Try:** plot `pt` instead of `m_inv`.

# ## 5. The toolbox (provided)
#
# The next steps need a bit of geometry — reconstructing the decay vertex from two
# tracks, and a smearing function. Here they are as small C++ helpers, declared once so
# RDataFrame can use them by name. **You don't need to read this cell** — just run it.

if not hasattr(ROOT, "docaVertex"):
    ROOT.gInterpreter.Declare(r'''
    #include "Math/Vector3D.h"
    #include <ROOT/RVec.hxx>
    #include "TRandom.h"
    using ROOT::Math::XYZVector;

    // unit direction of daughter 0 / daughter 1
    XYZVector dir0(const ROOT::RVec<float>&px,const ROOT::RVec<float>&py,const ROOT::RVec<float>&pz){
        return XYZVector(px[0],py[0],pz[0]).Unit(); }
    XYZVector dir1(const ROOT::RVec<float>&px,const ROOT::RVec<float>&py,const ROOT::RVec<float>&pz){
        return XYZVector(px[1],py[1],pz[1]).Unit(); }

    // anchor point of a track: foot of the perpendicular from the beam origin
    XYZVector perigee(const XYZVector&V,const XYZVector&d){
        XYZVector u=d.Unit(); return V-u*(V.Dot(u)); }

    // decay vertex = midpoint of the closest approach of the two tracks (DOCA)
    // returns {doca, x, y, z}
    ROOT::RVec<double> docaVertex(const XYZVector&A1,const XYZVector&u1in,
                                  const XYZVector&A2,const XYZVector&u2in){
        XYZVector u1=u1in.Unit(),u2=u2in.Unit(),w0=A1-A2;
        double b=u1.Dot(u2),d=u1.Dot(w0),e=u2.Dot(w0),den=1.0-b*b,sc,tc;
        if(den<1e-9){sc=0.0;tc=e;}else{sc=(b*e-d)/den;tc=(e-b*d)/den;}
        XYZVector P1=A1+u1*sc,P2=A2+u2*tc,mid=(P1+P2)*0.5;
        return {(P1-P2).R(),mid.X(),mid.Y(),mid.Z()}; }

    // toy momentum smearing: scale each component by (1 + sigma*Gaussian)
    ROOT::RVec<float> smear(const ROOT::RVec<float>&v,double sigma){
        ROOT::RVec<float> o(v.size());
        for(size_t i=0;i<v.size();++i) o[i]=v[i]*(1.0+sigma*gRandom->Gaus());
        return o; }
    ''')
ROOT.gRandom.SetSeed(123)
print("toolbox ready")

# ## 6. Add detector resolution to the mass
#
# A real detector doesn't measure momenta perfectly. Smear each daughter's $p_x,p_y$ by
# $(1+\sigma\,\mathcal{N}(0,1))$ and rebuild the mass. This wraps the mass reconstruction
# from step 4 so we can try several $\sigma$:

def smeared_mass(sigma, name):
    return (base
        .Define("px_s", f"smear(d_px, {sigma})")
        .Define("py_s", f"smear(d_py, {sigma})")
        .Define("E_s",  "sqrt(px_s*px_s + py_s*py_s + d_pz*d_pz + d_m*d_m)")
        .Define("m_s",  "sqrt(pow(Sum(E_s),2)-pow(Sum(px_s),2)-pow(Sum(py_s),2)-pow(Sum(d_pz),2))")
        .Histo1D((name, "mass vs resolution; m_{inv} [GeV]; events", 60, 0.6, 1.4), "m_s"))

sigmas = [0.00, 0.01, 0.05, 0.10]
colors = [ROOT.kCyan+1, ROOT.kBlue, ROOT.kViolet, ROOT.kMagenta]
c3 = ROOT.TCanvas("c3", "", 800, 600)
leg = ROOT.TLegend(0.62, 0.66, 0.88, 0.88); leg.SetBorderSize(0)
hists = []
for i,(s,col) in enumerate(zip(sigmas, colors)):
    h = smeared_mass(s, f"hm_{int(s*100)}").GetValue()
    h.SetLineColor(col); h.SetLineWidth(2); h.SetDirectory(0); hists.append(h)
    h.Draw("hist" if i==0 else "hist same")
    leg.AddEntry(h, f"#sigma = {int(s*100)}%", "l")
leg.Draw(); c3.Draw()

# The peak **broadens** as resolution worsens, while its centre stays at 1 GeV.

# ## 7. Reconstruct the decay vertex (DOCA)
#
# Now the vertex. Each daughter is a straight **track**: a line through an anchor point
# (`perigee`) pointing along its momentum. The decay vertex is where the two lines come
# **closest** — the `docaVertex` helper returns that midpoint.
#
# The function below builds the vertex from the *smeared* directions and compares it to
# the truth vertex, giving the residual in $z$: `dz_res = z_reco − z_true`.

def vertex_residuals(sigma):
    return (base
        .Define("vtx",  "ROOT::Math::XYZVector(vtx_x, vtx_y, vtx_z)")
        .Define("u0",   "dir0(d_px, d_py, d_pz)")     # true directions
        .Define("u1",   "dir1(d_px, d_py, d_pz)")
        .Define("A1",   "perigee(vtx, u0)")           # track anchors
        .Define("A2",   "perigee(vtx, u1)")
        .Define("px_s", f"smear(d_px, {sigma})")      # smear, then re-take directions
        .Define("py_s", f"smear(d_py, {sigma})")
        .Define("u0_s", "dir0(px_s, py_s, d_pz)")
        .Define("u1_s", "dir1(px_s, py_s, d_pz)")
        .Define("doca", "docaVertex(A1, u0_s, A2, u1_s)")   # {doca, x, y, z}
        .Define("dz_res", "doca[3] - vtx_z")          # longitudinal vertex residual
    )
print("built vertex_residuals()")

# ### Closure check
#
# At $\sigma=0$ the smeared directions equal the true ones, so the reconstructed vertex
# **is** the truth and `dz_res` is essentially zero. Turn on smearing and it spreads out.
# This zero-at-zero behaviour is a built-in sanity check that the reconstruction is
# consistent.

c4 = ROOT.TCanvas("c4", "", 800, 600)
leg = ROOT.TLegend(0.62, 0.66, 0.88, 0.88); leg.SetBorderSize(0)
hists = []
for i,(s,col) in enumerate(zip(sigmas, colors)):
    h = vertex_residuals(s).Histo1D(
        (f"hz_{int(s*100)}", "vertex residual; #Deltaz = z_{reco}-z_{true} [m]; events", 80, -20, 20),
        "dz_res").GetValue()
    h.SetLineColor(col); h.SetLineWidth(2); h.SetDirectory(0); hists.append(h)
    h.Draw("hist" if i==0 else "hist same")
    leg.AddEntry(h, f"#sigma = {int(s*100)}%", "l")
leg.Draw(); c4.Draw()

# ## 8. Put a number on it: vertex resolution vs smearing
#
# The width of that residual **is** the vertex resolution. Print the standard deviation
# of `dz_res` at each smearing level — ~0 at $\sigma=0$ (closure), growing with $\sigma$.

print("sigma   vertex resolution sigma_z")
for s in sigmas:
    sz = vertex_residuals(s).StdDev("dz_res").GetValue()
    print(f"{s:>5}   {sz:.3f} m")

# ## That's the analysis
#
# You reconstructed the HNL **mass** and its **decay vertex**, then measured how detector
# **momentum resolution** blurs each. That's exactly what `reference/llp_simple_analysis.py`
# does — plus two more observables ($p_T$, the transverse residual $\Delta r$, the vertex
# impact parameter) and a nicer `THStack` layout. Open it and compare.
#
# **▸ Try:** in `vertex_residuals`, also define a transverse residual
# `dr_res = sqrt(doca[1]^2 + doca[2]^2) - sqrt(vtx_x^2 + vtx_y^2)` and histogram it.
