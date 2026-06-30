"""
llp_simple_analysis.py
======================

A PyROOT translation of ``llp_simple_analysis.C``: a small RDataFrame
analysis that studies how detector resolution (a momentum-smearing
parameter ``sigma``) affects a few reconstructed observables for a
long-lived particle (LLP) decaying to two daughters:

    1. The reconstructed invariant mass of the LLP
       (from the smeared daughter momenta).
    2. The reconstructed transverse momentum p_T of the LLP.
    3. The transverse vertex residual  Delta r = r_reco - r_true.
    4. The longitudinal vertex residual Delta z = z_reco - z_true.
    5. The 3-D impact parameter of the reconstructed LLP direction
       with respect to the reconstructed decay vertex.

The reconstructed decay vertex is obtained from the two smeared
daughter directions via the **DOCA** method (Distance Of Closest
Approach between two skew lines).  At sigma = 0 (perfect knowledge
of the momenta) the reconstructed vertex equals the true vertex,
so all residuals collapse onto zero -- a sanity check that the
arithmetic and the smearing convention are consistent.

How it is structured (and why)
------------------------------

The C++ original uses ``ROOT::RDataFrame``, which is the high-level
ROOT analysis framework: you build a *graph* of operations
(``Filter``, ``Define``, ``Histo1D``, ...) and ROOT compiles and
runs the whole thing in a single multi-threaded pass over the
input file.  Python is a poor fit for the inner loop -- a Python
callback inside ``Define`` runs once per event in interpreted
Python, around two orders of magnitude slower than the
equivalent C++.

The idiomatic PyROOT solution, which we use here, is to declare
the small C++ helper functions to ROOT's interpreter via
``ROOT.gInterpreter.Declare(...)``.  RDataFrame then resolves the
``Define`` arguments to those compiled C++ symbols, exactly as in
the C++ source.  The result is identical performance to the C++
version, while every line that drives the analysis sits in
Python.

Students who want to inspect or modify the per-event logic can
simply edit the embedded C++ string at the top of this file.

Output
------

``llp_smearing.png``  -- a 3x2 grid of TH1D distributions, one
per observable, with the four sigma values overlaid as
non-stacked coloured curves.
"""

import ROOT


# ---------------------------------------------------------------
# 1.  Helper functions in C++, declared once at module load.
# ---------------------------------------------------------------
# RDataFrame's Define is happiest with compiled C++ callables --
# Python lambdas there force a per-event interpreter callback.
# We declare the helpers to ROOT's C++ JIT here, then reference
# them by *name* in the RDataFrame graph below.
#
# This is the only chunk of C++ in the file.  Everything that
# drives the analysis -- the graph construction, the histogram
# parameters, the canvas layout -- is plain Python.

ROOT.gInterpreter.Declare(r"""
#include <cmath>
#include <ROOT/RVec.hxx>
#include "Math/Vector4D.h"     // ROOT::Math::PxPyPzEVector
#include "Math/Vector3D.h"     // ROOT::Math::XYZVector
#include "TRandom3.h"

using ROOT::Math::PxPyPzEVector;
using ROOT::Math::XYZVector;

// LLP 4-momentum: vector sum of the daughter 4-momenta.
inline PxPyPzEVector llpFourMomentum(const ROOT::RVec<float>& px,
                                     const ROOT::RVec<float>& py,
                                     const ROOT::RVec<float>& pz,
                                     const ROOT::RVec<float>& e)
{
    PxPyPzEVector sum(0., 0., 0., 0.);
    for (std::size_t i = 0; i < px.size(); ++i)
        sum += PxPyPzEVector(px[i], py[i], pz[i], e[i]);
    return sum;
}

// 3-D impact parameter of `point' w.r.t. a line through the
// origin with direction `dir' (length of the rejection of point
// onto dir).  Returns 0 if dir is degenerate.
inline double impactParameter3D(const XYZVector& point,
                                const XYZVector& dir)
{
    const double n = dir.R();
    if (n == 0.) return 0.;
    return point.Cross(dir / n).R();
}

// Reconstruct daughter energies from (smeared) 3-momenta + masses
// using E = sqrt(p^2 + m^2). The daughter masses are taken as
// truth -- only the momentum is smeared, not the mass.
inline ROOT::RVec<float> recoEnergy(const ROOT::RVec<float>& px,
                                    const ROOT::RVec<float>& py,
                                    const ROOT::RVec<float>& pz,
                                    const ROOT::RVec<float>& m)
{
    ROOT::RVec<float> e(px.size());
    for (std::size_t i = 0; i < px.size(); ++i)
        e[i] = std::sqrt(px[i]*px[i] + py[i]*py[i]
                        + pz[i]*pz[i] + m[i]*m[i]);
    return e;
}

// Unit direction of daughter 0 / daughter 1.
inline XYZVector dir0(const ROOT::RVec<float>& px,
                      const ROOT::RVec<float>& py,
                      const ROOT::RVec<float>& pz)
{ return XYZVector(px[0], py[0], pz[0]).Unit(); }

inline XYZVector dir1(const ROOT::RVec<float>& px,
                      const ROOT::RVec<float>& py,
                      const ROOT::RVec<float>& pz)
{ return XYZVector(px[1], py[1], pz[1]).Unit(); }

// Perigee of the line {V + t * dir} relative to the origin --
// i.e. the point on that line closest to (0,0,0). Used as a
// stable "anchor" A on each track so DOCA can solve in
// closed form (see docaVertex).
inline XYZVector perigee(const XYZVector& V, const XYZVector& dir)
{
    const XYZVector u = dir.Unit();
    return V - u * (V.Dot(u));
}

// DOCA between two skew lines and the midpoint of the segment
// joining their closest approach. Returns {doca, mid.x, mid.y, mid.z}.
//
// Math: with line 1 = A1 + sc*u1 and line 2 = A2 + tc*u2,
// the closest-approach parameters (sc, tc) satisfy the linear
// system from minimising |line1(sc) - line2(tc)|^2. The denom
// 1 - (u1.u2)^2 vanishes when the lines are parallel; we then
// fall back to the parameter pair that handles that limit
// gracefully (sc=0, tc=u2.(A1-A2)).
inline ROOT::RVec<double> docaVertex(const XYZVector& A1,
                                     const XYZVector& u1in,
                                     const XYZVector& A2,
                                     const XYZVector& u2in)
{
    const XYZVector u1 = u1in.Unit();
    const XYZVector u2 = u2in.Unit();
    const XYZVector w0 = A1 - A2;
    const double b = u1.Dot(u2);
    const double d = u1.Dot(w0);
    const double e = u2.Dot(w0);
    const double denom = 1.0 - b * b;
    double sc, tc;
    if (denom < 1e-9) { sc = 0.0;  tc = e; }              // parallel
    else              { sc = (b*e - d) / denom;
                        tc = (e - b*d) / denom; }
    const XYZVector P1  = A1 + u1 * sc;
    const XYZVector P2  = A2 + u2 * tc;
    const XYZVector mid = (P1 + P2) * 0.5;
    const double doca = (P1 - P2).R();
    return {doca, mid.X(), mid.Y(), mid.Z()};
}

// Smear a single (px or py) column by a multiplicative Gaussian
// factor (1 + sigma * N(0,1)). Deterministic given `entry' and
// `salt' -- different (sigma, axis) combinations use different
// salt values so they don't collide.
//
// IMPORTANT: this is NOT a physically realistic detector model
// -- it is a simple multiplicative smearing used to study how
// the downstream observables respond to per-component momentum
// uncertainty. A real detector smears momentum magnitude with a
// p_T-dependent resolution, smears angles separately, etc.
inline ROOT::RVec<float> smearVec(ULong64_t entry, unsigned salt,
                                  double sigma,
                                  const ROOT::RVec<float>& v)
{
    const UInt_t seed =
       (UInt_t)((entry*2654435761ULL + salt*40503ULL + 12345ULL)
                & 0xffffffffULL);
    TRandom3 rng(seed);
    ROOT::RVec<float> out(v.size());
    for (std::size_t i = 0; i < v.size(); ++i)
        out[i] = v[i] * (1.f + (float)(sigma * rng.Gaus()));
    return out;
}
""")


# ---------------------------------------------------------------
# 2.  Build the smeared analysis sub-graph.
# ---------------------------------------------------------------
def make_smeared(base, sigma, lvl):
    """Attach a chain of ``Define`` columns that compute the
    smeared analysis for one value of sigma.

    Parameters
    ----------
    base   : ROOT.RDF.RNode
        The pre-built RDataFrame node carrying the truth columns
        (``vtx``, ``A1``, ``A2``, ``d_px``, ...).
    sigma  : float
        Relative Gaussian width of the multiplicative momentum
        smearing.  ``sigma = 0`` reproduces the truth quantities.
    lvl    : int
        Salt to keep different sigma levels' random sequences
        independent (px gets ``lvl*10 + 1``, py gets ``lvl*10 + 2``).

    Returns
    -------
    ROOT.RDF.RNode
        A new node with the smeared observables defined:
        ``m_inv_s, pt_s, dr_res, dz_res, ip_vtx``.

    Notes
    -----
    Only ``d_px`` and ``d_py`` are smeared -- ``d_pz`` is left
    untouched.  This matches the C++ source and is a deliberate
    pedagogical simplification: it isolates the effect of a
    transverse-momentum mis-measurement, the residuals in the
    bend plane being usually the dominant uncertainty.
    """
    # Smearing salt values: per-sigma-level identifiers for the
    # two smeared axes. Choosing distinct salts means px and py
    # smearings within one event are independent draws (and that
    # different sigma levels are independent of each other).
    salt_px = lvl * 10 + 1
    salt_py = lvl * 10 + 2

    # We use a small JIT-string per Define so we can embed the
    # sigma/salt values directly. (RDataFrame's Python overload
    # of Define accepts a callable name OR a free-form C++
    # expression string; the latter is what we lean on here.)
    smear_px = (
        f"smearVec(rdfentry_, {salt_px}u, {sigma}, d_px)"
    )
    smear_py = (
        f"smearVec(rdfentry_, {salt_py}u, {sigma}, d_py)"
    )

    return (base
        # ---- smear the two transverse-momentum components ----
        .Define("d_px_s", smear_px)
        .Define("d_py_s", smear_py)
        # ---- recompute daughter energies from smeared momenta ----
        # using the *true* masses (d_m).
        .Define("d_E_s", "recoEnergy(d_px_s, d_py_s, d_pz, d_m)")
        # ---- LLP 4-momentum from the smeared daughters ----
        .Define("p_sm", "llpFourMomentum(d_px_s, d_py_s, d_pz, d_E_s)")
        # ---- scalar observables: smeared invariant mass and pT ----
        .Define("m_inv_s", "p_sm.M()")
        .Define("pt_s",    "p_sm.Pt()")
        # ---- smeared daughter directions ----
        .Define("u0_s", "dir0(d_px_s, d_py_s, d_pz)")
        .Define("u1_s", "dir1(d_px_s, d_py_s, d_pz)")
        # ---- reconstructed vertex via DOCA of the two tracks ----
        # Each track is the line through its perigee anchor A_i
        # along the smeared direction u_i_s. docaVertex returns
        # {doca, mid.x, mid.y, mid.z}; we keep the midpoint as
        # the vertex estimate.
        .Define("doca_vec",      "docaVertex(A1, u0_s, A2, u1_s)")
        .Define("reco_vtx_doca", "ROOT::Math::XYZVector(doca_vec[1],"
                                  " doca_vec[2], doca_vec[3])")
        # ---- vertex residuals (signed, centred at 0) ----
        # transverse: rho_reco - rho_true; longitudinal: z_reco - z_true.
        # Both -> 0 at sigma=0, since reco_vtx_doca = vtx in that limit.
        .Define("dr_res", "reco_vtx_doca.Rho() - vtx.Rho()")
        .Define("dz_res", "reco_vtx_doca.Z()   - vtx.Z()")
        # ---- LLP 3-D impact parameter relative to its own
        # reconstructed direction -- expected to vanish if the
        # LLP momentum points back at the vertex, useful as a
        # closure check on the kinematic reconstruction.
        .Define("nhat_s",  "p_sm.Vect().Unit()")
        .Define("ip_vtx",  "reco_vtx_doca.Cross(nhat_s).R()")
    )


# ---------------------------------------------------------------
# 3.  Drawing helper: overlay a vector of histograms in one pad.
# ---------------------------------------------------------------
def draw_stack(canvas, pad, title, hv, labels, logy):
    """Overlay several histograms in a single canvas pad,
    coloured along the active palette.

    Parameters
    ----------
    canvas, pad : the TCanvas and pad number to draw into.
    title       : title string for the THStack (use
                  "title;xlabel;ylabel" semicolon syntax).
    hv          : list of RResultPtr<TH1D> (RDataFrame Histo1D
                  results).  Evaluating any of them triggers
                  the whole event-loop graph the first time.
    labels      : legend entries, one per histogram.
    logy        : set Y-axis logarithmic if True.
    """
    p = canvas.cd(pad)
    if logy:
        p.SetLogy()

    st = ROOT.THStack(f"st{pad}", title)
    lg = ROOT.TLegend(0.58, 0.64, 0.88, 0.88)
    lg.SetBorderSize(0)
    lg.SetFillStyle(0)

    n    = len(hv)
    ncol = ROOT.TColor.GetNumberOfColors()

    # Keep clones alive while the canvas references them by
    # appending to a list owned by the function caller's scope
    # (we attach to the pad through Draw, but THStack does not
    # own the source histograms).
    keepalive = []

    for i, hres in enumerate(hv):
        # Clone so we can recolour without mutating the
        # RDataFrame-owned result, and detach from any TDirectory.
        h = hres.GetValue().Clone(f"c_{pad}_{i}")
        h.SetDirectory(0)
        # Pick a palette index spanning the colour range:
        # 0 -> 0, n-1 -> ncol-1, with linear interpolation.
        if n <= 1:
            idx = ncol // 2
        else:
            idx = int(round(i / (n - 1) * (ncol - 1)))
        h.SetLineColor(ROOT.TColor.GetColorPalette(idx))
        h.SetLineWidth(2)
        st.Add(h, "hist")
        lg.AddEntry(h, labels[i], "l")
        keepalive.append(h)

    st.Draw("nostack")
    lg.Draw()
    # ROOT only keeps weak references to st/lg via the canvas;
    # return them so the caller can keep them alive until SaveAs.
    return st, lg, keepalive


# ---------------------------------------------------------------
# 4.  Main analysis.
# ---------------------------------------------------------------
def llp_simple_analysis():
    """Run the full LLP smearing study and write llp_smearing.png."""

    ROOT.EnableImplicitMT()

    # Input tree -- same hardcoded file as the C++ version.
    df = ROOT.RDataFrame(
        "Events",
        "./HNL_1.000e+00_1.000e+00_0.000e+00_1.000e+00_0.000e+00_data.root"
    )

    # The primary vertex (PV) is the origin in this MC.  We pass
    # it into perigee() implicitly: the perigee of a line
    # relative to the origin is V - u * (V . u).  The PV here is
    # NOT used directly -- it is documented for completeness.
    # PV = ROOT.Math.XYZVector(0., 0., 0.)

    # ---- truth-level columns ----
    # We require exactly two daughters and build the truth
    # quantities needed downstream by all sigma-levels:
    #   vtx, Lxy, L3D                 -- true decay vertex
    #   u0_true, u1_true              -- true daughter directions
    #   A1, A2                        -- perigee anchors of the
    #                                    two true tracks (used
    #                                    by docaVertex as the
    #                                    reference points)
    base = (df
        .Filter("d_px.size() == 2u", "exactly two daughters")
        .Define("vtx", "ROOT::Math::XYZVector(vtx_x, vtx_y, vtx_z)")
        .Define("Lxy", "vtx.Rho()")
        .Define("L3D", "vtx.R()")
        .Define("u0_true", "dir0(d_px, d_py, d_pz)")
        .Define("u1_true", "dir1(d_px, d_py, d_pz)")
        .Define("A1", "perigee(vtx, u0_true)")
        .Define("A2", "perigee(vtx, u1_true)")
    )

    # The smearing levels we sweep over and matching legend labels.
    sig = [0.00, 0.01, 0.05, 0.10]
    lab = [
        "true (#sigma=0)",
        "#sigma=1%",
        "#sigma=5%",
        "#sigma=10%",
    ]
    NB = 60  # number of histogram bins (same for every observable)

    # ---- adaptive histogram ranges ----
    # We compute the min/max of each smeared observable at the
    # WORST smearing (10%) to size the axis ranges sensibly.
    # RDataFrame's Min/Max return lazy result pointers --
    # dereferencing one triggers the event loop the first time.
    n10  = make_smeared(base, 0.10, 99)
    Mmin = n10.Min("m_inv_s");   Mmax = n10.Max("m_inv_s")
    Pmax = n10.Max("pt_s")
    DRlo = n10.Min("dr_res");    DRhi = n10.Max("dr_res")
    DZlo = n10.Min("dz_res");    DZhi = n10.Max("dz_res")
    IVmx = n10.Max("ip_vtx")

    mlo, mhi = Mmin.GetValue(), Mmax.GetValue()
    mpad = 0.1 * (mhi - mlo) + 1e-6
    MLO, MHI = mlo - mpad, mhi + mpad
    PHI = Pmax.GetValue() * 1.1 + 1e-6
    DRA = max(abs(DRlo.GetValue()), abs(DRhi.GetValue())) * 1.1 + 1e-9
    DZA = max(abs(DZlo.GetValue()), abs(DZhi.GetValue())) * 1.1 + 1e-9
    IVH = IVmx.GetValue() * 1.1 + 1e-9

    # ---- book one histogram per observable, per sigma level ----
    # We collect the result pointers into lists and let RDataFrame
    # batch them into a single multi-threaded pass over the
    # input tree on first dereference (done implicitly inside
    # draw_stack via Histo1D->Clone()).
    HM, HP, HDR, HDZ, HIV = [], [], [], [], []

    # Mean / std collectors for the quick console summary at the end.
    mass_mean, dr_std, dz_std = [], [], []

    for i, s in enumerate(sig):
        nd = make_smeared(base, s, i)
        HM.append(nd.Histo1D(
            (f"hM{i}",  "invariant mass;m_{inv} [GeV];events",
             NB, MLO, MHI), "m_inv_s"))
        HP.append(nd.Histo1D(
            (f"hP{i}",  "reconstructed p_{T};p_{T} [GeV];events",
             NB, 0., PHI), "pt_s"))
        HDR.append(nd.Histo1D(
            (f"hDR{i}", "transverse vertex residual;"
             "#Deltar = r_{reco}-r_{true} [m];events",
             NB, -DRA, DRA), "dr_res"))
        HDZ.append(nd.Histo1D(
            (f"hDZ{i}", "longitudinal vertex residual;"
             "#Deltaz = z_{reco}-z_{true} [m];events",
             NB, -DZA, DZA), "dz_res"))
        HIV.append(nd.Histo1D(
            (f"hIV{i}", "vertex IP;d^{3D}_{LLP} [m];events",
             NB, 0., IVH), "ip_vtx"))
        mass_mean.append(nd.Mean("m_inv_s"))
        dr_std.append(nd.StdDev("dr_res"))
        dz_std.append(nd.StdDev("dz_res"))

    # ---- draw ----
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kCool)

    c = ROOT.TCanvas("c", "LLP momentum-smearing study", 1500, 900)
    c.Divide(3, 2)

    # The draw_stack helper returns the THStack/TLegend it makes
    # so Python's GC doesn't collect them before SaveAs fires.
    # We collect them all into one list.
    keepalive = []
    keepalive.append(draw_stack(c, 1, "invariant mass;m_{inv} [GeV];events",
                                HM,  lab, False))
    keepalive.append(draw_stack(c, 2, "reconstructed p_{T};p_{T} [GeV];events",
                                HP,  lab, False))
    keepalive.append(draw_stack(c, 3, "transverse vertex residual;"
                                "#Deltar [m];events", HDR, lab, False))
    keepalive.append(draw_stack(c, 4, "longitudinal vertex residual;"
                                "#Deltaz [m];events", HDZ, lab, False))
    keepalive.append(draw_stack(c, 5, "vertex impact parameter;"
                                "d^{3D}_{LLP} [m];events", HIV, lab, True))

    c.SaveAs("llp_smearing.png")

    # ---- short console summary ----
    # Useful sanity check: <m_inv> should hardly move with sigma
    # (the smearing is multiplicative and unbiased on average),
    # while StdDev(dr_res) and StdDev(dz_res) should grow with
    # sigma roughly linearly at small sigma.
    print("\n  sigma     <m_inv>     std(dr)     std(dz)")
    print(  "  ------  ----------  ----------  ----------")
    for s, mm, dr, dz in zip(sig, mass_mean, dr_std, dz_std):
        print(f"  {s:5.3f}  {mm.GetValue():10.4f}"
              f"  {dr.GetValue():10.4e}  {dz.GetValue():10.4e}")


# Run when executed as a script: ``python llp_simple_analysis.py``
if __name__ == "__main__":
    llp_simple_analysis()
