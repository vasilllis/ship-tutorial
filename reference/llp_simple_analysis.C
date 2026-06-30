using ROOT::Math::PxPyPzEVector;
using ROOT::Math::XYZVector;

// LLP 4-momentum = vector sum of the daughter 4-momenta.
PxPyPzEVector llpFourMomentum(const ROOT::RVec<float> &px,const ROOT::RVec<float> &py,const ROOT::RVec<float> &pz,const ROOT::RVec<float> &e){

  PxPyPzEVector sum(0., 0., 0., 0.);
   for (std::size_t i = 0; i < px.size(); ++i)
      sum += PxPyPzEVector(px[i], py[i], pz[i], e[i]);
   return sum;
}

// 3-D impact parameter 
double impactParameter3D(const XYZVector &point, const XYZVector &dir){
   
   const double n = dir.R();
   if (n == 0.) return 0.;
   return point.Cross(dir / n).R();
}

// Reconstruct daughter energies from (smeared) momenta and the true masses.
ROOT::RVec<float> recoEnergy(const ROOT::RVec<float> &px,const ROOT::RVec<float> &py,const ROOT::RVec<float> &pz,const ROOT::RVec<float> &m){
   
  ROOT::RVec<float> e(px.size());
   for (std::size_t i = 0; i < px.size(); ++i)
      e[i] = std::sqrt(px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i] + m[i] * m[i]);
   return e;
}

// Unit direction of daughter 0 / daughter 1.
XYZVector dir0(const ROOT::RVec<float> &px, const ROOT::RVec<float> &py,const ROOT::RVec<float> &pz){
  return XYZVector(px[0], py[0], pz[0]).Unit(); 
}

XYZVector dir1(const ROOT::RVec<float> &px, const ROOT::RVec<float> &py,const ROOT::RVec<float> &pz){ 
  return XYZVector(px[1], py[1], pz[1]).Unit(); 
}

XYZVector perigee(const XYZVector &V, const XYZVector &dir){
   const XYZVector u = dir.Unit();
   return V - u * (V.Dot(u));
}

ROOT::RVec<double> docaVertex(const XYZVector &A1, const XYZVector &u1in,const XYZVector &A2, const XYZVector &u2in){

   const XYZVector u1 = u1in.Unit();
   const XYZVector u2 = u2in.Unit();
   const XYZVector w0 = A1 - A2;
   const double b = u1.Dot(u2);          
   const double d = u1.Dot(w0);
   const double e = u2.Dot(w0);
   const double denom = 1.0 - b * b;
   double sc, tc;
   if (denom < 1e-9) {                  
      sc = 0.0;
      tc = e;
   } else {
      sc = (b * e - d) / denom;
      tc = (e - b * d) / denom;
   }
   const XYZVector P1  = A1 + u1 * sc;
   const XYZVector P2  = A2 + u2 * tc;
   const XYZVector mid = (P1 + P2) * 0.5;
   const double doca = (P1 - P2).R();
   return {doca, mid.X(), mid.Y(), mid.Z()};
}

ROOT::RDF::RNode makeSmeared(ROOT::RDF::RNode base, double sigma, unsigned lvl){

    auto smear = [sigma](unsigned salt) {
      return [sigma, salt](unsigned /*slot*/, ULong64_t entry,
                           const ROOT::RVec<float> &v) {
         const UInt_t seed =
            (UInt_t)((entry * 2654435761ULL + salt * 40503ULL + 12345ULL) & 0xffffffffULL);
         TRandom3 rng(seed);
         ROOT::RVec<float> out(v.size());
         for (std::size_t i = 0; i < v.size(); ++i)
            out[i] = v[i] * (1.f + (float)(sigma * rng.Gaus()));
         return out;
      };
   };

   return base
      .DefineSlotEntry("d_px_s", smear(lvl * 10u + 1u), {"d_px"})
      .DefineSlotEntry("d_py_s", smear(lvl * 10u + 2u), {"d_py"})
      .Define("d_E_s", recoEnergy, {"d_px_s", "d_py_s", "d_pz", "d_m"})
      .Define("p_sm",  llpFourMomentum, {"d_px_s", "d_py_s", "d_pz", "d_E_s"})
      // invariant mass / pT from the smeared momenta
      .Define("m_inv_s", [](const PxPyPzEVector &p) { return p.M();  }, {"p_sm"})
      .Define("pt_s",    [](const PxPyPzEVector &p) { return p.Pt(); }, {"p_sm"})
      // smeared daughter directions and their impact parameters
      .Define("u0_s", dir0, {"d_px_s", "d_py_s", "d_pz"})
      .Define("u1_s", dir1, {"d_px_s", "d_py_s", "d_pz"})
      // DOCA-based reconstructed vertex 
      .Define("doca_vec", docaVertex, {"A1", "u0_s", "A2", "u1_s"})
      .Define("reco_vtx_doca", [](const ROOT::RVec<double> &r) {
                                  return XYZVector(r[1], r[2], r[3]); }, {"doca_vec"})
      // vertex position residuals: transverse and longitudinal, signed
      // and centred at 0 (both -> 0 at sigma=0, since reco vertex == truth).
      .Define("dr_res", [](const XYZVector &v, const XYZVector &r) {
                           return r.Rho() - v.Rho(); }, {"vtx", "reco_vtx_doca"})
      .Define("dz_res", [](const XYZVector &v, const XYZVector &r) {
                           return r.Z() - v.Z(); }, {"vtx", "reco_vtx_doca"})
      .Define("nhat_s", [](const PxPyPzEVector &p) { return p.Vect().Unit(); }, {"p_sm"})
      .Define("ip_vtx", [](const XYZVector &rv, const XYZVector &n) {
                           return rv.Cross(n).R(); }, {"reco_vtx_doca", "nhat_s"});
}

void drawStack(TCanvas &c, int pad, const char *title,std::vector<ROOT::RDF::RResultPtr<TH1D>> &hv,const std::vector<std::string> &labels, bool logy){
   TVirtualPad *p = c.cd(pad);
   if (logy) p->SetLogy();

   auto *st = new THStack(Form("st%d", pad), title);
   auto *lg = new TLegend(0.58, 0.64, 0.88, 0.88);
   lg->SetBorderSize(0);
   lg->SetFillStyle(0);

   const int n    = (int)hv.size();
   const int ncol = TColor::GetNumberOfColors();
   for (int i = 0; i < n; ++i) {
      auto *h = (TH1D *)hv[i]->Clone(Form("c_%d_%d", pad, i));
      h->SetDirectory(nullptr);
      const int idx = (n <= 1) ? ncol / 2
                               : (int)std::lround((double)i / (n - 1) * (ncol - 1));
      h->SetLineColor(TColor::GetColorPalette(idx));
      h->SetLineWidth(2);
      st->Add(h, "hist");
      lg->AddEntry(h, labels[i].c_str(), "l");
   }
   st->Draw("nostack");
   lg->Draw();
}

void llp_simple_analysis(){

   ROOT::EnableImplicitMT();

   ROOT::RDataFrame df("Events", "./HNL_1.000e+00_1.000e+00_0.000e+00_1.000e+00_0.000e+00_data.root");
   const XYZVector PV(0., 0., 0.);

   auto base =
      df
         .Filter([](const ROOT::RVec<float> &px) { return px.size() == 2u; },
                 {"d_px"}, "exactly two daughters")
         .Define("vtx", [](float x, float y, float z) { return XYZVector(x, y, z); },
                 {"vtx_x", "vtx_y", "vtx_z"})
         .Define("Lxy", [](const XYZVector &v) { return v.Rho(); }, {"vtx"})
         .Define("L3D", [](const XYZVector &v) { return v.R();   }, {"vtx"})
         .Define("u0_true", dir0, {"d_px", "d_py", "d_pz"})
         .Define("u1_true", dir1, {"d_px", "d_py", "d_pz"})
         .Define("A1", [](const XYZVector &v, const XYZVector &u) {
                          return perigee(v, u); }, {"vtx", "u0_true"})
         .Define("A2", [](const XYZVector &v, const XYZVector &u) {
                          return perigee(v, u); }, {"vtx", "u1_true"});

   ROOT::RDF::RNode root = base;

   const std::vector<double>      sig = {0.00, 0.01, 0.05, 0.10};
   const std::vector<std::string> lab = {"true (#sigma=0)", "#sigma=1%",
                                         "#sigma=5%", "#sigma=10%"};
   const int NB = 60;

   auto n10  = makeSmeared(root, 0.10, 99);
   auto Mmin = n10.Min("m_inv_s");  auto Mmax = n10.Max("m_inv_s");
   auto Pmax = n10.Max("pt_s");
   auto DRlo = n10.Min("dr_res");   auto DRhi = n10.Max("dr_res");
   auto DZlo = n10.Min("dz_res");   auto DZhi = n10.Max("dz_res");
   auto IVmx = n10.Max("ip_vtx");

   const double mlo = *Mmin, mhi = *Mmax;
   const double mpad = 0.1 * (mhi - mlo) + 1e-6;
   const double MLO = mlo - mpad,  MHI = mhi + mpad;
   const double PHI = (*Pmax) * 1.1 + 1e-6;
   const double DRA = std::max(std::fabs(*DRlo), std::fabs(*DRhi)) * 1.1 + 1e-9;
   const double DZA = std::max(std::fabs(*DZlo), std::fabs(*DZhi)) * 1.1 + 1e-9;
   const double IVH = (*IVmx) * 1.1 + 1e-9;

   std::vector<ROOT::RDF::RResultPtr<TH1D>> HM, HP, HDR, HDZ, HIV;
   std::vector<ROOT::RDF::RResultPtr<double>> massMean, drStd, dzStd;

   for (std::size_t i = 0; i < sig.size(); ++i) {
      auto nd = makeSmeared(root, sig[i], (unsigned)i);
      HM .push_back(nd.Histo1D({Form("hM%zu", i),  "invariant mass;m_{inv} [GeV];events",                NB, MLO, MHI}, "m_inv_s"));
      HP .push_back(nd.Histo1D({Form("hP%zu", i),  "reconstructed p_{T};p_{T} [GeV];events",             NB, 0.,  PHI}, "pt_s"));
      HDR.push_back(nd.Histo1D({Form("hDR%zu", i), "transverse vertex residual;#Deltar = r_{reco}-r_{true} [m];events", NB, -DRA, DRA}, "dr_res"));
      HDZ.push_back(nd.Histo1D({Form("hDZ%zu", i), "longitudinal vertex residual;#Deltaz = z_{reco}-z_{true} [m];events", NB, -DZA, DZA}, "dz_res"));
      HIV.push_back(nd.Histo1D({Form("hIV%zu", i), "vertex IP;d^{3D}_{LLP} [m];events",                  NB, 0.,  IVH}, "ip_vtx"));
      massMean.push_back(nd.Mean("m_inv_s"));
      drStd.push_back(nd.StdDev("dr_res"));
      dzStd.push_back(nd.StdDev("dz_res"));
   }

   gStyle->SetOptStat(0);
   gStyle->SetPalette(kCool);

   TCanvas c("c", "LLP momentum-smearing study", 1500, 900);
   c.Divide(3, 2);
   drawStack(c, 1, "invariant mass;m_{inv} [GeV];events",                       HM,  lab, false);
   drawStack(c, 2, "reconstructed p_{T};p_{T} [GeV];events",                    HP,  lab, false);
   drawStack(c, 3, "transverse vertex residual;#Deltar [m];events",            HDR, lab, false);
   drawStack(c, 4, "longitudinal vertex residual;#Deltaz [m];events",          HDZ, lab, false);
   drawStack(c, 5, "vertex impact parameter;d^{3D}_{LLP} [m];events",          HIV, lab, true);
   c.SaveAs("llp_smearing.png");

}
