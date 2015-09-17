

void total_xsections()
{
  const int n_scen   = 3;
  const int n_scales = 3;

  
  // ================ 8 TeV ================ 
  double lo_qcd_8[n_scen*n_scales] = {
    // scenario 1: mu = 265, 530 , 132.5 
    121.59702301025391,
    94.315895080566406,
    160.01971435546875,
    // scenario 2: mu = 312.5, 625 , 156.25
    114.26799774169922,
    89.024284362792969,
    149.59356689453125,
    // scenario 3: mu = 325, 650, 162.5
    112.60903930664062,
    87.822532653808594,
    147.24418640136719
  };
  double delta_nlow_qcd_8[n_scen*n_scales] = {
    // scenario 1: mu = 265, 530 , 132.5 
    70.259955644607544,
    72.170376062393188,
    57.502393722534180,
    // scenario 2: mu = 312.5, 625, 156.25
    71.419341802597046,
    71.737211942672729,
    62.003070354461670,
    // scenario 3: mu = 325, 650, 162.5
    71.613635063171387,
    71.599982380867004,
    62.913398981094360
  };
  double delta_nlo_phi_8[n_scen*n_scales] = {
    // scenario 1: mu = 265, 530 , 132.5
    4.718085,
    4.2541863,
    5.0965632,
    // scenario 2: mu = 312.5, 625, 156.25
    2.46246962,
    2.25803274,
    2.5761242,
    // scenario 3: mu = 325, 650, 162.5
    3.4603297,
    3.0912555,
    3.7933129
    // // without eff. K-factor
    // // scenario 1: mu = 265, 530 , 132.5
    // 2.3545837,
    // 2.0626584,
    // 2.6682418,
    // // scenario 2: mu = 312.5, 625, 156.25
    // 1.53692138,
    // 1.37847951,
    // 1.6666372,
    // // scenario 3: mu = 325, 650, 162.5
    // 1.8945024,
    // 1.6515514,
    // 2.1529071       
  };

  
  // ================ 13 TeV ================ 
  double lo_qcd_13[n_scen*n_scales] = {
    // scenario 1: mu = 265, 530 , 132.5
    404.50357055664062,
    322.00964355468750,
    516.58612060546875,
    // scenario 2: mu = 312.5, 625 , 156.25
    382.60446166992188,
    305.67599487304688,
    486.59451293945312,
    // scenario 3: mu = 325, 650, 162.5
    377.62216186523438,
    301.94985961914062,
    479.79479980468750
  };
  double delta_nlow_qcd_13[n_scen*n_scales] = {
    // scenario 1: mu = 265, 530 , 132.5
    233.35254096984863,
    238.21308898925781,
    201.90867900848389,
    // scenario 2: mu = 312.5, 625, 156.25
    236.30378627777100,
    237.04260063171387,
    212.84232044219971,
    // scenario 3: mu = 325, 650, 162.5
    236.81546401977539,
    236.67661857604980,
    215.03704643249512
  };
  double delta_nlo_phi_13[n_scen*n_scales] = {
    // scenario 1: mu = 265, 530 , 132.5
    15.208241,
    13.825387,
    16.545652,
    // scenario 2: mu = 312.5, 625, 156.25
    8.1784890,
    7.5939478,
    8.567177,
    // scenario 3: mu = 325, 650, 162.5
    11.256460,
    10.218260,
    12.337870    
  };


  const double eps_ll = 0.22;
  const double br_ll  = 4.0/81.0;
  const double eps_lj = 0.12;
  const double br_lj  = 24.0/81.0;
  double L = 20.0*1000.0; // integrated lumi of 8TeV run in [pb^-1]


  cout << endl << " ========== sqrt(s) = 8 TeV ========== " << endl;
  for (int i = 0; i<n_scen; ++i)
    {
      cout << endl << " Scenario " << i+1 << ":" << endl;
      cout << setw(25) << " " << setw(15) << " total xsection"  << setw(15) << "mu_hi" << setw(15) << " mu_lo " << endl;
      double val = lo_qcd_8[n_scen*i+0];
      double vhi = lo_qcd_8[n_scen*i+1];
      double vlo = lo_qcd_8[n_scen*i+2];
      cout << setw(25) << " QCD LO "   << setw(15) << val << setw(15) << vhi - val << setw(15) << vlo - val << endl;
      val += delta_nlow_qcd_8[n_scen*i+0];
      vhi += delta_nlow_qcd_8[n_scen*i+1];
      vlo += delta_nlow_qcd_8[n_scen*i+2];
      cout << setw(25) << " QCD NLOW "   << setw(15) << val << setw(15) << vhi - val << setw(15) << vlo - val << endl;
      val += delta_nlo_phi_8[n_scen*i+0];
      vhi += delta_nlo_phi_8[n_scen*i+1];
      vlo += delta_nlo_phi_8[n_scen*i+2];
      cout << setw(25) << " QCD NLOW + NLO PHI "   << setw(15) << val << setw(15) << vhi - val << setw(15) << vlo - val << endl;

      cout << setw(25) << "                  R "   << val/(lo_qcd_8[n_scen*i+0]+delta_nlow_qcd_8[n_scen*i+0]) << endl;
      
      // additional events from phi production assuming efficiencies, luminosity given above
      // dilepton channel
      double N   = delta_nlo_phi_8[n_scen*i+0]*L*eps_ll*br_ll;
      double Nhi = delta_nlo_phi_8[n_scen*i+1]*L*eps_ll*br_ll;
      double Nlo = delta_nlo_phi_8[n_scen*i+2]*L*eps_ll*br_ll;
      cout << setw(25) << " delta N NLO PHI (ll) "   << setw(15) << N << setw(15) << Nhi - N << setw(15) << Nlo - N << endl;
      // lepton+jets channel
      N   = delta_nlo_phi_8[n_scen*i+0]*L*eps_lj*br_lj;
      Nhi = delta_nlo_phi_8[n_scen*i+1]*L*eps_lj*br_lj;
      Nlo = delta_nlo_phi_8[n_scen*i+2]*L*eps_lj*br_lj;
      cout << setw(25) << " delta N NLO PHI (lj) "   << setw(15) << N << setw(15) << Nhi - N << setw(15) << Nlo - N << endl;
      cout << endl;
    }
  
  cout << endl << " ========== sqrt(s) = 13 TeV ========== " << endl;
  for (int i = 0; i<n_scen; ++i)
    {
      cout << endl << " Scenario " << i+1 << ":" << endl;
      cout << setw(25) << " " << setw(15) << " total xsection"  << setw(15) << "mu_hi" << setw(15) << " mu_lo " << endl;
      double val = lo_qcd_13[n_scen*i+0];
      double vhi = lo_qcd_13[n_scen*i+1];
      double vlo = lo_qcd_13[n_scen*i+2];
      cout << setw(25) << " QCD LO "   << setw(15) << val << setw(15) << vhi - val << setw(15) << vlo - val << endl;
      val += delta_nlow_qcd_13[n_scen*i+0];
      vhi += delta_nlow_qcd_13[n_scen*i+1];
      vlo += delta_nlow_qcd_13[n_scen*i+2];
      cout << setw(25) << " QCD NLOW "   << setw(15) << val << setw(15) << vhi - val << setw(15) << vlo - val << endl;
      val += delta_nlo_phi_13[n_scen*i+0];
      vhi += delta_nlo_phi_13[n_scen*i+1];
      vlo += delta_nlo_phi_13[n_scen*i+2];
      cout << setw(25) << " QCD NLOW + NLO PHI "   << setw(15) << val << setw(15) << vhi - val << setw(15) << vlo - val << endl;
      cout << endl;

      cout << setw(25) << "                  R "   << val/(lo_qcd_13[n_scen*i+0]+delta_nlow_qcd_13[n_scen*i+0]) << endl;
      
      // assume lumi 100 fm^-1 for 13 TeV run
      L = 100.0*1000.0;
      // additional events from phi production assuming efficiencies, luminosity given above
      // dilepton channel
      double N   = delta_nlo_phi_13[n_scen*i+0]*L*eps_ll*br_ll;
      double Nhi = delta_nlo_phi_13[n_scen*i+1]*L*eps_ll*br_ll;
      double Nlo = delta_nlo_phi_13[n_scen*i+2]*L*eps_ll*br_ll;
      cout << setw(25) << " delta N NLO PHI (ll) "   << setw(15) << N << setw(15) << Nhi - N << setw(15) << Nlo - N << endl;
      // lepton+jets channel
      N   = delta_nlo_phi_13[n_scen*i+0]*L*eps_lj*br_lj;
      Nhi = delta_nlo_phi_13[n_scen*i+1]*L*eps_lj*br_lj;
      Nlo = delta_nlo_phi_13[n_scen*i+2]*L*eps_lj*br_lj;
      cout << setw(25) << " delta N NLO PHI (lj) "   << setw(15) << N << setw(15) << Nhi - N << setw(15) << Nlo - N << endl;
      cout << endl;
    }
  
}
