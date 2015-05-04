

int plot_NLO_corrections_B()
{
  const int num_vals = 19;
  double vals_mH[] = {
    100 ,
    125 ,
    150 ,
    200 ,
    250 ,
    300 ,
    320 ,
    340 ,
    360 ,
    380 ,
    400 ,
    450 ,
    500 ,
    550 ,
    600 ,
    700 ,
    800 ,
    900 ,
    1000
  };
  //////////////////////////////////////
  //////////////////////////////////////  
  double vals_LO[] = {
    // mu = 0.5*mH
    23.5263, 0.138981,
    16.9281, 0.114975,
    12.7646, 0.0960846,
    8.08539, 0.0746458,
    5.77303, 0.0626781,
    4.6548, 0.0570078,
    4.47128, 0.0549722,
    4.51413, 0.05414,
    5.07187, 0.0575082,
    5.0987, 0.0574883,
    4.75276, 0.0536175,
    3.51263, 0.0474603,
    2.41664, 0.0382267,
    1.63656, 0.0315347,
    1.10987, 0.0251373,
    0.526043, 0.0176245,
    0.261861, 0.0121219,
    0.136988, 0.00890458,
    0.074827, 0.00634679
  };
  //////////////////////////////////////
  //////////////////////////////////////  
  double vals_NLO_gg[] = {
    28.1233, 0.191763,
    19.7934, 0.158016,
    14.6978, 0.133611,
    9.04298, 0.100554,
    6.32715, 0.0838813,
    5.03356, 0.0728151,
    4.81211, 0.0717824,
    4.83563, 0.0728466,
    5.42096, 0.0760111,
    5.41542, 0.0748692,
    5.0307, 0.0714859,
    3.69665, 0.061497,
    2.52271, 0.0509918,
    1.7001, 0.0408935,
    1.14441, 0.033019,
    0.538549, 0.0228809,
    0.266949, 0.0161586,
    0.138746, 0.0116376,
    0.0754383, 0.00854131
  };
  //////////////////////////////////////
  //////////////////////////////////////  
  double vals_NLO_qg[] = {
    // mu = 0.5*mH
    2.48101, 0.116617,
    1.54696, 0.0984597,
    1.06522, 0.0834801,
    0.557342, 0.0658423,
    0.33309, 0.0553358,
    0.232833, 0.0493067,
    0.210152, 0.0501656,
    0.202903, 0.0486591,
    0.215392, 0.0521189,
    0.204025, 0.0526732,
    0.17835, 0.0507442,
    0.114558, 0.0433211,
    0.0665019, 0.0360875,
    0.040681, 0.0295398,
    0.0224442, 0.0249562,
    0.00620316, 0.0172435,
    0.00176927, 0.0123582,
    4.52996e-05, 0.00911115,
    -0.000400108, 0.0066855
  };
//////////////////////////////////////
//////////////////////////////////////
  double vals_NLO_qq[] = {
    // mu = 0.5*mH
    0.0343459, 0.0065151,
    0.024167, 0.00519118,
    0.0179278, 0.00383607,
    0.0112211, 0.00345084,
    0.00797839, 0.00296754,
    0.00642101, 0.00296468,
    0.00616298, 0.00257555,
    0.00622372, 0.00292873,
    0.0069806, 0.00276381,
    0.0070135, 0.00231969,
    0.00654733, 0.0022522,
    0.00483518, 0.00254331,
    0.00333124, 0.00212272,
    0.00224432, 0.001296,
    0.00151701, 0.00105803,
    0.000712084, 0.000720863,
    0.00035028, 0.000654966,
    0.000180437, 0.00042493,
    9.6742e-05, 0.000331039
  };
//////////////////////////////////////
//////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////
  graph_LO  = new TGraph(num_vals);
  graph_LO->SetLineColor(1);
  graph_LO->SetLineWidth(1);
  graph_LO->SetLineStyle(2);
  
  graph_LO_scale = new TGraphAsymmErrors(num_vals);
  graph_LO_scale->SetFillColor(8);
  graph_LO_scale->SetFillStyle(3002);
  graph_LO_scale->SetTitle("");
  graph_LO_scale->GetXaxis()->SetTitle("m_{H} [GeV]");
  graph_LO_scale->GetYaxis()->SetTitle("#sigma_{pp#rightarrow H+X} [pb]");
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////  
  graph_NLO_gg  = new TGraph(num_vals);
  graph_NLO_gg->SetLineColor(1);
  graph_NLO_gg->SetLineWidth(1);
  graph_NLO_gg->SetLineStyle(3);
  
  graph_NLO_gg_scale = new TGraphAsymmErrors(num_vals);
  graph_NLO_gg_scale->SetFillColor(9);
  graph_NLO_gg_scale->SetFillStyle(3002);
  graph_NLO_gg_scale->GetXaxis()->SetTitle("m_{H} [GeV]");
  graph_NLO_gg_scale->GetYaxis()->SetTitle("#sigma_{pp#rightarrow H+X} [pb]");
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  graph_NLO_qg  = new TGraph(num_vals,vals_mH,vals_NLO_qg);
  graph_NLO_qg->SetLineColor(1);
  graph_NLO_qg->SetLineWidth(1);
  graph_NLO_qg->SetLineStyle(4);
  ////////////////////////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////////////////////////////// 
  graph_NLO_qq  = new TGraph(num_vals,vals_mH,vals_NLO_qq);
  graph_NLO_qq->SetLineColor(1);
  graph_NLO_qq->SetLineWidth(1);
  graph_NLO_qq->SetLineStyle(5);
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  graph_NLO  = new TGraph(num_vals);
  graph_NLO->SetLineColor(1);
  graph_NLO->SetLineWidth(1);
  graph_NLO->SetLineStyle(1);
  
  graph_NLO_scale  = new TGraphAsymmErrors(num_vals);
  graph_NLO_scale->SetFillColor(9);
  graph_NLO_scale->SetFillStyle(3002);
  graph_NLO_scale->SetTitle("");
  graph_NLO_scale->GetXaxis()->SetTitle("m_{H} [GeV]");
  graph_NLO_scale->GetYaxis()->SetTitle("#sigma_{pp#rightarrow H+X} [pb]");
  ////////////////////////////////////////////////////////////////////////////////


  int M = 0; // 0: mu=0.5, 2: mu=1.0, 4: mu=2.0
  cout << setprecision(4) << endl;
  cout << setw(10) <<  "nlo"  << setw(10) <<  "mu=0.5" << setw(10) << "mu=2" << endl;   
  for (int i=0;i<num_vals;++i)
    {
      // LO
      graph_LO->SetPoint(i,vals_mH[i],vals_LO[M*num_vals+2*i]);

      // NLO - GG
      graph_NLO_gg->SetPoint(i,vals_mH[i],(vals_NLO_gg[M*num_vals+2*i]));

      // NLO - QG
      graph_NLO_qg->SetPoint(i,vals_mH[i],fabs(vals_NLO_qg[M*num_vals+2*i])*10.0);

      // NLO - QQ      
      graph_NLO_qq->SetPoint(i,vals_mH[i],(vals_NLO_qq[M*num_vals+2*i])*100.0);

      // complete NLO
      double nlo = vals_LO[M*num_vals+2*i]+vals_NLO_gg[M*num_vals+2*i]+(vals_NLO_qq[M*num_vals+2*i]+vals_NLO_qg[M*num_vals+2*i]);
      graph_NLO->SetPoint(i,vals_mH[i],nlo);
      cout << setw(15) << setprecision(10) << nlo << endl;
    }
  cout << endl;

  graph_NLO->SetTitle("");
  graph_NLO->GetXaxis()->SetTitle("m_{H} [GeV]");
  graph_NLO->GetYaxis()->SetTitle("#sigma_{pp#rightarrow H+X} [pb]");
  
  // graph_NLO_scale->SetFillColor(9);
  // graph_NLO_scale->SetFillStyle(3002);
  // graph_NLO_scale->GetYaxis()->SetRangeUser(0.03,70.0);
  // graph_NLO_scale->GetXaxis()->SetRangeUser(100,1000);
  
  c1 = new TCanvas("c1","Comparison of pp #rightarrow H+X cross-sections",0,0,1800,1200);
  c1->cd(1);
  gPad->SetLogy();
  gPad->SetTopMargin(0.025);	     
  gPad->SetBottomMargin(0.075);
  gPad->SetRightMargin(0.025);
  gPad->SetLeftMargin(0.075);
	
  leg1 = new TLegend(
		     0.6565 ,
		     0.8063 ,
		     0.9565 ,
		     0.9540 
		     );

  leg1->AddEntry(graph_LO    ,"LO","L");
  leg1->AddEntry(graph_NLO   ,"NLO (method B)","L");
  leg1->AddEntry(graph_NLO_gg,"K #times #delta_{NLO} (gg) ","L");
  leg1->AddEntry(graph_NLO_qg,"K #times #delta_{NLO} (qg)  #times 10 ","L");
  leg1->AddEntry(graph_NLO_qq,"K #times #delta_{NLO} (q#bar{q})  #times 100","L");
  leg1->AddEntry(graph_NLO,"","");

  //graph_LO_scale->Draw("A3");
  // graph_NLO_scale->Draw("A3");

  graph_NLO->GetXaxis()->SetRangeUser(100,1000);
  graph_NLO->Draw("AC");
  graph_LO->Draw("Csame");

  graph_NLO_gg->Draw("Csame");
  graph_NLO_qg->Draw("Csame");
  graph_NLO_qq->Draw("Csame");
      
  leg1->Draw("same");

  return 1;
}
