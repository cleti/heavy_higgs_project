

int plot_NLO_corrections_A()
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
  double vals_LO[] = {
    // mu=0.5 //////////
23.5595, 0.14025,
16.9304, 0.116149,
12.7476, 0.0988597,
8.08533, 0.0751173,
5.77408, 0.0625984,
4.6581, 0.0565473,
4.46915, 0.053997,
4.52117, 0.0544685,
5.07072, 0.057133,
5.10293, 0.0582127,
4.75696, 0.0550162,
3.51431, 0.0466831,
2.41666, 0.0375046,
1.6379, 0.0306601,
1.10852, 0.0247366,
0.525568, 0.016889,
0.261853, 0.012077,
0.136924, 0.00871408,
0.0747957, 0.00635881
    ////////////////////
  };
  double vals_NLO_gg[] = {
    // mu=0.5 /////////////
29.7901, 0.197721,
19.9811, 0.153935,
14.1229, 0.130467,
7.87722, 0.0933073,
4.83588, 0.072974,
3.17346, 0.0590687,
2.72678, 0.0528834,
2.35468, 0.0503759,
2.04501, 0.0468736,
1.7829, 0.0436445,
1.56379, 0.0406828,
1.14714, 0.0346467,
0.861576, 0.0294562,
0.662022, 0.0253587,
0.515289, 0.0223819,
0.324731, 0.0178944,
0.213494, 0.0142573,
0.145037, 0.0116395,
0.101156, 0.00985452
    ////////////////////
  };
  double vals_NLO_qg[] = {
    // mu=0.5 //////////
2.6433, 0.117306,
1.59349, 0.0958385,
1.01339, 0.0804794,
0.476835, 0.061565,
0.260868, 0.0480207,
0.145607, 0.0397546,
0.118845, 0.0373139,
0.0989376, 0.0339203,
0.0804567, 0.0316876,
0.0673107, 0.0296641,
0.0558977, 0.0294064,
0.0361949, 0.0245614,
0.0233391, 0.0211623,
0.0145086, 0.0193697,
0.0100373, 0.0167665,
0.00382183, 0.0143816,
0.00114831, 0.0110823,
-0.000130658, 0.00924216,
-0.000467171, 0.00785675
    ////////////////////
  };
  double vals_NLO_qq[] = {
    // mu=0.5 //////////
0.036402, 0.0063341,
0.0243833, 0.00529265,
0.0173107, 0.00481133,
0.00977413, 0.00285643,
0.00609285, 0.00220486,
0.00405654, 0.00217788,
0.0034911, 0.00199368,
0.00302697, 0.00152279,
0.0026363, 0.00187778,
0.00230559, 0.00153254,
0.00203601, 0.00145984,
0.00150277, 0.00146958,
0.00113223, 0.00130485,
0.000873668, 0.00081709,
0.000680861, 0.000711511,
0.000428969, 0.000561622,
0.000280121, 0.000550684,
0.000188243, 0.000504688,
0.000129574, 0.000403794
    ////////////////////
  };


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
  leg1->AddEntry(graph_NLO   ,"NLO (eff.)","L");
  leg1->AddEntry(graph_NLO_gg," #delta_{NLO} (gg) ","L");
  leg1->AddEntry(graph_NLO_qg," #delta_{NLO} (qg)  #times 10 ","L");
  leg1->AddEntry(graph_NLO_qq," #delta_{NLO} (q#bar{q})  #times 100","L");
  //leg1->AddEntry(graph_NLO,"","");

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
