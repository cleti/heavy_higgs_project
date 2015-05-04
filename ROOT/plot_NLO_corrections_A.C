

int plot_NLO_corrections_A()
{
  const int num_vals = 18;
  double vals_mH[] = {
    100 ,
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
37.4735 ,
20.3146 ,
12.8671 ,
9.19128 ,
7.41084 ,
7.11907 ,
7.19512 ,
8.0774 ,
8.11055 ,
7.57484 ,
5.58801 ,
3.8448 ,
2.6026 ,
1.76551 ,
0.836839 ,
0.416882 ,
0.217955 ,
0.118973 ,
    // mu=0.25 //////////
42.0983 ,
23.8326 ,
15.4917 ,
11.2634 ,
9.19282 ,
8.86754 ,
8.99957 ,
10.1428 ,
10.2216 ,
9.57524 ,
7.11185 ,
4.92291 ,
3.34996 ,
2.28252 ,
1.09011 ,
0.546522,
0.287284,
0.157542,
// mu=1 ////////////    
33.0326 ,
17.3187 ,
10.745 ,
7.56603 ,
6.03535 ,
5.77741 ,
5.82021 ,
6.51356 ,
6.52273 ,
6.07572 ,
4.4558 ,
3.051 ,
2.05611 ,
1.38922 ,
0.654042,
0.323976,
0.168548,
0.091597
// mu=2 ////////////
// 29.0249 ,
// 14.8032 ,
// 9.03244 ,
// 6.28542 ,
// 4.9684 ,
// 4.74123 ,
// 4.76269 ,
// 5.31614 ,
// 5.31077 ,
// 4.93538 ,
// 3.60184 ,
// 2.4551 ,
// 1.64821 ,
// 1.10967 ,
// 0.519336,
// 0.255932,
// 0.132562,
// 0.071756
////////////////////
  };
  double vals_NLO_gg[] = {
    // mu=1 /////////////
14.0416,
7.77374,
4.70572,
3.04993,
2.07909,
1.80434,
1.57352,
1.37957,
1.21494,
1.07449,
0.80341,
0.61284,
0.47526,
0.37396,
0.23981,
0.15967,
0.10959,
0.07709,
// mu=0.5 //////////
-9.75698,
-2.56203,
-0.73527 ,
-0.15537 ,
0.042621 ,
0.077963 ,
0.100083,
0.113045,
0.119673,
0.12206 ,
0.117256,
0.105658,
0.092546 ,
0.079931 ,
0.058647 ,
0.042926 ,
0.031594 ,
0.023483,
// mu=2 ////////////
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
////////////////////
  };
  double vals_NLO_qg[] = {
// mu=1 /////////////
-1.50575,
-2.03647,
-1.64306,
-1.26424,
-0.97931,
-0.88563,
-0.79890,
-0.72239,
-0.65882,
-0.60109,
-0.48108,
-0.38975 ,
-0.32013 ,
-0.26402 ,
-0.18367 ,
-0.13184 ,
-0.09646,
-0.07203,
    // mu=0.5 //////////
25.409 ,
9.87016,
4.69347,
2.51378,
1.45095,
1.18706,
0.96502,
0.79777,
0.66018,
0.55591,
0.35129,
0.23226,
0.1502 ,
0.09965,
0.04090,
0.01398,
0.00115 ,
-0.0052,
    // mu=2 ////////////
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
    ////////////////////
  };
  double vals_NLO_qq[] = {
    // mu=1 /////////////
0.270247 ,
0.127668 ,
0.071877 ,
0.0447321,
0.0297195,
0.0255854,
0.0221677,
0.0193112,
0.0169147,
0.0148812,
0.0109993,
0.0083049,
0.0063810,
0.0049732,
0.0031279,
0.0020408,
0.0013699,
0.0009403,
    // mu=0.5 //////////
0.36009 ,
0.171303,
0.096870,
0.060481 ,
0.040275 ,
0.034690 ,
0.030076 ,
0.026219 ,
0.022979 ,
0.020228,
0.014970 ,
0.011315 ,
0.008703,
0.006788,
0.004278,
0.002795,
0.001879,
0.001292,
    // mu=2 ////////////
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
    ////////////////////
  };
  double vals_NLO[] = {
        // mu=0.5 //////////
53.5143,
27.7944,
16.93 ,
11.6029,
8.94029,
8.40973,
8.29507,
9.01455,
8.91988,
8.26696,
6.0771 ,
4.19162,
2.85587,
1.95219,
0.94107,
0.47594,
0.25233,
0.13874,
    // mu=0.25 /////////////
65.5188,
34.4492,
21.2425,
14.7101,
11.392 ,
10.7357,
10.5965,
11.516 ,
11.4055,
10.6021,
7.84651,
5.45914,
3.74693,
2.58234,
1.26497,
0.65389,
0.35381,
0.19952,
    // mu=1 /////////////
45.8852,
23.2104,
13.878 ,
9.39214,
7.17065,
6.72389,
6.62004,
7.18617,
7.09239,
6.55731,
4.78739,
3.28139,
2.21905,
1.50285,
0.71293,
0.35399,
0.18323,
0.09742
    // mu=2 ////////////
// 40.1579,
// 19.7754,
// 11.5736,
// 7.7326 ,
// 5.83346,
// 5.45961,
// 5.35852,
// 5.81726,
// 5.73007,
// 5.28537,
// 3.8344 ,
// 2.60866,
// 1.75212,
// 1.17979,
// 0.54905,
// 0.26703,
// 0.13434,
// 0.06916
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


  
  cout << setprecision(4) << endl;
  cout << setw(10) <<  "nlo"  << setw(10) <<  "mu=0.5" << setw(10) << "mu=2" << endl;   
  for (int i=0;i<num_vals;++i)
    {
      // LO
      graph_LO->SetPoint(i,vals_mH[i],vals_LO[i]);
      graph_LO_scale->SetPoint(i,vals_mH[i],vals_LO[i]);
      graph_LO_scale->SetPointError(i,0,0,vals_LO[i]-vals_LO[2*num_vals+i],vals_LO[num_vals+i]-vals_LO[i]);

      // NLO - GG
      graph_NLO_gg->SetPoint(i,vals_mH[i],fabs(vals_NLO_gg[num_vals+i]));

      // NLO - QG
      graph_NLO_qg->SetPoint(i,vals_mH[i],(vals_NLO_qg[num_vals+i]));

      // NLO - QQ      
      graph_NLO_qq->SetPoint(i,vals_mH[i],(vals_NLO_qq[num_vals+i])*10.0);

      // complete NLO
      double nlo = vals_NLO[i];//vals_LO[i]+vals_NLO_gg[i]+0.1*(vals_NLO_qq[i]+vals_NLO_qg[i]);
      graph_NLO->SetPoint(i,
			  vals_mH[i],
			  nlo
			  );
      graph_NLO_scale->SetPoint(i,
				vals_mH[i],
				nlo
				);
      double s_lo = vals_NLO[num_vals+i];//vals_LO[num_vals+i]+vals_NLO_gg[num_vals+i];
      double s_hi = vals_NLO[2*num_vals+i];//vals_NLO_gg[2*num_vals+i]+vals_LO[2*num_vals+i];


      cout << setw(10) <<  nlo    << setw(10) <<  s_lo     << setw(10) <<  s_hi  << endl;
      
      double nlo_lo = s_lo;
      double nlo_hi = s_hi;
      if (s_lo<nlo && s_hi>nlo)
      	{
      	  nlo_lo = nlo-s_lo;
      	  nlo_hi = s_hi-nlo;
      	}
      else (s_lo>nlo && s_hi<nlo)
      	{
      	  nlo_lo =nlo-s_hi;
      	  nlo_hi = s_lo-nlo;
      	}
      // else
      // 	{
      // 	  cout << endl << " invalid scale dependence! " << endl;
      // 	  exit(1);
      // 	}
      graph_NLO_scale->SetPointError(i,0,0,
				     nlo_lo, // error x-low
				     nlo_hi  // error x-high
				     );
    }
  cout << endl;


  graph_NLO_scale->SetFillColor(9);
  graph_NLO_scale->SetFillStyle(3002);
  graph_NLO_scale->SetTitle("");
  graph_NLO_scale->GetXaxis()->SetTitle("m_{H} [GeV]");
  graph_NLO_scale->GetYaxis()->SetTitle("#sigma_{pp#rightarrow H+X} [pb]");
  graph_NLO_scale->GetYaxis()->SetRangeUser(0.03,70.0);
  graph_NLO_scale->GetXaxis()->SetRangeUser(100,1000);
    
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
  leg1->AddEntry(graph_NLO   ,"NLO (method A)","L");
  leg1->AddEntry(graph_NLO_gg,"| #delta_{NLO} (gg) |","L");
  leg1->AddEntry(graph_NLO_qg,"#delta_{NLO} (qg)","L");
  leg1->AddEntry(graph_NLO_qq,"#delta_{NLO} (q#bar{q}) #times 10","L");

  //graph_LO_scale->Draw("A3");
  graph_NLO_scale->Draw("A3");

  
  graph_LO->Draw("Csame");
  graph_NLO->Draw("Csame");
  graph_NLO_gg->Draw("Csame");
  graph_NLO_qg->Draw("Csame");
  graph_NLO_qq->Draw("Csame");
      
  leg1->Draw("same");

  return 1;
}
