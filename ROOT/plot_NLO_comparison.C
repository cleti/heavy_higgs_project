

#include "TMath.h"

int plot_NLO_comparison()
{
  const int num_vals = 19;
  double vals_mH[] = {
100,
125,
150,
200,
250,
300,
320,
340,
360,
380,
400,
450,
500,
550,
600,
700,
800,
900,
1000
  };


  
  double vals_dFG[] = {
    73.38 ,// 100
    49.97 ,// 125    
    36.38 ,
    20.18 ,// 200
    13.92 ,
    11.07 ,// 300
    10.60 ,
    10.69 ,
    11.81 ,
    11.66 ,
    10.76 ,// 400
    7.80 ,
    5.31 ,// 500
    3.54 ,
    2.37 ,// 600
    1.10 ,// 700
    0.539 ,// 800
    0.279 ,// 900
    0.151  // 1000
  };
  double vals_ABPS[] = {
    76.07 ,// 100
    51.47 ,// 125    
    37.43 ,
    20.64 ,// 200
    14.23 ,
    11.28 ,// 300
    10.81 ,
    11.00 ,
    12.30 ,
    12.01 ,
    10.98 ,// 400
    7.81 ,
    5.24 ,// 500
    3.48 ,
    2.32 ,// 600
    1.07 ,// 700
    0.525 ,// 800
    0.270 ,// 900
    0.146  // 1000
  };


  double vals_NLO_A[] = {
        // mu=0.5 //////////
56.0927, 0.261348,
38.6029, 0.207474,
27.917, 0.172397,
16.4494, 0.13109,
10.9076, 0.104772,
8.00263, 0.0879754,
7.32475, 0.0815586,
6.96302, 0.0782515,
7.19512, 0.0787171,
6.95122, 0.0775717,
6.38794, 0.0721663,
4.69638, 0.0622772,
3.30517, 0.0514557,
2.31486, 0.0433776,
1.633, 0.0372043,
0.85489, 0.0282657,
0.477095, 0.0211161,
0.28199, 0.0168041,
0.175409, 0.013945,
        // mu=0.25 /////////////
50.9357, 0.245922,
35.9863, 0.20525,
26.6615, 0.173665,
16.2743, 0.132086,
11.1184, 0.107991,
8.43122, 0.0901604,
7.83923, 0.0851644,
7.61101, 0.0831326,
8.07374, 0.0839911,
7.90724, 0.0793804,
7.3247, 0.0759062,
5.42445, 0.0657698,
3.81441, 0.0559995,
2.65469, 0.0457416,
1.86243, 0.0395563,
0.953359, 0.0294855,
0.519111, 0.0226278,
0.299424, 0.0173962,
0.18162, 0.0145131,
        // mu=1 /////////////
57.0784, 0.261795,
38.4736, 0.212182,
27.4359, 0.175874,
15.7367, 0.126832,
10.2033, 0.100506,
7.32788, 0.0847025,
6.63395, 0.0799506,
6.23436, 0.0758508,
6.32926, 0.0742572,
6.04428, 0.0703387,
5.52149, 0.0687406,
4.03134, 0.0561404,
2.84037, 0.0485205,
1.9936, 0.0402929,
1.41652, 0.0342772,
0.750003, 0.0255001,
0.423986, 0.0202694,
0.254391, 0.0159937,
0.160187, 0.013342
////////////////////
  };
  
  double vals_NLO_B[] = {
        // mu=0.5 //////////
54.2201, 0.0904763,
38.3308, 0.0748134,
28.4887, 0.0637859,
17.69, 0.0495137,
12.4556, 0.0415652,
9.93146, 0.0369377,
9.50261, 0.0363223,
9.56711, 0.0360376,
10.7022, 0.0381522,
10.7205, 0.0381434,
9.98485, 0.0369514,
7.32117, 0.0315757,
5.01033, 0.0262592,
3.37463, 0.0214178,
2.27896, 0.0177098,
1.07229, 0.0121616,
0.530873, 0.0086886,
0.276036, 0.00621192,
0.149978, 0.00464185,    
// 54.1948, 0.246563,
// 38.3424, 0.208109,
// 28.4517, 0.178167,
// 17.7003, 0.137831,
// 12.466, 0.114418,
// 9.92046, 0.102764,
// 9.50947, 0.102717,
// 9.57937, 0.0978227,
// 10.7031, 0.102923,
// 10.7329, 0.104088,
// 9.98647, 0.0972255,
// 7.31499, 0.0867369,
// 5.00986, 0.072957,
// 3.37745, 0.0582054,
// 2.27991, 0.0468277,
// 1.07226, 0.0332855,
// 0.530876, 0.0230536,
// 0.276029, 0.0165435,
// 0.149964, 0.0125778,
        // mu=0.25 //////////
49.3993, 0.0878093,
35.7473, 0.073906,
27.0585, 0.0640479,
17.2318, 0.0500637,
12.3458, 0.0422136,
9.96436, 0.0378336,
9.57222, 0.0371672,
9.67858, 0.037091,
10.8735, 0.0393718,
10.9298, 0.0398647,
10.209, 0.0386041,
7.53717, 0.0333519,
5.19045, 0.0271687,
3.51507, 0.0226384,
2.38485, 0.0188045,
1.13074, 0.0129687,
0.563396, 0.00911744,
0.294582, 0.00664473,
0.160828, 0.00495837,
// 49.4305, 0.244441,
// 35.7194, 0.200769,
// 27.0518, 0.175208,
// 17.2277, 0.135164,
// 12.3443, 0.112066,
// 9.98405, 0.101634,
// 9.58295, 0.100093,
// 9.67746, 0.100661,
// 10.8849, 0.10763,
// 10.9448, 0.108306,
// 10.211, 0.1021,
// 7.53683, 0.0873907,
// 5.19559, 0.0722312,
// 3.51404, 0.0598913,
// 2.38501, 0.0491233,
// 1.13031, 0.0346162,
// 0.563349, 0.0241342,
// 0.294624, 0.0170053,
// 0.160877, 0.012863,
        // mu=1 /////////////
55.2346, 0.0892731,
38.3402, 0.0742535,
28.1357, 0.0626167,
17.127, 0.0482533,
11.8905, 0.0404941,
9.38372, 0.0355531,
8.94607, 0.0349734,
8.9792, 0.0345245,
10.0155, 0.0365991,
10.0049, 0.0366973,
9.29463, 0.0351859,
6.77708, 0.0301125,
4.61532, 0.02505,
3.0953, 0.0204381,
2.08239, 0.0167798,
0.97343, 0.0116087,
0.479272, 0.00812208,
0.248043, 0.00590364,
0.134208, 0.00436247
// 55.2885, 0.25267,
// 38.3024, 0.209036,
// 28.1438, 0.175186,
// 17.1055, 0.135011,
// 11.8873, 0.111726,
// 9.38696, 0.098831,
// 8.95186, 0.0971281,
// 8.98173, 0.0941189,
// 10.0176, 0.101111,
// 10.0077, 0.10065,
// 9.30518, 0.0966328,
// 6.78018, 0.0813513,
// 4.61506, 0.0693326,
// 3.09406, 0.0553635,
// 2.0815, 0.045501,
// 0.974787, 0.0312068,
// 0.479211, 0.0219062,
// 0.248126, 0.0155687,
// 0.134238, 0.0115568
////////////////////
  };


  ////////////////////////////////////////////////////////////////////////////////
  graph_dFG  = new TGraph(num_vals);
  graph_dFG->SetLineColor(1);
  graph_dFG->SetLineWidth(1);
  graph_dFG->SetLineStyle(3);
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////  
  graph_ABPS  = new TGraph(num_vals);
  graph_ABPS->SetLineColor(1);
  graph_ABPS->SetLineWidth(1);
  graph_ABPS->SetLineStyle(4);
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  graph_NLO_A  = new TGraph(num_vals);
  graph_NLO_A->SetLineColor(1);
  graph_NLO_A->SetLineWidth(1);
  graph_NLO_A->SetLineStyle(2);
  
  graph_NLO_A_scale  = new TGraphAsymmErrors(num_vals);
  graph_NLO_A_scale->SetFillColor(9);
  graph_NLO_A_scale->SetFillStyle(3002);
  graph_NLO_A_scale->SetTitle("");
  graph_NLO_A_scale->GetXaxis()->SetTitle("m_{H} [GeV]");
  graph_NLO_A_scale->GetYaxis()->SetTitle("#sigma_{pp#rightarrow H+X} [pb]");
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  graph_NLO_B  = new TGraph(num_vals);
  graph_NLO_B->SetLineColor(1);
  graph_NLO_B->SetLineWidth(1);
  graph_NLO_B->SetLineStyle(1);
  
  graph_NLO_B_scale  = new TGraphAsymmErrors(num_vals);
  graph_NLO_B_scale->SetFillColor(9);
  graph_NLO_B_scale->SetFillStyle(3002);
  graph_NLO_B_scale->SetTitle("");
  graph_NLO_B_scale->GetXaxis()->SetTitle("m_{H} [GeV]");
  graph_NLO_B_scale->GetYaxis()->SetTitle("#sigma_{pp#rightarrow H+X} [pb]");
  ////////////////////////////////////////////////////////////////////////////////



  
  cout << setprecision(4) << endl;
  cout << setw(10) <<  "nlo"  << setw(10) <<  "mu=0.5" << setw(10) << "mu=2" << endl;   
  for (int i=0;i<num_vals;++i)
    {
      // results from literature
      graph_dFG->SetPoint(i,vals_mH[i],vals_dFG[i]);
      graph_ABPS->SetPoint(i,vals_mH[i],vals_ABPS[i]);

      // complete NLO (method A)
      double nlo_A = vals_NLO_A[2*i];//vals_LO[i]+vals_NLO_gg[i]+0.1*(vals_NLO_qq[i]+vals_NLO_qg[i]);
      graph_NLO_A->SetPoint(i,
			  vals_mH[i],
			  nlo_A
			  );
      graph_NLO_A_scale->SetPoint(i,
				vals_mH[i],
				nlo_A
				);
      double s_lo_A = vals_NLO_A[2*(num_vals+i)];//vals_LO[num_vals+i]+vals_NLO_gg[num_vals+i];
      double s_hi_A = vals_NLO_A[2*(2*num_vals+i)];//vals_NLO_gg[2*num_vals+i]+vals_LO[2*num_vals+i];


      // cout << setw(10) <<  nlo_A    << setw(10) <<  s_lo_A     << setw(10) <<  s_hi_A  << endl;
      
      double nlo_lo_A = s_lo_A;
      double nlo_hi_A = s_hi_A;
      if (s_lo_A<nlo_A && s_hi_A>nlo_A)
      	{
      	  nlo_lo_A = nlo_A-s_lo_A;
      	  nlo_hi_A = s_hi_A-nlo_A;
      	}
      else (s_lo_A>nlo_A && s_hi_A<nlo_A)
      	{
      	  nlo_lo_A = nlo_A-s_hi_A;
      	  nlo_hi_A = s_lo_A-nlo_A;
      	}
      // else
      // 	{
      // 	  cout << endl << " invalid scale dependence! " << endl;
      // 	  exit(1);
      // 	}
      graph_NLO_A_scale->SetPointError(i,0,0,
				     nlo_lo_A, // error x-low
				     nlo_hi_A  // error x-high
				     );

      // complete NLO (method B)
      double nlo_B = vals_NLO_B[2*i];//vals_LO[i]+vals_NLO_gg[i]+0.1*(vals_NLO_qq[i]+vals_NLO_qg[i]);
      graph_NLO_B->SetPoint(i,
			  vals_mH[i],
			  nlo_B
			  );
      graph_NLO_B_scale->SetPoint(i,
				vals_mH[i],
				nlo_B
				);

      // scale dependence
      double s_lo = vals_NLO_B[2*(num_vals+i)];
      double s_hi = vals_NLO_B[2*(2*num_vals+i)];
      double nlo_lo_B,nlo_hi_B;

      cout << setw(4) << i << setw(10) <<  nlo_B    << setw(10) <<  s_lo     << setw(10) <<  s_hi;

      
      if (s_lo<nlo_B && s_hi>nlo_B)
      	{
      	  nlo_lo_B = nlo_B-s_lo;
      	  nlo_hi_B = s_hi-nlo_B;
      	}
      else if (s_lo>nlo_B && s_hi<nlo_B)
      	{
      	  nlo_lo_B = nlo_B-s_hi;
      	  nlo_hi_B = s_lo-nlo_B;
      	}
      else if (s_lo>nlo_B && s_hi>nlo_B)
      	{
      	  nlo_lo_B = 0.0;
      	  nlo_hi_B = TMath::Max(s_lo,s_hi)-nlo_B;
      	}
      else if (s_lo<nlo_B && s_hi<nlo_B)
      	{
      	  nlo_lo_B = nlo_B-TMath::Min(s_lo,s_hi);
      	  nlo_hi_B = 0.0;
      	}
      else
	{
      	  cout << endl << " invalid scale dependence! " << endl;
      	  exit(1);
	}
      cout << setw(10) <<  nlo_lo_B/nlo_B    << setw(10) <<  nlo_hi_B/nlo_B  << endl;
      
      graph_NLO_B_scale->SetPointError(i,0,0,
				       nlo_lo_B, // error x-low
				       nlo_hi_B  // error x-high
				       );
      
    }
  cout << endl;


  graph_NLO_A_scale->SetFillColor(9);
  graph_NLO_A_scale->SetFillStyle(3002);
  graph_NLO_A_scale->SetTitle("");
  graph_NLO_A_scale->GetXaxis()->SetTitle("m_{H} [GeV]");
  graph_NLO_A_scale->GetYaxis()->SetTitle("#sigma_{pp#rightarrow H+X} [pb]");
  graph_NLO_A_scale->GetYaxis()->SetRangeUser(0.1,80.0);
  graph_NLO_A_scale->GetXaxis()->SetRangeUser(100,1000);
  
  graph_NLO_B_scale->SetFillColor(9);
  graph_NLO_B_scale->SetFillStyle(3002);
  graph_NLO_B_scale->SetTitle("");
  graph_NLO_B_scale->GetXaxis()->SetTitle("m_{H} [GeV]");
  graph_NLO_B_scale->GetYaxis()->SetTitle("#sigma_{pp#rightarrow H+X} [pb]");
  graph_NLO_B_scale->GetYaxis()->SetRangeUser(0.1,80.0);
  graph_NLO_B_scale->GetXaxis()->SetRangeUser(100,1000);
  
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

  // leg1->AddEntry(graph_dFG    ,"dFG: NLO + NNLO (eff.) + EW, #mu=0.5 m_{h},","L");
  leg1->AddEntry(graph_ABPS   ,"ABPS:  NLO + NNLO (eff.) + EW","L");
  leg1->AddEntry(graph_NLO_A   ,"NLO (eff.)","L"); 
  leg1->AddEntry(graph_NLO_B   ,"NLO (eff.), rescaled #delta_{NLO}","L");


  graph_NLO_B_scale->Draw("AC3");
  //graph_NLO_A_scale->Draw("C3same");
  
  // graph_dFG->Draw("Csame");
  graph_ABPS->Draw("Csame");
  graph_NLO_A->Draw("Csame");
  graph_NLO_B->Draw("Csame");

      
  leg1->Draw("same");

  return 1;
}
