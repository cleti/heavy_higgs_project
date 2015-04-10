

int test()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleAlign(12);
  gStyle->SetPadTopMargin(0.03);
  
  h = new TH1D("h1","test histogram",20,0.0,1.0);
  h->SetTitle("something here");
  h->Draw("");

}
