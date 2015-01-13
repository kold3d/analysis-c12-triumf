{
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TFile* file = TFile::Open("cluster_out.root");
TCanvas* c1 = file->Get("c1");
c1->Draw("goff");
TH1F* s1 = (TH1F*) c1->cd(1)->GetPrimitive("s1")->Clone("s1");
TH1F* s2 = (TH1F*) c1->cd(2)->GetPrimitive("s2")->Clone("s2");
TH1F* s3 = (TH1F*) c1->cd(3)->GetPrimitive("s3")->Clone("s3");

TGraph* eff1 = new TGraph();
TGraph* eff2 = new TGraph();
TGraph* eff3 = new TGraph();
std::ifstream in("efficiencies.out");
while(!in.eof()) {
  int region;
  float x,y;
  in >> region >> x >> y;
  if(!in.eof()) {
    if(region==1) eff1->SetPoint(eff1->GetN(),x,y);
    else if(region==2) eff2->SetPoint(eff2->GetN(),x,y);
    else if(region==3) eff3->SetPoint(eff3->GetN(),x,y);
  }
}
in.close();

s1->Scale(1000.);
s2->Scale(1000.);
s3->Scale(1000.);

s1->GetXaxis()->SetRangeUser(0.5,3.4);
s1->GetXaxis()->SetTitle("Center of Mass Energy [MeV]");
s1->GetXaxis()->CenterTitle();
s1->GetYaxis()->SetTitle("Differential Cross Section [mb/sr]");
s1->GetYaxis()->CenterTitle();

s2->SetBinContent(16,0);
s2->SetBinContent(17,0);
s2->SetBinError(16,0);
s2->SetBinError(17,0);

s1->SetMarkerStyle(kFullCircle);
s2->SetMarkerStyle(kFullCircle);
s3->SetMarkerStyle(kFullCircle);
s1->SetLineWidth(2);
s2->SetLineWidth(2);
s3->SetLineWidth(2);
s1->SetMarkerColor(kBlue);
s2->SetMarkerColor(kRed);
s3->SetMarkerColor(kGreen);
s1->SetLineColor(kBlue);
s2->SetLineColor(kRed);
s3->SetLineColor(kGreen);
	


 FILE* file1 = fopen("region1.dat","w");
 TFile* afile1 = TFile::Open("angle_dists/region_1.root");
 for(int i = 1;i<=s1->GetNbinsX();i++) {
   TH1F* hist = afile1->Get(Form("bin_%d_cm_fk",i+1));
   if(s1->GetBinContent(i)>0) {
     double eff = eff1->Eval(s1->GetBinCenter(i));
     double binContent = (eff < 1e-5) ? 0. : s1->GetBinContent(i)/eff;
     double binError = (eff < 1e-5) ? 0. : s1->GetBinError(i)/eff;
     s1->SetBinContent(i,binContent);
     s1->SetBinError(i,binError);
   }
   float angle = (hist) ? hist->GetMean() : 0.;
   fprintf(file1,"%f %f %f %f\n",s1->GetBinCenter(i),s1->GetBinContent(i),s1->GetBinError(i),angle);
 }
 afile1->Close();
 fclose(file1);
 FILE* file2 = fopen("region2.dat","w");
 TFile* afile2 = TFile::Open("angle_dists/region_2.root");
 for(int i = 1;i<=s2->GetNbinsX();i++) {
   TH1F* hist = afile2->Get(Form("bin_%d_cm_fk",i+1));
   if(s2->GetBinContent(i)>0) {
     double eff = eff2->Eval(s2->GetBinCenter(i));
     double binContent = (eff < 1e-5) ? 0. : s2->GetBinContent(i)/eff;
     double binError = (eff < 1e-5) ? 0. : s2->GetBinError(i)/eff;
     s2->SetBinContent(i,binContent);
     s2->SetBinError(i,binError);
   }
   float angle = (hist) ? hist->GetMean() : 0.;
   fprintf(file2,"%f %f %f %f\n",s2->GetBinCenter(i),s2->GetBinContent(i),s2->GetBinError(i),angle);
 }
 afile2->Close();
 fclose(file2);
 FILE* file3 = fopen("region3.dat","w");
 TFile* afile3 = TFile::Open("angle_dists/region_3.root");
 for(int i = 1;i<=s3->GetNbinsX();i++) {
   TH1F* hist = afile3->Get(Form("bin_%d_cm_fk",i+1));
   if(s3->GetBinContent(i)>0) {
     double eff = eff3->Eval(s3->GetBinCenter(i));
     double binContent = (eff < 1e-5) ? 0. : s3->GetBinContent(i)/eff;
     double binError = (eff < 1e-5) ? 0. : s3->GetBinError(i)/eff;
     s3->SetBinContent(i,binContent);
     s3->SetBinError(i,binError);
   }
   float angle = (hist) ? hist->GetMean() : 0.;
   fprintf(file3,"%f %f %f %f\n",s3->GetBinCenter(i),s3->GetBinContent(i),s3->GetBinError(i),angle);
 }
 afile3->Close();
 fclose(file3);

 TCanvas* c2 = new TCanvas;
 s1->Draw();
 s2->Draw("same");
 s3->Draw("same");
 c2->Print("He8pp_final.pdf","pdf");
}
