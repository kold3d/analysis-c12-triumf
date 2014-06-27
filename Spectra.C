#define Spectra_cxx
#include "Spectra.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include "EnergyLoss.h"

void Spectra::InitParameters() {
  pressure         = 397.;   //In Torr
  temperature      = 290.;   //In Kelvin
}

void Spectra::Loop(Int_t incoming)
{
//   In a ROOT session, you can do:
//      Root > .L Spectra.C
//      Root > Spectra t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;


  TFile* cutFile = TFile::Open("proton_cut.root");
  TCutG* proton_cut = (TCutG*)cutFile->Get("proton_cut");
  cutFile->Close();

  Int_t numBins = 50;

  TH1F* s1 = new TH1F("s1","Forward Angles",numBins,0,3.);
  s1->Sumw2();
  s1->GetXaxis()->SetTitle("Center of Mass Energy [MeV]");
  s1->GetYaxis()->SetTitle("Yield [arb. units]");
  s1->GetYaxis()->SetTitleOffset(1.6);
  TH1F* s2 = new TH1F("s2","First Ring",numBins,0,3.);
  s2->Sumw2();
  s2->GetXaxis()->SetTitle("Center of Mass Energy [MeV]");
  s2->GetYaxis()->SetTitle("Yield [arb. units]");
  s2->GetYaxis()->SetTitleOffset(1.6);
  TH1F* s3 = new TH1F("s3","Second Ring",numBins,0,3.);
  s3->Sumw2();
  s3->GetXaxis()->SetTitle("Center of Mass Energy [MeV]");
  s3->GetYaxis()->SetTitle("Yield [arb. units]");
  s3->GetYaxis()->SetTitleOffset(1.6);
  TH1F* s4 = new TH1F("s4","Third Ring",numBins,0,3.);
  s4->Sumw2();
  s4->GetXaxis()->SetTitle("Center of Mass Energy [MeV]");
  s4->GetYaxis()->SetTitle("Yield [arb. units]");
  s4->GetYaxis()->SetTitleOffset(1.6);
  TH1F* s5 = new TH1F("s5","Fourth Ring",numBins,0,3.);
  s5->Sumw2();
  s5->GetXaxis()->SetTitle("Center of Mass Energy [MeV]");
  s5->GetYaxis()->SetTitle("Yield [arb. units]");
  s5->GetYaxis()->SetTitleOffset(1.6);
  TH1F* s6 = new TH1F("s6","Fifth Ring",numBins,0,3.);
  s6->Sumw2();
  s6->GetXaxis()->SetTitle("Center of Mass Energy [MeV]");
  s6->GetYaxis()->SetTitle("Yield [arb. units]");
  s6->GetYaxis()->SetTitleOffset(1.6);

  
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if(!proton_cut->IsInside(measured_energy,sum_dE[0])) continue;

      if(detector == 2 && 
	 (wire[0] == 2  || wire[0] == 3 || wire[0] == 4) && 
	 fabs(position[0]) < 15) s1->Fill(cm_energy[0]);

      if(detector == 2 && 
	 (wire[0] == 1  || wire[0] == 5 || 
	  ((wire[0] == 2  || wire[0] == 3 || wire[0] == 4) && fabs(position[0]) > 15)))
	 s2->Fill(cm_energy[0]);

      if((detector == 3 || detector == 1) && 
	 fabs(position[0]) > 35.96 && fabs(position[0]) < 48.46) 
	s3->Fill(cm_energy[0]);

      if((detector == 3 || detector == 1) && 
	 fabs(position[0]) > 48.46 && fabs(position[0]) < 60.96) 
	s4->Fill(cm_energy[0]);

      if((detector == 3 || detector == 1) && 
	 fabs(position[0]) > 60.96 && fabs(position[0]) < 73.46) 
	s5->Fill(cm_energy[0]);

     if((detector == 3 || detector == 1) && 
	 fabs(position[0]) > 73.46) 
	s6->Fill(cm_energy[0]);
   }
   
   TCanvas* c1 = new TCanvas();
   c1->Divide(3,2);
   c1->cd(1);
   DivideTargetThickness(s1);
   EstimateSolidAngleNorm(s1,1);
   s1->Draw();
   c1->cd(2);
   DivideTargetThickness(s2);
   EstimateSolidAngleNorm(s2,2);
   s2->Draw();
   c1->cd(3);
   DivideTargetThickness(s3);
   EstimateSolidAngleNorm(s3,3);
   s3->Draw();
   c1->cd(4);
   DivideTargetThickness(s4);
   EstimateSolidAngleNorm(s4,4);
   s4->Draw();
   c1->cd(5);
   DivideTargetThickness(s5);
   EstimateSolidAngleNorm(s5,5);
   s5->Draw();
   c1->cd(6);
   DivideTargetThickness(s6);
   EstimateSolidAngleNorm(s6,6);
   s6->Draw();
   
   s1->Scale(1./incoming);
   s2->Scale(1./incoming);
   s3->Scale(1./incoming);
   s4->Scale(1./incoming);
   s5->Scale(1./incoming);
   s6->Scale(1./incoming);
}

void Spectra::DivideTargetThickness(TH1F *f){
  Float_t gasConstant = 8.3144621;
  Float_t torrInPa = 133.322368;
  Float_t molarMassMethane = 0.01604;
  Float_t density = pressure*torrInPa*molarMassMethane/
    gasConstant/temperature*0.001;
  
  EnergyLoss methane("dEdx_carbon_methane_290K_400torr.dat",
		     density*100);
    
  Int_t i_size = f->GetSize();
  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1;i<i_size-1;i++){
    Double_t binLowEdge = xaxis->GetBinLowEdge(i);
    Double_t binUpEdge = xaxis->GetBinUpEdge(i);
    binLowEdge *= 13.0;
    binUpEdge *= 13.0; // From C.M. to Lab Frame
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);
    Double_t delta_x = methane.CalcRange(binUpEdge,binLowEdge);
    delta_x /= 10.0;
    Double_t factor = 4.e-27*density*delta_x*TMath::Na()/molarMassMethane;
    binContent /= factor;
    binError /= factor;
    f->SetBinContent(i,binContent);
    f->SetBinError(i,binError);
  }
}

void Spectra::EstimateSolidAngleNorm(TH1F* f, Int_t region) {
  Float_t gasConstant = 8.3144621;
  Float_t torrInPa = 133.322368;
  Float_t molarMassMethane = 0.01604;
  Float_t density = pressure*torrInPa*molarMassMethane/
    gasConstant/temperature*0.001;

  EnergyLoss methane("dEdx_carbon_methane_290K_400torr.dat",
		     density*100);
    
  Int_t i_size = f->GetSize();
  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1;i<i_size-1;i++){
    Double_t binCenter = xaxis->GetBinCenter(i);
    binCenter *= 13.0;
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);
    Double_t delta_x = methane.CalcRange(79.76,binCenter);
    Double_t x;
    if(region == 1) x = 0;
    else if(region == 2) x = 20.;
    else if(region == 3) x = 42.21;
    else if(region == 4) x = 54.71;
    else if(region == 5) x = 67.21;
    else if(region == 6) x = 79.71;
    Double_t r2 = x*x+(484.5-delta_x)*(484.5-delta_x);
    binContent *= r2;
    binError *= r2;
    f->SetBinContent(i,binContent);
    f->SetBinError(i,binError);
  }
}

Float_t Spectra::pressure;
Float_t Spectra::temperature;

