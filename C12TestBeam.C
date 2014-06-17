#define C12TestBeam_cxx
#include "C12TestBeam.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void C12TestBeam::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L C12TestBeam.C
//      Root > C12TestBeam t
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


   TH1F* h_si[3][4];
   for(Int_t j =0;j<3;j++) {
     for(Int_t k = 0;k<4;k++) {
       TString name = Form("si_%d_%d",j+1,k+1);
       h_si[j][k] = new TH1F(name,name,400,0,12000);
     }
   }

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%1000==0&&jentry!=0)
	printf("Processing event %d of %d...\n",jentry,nentries);
      
      for(Int_t i =0;i<si_mul;i++) {
	h_si[si_det[i]-1][si_quad[i]-1]->Fill(si_cal[i]);
      }
   }
   
}
