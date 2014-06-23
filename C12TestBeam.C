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

   Int_t positionThreshold=120;

   TH1F* h_si[3][4];
   TH1F* h_pc[8][2];
   TH1F* h_pos[8];

   for(Int_t j =0;j<3;j++) {
     for(Int_t k = 0;k<4;k++) {
       TString name = Form("si_%d_%d",j+1,k+1);
       h_si[j][k] = new TH1F(name,name,400,0,12000);
     }
   }

   for(Int_t i = 0;i<8;i++) {
     TString name = Form("pos_%d",i+1);
     h_pos[i] = new TH1F(name,name,200,-1,1);
     name = Form("pc_%d_left",i+1);
     h_pc[i][0]= new TH1F(name,name,4000,0,8000);
     name = Form("pc_%d_right",i+1);
     h_pc[i][1]= new TH1F(name,name,4000,0,8000);
   }

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%1000==0&&jentry!=0){
	      printf("Processing event %lld of %lld...\n",jentry,nentries);
      }
      
      Int_t highEDet = -1;
      Int_t highEQuad = -1;
      Float_t highE = 0;

      for(Int_t i =0;i<si_mul;i++) {
	      h_si[si_det[i]-1][si_quad[i]-1]->Fill(si_cal[i]);
	      if(si_cal[i]>highE) {
	        highE = si_cal[i];
	        highEDet = si_det[i];
	        highEQuad = si_quad[i];
	      }
      }

      for(Int_t i=0;i<pc_mul;i++)
      {
	      Float_t left = pc_ch_left[i];
        Float_t right = pc_ch_right[i];
        h_pc[pc_wire[i]-1][0]->Fill(left);
        h_pc[pc_wire[i]-1][1]->Fill(right);	
        if(left<positionThreshold || right<positionThreshold) continue;
        Float_t pos = (right-left)/(left+right);
        if(highEDet == 3 && (highEQuad == 1 || highEQuad == 3) &&
           highE>800 && highE<1200){
          h_pos[pc_wire[i]-1]->Fill(pos);
        }
      }
   }
   
   TCanvas* c = new TCanvas();
   c->Divide(4,2);
   for(Int_t i = 0;i<8;i++) {
     c->cd(i+1);
     h_pos[i]->Draw();
   }
}
