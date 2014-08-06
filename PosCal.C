#define PosCal_cxx
#include "PosCal.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <cstdio>

void PosCal::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L PosCal.C
//      Root > PosCal t
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

   std::map<int,std::pair<float,float> > wire_offset;
   std::map<int,std::pair<float,float> > wire_gain_diff;

     //Noise peaks
  wire_offset[7] = std::pair<float,float>(61.39516087,63.36229816);
  wire_offset[6] = std::pair<float,float>(62.24532941,77.76395862);
  wire_offset[5] = std::pair<float,float>(62.72446862,77.51005426);
  wire_offset[4] = std::pair<float,float>(67.83539408,62.59553001);
  wire_offset[3] = std::pair<float,float>(47.43670403,47.48292601);
  wire_offset[2] = std::pair<float,float>(50.65317676,49.86160510);
  wire_offset[1] = std::pair<float,float>(51.20105207,43.77215703);
  wire_offset[0] = std::pair<float,float>(54.65770540,65.99354481);

  //Gain matches from pulser (ratio left to right);
  wire_gain_diff[0] = std::pair<float,float>(2.185976587,2.139817530);
  wire_gain_diff[1] = std::pair<float,float>(2.055031890,1.977072138);
  wire_gain_diff[2] = std::pair<float,float>(2.104007381,2.033647438);
  wire_gain_diff[3] = std::pair<float,float>(2.000304695,2.037351174);
  wire_gain_diff[4] = std::pair<float,float>(3.090042761,3.101041098);
  wire_gain_diff[5] = std::pair<float,float>(3.040920343,3.167046917);
  wire_gain_diff[6] = std::pair<float,float>(2.976782799,3.091595533);
  wire_gain_diff[7] = std::pair<float,float>(3.051640102,3.045576998);

  std::map<int,float> energy_cuts;
  energy_cuts[0]=1.97635e+03;
  energy_cuts[1]=2.19773e+03;
  energy_cuts[2]=2.35723e+03;
  energy_cuts[3]=2.38259e+03;
  energy_cuts[4]=2.26268e+03;
  energy_cuts[5]=2.05680e+03;

  std::map<int,float> centers;
  centers[0]=-68.88180188;
  centers[1]=-45.43986005;
  centers[2]=-11.72097092;
  centers[3]=11.72097092;
  centers[4]=45.43986005;
  centers[5]=68.88180188;

   Long64_t nentries = fChain->GetEntriesFast();
   TH1F* h[8][6];

   for(int i = 0;i<8;i++) {
     for(int k = 0;k<6;k++) {
       TString name = Form("h_%d_%d",i,k);
       h[i][k]= new TH1F(name,name,100,-1,1);
     }
   }
      
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      for(int i = 0;i<si_mul;i++) {
	if(si_det[i]==3 && (si_quad[i]==1 || si_quad[i]==3)) {
	  if(si_cal_e[i]>energy_cuts[0]+100||si_cal_e[i]<energy_cuts[0]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left =(pc_ch_left_e[j]-wire_offset[pc_wire[j]-1].first)/wire_gain_diff[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]-wire_offset[pc_wire[j]-1].second)/wire_gain_diff[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][0]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==3 && (si_quad[i]==2 || si_quad[i]==4)) {
	  if(si_cal_e[i]>energy_cuts[1]+100||si_cal_e[i]<energy_cuts[1]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left =(pc_ch_left_e[j]-wire_offset[pc_wire[j]-1].first)/wire_gain_diff[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]-wire_offset[pc_wire[j]-1].second)/wire_gain_diff[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][1]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==2 && (si_quad[i]==1 || si_quad[i]==3)) {
	  if(si_cal_e[i]>energy_cuts[2]+100||si_cal_e[i]<energy_cuts[2]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left =(pc_ch_left_e[j]-wire_offset[pc_wire[j]-1].first)/wire_gain_diff[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]-wire_offset[pc_wire[j]-1].second)/wire_gain_diff[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][2]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==2 && (si_quad[i]==2 || si_quad[i]==4)) {
	  if(si_cal_e[i]>energy_cuts[3]+100||si_cal_e[i]<energy_cuts[3]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left =(pc_ch_left_e[j]-wire_offset[pc_wire[j]-1].first)/wire_gain_diff[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]-wire_offset[pc_wire[j]-1].second)/wire_gain_diff[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][3]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==1 && (si_quad[i]==1 || si_quad[i]==3)) {
	  if(si_cal_e[i]>energy_cuts[4]+100||si_cal_e[i]<energy_cuts[4]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left =(pc_ch_left_e[j]-wire_offset[pc_wire[j]-1].first)/wire_gain_diff[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]-wire_offset[pc_wire[j]-1].second)/wire_gain_diff[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][4]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==1 && (si_quad[i]==2 || si_quad[i]==4)) {
	  if(si_cal_e[i]>energy_cuts[5]+100||si_cal_e[i]<energy_cuts[5]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left =(pc_ch_left_e[j]-wire_offset[pc_wire[j]-1].first)/wire_gain_diff[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]-wire_offset[pc_wire[j]-1].second)/wire_gain_diff[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][5]->Fill((right-left)/(left+right));
	  }
	}
     }
   }
   
   float mean[8][6];
   for(int i = 0;i<8;i++) {
     for(int k = 0;k<6;k++) {
       h[i][k]->Fit("gaus");
       mean[i][k] = h[i][k]->GetFunction("gaus")->GetParameter(1);
       //mean[i][k] = h[i][k]->GetMean();
     }
   }

   TGraph* line[8];
   float slope[8],offset[8];
   TCanvas *c1 = new TCanvas;
   c1->Divide(4,2);
   for(int i = 0;i<8;i++) {
     line[i] = new TGraph(6);
     line[i]->SetTitle(Form("Wire %d",i+1));
     c1->cd(i+1);
     for(int k = 0;k<6;k++) {
       line[i]->SetPoint(k,mean[i][k],centers[k]);
     }
     line[i]->Draw("ap");
     line[i]->SetMarkerStyle(kFullCircle);
     line[i]->Fit("pol1");
     slope[i]=line[i]->GetFunction("pol1")->GetParameter(1);
     offset[i]=line[i]->GetFunction("pol1")->GetParameter(0);
   }
   
   FILE *out;
   out = fopen("position_cal.txt","w");
   for(int i = 0;i<8;i++) {
     fprintf(out,"%d   %f   %f\n",i+1,slope[i],offset[i]);
   }
   fclose(out);
   
}
