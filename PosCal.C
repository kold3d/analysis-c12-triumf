#define PosCal_cxx
#include "PosCal.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <cstdio>
#include "Calibrations.h"

void PosCal::Loop()
{
  if (fChain == 0) return;

  Float_t resonanceToSi = 459.35;

  std::map<int,float> energy_cuts;
  energy_cuts[0]=1.97635e+03;
  energy_cuts[1]=2.19773e+03;
  energy_cuts[2]=2.35723e+03;
  energy_cuts[3]=2.38259e+03;
  energy_cuts[4]=2.26268e+03;
  energy_cuts[5]=2.05680e+03;

  std::map<int, std::map<int,float> > centers;
  centers[0][0]=-73.46*(resonanceToSi-Calibrations::anode_to_si)/resonanceToSi;
  centers[0][1]=-48.46*(resonanceToSi-Calibrations::anode_to_si)/resonanceToSi;
  centers[0][2]=-12.5*(resonanceToSi-Calibrations::anode_to_si)/resonanceToSi;
  centers[0][3]=-1*centers[0][2];
  centers[0][4]=-1*centers[0][1];
  centers[0][5]=-1*centers[0][0];
  for(int i = 1; i<5; i++) {
    for(int j = 0;j<6;j++) {
      centers[i][j]=centers[0][j];
    }
  }
  centers[5][0]=-73.46*(resonanceToSi-Calibrations::anode_to_si-Calibrations::anode_sep)/resonanceToSi;
  centers[5][1]=-48.46*(resonanceToSi-Calibrations::anode_to_si-Calibrations::anode_sep)/resonanceToSi;
  centers[5][2]=-12.5*(resonanceToSi-Calibrations::anode_to_si-Calibrations::anode_sep)/resonanceToSi;
  centers[5][3]=-1*centers[0][2];
  centers[5][4]=-1*centers[0][1];
  centers[5][5]=-1*centers[0][0];
  for(int i = 6; i<8; i++) {
    for(int j = 0;j<6;j++) {
      centers[i][j]=centers[5][j];
    }
  }

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
	    Float_t left = pc_ch_left_e[j];
	    Float_t right = pc_ch_right_e[j];
	    Calibrations::MatchPC(left,right,pc_wire[j]-1);
	    h[pc_wire[j]-1][0]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==3 && (si_quad[i]==2 || si_quad[i]==4)) {
	  if(si_cal_e[i]>energy_cuts[1]+100||si_cal_e[i]<energy_cuts[1]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    Float_t left = pc_ch_left_e[j];
	    Float_t right = pc_ch_right_e[j];
	    Calibrations::MatchPC(left,right,pc_wire[j]-1);
	    h[pc_wire[j]-1][1]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==2 && (si_quad[i]==1 || si_quad[i]==3)) {
	  if(si_cal_e[i]>energy_cuts[2]+100||si_cal_e[i]<energy_cuts[2]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    Float_t left = pc_ch_left_e[j];
	    Float_t right = pc_ch_right_e[j];
	    Calibrations::MatchPC(left,right,pc_wire[j]-1);
	    h[pc_wire[j]-1][2]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==2 && (si_quad[i]==2 || si_quad[i]==4)) {
	  if(si_cal_e[i]>energy_cuts[3]+100||si_cal_e[i]<energy_cuts[3]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    Float_t left = pc_ch_left_e[j];
	    Float_t right = pc_ch_right_e[j];
	    Calibrations::MatchPC(left,right,pc_wire[j]-1);
	    h[pc_wire[j]-1][3]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==1 && (si_quad[i]==1 || si_quad[i]==3)) {
	  if(si_cal_e[i]>energy_cuts[4]+100||si_cal_e[i]<energy_cuts[4]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    Float_t left = pc_ch_left_e[j];
	    Float_t right = pc_ch_right_e[j];
	    Calibrations::MatchPC(left,right,pc_wire[j]-1);
	    h[pc_wire[j]-1][4]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==1 && (si_quad[i]==2 || si_quad[i]==4)) {
	  if(si_cal_e[i]>energy_cuts[5]+100||si_cal_e[i]<energy_cuts[5]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    Float_t left = pc_ch_left_e[j];
	    Float_t right = pc_ch_right_e[j];
	    Calibrations::MatchPC(left,right,pc_wire[j]-1);
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
       line[i]->SetPoint(k,mean[i][k],centers[i][k]);
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

