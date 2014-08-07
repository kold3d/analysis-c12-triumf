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

   std::map<int,std::pair<float,float> > wire_matching;
   std::map<int,std::pair<float,float> > wire_offset_low;
   std::map<int,std::pair<float,float> > wire_gain_diff_low;
   std::map<int,std::pair<float,float> > wire_offset_high;
   std::map<int,std::pair<float,float> > wire_gain_diff_high;
   
  //Matching point
  wire_matching[0] = std::pair<float,float>(115.205742,111.585289);
  wire_matching[1] = std::pair<float,float>(115.205742,119.102402);
  wire_matching[2] = std::pair<float,float>(115.205742,111.585289);
  wire_matching[3] = std::pair<float,float>(115.205742,111.585289);
  wire_matching[4] = std::pair<float,float>(121.651398,120.502197);
  wire_matching[5] = std::pair<float,float>(116.533241,112.460693);
  wire_matching[6] = std::pair<float,float>(115.205742,112.419273);
  wire_matching[7] = std::pair<float,float>(116.535904,119.683525);

  //Gain intercept from channel to mV for channels < 120
  wire_offset_low[0] = std::pair<float,float>(-16.597237,-16.564661);
  wire_offset_low[1] = std::pair<float,float>(-19.425051,-17.123547);
  wire_offset_low[2] = std::pair<float,float>(-18.949633,-20.48386);
  wire_offset_low[3] = std::pair<float,float>(-20.112179,-20.631563);
  wire_offset_low[4] = std::pair<float,float>(-11.075700,-11.315839);
  wire_offset_low[5] = std::pair<float,float>(-12.573323,-12.508942);
  wire_offset_low[6] = std::pair<float,float>(-13.271634,-12.684995);
  wire_offset_low[7] = std::pair<float,float>(-11.083432,-12.783762);

  //Gain slope from channel to mV for channels < 120
  wire_gain_diff_low[0] = std::pair<float,float>(0.584295,0.559508);
  wire_gain_diff_low[1] = std::pair<float,float>(0.623926,0.618987);
  wire_gain_diff_low[2] = std::pair<float,float>(0.620097,0.639888);
  wire_gain_diff_low[3] = std::pair<float,float>(0.645479,0.652981);
  wire_gain_diff_low[4] = std::pair<float,float>(0.395027,0.400060);
  wire_gain_diff_low[5] = std::pair<float,float>(0.416100,0.401508);
  wire_gain_diff_low[6] = std::pair<float,float>(0.437602,0.402177);
  wire_gain_diff_low[7] = std::pair<float,float>(0.405029,0.410348);

  //Gain intercept from channel to mV for channels > 120
  wire_offset_high[0] = std::pair<float,float>(-2.939841,-5.237459);
  wire_offset_high[1] = std::pair<float,float>(-4.236787,-5.647059);
  wire_offset_high[2] = std::pair<float,float>(-2.776138,-2.112775);
  wire_offset_high[3] = std::pair<float,float>(-3.376969,-3.132434);
  wire_offset_high[4] = std::pair<float,float>(-3.189401,-3.161160);
  wire_offset_high[5] = std::pair<float,float>(-3.138486,-4.980239);
  wire_offset_high[6] = std::pair<float,float>(-2.918014,-3.946469);
  wire_offset_high[7] = std::pair<float,float>(-2.594841,-4.308412);

  //Gain slope from channel to mV for channels > 120
  wire_gain_diff_high[0] = std::pair<float,float>(0.413765,0.415051);
  wire_gain_diff_high[1] = std::pair<float,float>(0.444331,0.470588);
  wire_gain_diff_high[2] = std::pair<float,float>(0.431105,0.443483);
  wire_gain_diff_high[3] = std::pair<float,float>(0.455244,0.447787);
  wire_gain_diff_high[4] = std::pair<float,float>(0.297162,0.298397);
  wire_gain_diff_high[5] = std::pair<float,float>(0.304225,0.287005);
  wire_gain_diff_high[6] = std::pair<float,float>(0.309277,0.292386);
  wire_gain_diff_high[7] = std::pair<float,float>(0.302560,0.304448);


  std::map<int,float> energy_cuts;
  energy_cuts[0]=1.97635e+03;
  energy_cuts[1]=2.19773e+03;
  energy_cuts[2]=2.35723e+03;
  energy_cuts[3]=2.38259e+03;
  energy_cuts[4]=2.26268e+03;
  energy_cuts[5]=2.05680e+03;

  std::map<int, std::map<int,float> > centers;
  centers[0][0]=-68.91;
  centers[0][1]=-45.46;
  centers[0][2]=-11.73;
  centers[0][3]=11.73;
  centers[0][4]=45.46;
  centers[0][5]=68.91;
  for(int i = 1; i<5; i++) {
    for(int j = 0;j<6;j++) {
      centers[i][j]=centers[0][j];
    }
  }
  centers[5][0]=-66.91;
  centers[5][1]=-44.14;
  centers[5][2]=-11.385;
  centers[5][3]=11.385;
  centers[5][4]=44.14;
  centers[5][5]=66.91;
  for(int i = 6; i<8; i++) {
    for(int j = 0;j<6;j++) {
      centers[i][j]=centers[0][j];
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
	    float left = (pc_ch_left_e[j]<120) ? pc_ch_left_e[j]*wire_gain_diff_low[pc_wire[j]-1].first+wire_offset_low[pc_wire[j]-1].first :
	      pc_ch_left_e[j]*wire_gain_diff_high[pc_wire[j]-1].first+wire_offset_high[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]<120) ? pc_ch_right_e[j]*wire_gain_diff_low[pc_wire[j]-1].second+wire_offset_low[pc_wire[j]-1].second :
	      pc_ch_right_e[j]*wire_gain_diff_high[pc_wire[j]-1].second+wire_offset_high[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][0]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==3 && (si_quad[i]==2 || si_quad[i]==4)) {
	  if(si_cal_e[i]>energy_cuts[1]+100||si_cal_e[i]<energy_cuts[1]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left = (pc_ch_left_e[j]<120) ? pc_ch_left_e[j]*wire_gain_diff_low[pc_wire[j]-1].first+wire_offset_low[pc_wire[j]-1].first :
	      pc_ch_left_e[j]*wire_gain_diff_high[pc_wire[j]-1].first+wire_offset_high[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]<120) ? pc_ch_right_e[j]*wire_gain_diff_low[pc_wire[j]-1].second+wire_offset_low[pc_wire[j]-1].second :
	      pc_ch_right_e[j]*wire_gain_diff_high[pc_wire[j]-1].second+wire_offset_high[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][1]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==2 && (si_quad[i]==1 || si_quad[i]==3)) {
	  if(si_cal_e[i]>energy_cuts[2]+100||si_cal_e[i]<energy_cuts[2]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left = (pc_ch_left_e[j]<120) ? pc_ch_left_e[j]*wire_gain_diff_low[pc_wire[j]-1].first+wire_offset_low[pc_wire[j]-1].first :
	      pc_ch_left_e[j]*wire_gain_diff_high[pc_wire[j]-1].first+wire_offset_high[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]<120) ? pc_ch_right_e[j]*wire_gain_diff_low[pc_wire[j]-1].second+wire_offset_low[pc_wire[j]-1].second :
	      pc_ch_right_e[j]*wire_gain_diff_high[pc_wire[j]-1].second+wire_offset_high[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][2]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==2 && (si_quad[i]==2 || si_quad[i]==4)) {
	  if(si_cal_e[i]>energy_cuts[3]+100||si_cal_e[i]<energy_cuts[3]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left = (pc_ch_left_e[j]<120) ? pc_ch_left_e[j]*wire_gain_diff_low[pc_wire[j]-1].first+wire_offset_low[pc_wire[j]-1].first :
	      pc_ch_left_e[j]*wire_gain_diff_high[pc_wire[j]-1].first+wire_offset_high[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]<120) ? pc_ch_right_e[j]*wire_gain_diff_low[pc_wire[j]-1].second+wire_offset_low[pc_wire[j]-1].second :
	      pc_ch_right_e[j]*wire_gain_diff_high[pc_wire[j]-1].second+wire_offset_high[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][3]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==1 && (si_quad[i]==1 || si_quad[i]==3)) {
	  if(si_cal_e[i]>energy_cuts[4]+100||si_cal_e[i]<energy_cuts[4]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left = (pc_ch_left_e[j]<120) ? pc_ch_left_e[j]*wire_gain_diff_low[pc_wire[j]-1].first+wire_offset_low[pc_wire[j]-1].first :
	      pc_ch_left_e[j]*wire_gain_diff_high[pc_wire[j]-1].first+wire_offset_high[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]<120) ? pc_ch_right_e[j]*wire_gain_diff_low[pc_wire[j]-1].second+wire_offset_low[pc_wire[j]-1].second :
	      pc_ch_right_e[j]*wire_gain_diff_high[pc_wire[j]-1].second+wire_offset_high[pc_wire[j]-1].second;
	    h[pc_wire[j]-1][4]->Fill((right-left)/(left+right));
	  }
	}
	if(si_det[i]==1 && (si_quad[i]==2 || si_quad[i]==4)) {
	  if(si_cal_e[i]>energy_cuts[5]+100||si_cal_e[i]<energy_cuts[5]-100) continue;
	  for(int j = 0;j<pc_mul;j++) {
	    if(pc_ch_left_e[j]<80 || pc_ch_right_e[j]< 80) continue;
	    float left = (pc_ch_left_e[j]<120) ? pc_ch_left_e[j]*wire_gain_diff_low[pc_wire[j]-1].first+wire_offset_low[pc_wire[j]-1].first :
	      pc_ch_left_e[j]*wire_gain_diff_high[pc_wire[j]-1].first+wire_offset_high[pc_wire[j]-1].first;
	    float right = (pc_ch_right_e[j]<120) ? pc_ch_right_e[j]*wire_gain_diff_low[pc_wire[j]-1].second+wire_offset_low[pc_wire[j]-1].second :
	      pc_ch_right_e[j]*wire_gain_diff_high[pc_wire[j]-1].second+wire_offset_high[pc_wire[j]-1].second;
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
