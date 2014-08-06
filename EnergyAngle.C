#define EnergyAngle_cxx
#include "EnergyAngle.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "EnergyLoss.h"
#include "TMath.h"
#include <fstream>
#include <sstream>

LookupTable EnergyAngle::table;

void EnergyAngle::InitParameters() {

  sum_pc_threshold = 0.;   //In Channels
  si_threshold     = 350.;   //In KeV
  beam_energy      = 41.62; //In MeV, after window
  pressure         = 785.;   //In Torr
  temperature      = 295.;   //In Kelvin
  m1 = 12.;
  m2 = 1.;

  Float_t gasConstant = 8.3144621;
  Float_t torrInPa = 133.322368;
  Float_t molarMassMethane = 0.01604;

  Float_t density = 0.000622*pressure/760/temperature*293;

  projectile = new EnergyLoss("dEdx_carbon_methane_290K_400torr.dat",density*100.);
  proton = new EnergyLoss("dEdx_proton_methane_290K_400torr.dat",density*100.);

  
  //Noise peaks
  wire_offset[7] = std::pair<float,float>(0.0,0.0);
  wire_offset[6] = std::pair<float,float>(38.5129,39.6507);
  wire_offset[5] = std::pair<float,float>(39.8341,39.0006);
  wire_offset[4] = std::pair<float,float>(38.8358,39.2025);
  wire_offset[3] = std::pair<float,float>(41.139,39.7448);
  wire_offset[2] = std::pair<float,float>(39.1638,39.6441);
  wire_offset[1] = std::pair<float,float>(40.8464,39.9842);
  wire_offset[0] = std::pair<float,float>(38.1933,39.6582);

  //Gain matches from pulser (ratio left to right);
  wire_gain_diff[7] = std::pair<float,float>(1.00636,1.00000);
  wire_gain_diff[6] = std::pair<float,float>(0.97957,1.02192);
  wire_gain_diff[5] = std::pair<float,float>(1.00233,1.05483);
  wire_gain_diff[4] = std::pair<float,float>(1.02017,1.02399);
  wire_gain_diff[3] = std::pair<float,float>(0.97122,0.99870);
  wire_gain_diff[2] = std::pair<float,float>(1.02667,0.98590);
  wire_gain_diff[1] = std::pair<float,float>(0.99821,0.95724);
  wire_gain_diff[0] = std::pair<float,float>(1.0,1.0);

   //Position calibration (slope,intercept)
  std::ifstream in("position_cal.txt");
  while(!in.eof()) {
    std::string line;
    getline(in,line);
    if(in.eof()) continue;
    std::stringstream stm;
    stm.str(line);
    int wire;
    float m,b;
    if(stm >> wire >> m >> b) wire_pos_cal[wire-1] =  std::pair<Float_t,Float_t>(m,b);
    else printf("WHY DIDN'T THAT WORK?!?\n");
  }
  in.close();
}

void EnergyAngle::Loop()
{
  if (fChain == 0) return;

  Float_t measured_energy,cm_energy[2],lab_angle[2],position[2],sum_dE[2];
  Int_t detector,quadrant,wire[2];

  TFile* file = new TFile("energy_angle.root","recreate");

  TH1F* h_no_pos = new TH1F("h_no_pos","h_no_pos",400,0,12000);
  TH1F* h_all_si = new TH1F("h_all_si","h_all_si",400,0,12000);

  TH2F* dE_E[8];
  for(Int_t i = 0;i<8;i++) {
    TString name = Form("dE_E_%d",i+1);
    dE_E[i] = new TH2F(name,name,400,0,12000,200,0,8000);
  }

  TTree* outTree = new TTree("energyAngle","Events tagged with reconstructed energy and angle.");
  outTree->Branch("measured_energy",&measured_energy,"measured_energy/F");
  outTree->Branch("quadrant",&quadrant,"quadrant/I");
  outTree->Branch("detector",&detector,"detector/I");
  outTree->Branch("cm_energy",cm_energy,"cm_energy[2]/F");
  outTree->Branch("lab_angle",lab_angle,"lab_angle[2]/F");
  outTree->Branch("position",position,"position[2]/F");
  outTree->Branch("sum_dE",sum_dE,"sum_dE[2]/F");
  outTree->Branch("wire",wire,"wire[2]/I");

  Int_t goodSi = 0, noPosition=0, noCMEnergy=0;
  Int_t numBelowPCThreshold= 0.;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      //Find legit hit in si_detector
      Int_t highSiEDet = -1;
      Int_t highSiEQuad = -1;
      Float_t highSiE = 0;
      for(Int_t i =0;i<si_mul;i++) {
	if(si_cal_e[i]>highSiE) {
	  highSiE = si_cal_e[i];
	  highSiEDet = si_det[i];
	  highSiEQuad = si_quad[i];
	}
      }

      //If Si energy is less than threshold, next event
      if(highSiE < si_threshold) continue;
      goodSi++;
      
      if(highSiEDet==2 && highSiEQuad ==1) h_all_si->Fill(highSiE);

      //Add proton information to tree
      measured_energy = highSiE;
      quadrant = highSiEQuad;
      detector = highSiEDet;

      //Find position in front and back wire
      Float_t positionFront = 0.;
      Float_t positionRear  = 0.;
      Bool_t goodRearPosition = false;
      Bool_t goodFrontPosition = false;
      Float_t highEPCFront  = 0.;
      Float_t highEPCRear   = 0.;
      UChar_t wireFront     = 0;
      UChar_t wireRear      = 0;
      Bool_t wireAboveThreshold = false;
      for(Int_t i = 0;i<pc_mul;i++) {
	Float_t left = pc_ch_left_e[i];
	Float_t right = pc_ch_right_e[i];
	MatchPC(pc_wire[i]-1,left,right);
	Float_t sum = left+right;
	if( sum < sum_pc_threshold ) {
	  continue;
	}
	wireAboveThreshold = true;
	if(pc_wire[i] < 6) {
	  if(sum > highEPCRear) {
	    highEPCRear = sum;
	    positionRear = CalcPosition(pc_wire[i]-1,left,right);
	    wireRear = pc_wire[i];
	    goodRearPosition = true;
	  }
	} else {
	  if(sum > highEPCFront) {
	    highEPCFront = sum;
	    positionFront = CalcPosition(pc_wire[i]-1,left,right);
	    wireFront = pc_wire[i];
	    goodFrontPosition = true;
	  }
	}
      }
      
      if(!wireAboveThreshold) numBelowPCThreshold++;

      if(!goodRearPosition&&!goodFrontPosition) {
	noPosition++;
	if(highSiEDet==2 && highSiEQuad ==1) h_no_pos->Fill(highSiE);      
	continue;
      }

      if(wireRear==1) dE_E[0]->Fill(highSiE,highEPCRear);
      else if(wireRear==2) dE_E[1]->Fill(highSiE,highEPCRear);
      else if(wireRear==3) dE_E[2]->Fill(highSiE,highEPCRear);
      else if(wireRear==4) dE_E[3]->Fill(highSiE,highEPCRear);
      else if(wireRear==5) dE_E[4]->Fill(highSiE,highEPCRear);

      if(wireFront == 6) dE_E[5]->Fill(highSiE,highEPCFront);
      else if(wireFront == 7) dE_E[6]->Fill(highSiE,highEPCFront);
      else if(wireFront == 8) dE_E[7]->Fill(highSiE,highEPCFront);

      Float_t angleRear=0.,cmEnergyRear = 0.;
      if(positionRear != 0. && fabs(positionRear)<=80.) {
	std::pair<Float_t,Float_t> entry = LookupCMEnergyAngle(wireRear-1,positionRear,highSiE/1000.);
	if(entry.first!=0.) {
	  angleRear = entry.second;
	  cmEnergyRear = entry.first;
	}
      }
      Float_t angleFront=0.,cmEnergyFront = 0.;
      if(positionFront != 0. && fabs(positionFront)<=80.) {
	std::pair<Float_t,Float_t> entry = LookupCMEnergyAngle(wireFront-1,positionFront,highSiE/1000.);
	if(entry.first!=0.) {
	  angleFront = entry.second;
	  cmEnergyFront = entry.first;
	}
      }
      if(cmEnergyFront == 0. && cmEnergyRear == 0.) {
	noCMEnergy++;
	continue;
      }
      
      wire[0] = wireRear;
      wire[1] = wireFront;
      cm_energy[0] = cmEnergyRear;
      cm_energy[1] = cmEnergyFront;
      lab_angle[0] = angleRear;
      lab_angle[1] = angleFront;
      sum_dE[0] = highEPCRear;
      sum_dE[1] = highEPCFront;
      position[0] = positionRear;
      position[1] = positionFront;

      outTree->Fill();
   }

   for(Int_t i = 0;i<8;i++) dE_E[i]->Write();
   outTree->Write();
   h_no_pos->Write();
   h_all_si->Write();

   file->Close();

   printf("Total Events: %d Below PC Threshold: %d No Position: %d No CM Energy %d\n",goodSi,numBelowPCThreshold,noPosition,noCMEnergy);   
}

void EnergyAngle::MatchPC(UChar_t wire, Float_t& left_ch,Float_t& right_ch) {
  Float_t left = left_ch-wire_offset[wire].first;
  Float_t right = right_ch-wire_offset[wire].second;
  left *= wire_gain_diff[wire].first;
  right *= wire_gain_diff[wire].second;

  left_ch = left;
  right_ch = right;
}

Float_t EnergyAngle::CalcPosition(UChar_t wire, Float_t left_ch, Float_t right_ch) {
  Float_t x = (right_ch-left_ch)/(right_ch+left_ch);
  
  Float_t pos = wire_pos_cal[wire].first*x+wire_pos_cal[wire].second;
  
  return pos;
}

std::pair<Float_t,Float_t> EnergyAngle::LookupCMEnergyAngle(UChar_t wire, Float_t x, Float_t protonEnergy) {
  Float_t newX = floor(10.*fabs(x)/5.)*5./10.;
  if(x-newX>=0.25) newX += 0.5;

  if(wire > 7 || wire < 0 || newX < 0. || newX > 80. ) {
    printf("wire=%d x=%f out of Range.\n",wire,newX);
    return std::pair<Float_t,Float_t>(0.,0.);
  }

  Int_t mappedWire = wire;
  if(wire == 3) mappedWire = 1;
  else if(wire == 4) mappedWire = 0;
  else if(wire == 7) mappedWire = 5;

  std::vector<LookupEntry> vec = table[mappedWire][newX];
  Int_t size = vec.size();

  Int_t start = floor(size/2);
  Int_t direction = 0;
  if(vec[start].protonEnergy > protonEnergy) direction = 1;
  else if(vec[start].protonEnergy < protonEnergy) direction  = -1;
  else return std::pair<Float_t,Float_t>(vec[start].cmEnergy,vec[start].angle);

  Bool_t done = false;
  Float_t previousAngle = vec[start].angle;
  Float_t previousCMEnergy = vec[start].cmEnergy;
  Float_t previousProtonEnergy = vec[start].protonEnergy;

  Float_t foundAngle =0;
  Float_t foundCMEnergy = 0.;

  Int_t i = start+direction;
  while(!done && i>0 && i<size) {
    Float_t thisAngle = vec[i].angle;
    Float_t thisCMEnergy = vec[i].cmEnergy;
    Float_t thisProtonEnergy = vec[i].protonEnergy;

    if((direction < 0 && thisProtonEnergy >=protonEnergy) ||
	      (direction >0 && thisProtonEnergy <= protonEnergy)) {
      Float_t energy_slope = (thisCMEnergy-previousCMEnergy)/(thisProtonEnergy-previousProtonEnergy);
      Float_t energy_intercept = thisCMEnergy-thisProtonEnergy*energy_slope;
      Float_t angle_slope = (thisAngle-previousAngle)/(thisProtonEnergy-previousProtonEnergy);
      Float_t angle_intercept = thisAngle-thisProtonEnergy*angle_slope;
      foundAngle = angle_slope*protonEnergy+angle_intercept;
      foundCMEnergy = energy_slope*protonEnergy+energy_intercept;
      done = true;
    } else {
      previousAngle = thisAngle;
      previousCMEnergy = thisCMEnergy;  
      previousProtonEnergy = thisProtonEnergy;
      i += direction;
    }
  }

  if(!done) {
    //printf("Entry wire=%d x=%f energy=%f not found in table!\n",wire,newX,protonEnergy);
  }
  
  return std::pair<Float_t,Float_t>(foundCMEnergy,foundAngle);
}

void EnergyAngle::ReadLookupTable() {
  table.clear();
  std::ifstream in("lookup_table.out");
  std::string line;
  while(!in.eof()) {
    getline(in,line);
    if(!in.eof()) {
      std::istringstream stm;
      stm.str(line);
      Int_t wire;
      Float_t position,range,r,angle,cmEnergy,protonEnergy;
      stm >> wire >> position >> range >> r >> angle >> cmEnergy >> protonEnergy;
      LookupEntry entry = {range,r,angle,cmEnergy,protonEnergy};
      table[wire][position].push_back(entry);
    }
  }
}

void EnergyAngle::CalcLookupTable() {
  InitParameters();

  Float_t deltaBeamE = 0.10;
  Float_t distanceToFirstWire = 472.;
  Float_t distanceToSecondWire = 484.5; 
  Float_t distanceToSiDets= 513.;

  Float_t wireHeight[8] = {
    -20.32,
    -10.16,
    0.0,
    10.16,
    20.32,
    -15.24,
    0.00,
    15.24
  };
  
  table.clear();
  FILE* out = fopen("lookup_table.out","w");
  for(UChar_t wire = 0;wire<8;wire++) {
    if(wire==7) {
      table[wire]=table[5];
      continue;
    } else if(wire == 4) {
      table[wire]=table[0];
      continue;
    } else if(wire == 3) {
      table[wire] = table[1];
      continue;
    }
    
    for(Float_t x=0.;x<=80.;x+=0.5) {
      for(Float_t E = beam_energy;E>0.;E-=deltaBeamE) {
	double range = projectile->CalcRange(beam_energy,E);
	double pointToWire = (wire<5) ? distanceToSecondWire-range :
	  distanceToFirstWire-range;
	double distance0 = sqrt(pointToWire*pointToWire+x*x);
	double distance1 = sqrt(wireHeight[wire]*wireHeight[wire]+
				distance0*distance0);

	double r = (distanceToSiDets-range)*distance1/pointToWire;
	double angle = acos((distanceToSiDets-range)/r)*180/M_PI;
	
	double protonEmissionEnergy = 4.*m1*m2/(m1+m2)/(m1+m2)*
	  (distanceToSiDets-range)/r*(distanceToSiDets-range)/r*E;
	double protonEnergy = proton->CalcRemainder(protonEmissionEnergy,r);

	if(protonEnergy<0.005) break;

	double cmEnergy = m2/(m1+m2)*E;

	LookupEntry entry = {range,r,angle,cmEnergy,protonEnergy};
	table[wire][x].push_back(entry);

	printf("wire: %d position: %f range: %f r: %f angle: %f CM Energy: %f Proton Energy: %f\n",wire,x,range,r,angle,cmEnergy,protonEnergy);
	fprintf(out,"%d %f %f %f %f %f %f\n",wire,x,range,r,angle,cmEnergy,protonEnergy);
      }
    }
  }
  fflush(out);
  fclose(out);
  
}

Float_t EnergyAngle::sum_pc_threshold;
Float_t EnergyAngle::si_threshold;
Float_t EnergyAngle::beam_energy;   
Float_t EnergyAngle::pressure;  
Float_t EnergyAngle::temperature;      
std::map<int,std::pair<float,float> > EnergyAngle::wire_offset;
std::map<int,std::pair<float,float> > EnergyAngle::wire_gain_diff;
std::map<int,std::pair<float,float> > EnergyAngle::wire_pos_cal;
EnergyLoss* EnergyAngle::projectile;
EnergyLoss* EnergyAngle::proton;
Float_t EnergyAngle::m1;
Float_t EnergyAngle::m2;
