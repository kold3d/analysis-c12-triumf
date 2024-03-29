#define EnergyAngle_cxx
#include "EnergyAngle.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "EnergyLoss.h"
#include "Calibrations.h"
#include "TMath.h"
#include <fstream>
#include <sstream>

LookupTable EnergyAngle::table;

void EnergyAngle::Loop(Int_t index_number)
{
  if (fChain == 0) return;

  Float_t measured_energy,measured_time,cm_energy[2],lab_angle[2],position[2],sum_dE[2],left_dE[2],right_dE[2];
  Int_t detector,quadrant,wire[2];

  TFile* file = new TFile(Form("energy_angle_%d.root",index_number),"recreate");

  TH2F* dE_E[8];
  for(Int_t i = 0;i<8;i++) {
    TString name = Form("dE_E_%d",i+1);
    dE_E[i] = new TH2F(name,name,400,0,12000,200,0,8000);
  }

  TTree* outTree = new TTree("energyAngle","Events tagged with reconstructed energy and angle.");
  outTree->Branch("measured_energy",&measured_energy,"measured_energy/F");
  outTree->Branch("measured_time",&measured_time,"measured_time/F");
  outTree->Branch("quadrant",&quadrant,"quadrant/I");
  outTree->Branch("detector",&detector,"detector/I");
  outTree->Branch("cm_energy",cm_energy,"cm_energy[2]/F");
  outTree->Branch("lab_angle",lab_angle,"lab_angle[2]/F");
  outTree->Branch("position",position,"position[2]/F");
  outTree->Branch("left_dE",left_dE,"left_dE[2]/F");
  outTree->Branch("right_dE",right_dE,"right_dE[2]/F");
  outTree->Branch("sum_dE",sum_dE,"sum_dE[2]/F");
  outTree->Branch("wire",wire,"wire[2]/I");
  outTree->Branch("rf_t",&rf_t,"rf_t/I");
  outTree->Branch("ic_ch_e",&ic_ch_e,"ic_ch_e/I");
  outTree->Branch("ic_trig",&ic_trig,"ic_trig/O");

  Int_t goodSi = 0, noPosition=0, noCMEnergy=0;
  Int_t numBelowPCThreshold= 0.;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   char thisFullPath[256] = "";
   int thisRunNumber = -1;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      TFile* currentFile = fChain->GetCurrentFile();
      if(strcmp(currentFile->GetName(),thisFullPath)!=0) {
	sprintf(thisFullPath,"%s",currentFile->GetName());
	std::string fullpath(thisFullPath);
	std::string runString = fullpath.substr(fullpath.length()-13,3);
	std::istringstream stm;
	stm.str(runString);
	stm >> thisRunNumber;
	printf("Processing run number %d...\n",thisRunNumber);	
      }      
      
      //Find legit hit in si_detector
      Int_t highSiEDet = -1;
      Int_t highSiEQuad = -1;
      Float_t highSiE = 0;
      Float_t highSiT = 0;
      Int_t mult = 0;
      for(Int_t i =0;i<si_mul;i++) {
	Float_t cal_e = Calibrations::CalibrateSi(si_ch_e[i],si_det[i]-1,si_quad[i]-1); 
	if(cal_e>Calibrations::si_threshold) mult++;
	if(cal_e>highSiE) {
	  highSiE = cal_e;
	  highSiT = si_ch_t[i];
	  highSiEDet = si_det[i];
	  highSiEQuad = si_quad[i];
	}
      }

      //if(mult>1) continue;

      //If Si energy is less than threshold, next event
      if(highSiE < Calibrations::si_threshold) continue;
      goodSi++;
      
      //Add proton information to tree
      measured_energy = highSiE;
      measured_time = highSiT;
      quadrant = highSiEQuad;
      detector = highSiEDet;

      //Find position in front and back wire
      Float_t positionFront = 0.;
      Float_t positionRear  = 0.;
      Bool_t goodRearPosition = false;
      Bool_t goodFrontPosition = false;
      Float_t highEPCFront  = 0.;
      Float_t highEPCRear   = 0.;
      Float_t highEPCRearLeft = 0.;
      Float_t highEPCRearRight = 0.;
      Float_t highEPCFrontLeft = 0.;
      Float_t highEPCFrontRight = 0.;
      UChar_t wireFront     = 0;
      UChar_t wireRear      = 0;
      Bool_t wireAboveThreshold = false;
      for(Int_t i = 0;i<pc_mul;i++) {
	Float_t left = pc_ch_left_e[i];
	Float_t right = pc_ch_right_e[i];
	left = Calibrations::MatchPCLeft(left,pc_wire[i]-1,thisRunNumber);
	right = Calibrations::MatchPCRight(right,pc_wire[i]-1,thisRunNumber);
	Float_t sum = left+right;
	if( sum < Calibrations::sum_pc_threshold ) {
	  continue;
	}
	wireAboveThreshold = true;
	if(pc_wire[i] < 6) {
	  if(sum > highEPCRear) {
	    highEPCRear = sum;
	    highEPCRearLeft = left;
	    highEPCRearRight = right;
	    positionRear = Calibrations::CalcPosition(pc_wire[i]-1,left,right);
	    wireRear = pc_wire[i];
	    goodRearPosition = true;
	  }
	} else {
	  if(sum > highEPCFront) {
	    highEPCFront = sum;
	    highEPCFrontLeft = left;
	    highEPCFrontRight = right;
	    positionFront = Calibrations::CalcPosition(pc_wire[i]-1,left,right);
	    wireFront = pc_wire[i];
	    goodFrontPosition = true;
	  }
	}
      }
      
      if(!wireAboveThreshold) numBelowPCThreshold++;

      if(!goodRearPosition&&!goodFrontPosition) {
	noPosition++;
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
      if(goodRearPosition && fabs(positionRear)<=90.) {
	LookupResult entry = LookupCMEnergyAngle(wireRear-1,positionRear,highSiE/1000.);
	if(entry.energy!=0.) {
	  angleRear = entry.angle;
	  cmEnergyRear = entry.energy;
	  Float_t norm = Calibrations::CalcPathNormalization(wireRear,positionRear,entry.range);
	  highEPCRear *= norm;
	}
      }
      Float_t angleFront=0.,cmEnergyFront = 0.;
      if(goodFrontPosition && fabs(positionFront)<=90.) {
	LookupResult entry = LookupCMEnergyAngle(wireFront-1,positionFront,highSiE/1000.);
	if(entry.energy!=0.) {
	  angleFront = entry.angle;
	  cmEnergyFront = entry.energy;
	  Float_t norm = Calibrations::CalcPathNormalization(wireFront,positionFront,entry.range);
	  highEPCFront *= norm;
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
      left_dE[0] = highEPCRearLeft;
      left_dE[1] = highEPCFrontLeft;
      right_dE[0] = highEPCRearRight;
      right_dE[1] = highEPCFrontRight;
      position[0] = positionRear;
      position[1] = positionFront;

      outTree->Fill();
   }

   for(Int_t i = 0;i<8;i++) dE_E[i]->Write();
   outTree->Write();

   file->Close();

   printf("Total Events: %d Below PC Threshold: %d No Position: %d No CM Energy %d\n",goodSi,numBelowPCThreshold,noPosition,noCMEnergy);   
}

LookupResult EnergyAngle::LookupCMEnergyAngle(UChar_t wire, Float_t x, Float_t protonEnergy) {
  LookupResult result;

  Float_t newX = floor(10.*fabs(x)/5.)*5./10.;
  if(x-newX>=0.25) newX += 0.5;

  if(wire > 7 || wire < 0 || newX < 0. || newX > 90. ) {
    //printf("wire=%d x=%f out of Range.\n",wire,newX);
    result.energy = 0.;
    result.angle = 0.;
    result.range = 0.;
    return result;
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
  else {
    result.energy = vec[start].cmEnergy;
    result.angle = vec[start].angle;
    result.range = vec[start].range;
    return result;

  }

  Bool_t done = false;
  Float_t previousAngle = vec[start].angle;
  Float_t previousCMEnergy = vec[start].cmEnergy;
  Float_t previousProtonEnergy = vec[start].protonEnergy;
  Float_t previousRange = vec[start].range;

  Float_t foundAngle =0;
  Float_t foundCMEnergy = 0.;
  Float_t foundRange = 0.;

  Int_t i = start+direction;
  while(!done && i>0 && i<size) {
    Float_t thisAngle = vec[i].angle;
    Float_t thisCMEnergy = vec[i].cmEnergy;
    Float_t thisProtonEnergy = vec[i].protonEnergy;
    Float_t thisRange = vec[i].range;

    if((direction < 0 && thisProtonEnergy >=protonEnergy) ||
	      (direction >0 && thisProtonEnergy <= protonEnergy)) {
      Float_t energy_slope = (thisCMEnergy-previousCMEnergy)/(thisProtonEnergy-previousProtonEnergy);
      Float_t energy_intercept = thisCMEnergy-thisProtonEnergy*energy_slope;
      Float_t angle_slope = (thisAngle-previousAngle)/(thisProtonEnergy-previousProtonEnergy);
      Float_t angle_intercept = thisAngle-thisProtonEnergy*angle_slope;
      Float_t range_slope = (thisRange-previousRange)/(thisProtonEnergy-previousProtonEnergy);
      Float_t range_intercept = thisRange-thisProtonEnergy*range_slope;
      foundAngle = angle_slope*protonEnergy+angle_intercept;
      foundCMEnergy = energy_slope*protonEnergy+energy_intercept;
      foundRange = range_slope*protonEnergy+range_intercept;
      done = true;
    } else {
      previousAngle = thisAngle;
      previousCMEnergy = thisCMEnergy;  
      previousProtonEnergy = thisProtonEnergy;
      previousRange = thisRange;
      i += direction;
    }
  }

  if(!done) {
    //printf("Entry wire=%d x=%f energy=%f not found in table!\n",wire,newX,protonEnergy);
  }
  
  result.energy = foundCMEnergy;
  result.angle = foundAngle;
  result.range = foundRange;
  return result;
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
  Float_t deltaBeamE = 0.10;
  Float_t distanceToSiDets= Calibrations::window_to_si;
  Float_t distanceToSecondWire = distanceToSiDets-Calibrations::anode_to_si; 
  Float_t distanceToFirstWire = distanceToSecondWire-Calibrations::anode_sep;

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
    
    for(Float_t x=0.;x<=90.;x+=0.5) {
      for(Float_t E = Calibrations::beam_energy;E>0.;E-=deltaBeamE) {
	Float_t range = Calibrations::projectile->CalcRange(Calibrations::beam_energy,E);
	Float_t pointToWire = (wire<5) ? distanceToSecondWire-range :
	  distanceToFirstWire-range;
	Float_t distance0 = sqrt(pointToWire*pointToWire+x*x);
	Float_t distance1 = sqrt(wireHeight[wire]*wireHeight[wire]+
				distance0*distance0);

	Float_t r = (distanceToSiDets-range)*distance1/pointToWire;
	Float_t angle = acos((distanceToSiDets-range)/r)*180/M_PI;
	
	Float_t protonEmissionEnergy = 4.*Calibrations::m1*Calibrations::m2/(Calibrations::m1+Calibrations::m2)/
	  (Calibrations::m1+Calibrations::m2)*(distanceToSiDets-range)/r*(distanceToSiDets-range)/r*E;
	Float_t protonEnergy = Calibrations::proton->CalcRemainder(protonEmissionEnergy,r);
        protonEnergy = Calibrations::proton_aluminum->CalcRemainder(protonEnergy,0.0005/cos(angle*M_PI/180.0));
        protonEnergy = Calibrations::proton_silicon->CalcRemainder(protonEnergy,0.0005/cos(angle*M_PI/180.0));

	if(protonEnergy<0.005) break;

	Float_t cmEnergy = Calibrations::m2/(Calibrations::m1+Calibrations::m2)*E;

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

