#define EnergyAngle_cxx
#include "EnergyAngle.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "EnergyLoss.h"
#include "TMath.h"

void EnergyAngle::InitParameters() {
  sum_pc_threshold = 110;  //In Channels
  si_threshold     = 300.; //In KeV
  beam_energy      = 79.76 ; //In MeV, after window
  pressure         = 397.; //In Torr
  temperature      = 290.; //In Kelvin

  //Noise peaks
  wire_offset[0] = std::pair<float,float>(0.0,0.0);
  wire_offset[1] = std::pair<float,float>(39.6507,38.5129);
  wire_offset[2] = std::pair<float,float>(39.0006,39.8341);
  wire_offset[3] = std::pair<float,float>(39.2025,38.8358);
  wire_offset[4] = std::pair<float,float>(39.7448,41.139);
  wire_offset[5] = std::pair<float,float>(39.6441,39.1638);
  wire_offset[6] = std::pair<float,float>(39.9842,40.8464);
  wire_offset[7] = std::pair<float,float>(39.6582,38.1933);

  //Gain matches from pulser (ratio left to right);
   wire_gain_diff[0] = 1.00636;
   wire_gain_diff[1] = 0.95856;
   wire_gain_diff[2] = 0.95023;
   wire_gain_diff[3] = 0.99627;
   wire_gain_diff[4] = 0.97248;
   wire_gain_diff[5] = 1.04135;
   wire_gain_diff[6] = 1.04280;
   wire_gain_diff[7] = 1.01840;

   //Position calibration (slope,intercept)
   wire_pos_cal[0] = std::pair<Float_t,Float_t>( 0.000, 0.000);
   wire_pos_cal[1] = std::pair<Float_t,Float_t>(90.150, 4.367);
   wire_pos_cal[2] = std::pair<Float_t,Float_t>(86.270, 5.192);
   wire_pos_cal[3] = std::pair<Float_t,Float_t>(87.699, 1.525);
   wire_pos_cal[4] = std::pair<Float_t,Float_t>(82.029, 0.407);
   wire_pos_cal[4] = std::pair<Float_t,Float_t>(82.657, 1.238);
   wire_pos_cal[5] = std::pair<Float_t,Float_t>(82.890,-3.905);
   wire_pos_cal[7] = std::pair<Float_t,Float_t>(87.144, 1.330);
}

void EnergyAngle::Loop()
{
   if (fChain == 0) return;

   Int_t totalCounts=0,numFrontPosition=0,numRearPosition=0;

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
	if(si_cal[i]>highSiE) {
	  highSiE = si_cal[i];
	  highSiEDet = si_det[i];
	  highSiEQuad = si_quad[i];
	}
      }

      //If Si energy is less than threshold, next event
      if(highSiE < si_threshold) continue;

      //Find position in front and back wire
      Float_t positionFront = 0.;
      Float_t positionRear  = 0.;
      Int_t highEPCFront    = 0.;
      Int_t highEPCRear     = 0.;
      UChar_t wireFront     = 0;
      UChar_t wireRear      = 0;
      for(Int_t i = 0;i<pc_mul;i++) {
	Int_t sum = pc_ch_left[i]+pc_ch_right[i];
	if( sum < sum_pc_threshold ) continue;
	if(pc_wire[i] < 6) {
	  if(sum > highEPCRear) {
	    highEPCRear = sum;
	    positionRear = CalcPosition(pc_wire[i]-1,pc_ch_left[i],pc_ch_right[i]);
	    wireRear = pc_wire[i];
	  }
	} else {
	  if(sum > highEPCFront) {
	    highEPCFront = sum;
	    positionFront = CalcPosition(pc_wire[i]-1,pc_ch_left[i],pc_ch_right[i]);
	    wireFront = pc_wire[i];
	  }
	}
      }
      
      if(positionFront != 0.) numFrontPosition++;
      if(positionRear != 0.) numRearPosition++;
      totalCounts++;
   }

   //printf("0:%d 1:%d 2:%d\n",numZeroPosition,numOnePosition,numTwoPosition);

   
}


Float_t EnergyAngle::CalcPosition(UChar_t wire, Int_t left_ch, Int_t right_ch) {
  Int_t left = left_ch-wire_offset[wire].first;
  Int_t right = right_ch-wire_offset[wire].second;
  right *= wire_gain_diff[wire];
  
  Float_t x = (right-left)/(right+left);
  
  Float_t pos = wire_pos_cal[wire].first*x+wire_pos_cal[wire].second;
  
  return pos;
}

void EnergyAngle::CalcLookupTable() {
  Float_t gasConstant = 8.3144621;
  Float_t torrInPa = 133.322368;
  Float_t molarMassMethane = 0.01604;

  Float_t density = pressure*torrInPa*molarMassMethane/
    gasConstant/temperature*0.001;

  printf("density: %f\n",density);

  EnergyLoss  carbon("dEdx_carbon_methane_290K_400torr.dat",density*100.);
  EnergyLoss  proton("dEdx_proton_methane_290K_400torr.dat",density*100.);

  Float_t deltaBeamE = 0.10;
  Float_t distanceToFirstWire = 472.;
  Float_t distanceToSecondWire = 484.5; 
  Float_t distanceToSiDets= 513.;
  Float_t m1 = 12.;
  Float_t m2 = 1.;

  Float_t wireHeight[8] = {
    -20.32,
    -10.16,
    0.0,
    10.16l,
    20.32,
    -15.24,
    0.00,
    15.24
  };
  
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
	double range = carbon.CalcRange(beam_energy,E);
	double pointToWire = (wire<5) ? distanceToSecondWire-range :
	  distanceToFirstWire-range;
	double distance0 = sqrt(pointToWire*pointToWire+x*x);
	double distance1 = sqrt(wireHeight[wire]*wireHeight[wire]+
				distance0*distance0);

	double r = (distanceToSiDets-range)*distance1/pointToWire;
	double angle = acos((distanceToSiDets-range)/r)*180/M_PI;
	
	double protonEmissionEnergy = 4.*m1*m2/(m1+m2)/(m1+m2)*
	  (distanceToSiDets-range)/r*(distanceToSiDets-range)/r*E;
	double protonEnergy = proton.CalcRemainder(protonEmissionEnergy,r);

	if(protonEnergy<0.005) break;

	double cmEnergy = m2/(m1+m2)*E;

	LookupEntry entry = {range,r,angle,cmEnergy,protonEnergy};
	table[wire][x].push_back(entry);

	printf("wire: %d position: %f range: %f r: %f angle: %f CM Energy: %f Proton Energy: %f\n",wire,x,range,r,angle,cmEnergy,protonEnergy);
	fprintf(out,"%d %f %f %f %f %f %f\n",wire,x,range,r,angle,cmEnergy,protonEnergy);
      }
    }
  }
  fclose(out);
  
}
