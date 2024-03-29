#define Spectra_cxx
#include "Spectra.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include "EnergyLoss.h"
#include <fstream>
#include <sstream>
#include <TRandom.h>
#include "Calibrations.h"

#ifdef __CINT__
#pragma link C++ class std::vector<GoodEvent>+;
#endif

std::vector<GoodEvent> Spectra::Loop()
{    
  TFile* cutFile = TFile::Open("cuts.root");
  TCutG* cuts[8];
  TCutG* rf_cut = (TCutG*) cutFile->Get("RF");
  TCutG* ic_cut = (TCutG*) cutFile->Get("IC_CUT");
  for(int i = 0;i<8;i++) {
    cuts[i] = (TCutG*)cutFile->Get(Form("PROTONS_%d",i+1));
    
  }
  cutFile->Close();
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  std::vector<GoodEvent> good_events;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    Int_t which = 0;
    if(sum_dE[0]>0) which = 0;
    else if(sum_dE[1]>0) which = 1;
    else continue;
    
    if(cm_energy[which] == 0.) continue;
    if(!cuts[wire[which]-1]->IsInside(measured_energy,sum_dE[which])) {
      continue;
    }
    //if(!rf_cut->IsInside(measured_energy,rf_t)) continue;
    //if(!ic_cut->IsInside(ic_ch_e,cm_energy[which])) continue;
    //if(ic_ch_e<120||ic_ch_e>300) continue;
    GoodEvent thisEvent = {wire[which],quadrant,detector,position[which],cm_energy[which]};
    good_events.push_back(thisEvent);
  }
  return good_events;
}

void Spectra::DivideTargetThickness(TH1F *f){    
  Int_t i_size = f->GetSize();
  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1;i<i_size-1;i++){
    Double_t binLowEdge = xaxis->GetBinLowEdge(i);
    Double_t binUpEdge = xaxis->GetBinUpEdge(i);
    if(binLowEdge==0.) binLowEdge+=0.001;
    binLowEdge *= (Calibrations::m1+Calibrations::m2)/Calibrations::m2;
    binUpEdge *= (Calibrations::m1+Calibrations::m2)/Calibrations::m2; // From C.M. to Lab Frame
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);
    Double_t delta_x = Calibrations::projectile->CalcRange(binUpEdge,binLowEdge);
    delta_x /= 10.0;
    Float_t molarMassMethane = 0.01604;
    Double_t factor = 4.e-27*Calibrations::density*delta_x*TMath::Na()/molarMassMethane;
    binContent /= factor;
    binError /= factor;
    f->SetBinContent(i,binContent);
    f->SetBinError(i,binError);
  }
}

void Spectra::CalcSolidAngleNorm(TH1F* f, Int_t region) {
  Float_t elementSize = 1.;

  TFile* file = new TFile(Form("angle_dists/region_%d.root",region),"recreate");
       
  Int_t i_size = f->GetSize();
  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1;i<i_size-1;i++){
    TH1F* h = new TH1F(Form("bin_%d_lab_ik",i+1),Form("bin_%d_lab_ik",i+1),360,0,180);
    TH1F* h2 = new TH1F(Form("bin_%d_cm_fk",i+1),Form("bin_%d_cm_fk",i+1),360,0,180);
    printf("Calculating solid angle for Region %d, Bin %d\n",region,i);

    Double_t binCenter = xaxis->GetBinCenter(i);
    Double_t cmEnergy = binCenter;
    binCenter *=  (Calibrations::m1+Calibrations::m2)/Calibrations::m2;
    Double_t binContent = f->GetBinContent(i);
    if(binContent == 0.) continue;
    Double_t binError = f->GetBinError(i);
    Double_t depth = Calibrations::projectile->CalcRange(Calibrations::beam_energy,binCenter);
    Double_t z = Calibrations::window_to_si-depth;
    
    Float_t sum = 0.;
    if(region == 1) {
      for(Float_t dx = 0.;dx<25.;dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	  Float_t protonEnergy = 4.*Calibrations::m1/(Calibrations::m1+Calibrations::m2)*z*z/(dx*dx+dy*dy+z*z)*cmEnergy;
	  Float_t range = Calibrations::proton->CalcRange(protonEnergy,0);
	  if(range>=sqrt(dx*dx+dy*dy+z*z)) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    h->Fill(angle);
	    h2->Fill(180-2.*angle);
	  }
	}
      }
      sum*=2.;
    } else if(region == 2) {
      for(Float_t dx = 35.96; dx< 60.96; dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	  Float_t protonEnergy = 4.*Calibrations::m1/(Calibrations::m1+Calibrations::m2)*z*z/(dx*dx+dy*dy+z*z)*cmEnergy;
	  Float_t range = Calibrations::proton->CalcRange(protonEnergy,0);
	  if(range>=sqrt(dx*dx+dy*dy+z*z)) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    h->Fill(angle);
	    h2->Fill(180-2.*angle);
	  }
	}
      }
      sum*=2.;
    } else if(region == 3) {
      for(Float_t dx = 60.96; dx < 85.96; dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	  Float_t protonEnergy = 4.*Calibrations::m1/(Calibrations::m1+Calibrations::m2)*z*z/(dx*dx+dy*dy+z*z)*cmEnergy;
	  Float_t range = Calibrations::proton->CalcRange(protonEnergy,0);
	  if(range>=sqrt(dx*dx+dy*dy+z*z)) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    h->Fill(angle);
	    h2->Fill(180-2.*angle);
	  }
        }	
      }
      sum*=2.;
    } 

    if(sum == 0.) {
      printf("Solid angle was zero!\n");
      f->SetBinContent(i,0.);
      f->SetBinError(i,0.);
    } else {
      Float_t change_bin_content = 4.*sum*cos(h->GetMean()*3.14159/180);
      printf("Region: %d CM Energy: %f Solid Angle: %f, change_bin_content:%f\n",region,binCenter*Calibrations::m2/(Calibrations::m1+Calibrations::m2),sum,change_bin_content);
      binContent /= 4.*sum*cos(h->GetMean()*3.14159/180);
      binError /= 4.*sum*cos(h->GetMean()*3.14159/180);
      f->SetBinContent(i,binContent);
      f->SetBinError(i,binError);
    }
    h->Write();
    h2->Write();
  }
  file->Close();
}

void Spectra::CalcSolidAngleFast(TH1F* f, Int_t region) {
  Int_t i_size = f->GetSize();
  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1;i<i_size-1;i++){
    //printf("Calculating solid angle for Region %d, Bin %d\n",region,i);

    Double_t binCenter = xaxis->GetBinCenter(i);
    Double_t binContent = f->GetBinContent(i);
    if(binContent == 0.) continue;
    Double_t binError = f->GetBinError(i);
    
    Double_t change_bin_content = LookupSolidAngle(region,binCenter); 

    if(change_bin_content == 0.) {
      printf("Solid angle was zero!\n");
      f->SetBinContent(i,0.);
      f->SetBinError(i,0.);
    } else {
      binContent /= change_bin_content;
      binError /= change_bin_content;
      f->SetBinContent(i,binContent);
      f->SetBinError(i,binError);
    }
  }
}

std::pair<Float_t,Float_t> Spectra::CalcPCBoundary(Int_t which, Int_t region, Float_t cmEnergy){
	cmEnergy *= (Calibrations::m1+Calibrations::m2)/Calibrations::m2;

	Double_t depth_bound = Calibrations::projectile->CalcRange(Calibrations::beam_energy,cmEnergy);
	Double_t z_bound = Calibrations::window_to_si-depth_bound;
	Double_t dist_wire = (which == 1) ? z_bound-Calibrations::anode_to_si-Calibrations::anode_sep : z_bound-Calibrations::anode_to_si;

	Boundary wire_fixed_values[6];
	wire_fixed_values[0] = Boundary(0.,6.);
	wire_fixed_values[1] = Boundary(6.,25.);
	wire_fixed_values[2] = Boundary(35.96,48.46);
	wire_fixed_values[3] = Boundary(48.46,60.96);
	wire_fixed_values[4] = Boundary(60.96,73.46);
	wire_fixed_values[5] = Boundary(73.46,85.96);

	std::pair<Float_t,Float_t> returnBounds;
	returnBounds.first = wire_fixed_values[region-1].first*dist_wire/z_bound;
	returnBounds.second = wire_fixed_values[region-1].second*dist_wire/z_bound;
	return returnBounds;
}

std::pair<Float_t,Float_t> Spectra::LookupPCBound(Int_t which, Int_t region,Float_t cmEnergy){
  cmEnergy = floor(cmEnergy*100.0+0.5)/100.0; 
  if(region > 6 || region < 1){
    printf("Region=%d out of Range\n",region);
    return std::pair<Float_t,Float_t>(0.,0.);
  }

  Int_t mappedRegion = region;
  std::vector<PCBoundEntry> vec = pctable[which][mappedRegion];
  Int_t size = vec.size();

  Int_t start = floor(size/2);
  Int_t direction = 0;
  if(vec[start].cmEnergy > cmEnergy) direction = -1;
  else if(vec[start].cmEnergy < cmEnergy) direction = 1;
  else return std::pair<Float_t,Float_t>(vec[start].LeftPCBound,vec[start].RightPCBound);

  Bool_t done = false;
  Float_t previousLeftPC = vec[start].LeftPCBound;
  Float_t previousRightPC = vec[start].RightPCBound;
  Float_t previousCMEnergy = vec[start].cmEnergy;

  Float_t foundLeftPC = 0.0;
  Float_t foundRightPC = 0.0;

  Int_t i = start+direction;
  while(!done && i>0 && i<size){
    Float_t thisLeftPC = vec[i].LeftPCBound;
    Float_t thisRightPC = vec[i].RightPCBound;
    Float_t thisCMEnergy = vec[i].cmEnergy;

    if((direction < 0 && thisCMEnergy <= cmEnergy) || 
      (direction > 0 && thisCMEnergy >= cmEnergy)){
      Float_t LeftPC_slope = (thisLeftPC-previousLeftPC)/(thisCMEnergy-previousCMEnergy);
      Float_t LeftPC_intercept = thisLeftPC-thisCMEnergy*LeftPC_slope;
      Float_t RightPC_slope = (thisRightPC-previousRightPC)/(thisCMEnergy-previousCMEnergy);
      Float_t RightPC_intercept = thisRightPC-thisCMEnergy*LeftPC_slope;
      foundLeftPC = LeftPC_slope*cmEnergy+LeftPC_intercept;
      foundRightPC = RightPC_slope*cmEnergy+RightPC_intercept;
      done = true;
    } else{
      previousLeftPC = thisLeftPC;
      previousRightPC = thisRightPC;
      previousCMEnergy = thisCMEnergy;
      i += direction;
    }

  }

  if(!done){
    printf("Region = %d CM Energy = %f Not found in boundary table!\n",region,cmEnergy);
  }

  return std::pair<Float_t,Float_t>(foundLeftPC,foundRightPC);
}

void Spectra::ReadPCBoundTable(){
  pctable.clear();
  std::ifstream in("pc_boundary_table.out");
  std::string line;
  while(!in.eof()){
    getline(in,line);
    if(!in.eof()){
      std::istringstream stm;
      stm.str(line);
      Int_t plane,region;
      Float_t cmEnergy,LeftPCBound,RightPCBound;
      stm >> plane >> region >> cmEnergy >> LeftPCBound >> RightPCBound;
      PCBoundEntry entry = {cmEnergy,LeftPCBound,RightPCBound};
      pctable[plane][region].push_back(entry);
    }
  }
}

void Spectra::CalcPCBoundTable(){
  Spectra sp;
  pctable.clear();
  FILE* out = fopen("pc_boundary_table.out","w");
  std::pair<Float_t,Float_t> pc_bound;
  for(Int_t plane = 0;plane<2;plane++) {
    for(Int_t region = 1; region <= 3; region++){
      for(Double_t cmEnergy = 0.001; cmEnergy <= 5.0; cmEnergy+=0.01){
	pc_bound = sp.CalcPCBoundary(plane,region,cmEnergy);
	printf("Plane: %d Region: %d, CM Energy: %f, Left Boundary: %f, Right Boundary: %f\n",plane,region, cmEnergy,pc_bound.first,pc_bound.second);
	fprintf(out, "%d %d %f %f %f\n",plane, region, cmEnergy,pc_bound.first,pc_bound.second);
      }
    }
  }
  fflush(out);
  fclose(out);
}

Float_t Spectra::LookupSolidAngle(Int_t region,Float_t cmEnergy){
  cmEnergy = floor(cmEnergy*10.0+0.5)/10.0; 
  if(region > 6 || region < 1){
    printf("Region=%d out of Range\n",region);
    return 0.;
  }

  Int_t mappedRegion = region;
  std::vector<SolidAngleEntry> vec = satable[mappedRegion];
  Int_t size = vec.size();

  Int_t start = floor(size/2);
  Int_t direction = 0;
  if(vec[start].cmEnergy > cmEnergy) direction = -1;
  else if(vec[start].cmEnergy < cmEnergy) direction = 1;
  else return vec[start].change_bin_content;

  Bool_t done = false;
  Float_t previouscontent = vec[start].change_bin_content;
  Float_t previousCMEnergy = vec[start].cmEnergy;

  Float_t foundcontent = 0.0;

  Int_t i = start+direction;
  while(!done && i>0 && i<size){
    Float_t thiscontent = vec[i].change_bin_content;
    Float_t thisCMEnergy = vec[i].cmEnergy;

    if((direction < 0 && thisCMEnergy <= cmEnergy) || 
      (direction > 0 && thisCMEnergy >= cmEnergy)){
      Float_t content_slope = (thiscontent-previouscontent)/(thisCMEnergy-previousCMEnergy);
      Float_t content_intercept = thiscontent-thisCMEnergy*content_slope;
      foundcontent = content_slope*cmEnergy+content_intercept;
      done = true;
    } else{
      previouscontent = thiscontent;
      previousCMEnergy = thisCMEnergy;
      i += direction;
    }

  }

  if(!done){
    printf("Region = %d CM Energy = %f Not found in boundary table!\n",region,cmEnergy);
  }

  return foundcontent;
}

void Spectra::ReadSolidAngleTable(){
  satable.clear();
  std::ifstream in("solid_angle_table.out");
  std::string line;
  while(!in.eof()){
    getline(in,line);
    if(!in.eof()){
      std::istringstream stm;
      stm.str(line);
      Int_t region;
      Float_t cmEnergy,solidangle,change_bin_content;
      stm >> region >> cmEnergy >> solidangle >> change_bin_content;
      SolidAngleEntry entry = {cmEnergy,solidangle,change_bin_content};
      satable[region].push_back(entry);
    }
  }
}

void Spectra::CalcSolidAngleTable(){
  Spectra sp;
  satable.clear();
  FILE* out = fopen("solid_angle_table.out","w");
  Float_t elementSize = 1.;
  for(Int_t region = 1; region <= 3; region++){
    for(Double_t cmEnergy = 0.001; cmEnergy <= 5.0; cmEnergy+=0.1){
      TH1F* h = new TH1F(Form("region_%d_cmEnergy_%f_lab_ik",region,cmEnergy),Form("region_%d_cmEnergy_%f_lab_ik",region,cmEnergy),360,0,180);
      Double_t depth =Calibrations::projectile->CalcRange(Calibrations::beam_energy,cmEnergy*(Calibrations::m1+Calibrations::m2)/Calibrations::m2);
      Double_t z = Calibrations::window_to_si-depth;

    Float_t sum = 0.;
    if(region == 1) {
      for(Float_t dx = 0.;dx<25.;dx+=elementSize) {
	for(Float_t dy = -25;dy<25.;dy+=elementSize) {
	  Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	  Float_t protonEnergy = 4.*Calibrations::m1/(Calibrations::m1+Calibrations::m2)*z*z/(dx*dx+dy*dy+z*z)*cmEnergy;
	  Float_t range = Calibrations::proton->CalcRange(protonEnergy,0);
	  if(range>=sqrt(dx*dx+dy*dy+z*z)) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    h->Fill(angle);
	  }
	}
      }
      sum*=2.;
    } else if(region == 2) {
      for(Float_t dx = 35.96;dx<60.96;dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	  Float_t protonEnergy = 4.*Calibrations::m1/(Calibrations::m1+Calibrations::m2)*z*z/(dx*dx+dy*dy+z*z)*cmEnergy;
	  Float_t range = Calibrations::proton->CalcRange(protonEnergy,0);
	  if(range>=sqrt(dx*dx+dy*dy+z*z)) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    h->Fill(angle);
	  }
	}
      }
      sum*=2.;
    } else if(region == 3) {
      for(Float_t dx = 60.96; dx< 85.96; dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	  Float_t protonEnergy = 4.*Calibrations::m1/(Calibrations::m1+Calibrations::m2)*z*z/(dx*dx+dy*dy+z*z)*cmEnergy;
	  Float_t range = Calibrations::proton->CalcRange(protonEnergy,0);
	  if(range>=sqrt(dx*dx+dy*dy+z*z)) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    h->Fill(angle);
	  }
	}
      }
      sum*=2.;
    }  

    Float_t change_bin_content = 4.*sum*cos(h->GetMean()*3.14159/180);
    printf("Region: %d, CM Energy: %f, SolidAngle: %f, change_bin_content:%f\n",region,cmEnergy,sum,change_bin_content);
    fprintf(out, "%d %f %f %f\n",region, cmEnergy,sum,change_bin_content);

    if(z>=Calibrations::window_to_si) break;

    }
  }
  fflush(out);
  fclose(out);
}

PCBoundTable Spectra::pctable;
SolidAngleTable Spectra::satable;
