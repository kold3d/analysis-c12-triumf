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

void Spectra::InitParameters() {
  beam_energy      = 79.76;  //MeV, after havar
  pressure         = 397.;   //In Torr
  temperature      = 290.;   //In Kelvin
  m1 = 12.;
  m2 =  1.;

  Float_t gasConstant = 8.3144621;
  Float_t torrInPa = 133.322368;
  Float_t molarMassMethane = 0.01604;

  density = pressure*torrInPa*molarMassMethane/
    gasConstant/temperature*0.001;
 
  proton = new EnergyLoss("dEdx_proton_methane_290K_400torr.dat",density*100.);
  carbon = new EnergyLoss("dEdx_carbon_methane_290K_400torr.dat",density*100.);
}

void Spectra::Loop(Int_t incoming)
{
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

  TH1F* h_no_pro = new TH1F("h_no_pro","h_no_pro",400,0,12000);
  
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(cm_energy[0] == 0.) continue;
      if(!proton_cut->IsInside(measured_energy,sum_dE[0])) {
	h_no_pro->Fill(measured_energy);
	continue;
      }

      std::pair<Float_t,Float_t> pc_bound = LookupPCBound(1,cm_energy[0]);
      if(detector == 2 && 
	 (wire[0] == 2  || wire[0] == 3 || wire[0] == 4) && 
	 fabs(position[0]) < pc_bound.second) s1->Fill(cm_energy[0]);

      pc_bound = LookupPCBound(2,cm_energy[0]);
      if(detector == 2 && 
	 (wire[0] == 1  || wire[0] == 5 || 
	  ((wire[0] == 2  || wire[0] == 3 || wire[0] == 4) && fabs(position[0]) > pc_bound.first)))
	 s2->Fill(cm_energy[0]);

      pc_bound = LookupPCBound(3,cm_energy[0]);
      if((detector == 3 || detector == 1) && 
	 fabs(position[0]) > pc_bound.first && fabs(position[0]) < pc_bound.second) 
	s3->Fill(cm_energy[0]);

      pc_bound = LookupPCBound(4,cm_energy[0]);
      if((detector == 3 || detector == 1) && 
	 fabs(position[0]) > pc_bound.first && fabs(position[0]) < pc_bound.second) 
	s4->Fill(cm_energy[0]);

      pc_bound = LookupPCBound(5,cm_energy[0]);
      if((detector == 3 || detector == 1) && 
	 fabs(position[0]) > pc_bound.first && fabs(position[0]) < pc_bound.second) 
	s5->Fill(cm_energy[0]);

      pc_bound = LookupPCBound(6,cm_energy[0]);
      if((detector == 3 || detector == 1) && 
	 fabs(position[0]) > pc_bound.first && fabs(position[0]) < pc_bound.second) 
	s6->Fill(cm_energy[0]);
   }
   
   TCanvas* c1 = new TCanvas();
   c1->Divide(3,2);
   c1->cd(1);
   DivideTargetThickness(s1);
   CalcSolidAngleFast(s1,1);
   s1->Draw();
   c1->cd(2);
   DivideTargetThickness(s2);
   CalcSolidAngleFast(s2,2);
   s2->Draw();
   c1->cd(3);
   DivideTargetThickness(s3);
   CalcSolidAngleFast(s3,3);
   s3->Draw();
   c1->cd(4);
   DivideTargetThickness(s4);
   CalcSolidAngleFast(s4,4);
   s4->Draw();
   c1->cd(5);
   DivideTargetThickness(s5);
   CalcSolidAngleFast(s5,5);
   s5->Draw();
   c1->cd(6);
   DivideTargetThickness(s6);
   CalcSolidAngleFast(s6,6);
   s6->Draw();
   
   s1->Scale(1./incoming);
   s2->Scale(1./incoming);
   s3->Scale(1./incoming);
   s4->Scale(1./incoming);
   s5->Scale(1./incoming);
   s6->Scale(1./incoming);

   TFile* spec_file = new TFile("spectra.root","recreate");
   s1->Write();
   s2->Write();
   s3->Write();
   s4->Write();
   s5->Write();
   s6->Write();
   spec_file->Close();
}

void Spectra::DivideTargetThickness(TH1F *f){    
  Int_t i_size = f->GetSize();
  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1;i<i_size-1;i++){
    Double_t binLowEdge = xaxis->GetBinLowEdge(i);
    Double_t binUpEdge = xaxis->GetBinUpEdge(i);
    binLowEdge *= (m1+m2)/m2;
    binUpEdge *= (m1+m2)/m2; // From C.M. to Lab Frame
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);
    Double_t delta_x = carbon->CalcRange(binUpEdge,binLowEdge);
    delta_x /= 10.0;
    Float_t molarMassMethane = 0.01604;
    Double_t factor = 4.e-27*density*delta_x*TMath::Na()/molarMassMethane;
    binContent /= factor;
    binError /= factor;
    f->SetBinContent(i,binContent);
    f->SetBinError(i,binError);
  }
}

void Spectra::EstimateSolidAngleNorm(TH1F* f, Int_t region) {
  Int_t i_size = f->GetSize();
  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1;i<i_size-1;i++){
    Double_t binCenter = xaxis->GetBinCenter(i);
    binCenter *= (m1+m2)/m2;
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);
    Double_t delta_x = carbon->CalcRange(beam_energy,binCenter);
    Double_t x=0.;
    if(region == 1) x = 0;
    else if(region == 2) x = 17.5;
    else if(region == 3) x = 42.21;
    else if(region == 4) x = 54.71;
    else if(region == 5) x = 67.21;
    else if(region == 6) x = 79.71;
    Double_t r2 = x*x+(513.-delta_x)*(513.-delta_x);
    binContent *= r2;
    binError *= r2;
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
    binCenter *=  (m1+m2)/m2;
    Double_t binContent = f->GetBinContent(i);
    if(binContent == 0.) continue;
    Double_t binError = f->GetBinError(i);
    Double_t depth = carbon->CalcRange(beam_energy,binCenter);
    Double_t z = 513.-depth;
    
    Float_t sum = 0.;
    if(region == 1) {
      for(Float_t dx = -25.;dx<25.;dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  std::pair<Int_t,Float_t> pc = CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,binCenter);
	  std::pair<Float_t,Float_t> pc_bound = LookupPCBound(1,binCenter*m2/(m1+m2));
	  //printf("%f %f %f %f\n", dx, dy, z, pc.second);
	  if((pc.first == 2 || pc.first ==3 || pc.first ==4) && fabs(pc.second) < pc_bound.second) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	    h->Fill(angle);
	    h2->Fill(180-2.*angle);
	  }
	}
      }
    } else if(region == 2) {
      for(Float_t dx = -25.;dx<25.;dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  std::pair<Int_t,Float_t> pc = CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,binCenter);
	  std::pair<Float_t,Float_t> pc_bound = LookupPCBound(2,binCenter*m2/(m1+m2));
	  if(pc.first == 1  || pc.first == 5 || 
	     ((pc.first == 2  || pc.first == 3 || pc.first == 4) && fabs(pc.second) > pc_bound.first)) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	    h->Fill(angle);
	    h2->Fill(180-2.*angle);
	  }
	}
      }
    } else if(region == 3) {
      for(Float_t dx = -85.96; dx < -35.96; dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  std::pair<Int_t,Float_t> pc = CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,binCenter);
	  std::pair<Float_t,Float_t> pc_bound = LookupPCBound(3,binCenter*m2/(m1+m2));
	  if(fabs(pc.second) > pc_bound.first && fabs(pc.second) < pc_bound.second) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	    h->Fill(angle);	  	    
	    h2->Fill(180-2.*angle);
	  }
	}
      }
      sum*=2.;
    } else if(region == 4) {
      for(Float_t dx = -85.96; dx < -35.96; dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  std::pair<Int_t,Float_t> pc = CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,binCenter);
	  std::pair<Float_t,Float_t> pc_bound = LookupPCBound(4,binCenter*m2/(m1+m2));
	  if(fabs(pc.second) > pc_bound.first && fabs(pc.second) < pc_bound.second) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	    h->Fill(angle); 
	    h2->Fill(180-2.*angle); 	    
 	  }
        }	
      }
      sum*=2.;
    } else if(region == 5) {
      for(Float_t dx = -85.96; dx < -35.96; dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  std::pair<Int_t,Float_t> pc = CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,binCenter);
	  std::pair<Float_t,Float_t> pc_bound = LookupPCBound(5,binCenter*m2/(m1+m2));
	  if(fabs(pc.second) > pc_bound.first && fabs(pc.second) < pc_bound.second) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
	    h->Fill(angle);
	    h2->Fill(180-2.*angle);
	  }
	}
      }
      sum*=2.;
    } else if(region == 6) {
      for(Float_t dx = -85.96; dx < -35.96; dx+=elementSize) {
	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
	  std::pair<Int_t,Float_t> pc = CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,binCenter);
	  std::pair<Float_t,Float_t> pc_bound = LookupPCBound(6,binCenter*m2/(m1+m2));
	  if(fabs(pc.second) > pc_bound.first && fabs(pc.second) < pc_bound.second) {
	    sum+=elementSize*elementSize*z/pow(dx*dx+dy*dy+z*z,1.5);
	    Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
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
      printf("Region: %d CM Energy: %f Solid Angle: %f, change_bin_content:%f\n",region,binCenter*m2/(m1+m2),sum,change_bin_content);
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
  TFile* file = new TFile(Form("angle_dists/region_%d.root",region),"recreate");
       
  Int_t i_size = f->GetSize();
  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1;i<i_size-1;i++){
    printf("Calculating solid angle for Region %d, Bin %d\n",region,i);

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
  file->Close();
}

std::pair<Int_t,Float_t> Spectra::CalcPCCell(Float_t x1, Float_t y1, Float_t depth, Float_t carbonEnergy) {		    
  Float_t SiToAnode = 28.5;
  Float_t cellWidth = 10.16;
  Float_t windowToSi = 513.;

  Float_t SiToBack  = SiToAnode-cellWidth/2.;
  Float_t SiToFront = SiToAnode+cellWidth/2.;
  Float_t z2 = windowToSi-depth;
  
  Boundary cell[5];
  cell[0] = Boundary( -25.4,-15.24);
  cell[1] = Boundary(-15.24,-5.08);
  cell[2] = Boundary( -5.08, 5.08);
  cell[3] = Boundary(  5.08,15.24);
  cell[4] = Boundary( 15.24,25.4 );

  Point3d siPoint = {x1,y1,0.};
  Point3d chamberPoint = {0.,0.,z2};
  Vec3d vec(siPoint,chamberPoint);
  
  Char_t entranceCell = -1;
  Point3d entrancePoint = vec.PointIntersectZPlane(SiToBack);
  for(UChar_t i = 0;i<5;i++) {
    if(cell[i].first < entrancePoint.y && entrancePoint.y <= cell[i].second) {
      entranceCell = i;
      break;
    }
  }

  Char_t exitCell = -1;
  Point3d exitPoint = vec.PointIntersectZPlane(SiToFront);
  for(UChar_t i = 0;i<5;i++) {
    if(cell[i].first < exitPoint.y && exitPoint.y <= cell[i].second) {
      exitCell = i;
      break;
    }
  }

  if(entranceCell < 0 || exitCell < 0 ) {
    printf("Entrance or exit cell not calculated.");
    return std::pair<Int_t,Float_t>(0.,0.);
  }

  Float_t protonInitialEnergy = 4.*m1*m2/(m1+m2)/(m1+m2)*vec.z*vec.z*carbonEnergy;
  Float_t energyAtPC = proton->CalcRemainder(protonInitialEnergy,Vec3d(chamberPoint,exitPoint).mag);

  std::pair<Int_t,Float_t> returnValues;
  if(entranceCell != exitCell) {
    Point3d p = entrancePoint;
    Char_t direction = (exitCell-entranceCell);
    std::vector<std::pair<UChar_t,Float_t> > lengths;
    for(UChar_t i = entranceCell;i!=exitCell;i+=direction) {
      Point3d boundaryPoint = (direction>0) ? vec.PointIntersectYPlane(cell[i].second) :
	vec.PointIntersectYPlane(cell[i].first);
      Vec3d tempVec(p,boundaryPoint);
      Float_t length = tempVec.mag;
      lengths.push_back(std::pair<UChar_t,Float_t>(i,length));
      p = boundaryPoint;
    }
    Vec3d tempVec(p,exitPoint);
    lengths.push_back(std::pair<UChar_t,Float_t>(exitCell,tempVec.mag));
    Float_t lastEnergy = energyAtPC;
    Float_t highestEnergy = 0.;
    UChar_t highestCell = 0;
    for(UChar_t i = lengths.size()-1;i>0;i--) {
      Float_t newEnergy = proton->CalcRemainder(lastEnergy,lengths[i].second);
      if((lastEnergy-newEnergy) > highestEnergy) {
	highestEnergy = lastEnergy - newEnergy;
	highestCell = lengths[i].first;
      }
      lastEnergy = newEnergy;
    }
    returnValues.first = highestCell+1;
  } else {
    returnValues.first = entranceCell+1;
  }
  Float_t crossing = vec.PointIntersectZPlane(SiToAnode).x;
  returnValues.second = crossing;
  
  return returnValues;
}

std::pair<Float_t,Float_t> Spectra::CalcPCBoundary(Int_t region, Float_t cmEnergy){
	cmEnergy *= (m1+m2)/m2;

	Double_t depth_bound = carbon->CalcRange(79.76,cmEnergy);
	Double_t z_bound = 513.-depth_bound;
	//Double_t dist_first_wire = 29.7;
	Double_t dist_second_wire = z_bound-28.5;

	Boundary wire_fixed_values[5];
	wire_fixed_values[0] = Boundary(0.,10.);
	wire_fixed_values[1] = Boundary(10.,25.);
	wire_fixed_values[2] = Boundary(35.96,48.46);
	wire_fixed_values[3] = Boundary(48.46,60.96);
	wire_fixed_values[4] = Boundary(60.96,73.46);
	wire_fixed_values[5] = Boundary(73.46,85.96);

	Boundary wire_float[5];
	for(Int_t i_wire = 0; i_wire<=5;i_wire++) {
	  wire_float[i_wire] = Boundary(wire_fixed_values[i_wire].first*dist_second_wire/z_bound,wire_fixed_values[i_wire].second*dist_second_wire/z_bound);
	}

	std::pair<Float_t,Float_t> returnBounds;
	returnBounds.first = wire_float[region-1].first;
	returnBounds.second = wire_float[region-1].second;

	return returnBounds;
}

std::pair<Float_t,Float_t> Spectra::LookupPCBound(Int_t region,Float_t cmEnergy){
  cmEnergy = floor(cmEnergy*100.0+0.5)/100.0; 
  if(region > 6 || region < 1){
    printf("Region=%d out of Range\n",region);
    return std::pair<Float_t,Float_t>(0.,0.);
  }

  Int_t mappedRegion = region;
  std::vector<PCBoundEntry> vec = pctable[mappedRegion];
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
      Int_t region;
      Float_t cmEnergy,LeftPCBound,RightPCBound;
      stm >> region >> cmEnergy >> LeftPCBound >> RightPCBound;
      PCBoundEntry entry = {cmEnergy,LeftPCBound,RightPCBound};
      pctable[region].push_back(entry);
    }
  }
}

void Spectra::CalcPCBoundTable(){
  InitParameters();
  Spectra sp;
  pctable.clear();
  FILE* out = fopen("pc_boundary_table.out","w");
  std::pair<Float_t,Float_t> pc_bound;
  for(Int_t region = 1; region <= 6; region++){
    for(Double_t cmEnergy = 0.0; cmEnergy <= 5.0; cmEnergy+=0.01){
      pc_bound = sp.CalcPCBoundary(region,cmEnergy);
      printf("Region: %d, CM Energy: %f, Left Boundary: %f, Right Boundary: %f\n",region, cmEnergy,pc_bound.first,pc_bound.second);
      fprintf(out, "%d %f %f %f\n",region, cmEnergy,pc_bound.first,pc_bound.second);
    }
  }
  fclose(out);
}

Float_t Spectra::LookupSolidAngle(Int_t region,Float_t cmEnergy){
  cmEnergy = floor(cmEnergy*100.0+0.5)/100.0; 
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
	InitParameters();
	Spectra sp;
	//satable.clear();
	FILE* out = fopen("solid_angle_table.out","w");
	Float_t elementSize = 1.;
	for(Int_t region = 1; region <= 6; region++){
		for(Double_t cmEnergy = 0.00; cmEnergy <= 5.0; cmEnergy+=0.1){
      TH1F* h = new TH1F(Form("region_%d_cmEnergy_%f_lab_ik",region,cmEnergy),Form("region_%d_cmEnergy_%f_lab_ik",region,cmEnergy),360,0,180);
      TH1F* h2 = new TH1F(Form("region_%d_cmEnergy_%f_lab_ik",region,cmEnergy),Form("region_%d_cmEnergy_%f_lab_ik",region,cmEnergy),360,0,180);
			Double_t depth = carbon->CalcRange(beam_energy,cmEnergy*(m1+m2)/m2);
		  Double_t z = 513.-depth;
		    
		  Float_t sum = 0.;
		  if(region == 1) {
        std::pair<Float_t,Float_t> pc_bound = sp.LookupPCBound(1,cmEnergy);
		  	for(Float_t dx = -25.;dx<25.;dx+=elementSize) {
			   	for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
			  		std::pair<Int_t,Float_t> pc = sp.CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,cmEnergy*(m1+m2)/m2);
						if((pc.first == 2 || pc.first ==3 || pc.first ==4) && fabs(pc.second) < pc_bound.second) {
              sum+=elementSize*elementSize*z/pow((dx*dx+dy*dy+z*z),1.5);
              Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
              h->Fill(angle);         
              h2->Fill(180-2.*angle);
			  		}
					}
		    	}
		    } else if(region == 2) {
          std::pair<Float_t,Float_t> pc_bound = sp.LookupPCBound(2,cmEnergy);
		      for(Float_t dx = -25.;dx<25.;dx+=elementSize) {
					 for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
			  		std::pair<Int_t,Float_t> pc = sp.CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,cmEnergy*(m1+m2)/m2);
			  		if(pc.first == 1  || pc.first == 5 || 
			     		((pc.first == 2  || pc.first == 3 || pc.first == 4) && fabs(pc.second) > pc_bound.first)) {
              sum+=elementSize*elementSize*z/pow((dx*dx+dy*dy+z*z),1.5);
              Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
              h->Fill(angle);         
              h2->Fill(180-2.*angle);
			  		}
					}
		      	}
		    } else if(region == 3) {
          std::pair<Float_t,Float_t> pc_bound = sp.LookupPCBound(3,cmEnergy);
		      for(Float_t dx = -85.96; dx < -35.96; dx+=elementSize) {
					 for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
			  			std::pair<Int_t,Float_t> pc = sp.CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,cmEnergy*(m1+m2)/m2);
			  			if(fabs(pc.second) > pc_bound.first && fabs(pc.second) < pc_bound.second) {
              sum+=elementSize*elementSize*z/pow((dx*dx+dy*dy+z*z),1.5);
              Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
              h->Fill(angle);         
              h2->Fill(180-2.*angle);
			  			}
						}
		      }
		      sum*=2.;
		    } else if(region == 4) {
          std::pair<Float_t,Float_t> pc_bound = sp.LookupPCBound(4,cmEnergy);
		      for(Float_t dx = -85.96; dx < -35.96; dx+=elementSize) {
					 for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
			  			std::pair<Int_t,Float_t> pc = sp.CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,cmEnergy*(m1+m2)/m2);
			  			if(fabs(pc.second) > pc_bound.first && fabs(pc.second) < pc_bound.second) {
			    		sum+=elementSize*elementSize*z/pow((dx*dx+dy*dy+z*z),1.5);
              Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
              h->Fill(angle);         
              h2->Fill(180-2.*angle);	    
		 	  			}
		        }	
		      }
		      sum*=2.;
		    } else if(region == 5) {
          std::pair<Float_t,Float_t> pc_bound = sp.LookupPCBound(5,cmEnergy);
		      for(Float_t dx = -85.96; dx < -35.96; dx+=elementSize) {
					 for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
			  		std::pair<Int_t,Float_t> pc = sp.CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,cmEnergy*(m1+m2)/m2);
			  		if(fabs(pc.second) > pc_bound.first && fabs(pc.second) < pc_bound.second) {
              sum+=elementSize*elementSize*z/pow((dx*dx+dy*dy+z*z),1.5);
              Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
              h->Fill(angle);         
              h2->Fill(180-2.*angle);
			  		}
					 }
		      }
		      sum*=2.;
		    } else if(region == 6) {
          std::pair<Float_t,Float_t> pc_bound = sp.LookupPCBound(6,cmEnergy);
		      for(Float_t dx = -85.96; dx < -35.96; dx+=elementSize) {
					 for(Float_t dy = -25.;dy<25.;dy+=elementSize) {
			  		std::pair<Int_t,Float_t> pc = sp.CalcPCCell(dx+elementSize/2.,dy+elementSize/2.,depth,cmEnergy*(m1+m2)/m2);
			  		if(fabs(pc.second) > pc_bound.first && fabs(pc.second) < pc_bound.second) {
              sum+=elementSize*elementSize*z/pow((dx*dx+dy*dy+z*z),1.5);
              Float_t angle = acos(z/sqrt(dx*dx+dy*dy+z*z))/3.14159*180;
              h->Fill(angle);         
              h2->Fill(180-2.*angle);
			  		}
				    }
		        }
		      sum*=2.;
		     }
         Float_t change_bin_content = 4.*sum*cos(h->GetMean()*3.14159/180);
		     printf("Region: %d, CM Energy: %f, SolidAngle: %f, change_bin_content:%f\n",region, cmEnergy,sum,change_bin_content);
      		fprintf(out, "%d %f %f %f\n",region, cmEnergy,sum,change_bin_content);
		}
	}
}

Float_t Spectra::pressure;
Float_t Spectra::temperature;
EnergyLoss* Spectra::proton;
EnergyLoss* Spectra::carbon;
Float_t Spectra::density;
Float_t Spectra::m1;
Float_t Spectra::m2;
Float_t Spectra::beam_energy;
PCBoundTable Spectra::pctable;
SolidAngleTable Spectra::satable;
