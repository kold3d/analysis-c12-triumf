#include "EnergyLoss.h"
#include <fstream>
#include <string>

ClassImp(EnergyLoss)

EnergyLoss::EnergyLoss(const char* energyFile, double conversion) :
TObject(){
  ifstream in(energyFile);

  double file_energy;
  double file_elec;
  double file_nucl;
  double file_range;
  double file_long;
  double file_lat;
  
  std::string energy_unit;
  std::string range_unit;
  std::string long_strag_unit;
  std::string lat_strag_unit;

  while(in >> file_energy >> energy_unit >> file_elec >> 
	file_nucl >> file_range >> range_unit >> file_long >> 
	long_strag_unit >> file_lat >> lat_strag_unit)  {
    if(energy_unit=="eV") {
      energy_.push_back(file_energy*(1.e-6));
    } else if(energy_unit=="keV"){
      energy_.push_back(file_energy*(1.e-3));
    } else if(energy_unit=="MeV"){
      energy_.push_back(file_energy);
    }
    dEdx_.push_back(file_elec*conversion + file_nucl*conversion);
    dEdx_no_norm_.push_back(file_elec+file_nucl);
  }
}

void EnergyLoss::SetConversion(double conversion) {
  for(int i = 0;i<dEdx_.size();i++) {
    dEdx_[i]=dEdx_no_norm_[i]*conversion;
  }
}

double EnergyLoss::CalcRemainder(double initialEnergy, double distance) {

  if(distance == 0.) return initialEnergy;
  if(initialEnergy<1e-6) return 0.;

  double beam_e = initialEnergy;
  double dist_init =0.0;
  tk::spline f;
  f.set_points(energy_,dEdx_);
  while(dist_init<=distance){
    if(f(beam_e)<=0.05){ // For Protons
      beam_e = CompSimpSub(f,beam_e,0.0,0.2,10);
      dist_init += 0.2;
    }
    else if(f(beam_e)<=0.10){ // For Carbon
      beam_e = CompSimpSub(f,beam_e,0.0,0.2,160);
      dist_init += 0.2;
    }
    else if(f(beam_e)<=0.2){
      beam_e = CompSimpSub(f,beam_e,0.0,0.15,160);
      dist_init += 0.15;
    }
    else if(f(beam_e)<=0.3){
      beam_e = CompSimpSub(f,beam_e,0.0,0.1,160);
      dist_init += 0.1;
    }
    else if(f(beam_e)<=1){
      beam_e = CompSimpSub(f,beam_e,0.0,0.05,160);
      dist_init += 0.05;
    }
    else{ // For the Havar
      beam_e = CompSimpSub(f,beam_e,0.0,0.00001,160);
      dist_init += 0.00001;
    }
  }

  return (beam_e < 0) ? 0. : beam_e;
}

double EnergyLoss::CalcLoss(double initialEnergy, double distance) {
  return initialEnergy-CalcRemainder(initialEnergy,distance);
}

double EnergyLoss::AddBack(double initialEnergy, double distance) {
  double dist_init =0.0;
  tk::spline f;
  f.set_points(energy_,dEdx_);
  while(dist_init<=distance){
    if(f(initialEnergy)<=0.05){ // For protons
      initialEnergy = CompSimpAdd(f,initialEnergy,0.0,0.2,10);
      dist_init += 0.2;
    }
    else if(f(initialEnergy)<=0.10){ // For Carbon
      initialEnergy = CompSimpAdd(f,initialEnergy,0.0,0.2,160);
      dist_init += 0.2;
    }
    else if(f(initialEnergy)<=0.2){
      initialEnergy = CompSimpAdd(f,initialEnergy,0.0,0.15,160);
      dist_init += 0.15;
    }
    else if(f(initialEnergy)<=0.3){
      initialEnergy = CompSimpAdd(f,initialEnergy,0.0,0.1,160);
      dist_init += 0.1;
    }
    else if(f(initialEnergy)<=1){
      initialEnergy = CompSimpAdd(f,initialEnergy,0.0,0.05,160);
      dist_init += 0.05;
    }
    else{ // For the Havar
      initialEnergy = CompSimpAdd(f,initialEnergy,0.0,0.00001,160);
      dist_init += 0.00001;
    }
  }

  return initialEnergy;
}

double EnergyLoss::CalcRange(double initialEnergy, double remainder) {

  if(initialEnergy == remainder) return 0;

  double beam_e = initialEnergy;
  double dist_init = 0.0;
  tk::spline f;
  f.set_points(energy_,dEdx_);

  while(beam_e>=remainder){
    if(f(beam_e)<=0.05){ // For protons
      beam_e = CompSimpSub(f,beam_e,0.0,0.2,10);
      dist_init += 0.2;
    }
    else if(f(beam_e)<=0.10){ // For Carbon
      beam_e = CompSimpSub(f,beam_e,0.0,0.2,160);
      dist_init += 0.2;
    }
    else if(f(beam_e)<=0.2){
      beam_e = CompSimpSub(f,beam_e,0.0,0.15,160);
      dist_init += 0.15;
    }
    else if(f(beam_e)<=0.3){
      beam_e = CompSimpSub(f,beam_e,0.0,0.1,160);
      dist_init += 0.1;
    }
    else if(f(beam_e)<=1){
      beam_e = CompSimpSub(f,beam_e,0.0,0.05,160);
      dist_init += 0.05;
    }
    else{ // For the Havar
      beam_e = CompSimpSub(f,beam_e,0.0,0.00001,160);
      dist_init += 0.00001;
    }
  }

  return dist_init;
}

double EnergyLoss::CompSimpSub(tk::spline f,double initialEnergy,double a, double b, int num_steps){
  double h = (b-a)/((double) num_steps);
  double beam_e = initialEnergy;
  //tk::spline f;
  //f.set_points(energy_,dEdx_);

  beam_e -= f(beam_e)*(h/3.0);
  for(int j=1; j<num_steps; j++ )
  {
    if(j%2==1){
      beam_e -= 4.0*f(beam_e)*(h/3.0);
    }
    else if(j%2==0){
      beam_e -= 2.0*f(beam_e)*(h/3.0);
    }
  }
  beam_e -= f(beam_e)*(h/3.0);
  return beam_e;
}

double EnergyLoss::CompSimpAdd(tk::spline f,double initialEnergy,double a, double b, int num_steps){
  double h = (b-a)/((double) num_steps);
  double beam_e = initialEnergy;
  //tk::spline f;
  //f.set_points(energy_,dEdx_);

  beam_e += f(beam_e)*(h/3.0);
  for(int j=1; j<num_steps; j++ )
  {
    if(j%2==1){
      beam_e += 4.0*f(beam_e)*(h/3.0);
    }
    else if(j%2==0){
      beam_e += 2.0*f(beam_e)*(h/3.0);
    }
  }
  beam_e += f(beam_e)*(h/3.0);
  return beam_e;
}
