#include "EnergyLoss.h"
#include <fstream>
#include <string>
#include "spline.h"

ClassImp(EnergyLoss)

EnergyLoss::EnergyLoss(const char* energyFile, double conversion, double xStepSize) :
TObject(), xStepSize_(xStepSize) {
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
  }
}

double EnergyLoss::CalcRemainder(double initialEnergy, double distance) {
  double h = xStepSize_;
  long long n = int(distance/xStepSize_);
  double beam_e = initialEnergy;
  tk::spline f;
  f.set_points(energy_,dEdx_);

  beam_e -= f(beam_e)*(h/3.0);
  for(long long j=1; j<n; j++ )
    {
      if(j%2==1){
	beam_e -= 4.0*f(beam_e)*(h/3.0);
      }
      else if(j%2==0){
	beam_e -= 2.0*f(beam_e)*(h/3.0);
      }
      if(beam_e - f(beam_e)*(h/3.0) <= 0.) break;
    }
  beam_e -= f(beam_e)*(h/3.0);

  return (beam_e < 0) ? 0. : beam_e;
}

double EnergyLoss::CalcLoss(double initialEnergy, double distance) {
  return initialEnergy-CalcRemainder(initialEnergy,distance);
}

double EnergyLoss::AddBack(double initialEnergy, double distance) {
  double h = xStepSize_;
  long long n = int(distance/xStepSize_);
  double beam_e = initialEnergy;
  tk::spline f;
  f.set_points(energy_,dEdx_);

  beam_e += f(beam_e)*(h/3.0);
  for(long long j=1; j<n; j++ )
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

double EnergyLoss::CalcRange(double initialEnergy, double remainder) {
  double h = xStepSize_;
  long long n = 1e9;
  double beam_e = initialEnergy;
  tk::spline f;
  f.set_points(energy_,dEdx_);

  beam_e -= f(beam_e)*(h/3.0); //step 1
  int num_steps = 1;
  for(long long j=1; j<n; j++ ) 
    {
      if(j%2==1){
	beam_e -= 4.0*f(beam_e)*(h/3.0);
      }
      else if(j%2==0){
	beam_e -= 2.0*f(beam_e)*(h/3.0);
      }
      num_steps++;
      if(beam_e - f(beam_e)*(h/3.0) <= remainder) 
	break;
    }
  beam_e -= f(beam_e)*(h/3.0);
  num_steps++;

  return xStepSize_*num_steps;
}
