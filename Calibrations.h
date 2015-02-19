#ifndef Calibrations_h
#define Calibrations_h

#include <TROOT.h>
#include <map>

class EnergyLoss;

class Calibrations {
 public:
  static void InitParameters();
  static void ScaleDensity(Float_t);
  static Float_t  MatchPCLeft(Float_t,Int_t,Int_t);
  static Float_t  MatchPCRight(Float_t,Int_t,Int_t);
  static Float_t  CalibrateSi(Float_t,Int_t,Int_t);
  static Float_t  CalcPosition(UChar_t,Float_t,Float_t);
  static Float_t  CalcPathNormalization(UChar_t,Float_t,Float_t);

  static Float_t m1;
  static Float_t m2;
  static Float_t beam_energy;
  static Float_t sum_pc_threshold;
  static Float_t si_threshold;
  static Float_t pressure;
  static Float_t density_offset;
  static Float_t density_slope;
  static Float_t density;
  static EnergyLoss* projectile;
  static EnergyLoss* proton;
  static EnergyLoss* proton_aluminum;
  static EnergyLoss* proton_silicon;
  static Float_t anode_to_si;
  static Float_t anode_sep;
  static Float_t window_to_si;
  
  static std::map<int,std::pair<float,float> > wire_offset_low;
  static std::map<int,std::pair<float,float> > wire_gain_diff_low;
  static std::map<int,std::pair<float,float> > wire_offset_high;
  static std::map<int,std::pair<float,float> > wire_gain_diff_high;
  static std::map<int,std::pair<float,float> > wire_pos_cal;
  static std::map<int,std::map<int,std::pair<float,float> > > si_cal;
  static std::map<int,std::map<int,std::pair<Float_t,Float_t> > > pc_run_gain;
  static std::vector<float> wireHeight;
};

#endif
