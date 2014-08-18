#include "Calibrations.h"
#include "EnergyLoss.h"
#include <fstream>
#include <sstream>

void Calibrations::InitParameters() {
  m1               = 12.  ;   //AMU of projectile
  m2               = 1.   ;   //AMU of target
  beam_energy      = 41.62;   //In MeV, after havar window
  sum_pc_threshold = 0.   ;   //In Channels
  si_threshold     = 350. ;   //In KeV

  anode_to_si      = 28.5 ;   //Distance from Si detectors to back anode plane
  anode_sep        = 12.5 ;   //Distance from one anode to the next
  window_to_si     = 513. ;   //Distance from entrance window to si detectors

  Float_t pressure    = 785;  //In torr
  Float_t temperature = 295;  //In Kelvin, not used

  density = 9.95784e-05+7.20831e-07*pressure;  //From linear fit to energy dep
  projectile = new EnergyLoss("dEdx_carbon_methane_290K_400torr.dat",density*100.);
  proton     = new EnergyLoss("dEdx_proton_methane_290K_400torr.dat",density*100.);

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
  }
  in.close();

  //New si Energy calibration
  si_cal[0][0] = std::pair<Float_t,Float_t>(1.85660280260919,145.765962067280);
  si_cal[0][1] = std::pair<Float_t,Float_t>(1.90379631725431,129.419656625576);
  si_cal[0][2] = std::pair<Float_t,Float_t>(1.87854910481639,135.518280441151);
  si_cal[0][3] = std::pair<Float_t,Float_t>(1.92133095877702,157.401201770548);
  si_cal[1][0] = std::pair<Float_t,Float_t>(1.81731779514392,120.443559621123);
  si_cal[1][1] = std::pair<Float_t,Float_t>(1.90855965596577,119.093782917867);
  si_cal[1][2] = std::pair<Float_t,Float_t>(1.87882230561214,135.559631886550);
  si_cal[1][3] = std::pair<Float_t,Float_t>(1.89393626000964,129.916019810436);
  si_cal[2][0] = std::pair<Float_t,Float_t>(1.81833034732690,145.356790789386);
  si_cal[2][1] = std::pair<Float_t,Float_t>(1.86496555662152,132.646141278295);
  si_cal[2][2] = std::pair<Float_t,Float_t>(1.87315964529692,142.573229493874);
  si_cal[2][3] = std::pair<Float_t,Float_t>(1.88828167671430,155.291004367064);
}

Float_t Calibrations::CalibrateSi(Float_t ch,Int_t detector, Int_t quadrant) {
  std::pair<Float_t,Float_t> calib = si_cal[detector][quadrant];
  return ch*calib.first+calib.second;
}

Float_t Calibrations::MatchPCLeft(Float_t left, Int_t wire) {
  Float_t localLeft = (left<120) ? left*wire_gain_diff_low[wire].first+wire_offset_low[wire].first :
    left*wire_gain_diff_high[wire].first+wire_offset_high[wire].first;
  return localLeft;
}

Float_t Calibrations::MatchPCRight(Float_t right, Int_t wire) {
  Float_t localRight = (right<120) ? right*wire_gain_diff_low[wire].second+wire_offset_low[wire].second :
    right*wire_gain_diff_high[wire].second+wire_offset_high[wire].second;
  return localRight;
}

Float_t Calibrations::CalcPosition(UChar_t wire, Float_t left_ch, Float_t right_ch) {
  Float_t x = (right_ch-left_ch)/(right_ch+left_ch);
  
  Float_t pos = wire_pos_cal[wire].first*x+wire_pos_cal[wire].second;
  
  return pos;
}

Float_t Calibrations::m1;
Float_t Calibrations::m2;
Float_t Calibrations::beam_energy;
Float_t Calibrations::sum_pc_threshold;
Float_t Calibrations::si_threshold;
Float_t Calibrations::density;
EnergyLoss* Calibrations::projectile;
EnergyLoss* Calibrations::proton;
Float_t Calibrations::anode_to_si;
Float_t Calibrations::anode_sep;
Float_t Calibrations::window_to_si;
std::map<int,std::pair<float,float> > Calibrations::wire_offset_low;
std::map<int,std::pair<float,float> > Calibrations::wire_gain_diff_low;
std::map<int,std::pair<float,float> > Calibrations::wire_offset_high;
std::map<int,std::pair<float,float> > Calibrations::wire_gain_diff_high;
std::map<int,std::pair<float,float> > Calibrations::wire_pos_cal;
std::map<int,std::map<int,std::pair<float,float> > > Calibrations::si_cal;

