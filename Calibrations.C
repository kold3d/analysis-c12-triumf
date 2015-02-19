#include "Calibrations.h"
#include "EnergyLoss.h"
#include <fstream>
#include <sstream>

void Calibrations::InitParameters() {
  m1               = 12.   ;   //AMU of projectile
  m2               = 1.    ;   //AMU of target
  beam_energy      = 41.62 ;   //In MeV, after havar window
  sum_pc_threshold = 0.    ;   //In Channels
  si_threshold     = 350.  ;   //In KeV

  anode_to_si      = 28.5 ;   //Distance from Si detectors to back anode plane
  anode_sep        = 12.5 ;   //Distance from one anode to the next
  window_to_si     = 513. ;   //Distance from entrance window to si detectors

  pressure            = 805; //In torr
  Float_t temperature = 295;  //In Kelvin, not used

  density_offset = 9.95784e-05;
  density_slope  = 7.20831e-07;

  density = density_offset+density_slope*pressure;  //From linear fit to energy dep
  projectile = new EnergyLoss("dedx_carbon_methane.dat",density*100.);
  proton     = new EnergyLoss("dedx_proton_methane.dat",density*100.);
  proton_aluminum = new EnergyLoss("dedx_proton_aluminum.dat",2.3211E2);
  proton_silicon = new EnergyLoss("dedx_proton_silicon.dat",2.7019E2);
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

  si_cal[0][0] = std::pair<Float_t,Float_t>(1.9012190714,-127.8103431249);
  si_cal[0][1] = std::pair<Float_t,Float_t>(1.9495421417,-144.5380724261);
  si_cal[0][2] = std::pair<Float_t,Float_t>(1.9236917330,-138.3016685913);
  si_cal[0][3] = std::pair<Float_t,Float_t>(1.9674987331,-115.8856796064);
  si_cal[1][0] = std::pair<Float_t,Float_t>(1.8609848713,-153.7278107085);
  si_cal[1][1] = std::pair<Float_t,Float_t>(1.9544200208,-155.1122721806);
  si_cal[1][2] = std::pair<Float_t,Float_t>(1.9239685413,-138.2518396519);
  si_cal[1][3] = std::pair<Float_t,Float_t>(1.9394490857,-144.0396486743);
  si_cal[2][0] = std::pair<Float_t,Float_t>(1.8620282481,-128.2329020587);
  si_cal[2][1] = std::pair<Float_t,Float_t>(1.9097857150,-141.2529042250);
  si_cal[2][2] = std::pair<Float_t,Float_t>(1.9181770423,-131.0880331551);
  si_cal[2][3] = std::pair<Float_t,Float_t>(1.9336643926,-118.0693601288);

  //Read run by run corrections for PC
  in.open("wires_scaled_table.out");
  std::pair<Float_t,Float_t> sum[8];
  std::pair<Int_t,Int_t> entries[8];
  for(int i = 0;i<8;i++) {
    sum[i].first=0.;
    sum[i].second=0.;
    entries[i].first=0;
    entries[i].second=0;
  }
  while(!in.eof()) {
    std::string line;
    getline(in,line);
    if(in.eof()) continue;
    std::stringstream stm;
    stm.str(line);
    int run,wire;
    float left,right;
    if(stm >> run >> wire >> left >> right) {
      pc_run_gain[run][wire-1] = std::pair<Float_t,Float_t>(left,right);
      if(left>0.) {
	sum[wire-1].first += left;
	entries[wire-1].first++;
      } 
      if(right>0.) {
	sum[wire-1].second += right;
	entries[wire-1].second++;
      }
    }
  }
  in.close();
  //Normalize to the mean value
  for(std::map<int,std::map<int,std::pair<Float_t,Float_t> > >::iterator it1 = pc_run_gain.begin();
      it1 != pc_run_gain.end(); it1++) {
    for(std::map<int,std::pair<Float_t,Float_t> >::iterator it2 = it1->second.begin();
	it2 != it1->second.end(); it2++) {
      it2->second.first /= sum[it2->first].first/entries[it2->first].first;
      it2->second.second /= sum[it2->first].second/entries[it2->first].second;
      //printf("Run: %d Wire: %d Left: %f Right: %f\n",it1->first,it2->first+1,it2->second.first,it2->second.second);
    }
  }

  wireHeight.push_back(-20.32); 
  wireHeight.push_back(-10.16); 
  wireHeight.push_back(  0.00); 
  wireHeight.push_back( 10.16); 
  wireHeight.push_back( 20.32); 
  wireHeight.push_back(-15.24); 
  wireHeight.push_back(  0.00); 
  wireHeight.push_back( 15.25); 
}

Float_t Calibrations::CalibrateSi(Float_t ch,Int_t detector, Int_t quadrant) {
  std::pair<Float_t,Float_t> calib = si_cal[detector][quadrant];
  return ch*calib.first+calib.second;
}

Float_t Calibrations::MatchPCLeft(Float_t left, Int_t wire, Int_t run) {
  Float_t localLeft = (left<120) ? left*wire_gain_diff_low[wire].first+wire_offset_low[wire].first :
    left*wire_gain_diff_high[wire].first+wire_offset_high[wire].first;
  if(run>0&&pc_run_gain[run][wire].first>0.) localLeft /= pc_run_gain[run][wire].first;
  return localLeft;
}

Float_t Calibrations::MatchPCRight(Float_t right, Int_t wire, Int_t run) {
  Float_t localRight = (right<120) ? right*wire_gain_diff_low[wire].second+wire_offset_low[wire].second :
    right*wire_gain_diff_high[wire].second+wire_offset_high[wire].second;
  if(run>0&&pc_run_gain[run][wire].second>0.) localRight /= pc_run_gain[run][wire].second;
  return localRight;
}

void Calibrations::ScaleDensity(Float_t scaleFactor) {
  density = density_offset+density_slope*pressure;
  density *= scaleFactor;
  proton->SetConversion(density*100);
  projectile->SetConversion(density*100);
}

Float_t Calibrations::CalcPosition(UChar_t wire, Float_t left_ch, Float_t right_ch) {
  Float_t x = (right_ch-left_ch)/(right_ch+left_ch);
  
  Float_t pos = wire_pos_cal[wire].first*x+wire_pos_cal[wire].second;
  
  return pos;
}

Float_t Calibrations::CalcPathNormalization(UChar_t wire,Float_t wirePosition,Float_t depth) {
  Float_t frontPlaneDepth = (wire>5) ? 464.18 : 479.42;
  Float_t wireDepth = (wire>5) ? 471.8 : 484.5;
  Float_t distanceToFront = frontPlaneDepth - depth;
  Float_t distanceToWire = wireDepth-depth;
  Float_t d = distanceToFront/distanceToWire;
  Float_t pathLength = 2.*(1.-d)*sqrt(wirePosition*wirePosition+
				      wireHeight[wire-1]*wireHeight[wire-1]+
				      distanceToWire*distanceToWire);
  Float_t cellLength = (wire>5) ? 15.24 : 10.16;
  Float_t normalization = cellLength/pathLength;
  //printf("%f %f %f %f %f %f %f\n",wirePosition,depth,distanceToFront,distanceToWire,pathLength,cellLength,normalization);
  return normalization;
}

Float_t Calibrations::m1;
Float_t Calibrations::m2;
Float_t Calibrations::beam_energy;
Float_t Calibrations::sum_pc_threshold;
Float_t Calibrations::si_threshold;
Float_t Calibrations::pressure;
Float_t Calibrations::density_offset;
Float_t Calibrations::density_slope;
Float_t Calibrations::density;
EnergyLoss* Calibrations::projectile;
EnergyLoss* Calibrations::proton;
EnergyLoss* Calibrations::proton_aluminum;
EnergyLoss* Calibrations::proton_silicon;
Float_t Calibrations::anode_to_si;
Float_t Calibrations::anode_sep;
Float_t Calibrations::window_to_si;
std::map<int,std::pair<float,float> > Calibrations::wire_offset_low;
std::map<int,std::pair<float,float> > Calibrations::wire_gain_diff_low;
std::map<int,std::pair<float,float> > Calibrations::wire_offset_high;
std::map<int,std::pair<float,float> > Calibrations::wire_gain_diff_high;
std::map<int,std::pair<float,float> > Calibrations::wire_pos_cal;
std::map<int,std::map<int,std::pair<float,float> > > Calibrations::si_cal;
std::map<int,std::map<int,std::pair<Float_t,Float_t> > > Calibrations::pc_run_gain;
std::vector<float> Calibrations::wireHeight;
