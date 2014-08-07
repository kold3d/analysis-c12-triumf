#include "AZURE2/Config.h"
#include "AZURE2/CNuc.h"
#include "AZURE2/EPoint.h"
#include "AZURE2/AZUREParams.h"
#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TFile.h>

int main(int argc, const char** argv) {
  std::cout << "here" << std::endl;
  //Make config structure
  Config configure(std::cout);
  
  //Setup path to config file, read config options, and set parameter file
  configure.configfile = "12C+p.azr";
  configure.ReadConfigFile();
  configure.paramfile = configure.outputdir+"param.sav";

  //Set runtime options
  configure.paramMask |= Config::USE_BRUNE_FORMALISM;
  configure.paramMask |= Config::USE_GSL_COULOMB_FUNC;
  configure.paramMask |= Config::IGNORE_ZERO_WIDTHS;

  //Create new compound nucleus, fill, and initialize
  CNuc* compound = new CNuc;
  compound->Fill(configure);
  compound->GetPair(compound->GetPairNumFromKey(1))->SetEntrance();
  compound->Initialize(configure);
    
  //Create parameters object, initialize from compound, fill from file, and return values to compound
  AZUREParams params;
  compound->FillMnParams(params.GetMinuitParams());
  params.ReadUserParameters(configure);
  compound->FillCompoundFromParams(params.GetMinuitParams().Params());
  if(configure.paramMask & Config::USE_BRUNE_FORMALISM) compound->CalcShiftFunctions(configure);

  TFile* spec_file = new TFile(argv[1],"update");
  TH1F* spec = (TH1F*) spec_file->Get(argv[2]);
  TH1F* r_matrix = (TH1F*) spec->Clone(Form("%s_r_matrix",spec->GetName()));
  r_matrix->Reset();
  TFile* dist_file = new TFile(argv[3],"read");
  
  double sigma = 0.05;
  for(int i = 1; i <= spec->GetNbinsX() ; i++) {
    if(spec->GetBinContent(i) == 0.) continue;
    TH1F* dist = (TH1F*) dist_file->Get(Form("bin_%d_cm_fk",i));
    if(!dist) continue;
    std::cout << "Calculating R-Matrix For Bin " << i << " of " << argv[2] << std::endl;
    double sumNum = 0.;
    double sumDenom = 0.;
    double energy = spec->GetBinCenter(i);
    for(double dE = energy-4.*sigma; dE<=energy+4.*sigma;dE+=8.*sigma/20) {
      double sumNum2 = 0.;
      double sumDenom2 = 0.;
      for(int j = 1; j <= dist->GetNbinsX() ; j++) {
	double angle = dist->GetBinCenter(j);
	if(angle == 0.) continue;
	double value = dist->GetBinContent(j);
	EPoint *point = new EPoint(angle,dE,1,1,true,false,false,0.0,0,0);
	point->Initialize(compound,configure);
	point->Calculate(compound,configure);
	double crossSection=point->GetFitCrossSection();
	sumNum2 += crossSection*value;
	sumDenom2 += value;
	delete point;
      }
      sumNum += sumNum2/sumDenom2*exp(-(dE-energy)*(dE-energy)/2./sigma/sigma);
      sumDenom += exp(-(dE-energy)*(dE-energy)/2./sigma/sigma);
    }
    r_matrix->SetBinContent(i,sumNum/sumDenom);    
  }
  spec_file->cd();
  r_matrix->Write();
  spec_file->Close();
  dist_file->Close();
}
