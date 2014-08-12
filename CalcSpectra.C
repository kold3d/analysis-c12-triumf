void CalcSpectra(TTree* tree,Float_t incoming) {
  gSystem->CompileMacro("EnergyLoss.C");
  gSystem->CompileMacro("Calibrations.C");
  Calibrations::InitParameters();
  gSystem->CompileMacro("EnergyAngle.C");
  EnergyAngle::ReadLookupTable();
  gSystem->CompileMacro("Spectra.C");
  Spectra::ReadPCBoundTable();
  Spectra::ReadSolidAngleTable();

  //EnergyAngle t1(tree);
  //t1.Loop();
  TFile * file = TFile::Open("energy_angle.root");
  TTree* energyAngle = (TTree*) file->Get("energyAngle");
  Spectra t2(energyAngle);
  t2.Loop(incoming,true,false);
}
