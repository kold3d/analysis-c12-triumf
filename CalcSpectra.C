void CalcSpectra(TTree* tree) {
  gSystem->CompileMacro("EnergyLoss.C");
  gSystem->CompileMacro("EnergyAngle.C");
  EnergyAngle::ReadLookupTable();
  gSystem->CompileMacro("Spectra.C");

  EnergyAngle t1(tree);
  t1.Loop();
  TFile * file = TFile::Open("energy_angle.root");
  TTree* energyAngle = (TTree*) file->Get("energyAngle");
  Spectra t2(energyAngle);
  t2.Loop();
}
