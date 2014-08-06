void CalcSpectraRMatrix(TTree* tree,Float_t incoming) {
  gSystem->CompileMacro("EnergyLoss.C");
  gSystem->CompileMacro("EnergyAngle.C");
  EnergyAngle::ReadLookupTable();
  gSystem->CompileMacro("Spectra.C");
  Spectra::ReadPCBoundTable();
  Spectra::ReadSolidAngleTable();

  EnergyAngle t1(tree);
  t1.Loop();
  TFile * file = TFile::Open("energy_angle.root");
  TTree* energyAngle = (TTree*) file->Get("energyAngle");
  Spectra t2(energyAngle);
  t2.Loop(incoming,false,false);

  gSystem->Exec("./runRMatrix");
  gROOT->ProcessLine(".X DrawSpectra.C");
}
