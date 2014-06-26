void CalcEnergyAngle(TTree* chain) {
  gSystem->CompileMacro("EnergyLoss.C");
  gSystem->CompileMacro("EnergyAngle.C");
  EnergyAngle::ReadLookupTable();

  EnergyAngle t(chain);
  t.Loop();
}
