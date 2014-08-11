{
  gSystem->CompileMacro("EnergyLoss.C");
  gSystem->CompileMacro("Calibrations.C");
  Calibrations::InitParameters();
  gSystem->CompileMacro("EnergyAngle.C");
  EnergyAngle::ReadLookupTable();
  
  EnergyAngle t(chain);
  t.Loop();
}
