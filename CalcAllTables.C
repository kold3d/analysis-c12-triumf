{
  gSystem->CompileMacro("EnergyLoss.C");
  gSystem->CompileMacro("EnergyAngle.C");
  EnergyAngle::CalcLookupTable();
  gSystem->CompileMacro("Spectra.C");
  Spectra::CalcPCBoundTable();
  Spectra::CalcSolidAngleTable();
}
