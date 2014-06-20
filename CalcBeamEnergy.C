double CalcBeamEnergy(double measured) {
  gSystem->CompileMacro("EnergyLoss.C");

  EnergyLoss methane("dEdx_carbon_methane_290K_350torr.dat",
     3.1042e-2);  //conversion to MeV/mm

  EnergyLoss havar("dEdx_carbon_havar.dat", 
     8.2997E+02);
  
  double afterMethane = methane.AddBack(measured,490);
  //printf("AfterMethane: %f\n",afterMethane);
  double afterHavar = havar.AddBack(afterMethane,4e-3);

  //printf("Calculated Beam Energy: %f MeV\n",afterHavar);
  return afterHavar;
}
