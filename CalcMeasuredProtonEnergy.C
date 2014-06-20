double CalcMeasuredProtonEnergy(double beamEnergy, double angle, double cmEnergy) {
  gSystem->CompileMacro("EnergyLoss.C");

  EnergyLoss havar("dEdx_carbon_havar.dat",
		   8.2997E+02,
		   1e-7);
  EnergyLoss methane_carbon("dEdx_carbon_methane_290K_400torr.dat",
		     3.5484e-2,  //conversion to MeV/mm
		     2e-5);      //Step size in mm

  double labEnergy = cmEnergy*13;

  double afterHavar = havar.CalcRemainder(beamEnergy,4e-3);
  double rangeInMethane = methane_carbon.CalcRange(afterHavar,labEnergy);
  
  printf("Ion Range for CM Energy: %.2f mm\n",rangeInMethane);

  double protonDistance = (490.-rangeInMethane)*cos(3.14159/180.*angle);

  EnergyLoss methane_proton("dEdx_proton_methane_290K_400torr.dat",
			    3.5484e-2,  //conversion to MeV/mm
			    2e-5);      //Step size in mm
  double labEnergyProton = labEnergy*4.*12/13./13.*cos(3.14159/180.*angle)*cos(3.14159/180.*angle);

  printf("Proton Energy in Lab: %2f MeV\n",labEnergyProton);

  double energyInDetector = methane_proton.CalcRemainder(labEnergyProton,protonDistance);
  
  return energyInDetector;
}
