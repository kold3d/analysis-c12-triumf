double CalcResonanceDistance(double beam_energy) {
	gSystem->CompileMacro("EnergyLoss.C");

	EnergyLoss methane("dEdx_carbon_methane_290K_400torr.dat",
     3.5484E-2);  //conversion to MeV/mm

	EnergyLoss havar("dEdx_carbon_havar.dat", 
     8.2997E+02);

	double afterHavar = havar.CalcRemainder(beam_energy,4e-3);

	double first_resonance = methane.CalcRange(afterHavar,15);

	double second_resonance = methane.CalcRange(afterHavar,20.26);

	double third_resonance = methane.CalcRange(afterHavar,20.85);

	return first_resonance;
}