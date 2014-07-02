double CalcProtonEnergyLoss(double energy,double distance){
	gSystem->CompileMacro("EnergyLoss.C");

	EnergyLoss proton("dEdx_proton_methane_290K_400torr.dat",
     3.5484E-2);  //conversion to MeV/mm

	double afterMethane = proton.CalcRemainder(energy,distance);

	return afterMethane;
}