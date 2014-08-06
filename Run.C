{
TFile *file = TFile::Open("carbon_triumf_09-13_t.root");
gSystem->CompileMacro("PosCal.C");
PosCal t(rawData);
t.Loop();
}
