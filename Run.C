{
TFile *file = TFile::Open("/data/lab03/grgroup/he8_triumf_0714/tree/carbon_triumf_09-13_t.root");
gSystem->CompileMacro("PosCal.C");
PosCal t(rawData);
t.Loop();
}
